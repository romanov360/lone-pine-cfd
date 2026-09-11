"""
Two-dimensional incompressible Navier-Stokes solver with Boussinesq buoyancy,
an immersed solid, and conjugate heat transfer.

Numerics
--------
* Staggered MAC grid: u on x-faces, v on y-faces, p and T at cell centres.
  Staggering is what keeps the pressure-velocity coupling tight; a collocated
  grid needs Rhie-Chow interpolation to avoid checkerboarding.
* Advection by a van Leer TVD scheme -- second order where the solution is
  smooth, monotone at fronts. Cell Reynolds numbers here reach O(50), far past
  where pure central differencing stops being stable.
* Explicit viscous and conductive terms (at the grid spacings used the
  diffusive step limit is well above the CFL limit, so implicit treatment
  would buy nothing).
* Chorin projection with the separable Poisson solver in poisson.py.
* The solid is imposed by Brinkman volume penalisation: a drag term -chi*u/eta
  is added inside the solid, so the velocity there decays towards zero over a
  time eta while the velocity field stays globally divergence-free. The
  alternative -- hard-zeroing the velocity on solid faces after the projection
  -- was tried first and is unusable here. The glass shell is only three cells
  thick, and re-imposing the mask after each projection puts divergence back
  in faster than a repeated projection can take it out; max|div| sat at O(1)
  instead of O(1e-14), which would have quietly destroyed energy conservation
  in the very region the whole study is about. Penalisation keeps the
  pressure operator a plain Laplacian, so the fast separable solver still
  applies, and it has a proper convergence theory (Angot, Bruneau & Fabrie,
  Numer. Math. 81, 1999) with an O(sqrt(eta)) boundary-layer error.
* Temperature is solved on a single domain spanning fluid and solid with
  spatially varying rho*cp and k, the interface conductivity being the
  harmonic mean of the two sides -- the standard conjugate treatment, which
  enforces flux continuity exactly.

Validation lives in validate.py: lid-driven cavity against Ghia et al. (1982),
the differentially heated cavity against de Vahl Davis (1983), and cylinder
cross-flow against the accepted Cd/St/Nu data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from poisson import FastPoisson

G = 9.80665


# ---------------------------------------------------------------------------
# TVD interpolation
# ---------------------------------------------------------------------------


def _van_leer(r: np.ndarray) -> np.ndarray:
    return (r + np.abs(r)) / (1.0 + np.abs(r))


def tvd_internal(q: np.ndarray, vel: np.ndarray, axis: int) -> np.ndarray:
    """
    Face values on the N-1 internal faces of an N-cell row, van Leer limited.

    `q` holds the cell values, `vel` the advecting velocity on those internal
    faces. Ghosts are filled by edge replication, which is first-order at the
    two cells adjacent to a boundary and irrelevant in the interior.
    """
    q = np.moveaxis(q, axis, 0)
    v = np.moveaxis(vel, axis, 0)
    n = q.shape[0]
    pad = [(2, 2)] + [(0, 0)] * (q.ndim - 1)
    qp = np.pad(q, pad, mode="edge")

    qm1, q0, q1, q2 = qp[1:n], qp[2:n + 1], qp[3:n + 2], qp[4:n + 3]
    d = q1 - q0
    ds = np.where(np.abs(d) < 1e-300, 1e-300, d)
    f_pos = q0 + 0.5 * _van_leer((q0 - qm1) / ds) * d
    f_neg = q1 - 0.5 * _van_leer((q2 - q1) / ds) * d
    out = np.where(v >= 0.0, f_pos, f_neg)
    return np.moveaxis(out, 0, axis)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


@dataclass
class Domain:
    nx: int
    ny: int
    Lx: float
    Ly: float

    @property
    def dx(self) -> float:
        return self.Lx / self.nx

    @property
    def dy(self) -> float:
        return self.Ly / self.ny

    def centers(self):
        x = (np.arange(self.nx) + 0.5) * self.dx
        y = (np.arange(self.ny) + 0.5) * self.dy
        return np.meshgrid(x, y, indexing="ij")


@dataclass
class BC:
    """
    Boundary conditions.

    left: 'wall' | 'inflow'
    right: 'wall' | 'outflow'
    top/bottom: 'wall' | 'slip'
    Thermal: 'adiabatic' | 'fixed'
    """

    left: str = "wall"
    right: str = "wall"
    top: str = "wall"          # 'wall' | 'slip' | 'lid'
    bottom: str = "wall"
    U_lid: float = 0.0
    U_in: float = 0.0
    T_in: float = 0.0
    T_left: str = "adiabatic"
    T_right: str = "adiabatic"
    T_top: str = "adiabatic"
    T_bottom: str = "adiabatic"
    T_wall_left: float = 0.0
    T_wall_right: float = 0.0
    inflow_profile: np.ndarray | None = None


@dataclass
class Materials:
    """Per-cell material fields (all shape (nx, ny))."""

    rho_cp: np.ndarray
    k: np.ndarray
    solid: np.ndarray          # bool, True where velocity is forced to zero
    nu: float                  # fluid kinematic viscosity
    rho_f: float               # fluid density (Boussinesq reference)
    beta: float                # fluid thermal expansion
    T_ref: float
    gravity: float = 0.0       # set to G to enable buoyancy


class Solver:
    def __init__(self, dom: Domain, mat: Materials, bc: BC, T0: np.ndarray,
                 cfl: float = 0.35, diff_safety: float = 0.25,
                 n_proj: int = 1, eta: float | None = None,
                 hard_mask: bool = False):
        self.d, self.m, self.bc = dom, mat, bc
        self.n_proj = n_proj
        self.eta = eta            # Brinkman penalisation time, s
        self.hard_mask = hard_mask
        nx, ny = dom.nx, dom.ny
        self.u = np.zeros((nx + 1, ny))
        self.v = np.zeros((nx, ny + 1))
        self.p = np.zeros((nx, ny))
        self.T = T0.astype(float).copy()
        self.t = 0.0
        self.step_count = 0
        self.cfl, self.diff_safety = cfl, diff_safety

        xbc = ("neumann", "dirichlet") if bc.right == "outflow" else ("neumann", "neumann")
        self.pois = FastPoisson(nx, ny, dom.dx, dom.dy, xbc)

        # Faces touching a solid cell.
        s = mat.solid
        self.u_solid = np.zeros((nx + 1, ny), dtype=bool)
        self.u_solid[1:-1] = s[:-1] | s[1:]
        self.u_solid[0] = s[0]
        self.u_solid[-1] = s[-1]
        self.v_solid = np.zeros((nx, ny + 1), dtype=bool)
        self.v_solid[:, 1:-1] = s[:, :-1] | s[:, 1:]
        self.v_solid[:, 0] = s[:, 0]
        self.v_solid[:, -1] = s[:, -1]

        # Harmonic-mean face conductivities (exact for a series interface).
        kx = 2.0 * mat.k[:-1] * mat.k[1:] / (mat.k[:-1] + mat.k[1:])
        ky = 2.0 * mat.k[:, :-1] * mat.k[:, 1:] / (mat.k[:, :-1] + mat.k[:, 1:])
        self.k_fx, self.k_fy = kx, ky

        if bc.left == "inflow":
            prof = bc.inflow_profile
            self.u_in = np.full(ny, bc.U_in) if prof is None else prof

        self._apply_velocity_bc()

    # -- boundaries ---------------------------------------------------------

    def _apply_velocity_bc(self):
        u, v, bc = self.u, self.v, self.bc
        if bc.left == "inflow":
            u[0, :] = self.u_in
        else:
            u[0, :] = 0.0
        if bc.right == "outflow":
            u[-1, :] = u[-2, :]
            # Enforce global mass conservation on the outflow plane.
            q_in = u[0, :].sum()
            q_out = u[-1, :].sum()
            if q_out != 0.0:
                u[-1, :] *= q_in / q_out
            else:
                u[-1, :] = q_in / u.shape[1]
        else:
            u[-1, :] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0
        if self.hard_mask:
            u[self.u_solid] = 0.0
            v[self.v_solid] = 0.0

    def _apply_domain_bc(self):
        """Outer-boundary conditions only; the solid is handled by the
        penalisation term, which must not be re-imposed after a projection."""
        self._apply_velocity_bc()

    def _u_ghost_y(self, u):
        """u padded by one cell in y with the wall condition applied."""
        g = np.empty((u.shape[0], u.shape[1] + 2))
        g[:, 1:-1] = u
        g[:, 0] = u[:, 0] if self.bc.bottom == "slip" else -u[:, 0]
        if self.bc.top == "slip":
            g[:, -1] = u[:, -1]
        elif self.bc.top == "lid":
            g[:, -1] = 2.0 * self.bc.U_lid - u[:, -1]
        else:
            g[:, -1] = -u[:, -1]
        return g

    def _v_ghost_x(self, v):
        g = np.empty((v.shape[0] + 2, v.shape[1]))
        g[1:-1] = v
        if self.bc.left == "inflow":
            g[0] = -v[0]               # no-slip on the inflow plane in v
        else:
            g[0] = -v[0]
        if self.bc.right == "outflow":
            g[-1] = v[-1]              # zero gradient
        else:
            g[-1] = -v[-1]
        return g

    # -- time step ----------------------------------------------------------

    def dt_max(self) -> float:
        d, m = self.d, self.m
        umax = max(np.abs(self.u).max(), 1e-9)
        vmax = max(np.abs(self.v).max(), 1e-9)
        dt_c = self.cfl / (umax / d.dx + vmax / d.dy)
        alpha_max = (m.k / m.rho_cp).max()
        dt_d = self.diff_safety / (max(m.nu, alpha_max) *
                                   (2.0 / d.dx**2 + 2.0 / d.dy**2))
        return min(dt_c, dt_d)

    def _momentum_rhs(self):
        d, m = self.d, self.m
        dx, dy = d.dx, d.dy
        u, v = self.u, self.v
        nx, ny = d.nx, d.ny

        # ---- u-momentum ----
        # d(uu)/dx : fluxes at cell centres (nx of them)
        u_adv_c = 0.5 * (u[:-1] + u[1:])                       # (nx, ny)
        uu_face = tvd_internal(u, u_adv_c, axis=0)             # (nx, ny)
        Fx = u_adv_c * uu_face
        duudx = np.zeros_like(u)
        duudx[1:-1] = (Fx[1:] - Fx[:-1]) / dx

        # d(uv)/dy : fluxes at corners
        vg = self._v_ghost_x(v)                                # (nx+2, ny+1)
        v_corner = 0.5 * (vg[:-1] + vg[1:])                    # (nx+1, ny+1)
        uv_int = tvd_internal(u, v_corner[:, 1:-1], axis=1)    # (nx+1, ny-1)
        Gy = np.zeros((nx + 1, ny + 1))
        Gy[:, 1:-1] = v_corner[:, 1:-1] * uv_int
        # Walls: v = 0 there so the flux vanishes; slip walls likewise.
        duvdy = (Gy[:, 1:] - Gy[:, :-1]) / dy

        ug = self._u_ghost_y(u)
        lap_u = np.zeros_like(u)
        lap_u[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dx**2
        lap_u += (ug[:, 2:] - 2 * ug[:, 1:-1] + ug[:, :-2]) / dy**2

        Ru = -(duudx + duvdy) + m.nu * lap_u

        # ---- v-momentum ----
        v_adv_c = 0.5 * (v[:, :-1] + v[:, 1:])                 # (nx, ny)
        vv_face = tvd_internal(v, v_adv_c, axis=1)
        Fy = v_adv_c * vv_face
        dvvdy = np.zeros_like(v)
        dvvdy[:, 1:-1] = (Fy[:, 1:] - Fy[:, :-1]) / dy

        ugy = self._u_ghost_y(u)                               # (nx+1, ny+2)
        u_corner = 0.5 * (ugy[:, :-1] + ugy[:, 1:])            # (nx+1, ny+1)
        vu_int = tvd_internal(v, u_corner[1:-1], axis=0)       # (nx-1, ny+1)
        Gx = np.zeros((nx + 1, ny + 1))
        Gx[1:-1] = u_corner[1:-1] * vu_int
        if self.bc.left == "inflow":
            Gx[0] = u_corner[0] * v[0]
        if self.bc.right == "outflow":
            Gx[-1] = u_corner[-1] * v[-1]
        dvudx = (Gx[1:] - Gx[:-1]) / dx

        vgx = self._v_ghost_x(v)
        lap_v = (vgx[2:] - 2 * vgx[1:-1] + vgx[:-2]) / dx**2
        lap_v_y = np.zeros_like(v)
        lap_v_y[:, 1:-1] = (v[:, 2:] - 2 * v[:, 1:-1] + v[:, :-2]) / dy**2
        lap_v = lap_v + lap_v_y

        Rv = -(dvvdy + dvudx) + m.nu * lap_v

        # ---- buoyancy (Boussinesq) ----
        if m.gravity != 0.0:
            T_face = np.zeros_like(v)
            T_face[:, 1:-1] = 0.5 * (self.T[:, :-1] + self.T[:, 1:])
            T_face[:, 0] = self.T[:, 0]
            T_face[:, -1] = self.T[:, -1]
            Rv = Rv + m.gravity * m.beta * (T_face - m.T_ref)

        return Ru, Rv

    def _advance_temperature(self, dt):
        d, m = self.d, self.m
        dx, dy = d.dx, d.dy
        T = self.T

        # Advection (face velocities already satisfy continuity).
        Tfx = np.zeros((d.nx + 1, d.ny))
        Tfx[1:-1] = tvd_internal(T, self.u[1:-1], axis=0)
        if self.bc.left == "inflow":
            Tfx[0] = self.bc.T_in
        Tfx[-1] = T[-1]                       # outflow / wall (u=0 there)
        Fx = self.u * Tfx

        Tfy = np.zeros((d.nx, d.ny + 1))
        Tfy[:, 1:-1] = tvd_internal(T, self.v[:, 1:-1], axis=1)
        Fy = self.v * Tfy

        adv = (Fx[1:] - Fx[:-1]) / dx + (Fy[:, 1:] - Fy[:, :-1]) / dy

        # Conduction with harmonic-mean face conductivity.
        qx = np.zeros((d.nx + 1, d.ny))
        qx[1:-1] = -self.k_fx * (T[1:] - T[:-1]) / dx
        if self.bc.T_left == "fixed":
            qx[0] = -m.k[0] * (T[0] - self.bc.T_wall_left) / (0.5 * dx)
        if self.bc.T_right == "fixed":
            qx[-1] = -m.k[-1] * (self.bc.T_wall_right - T[-1]) / (0.5 * dx)
        if self.bc.left == "inflow":
            qx[0] = 0.0                       # advection carries the inlet T
        qy = np.zeros((d.nx, d.ny + 1))
        qy[:, 1:-1] = -self.k_fy * (T[:, 1:] - T[:, :-1]) / dy

        cond = -((qx[1:] - qx[:-1]) / dx + (qy[:, 1:] - qy[:, :-1]) / dy)

        self.T = T + dt * (-adv * self._rho_cp_adv + cond) / m.rho_cp

    @property
    def _rho_cp_adv(self):
        # Advection only happens in the fluid, where rho*cp is the fluid value.
        return self.m.rho_cp

    def step(self, dt: float | None = None) -> float:
        dt = dt or self.dt_max()
        d, m = self.d, self.m

        Ru, Rv = self._momentum_rhs()
        us = self.u + dt * Ru
        vs = self.v + dt * Rv

        if self.bc.left == "inflow":
            us[0, :] = self.u_in
        else:
            us[0, :] = 0.0
        if self.bc.right == "outflow":
            us[-1, :] = us[-2, :]
            q_in, q_out = us[0, :].sum(), us[-1, :].sum()
            us[-1, :] = us[-1, :] * (q_in / q_out) if q_out != 0 else q_in / d.ny
        else:
            us[-1, :] = 0.0
        vs[:, 0] = 0.0
        vs[:, -1] = 0.0

        # Projection. Forcing the immersed-solid faces back to zero after a
        # projection reintroduces a little divergence next to the wall, so the
        # projection is repeated: each pass shrinks that residual while the
        # zero-flux condition on the solid stays exact. Two or three passes are
        # enough to push max|div| down by several orders of magnitude, which
        # matters here because the glass shell is only a few cells thick and a
        # leaky wall would short-circuit the whole conjugate problem.
        # Penalisation and projection do not commute, so they are iterated.
        # Damping the solid velocity BEFORE a single projection is not enough:
        # the projection then puts velocity straight back into the solid via
        # the pressure gradient, and at 2e-3 m/s that is 40 mm of drift in 20 s
        # through a 3 mm wall -- the bottle simply empties itself into the bath.
        # Alternating the two drives both residuals down together.
        self.u, self.v = us, vs
        phi_total = np.zeros((d.nx, d.ny))
        eta = self.eta if self.eta is not None else 0.02 * dt
        damp = 0.0 if self.hard_mask else 1.0 / (1.0 + dt / eta)
        for _ in range(self.n_proj):
            self.u[self.u_solid] *= damp
            self.v[self.v_solid] *= damp
            self._apply_domain_bc()
            div = (self.u[1:] - self.u[:-1]) / d.dx + \
                  (self.v[:, 1:] - self.v[:, :-1]) / d.dy
            phi = self.pois.solve(div / dt)
            phi_total += phi
            self.u[1:-1] -= dt * (phi[1:] - phi[:-1]) / d.dx
            self.v[:, 1:-1] -= dt * (phi[:, 1:] - phi[:, :-1]) / d.dy
            self._apply_domain_bc()
        self.p = phi_total
        self._advance_temperature(dt)

        self.t += dt
        self.step_count += 1
        return dt

    # -- diagnostics --------------------------------------------------------

    def divergence_norm(self) -> float:
        div = (self.u[1:] - self.u[:-1]) / self.d.dx + \
              (self.v[:, 1:] - self.v[:, :-1]) / self.d.dy
        return np.abs(div).max()

    def solid_energy(self) -> float:
        """Thermal energy of the immersed solid region, J/m (per unit depth)."""
        m, d = self.m, self.d
        cell = d.dx * d.dy
        return float((m.rho_cp * self.T * m.solid).sum() * cell)

    def interface_heat_rate(self, mask: np.ndarray | None = None):
        """
        Net conductive heat leaving `mask` through its boundary faces,
        W per metre of depth, plus the wetted perimeter.

        Defaults to the immersed solid. For the conjugate bottle the useful
        control volume is the WHOLE bottle -- glass plus the water inside it --
        because the flux across the glass's inner face is heat moving within
        the bottle, not heat lost to the bath.
        """
        d, m = self.d, self.m
        s = m.solid if mask is None else mask
        Q = 0.0
        per = 0.0
        # x-faces between a solid and a fluid cell
        fx = s[:-1] ^ s[1:]
        if fx.any():
            flux = -self.k_fx * (self.T[1:] - self.T[:-1]) / d.dx   # +x direction
            sgn = np.where(s[:-1], 1.0, -1.0)      # outward from the solid
            Q += float((flux * sgn * fx).sum() * d.dy)
            per += float(fx.sum() * d.dy)
        fy = s[:, :-1] ^ s[:, 1:]
        if fy.any():
            flux = -self.k_fy * (self.T[:, 1:] - self.T[:, :-1]) / d.dy
            sgn = np.where(s[:, :-1], 1.0, -1.0)
            Q += float((flux * sgn * fy).sum() * d.dx)
            per += float(fy.sum() * d.dx)
        return Q, per
