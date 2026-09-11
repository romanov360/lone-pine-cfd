"""
Verification and validation of the solver in nsolver.py.

Nothing downstream is worth reading unless these pass, so they are run first
and their output is kept in results/data/validation.json.

  1. Lid-driven cavity, Re = 100 / 400 / 1000, against Ghia, Ghia & Shin
     (J. Comput. Phys. 48, 1982) centreline velocities. Pure momentum.
  2. Differentially heated square cavity, Ra = 1e3..1e6, Pr = 0.71, against
     de Vahl Davis (Int. J. Numer. Meth. Fluids 3, 1983). Buoyancy + energy.
  3. Cylinder in cross-flow: steady recirculation length at Re = 40 against
     Coutanceau & Bouard (1977), vortex-shedding Strouhal number at Re = 100
     against Williamson (1989), and Nusselt number at Re = 20..200 against
     the Churchill-Bernstein correlation. Immersed boundary + convection.
"""

from __future__ import annotations

import json
import time

import numpy as np

from nsolver import BC, Domain, Materials, Solver

# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

GHIA_Y = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                   0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                   0.9688, 0.9766, 1.0000])
GHIA_U = {
    100: np.array([0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
                   -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151,
                   0.68717, 0.73722, 0.78871, 0.84123, 1.00000]),
    400: np.array([0.00000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299,
                   -0.32726, -0.17119, -0.11477, 0.02135, 0.16256, 0.29093,
                   0.55892, 0.61756, 0.68439, 0.75837, 1.00000]),
    1000: np.array([0.00000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289,
                    -0.27805, -0.10648, -0.06080, 0.05702, 0.18719, 0.33304,
                    0.46604, 0.51117, 0.57492, 0.65928, 1.00000]),
}

DVD_NU = {1e3: 1.118, 1e4: 2.243, 1e5: 4.519, 1e6: 8.800}


# ---------------------------------------------------------------------------
# 1. Lid-driven cavity
# ---------------------------------------------------------------------------


def lid_driven_cavity(Re: float, n: int = 128, t_end: float | None = None,
                      verbose: bool = True):
    L, U = 1.0, 1.0
    nu = U * L / Re
    dom = Domain(n, n, L, L)
    mat = Materials(rho_cp=np.ones((n, n)), k=np.full((n, n), 1e-9),
                    solid=np.zeros((n, n), bool), nu=nu, rho_f=1.0,
                    beta=0.0, T_ref=0.0, gravity=0.0)
    bc = BC(top="lid", U_lid=U)
    s = Solver(dom, mat, bc, np.zeros((n, n)), cfl=0.4)

    t_end = t_end or 40.0 * L / U
    last = None
    while s.t < t_end:
        s.step()
        if s.step_count % 2000 == 0:
            cur = s.u.copy()
            if last is not None:
                ch = np.abs(cur - last).max() / max(np.abs(cur).max(), 1e-12)
                if ch < 1e-6:
                    break
            last = cur

    # u on the vertical centreline: u lives on x-faces, take the face at x=L/2.
    i = n // 2
    y_u = (np.arange(n) + 0.5) / n
    u_line = s.u[i]
    u_interp = np.interp(GHIA_Y, np.r_[0.0, y_u, 1.0],
                         np.r_[0.0, u_line, U])
    err = np.abs(u_interp - GHIA_U[Re])
    rms = float(np.sqrt((err**2).mean()))
    if verbose:
        print(f"    Re={Re:5.0f}  n={n}  steps={s.step_count:6d}  "
              f"t={s.t:6.1f}  RMS(u - Ghia) = {rms:.4f}   "
              f"max div = {s.divergence_norm():.2e}")
    return {"Re": Re, "n": n, "rms": rms, "u": u_interp.tolist(),
            "ghia": GHIA_U[Re].tolist(), "y": GHIA_Y.tolist()}


# ---------------------------------------------------------------------------
# 2. Differentially heated cavity
# ---------------------------------------------------------------------------


def heated_cavity(Ra: float, n: int = 96, Pr: float = 0.71,
                  verbose: bool = True):
    L, dT = 1.0, 1.0
    alpha = 1.0
    nu = Pr * alpha
    g_beta = Ra * nu * alpha / L**3        # with dT = 1

    dom = Domain(n, n, L, L)
    mat = Materials(rho_cp=np.ones((n, n)), k=np.full((n, n), alpha),
                    solid=np.zeros((n, n), bool), nu=nu, rho_f=1.0,
                    beta=1.0, T_ref=0.5, gravity=g_beta)
    bc = BC(T_left="fixed", T_right="fixed", T_wall_left=1.0, T_wall_right=0.0)
    X, _ = dom.centers()
    s = Solver(dom, mat, bc, 1.0 - X / L, cfl=0.35)

    # Scale the run time by the convective turnover time.
    t_end = 6.0 * L / max(np.sqrt(g_beta * dT * L), 1e-9) * 40
    t_end = min(max(t_end, 2.0), 60.0)

    nu_hist = []
    while s.t < t_end:
        s.step()
        if s.step_count % 200 == 0:
            nu_hist.append(_cavity_nusselt(s, alpha, dT, L))
            if len(nu_hist) > 12:
                w = np.array(nu_hist[-12:])
                if (w.max() - w.min()) / w.mean() < 2e-4:
                    break

    Nu = _cavity_nusselt(s, alpha, dT, L)
    ref = DVD_NU[Ra]
    if verbose:
        print(f"    Ra={Ra:8.0e}  n={n}  steps={s.step_count:6d}  "
              f"Nu = {Nu:6.3f}   ref = {ref:6.3f}   "
              f"error = {(Nu/ref-1)*100:+6.2f} %")
    return {"Ra": Ra, "n": n, "Nu": float(Nu), "ref": ref,
            "error_pct": float((Nu / ref - 1) * 100)}


def _cavity_nusselt(s: Solver, alpha, dT, L):
    """Average Nusselt on the hot wall (k = alpha since rho*cp = 1)."""
    dx = s.d.dx
    q = alpha * (s.bc.T_wall_left - s.T[0]) / (0.5 * dx)
    return float(q.mean() * L / (alpha * dT))


# ---------------------------------------------------------------------------
# 3. Cylinder in cross-flow
# ---------------------------------------------------------------------------


def cylinder_mask(dom: Domain, xc, yc, D):
    X, Y = dom.centers()
    return ((X - xc) ** 2 + (Y - yc) ** 2) <= (D / 2.0) ** 2


def cylinder_flow(Re: float, D_cells: int = 40, Pr: float = 7.0,
                  heated: bool = True, t_end_D: float = 120.0,
                  Lx_D: float = 32.0, Ly_D: float = 16.0, verbose: bool = True):
    """
    Cylinder at Re in a channel. Returns recirculation length, Strouhal number
    and Nusselt number (isothermal cylinder, conjugate machinery switched off
    by making the "solid" a very high conductivity region held at T=1).
    """
    D = 1.0
    U = 1.0
    nu = U * D / Re
    dx = D / D_cells
    nx, ny = int(Lx_D * D / dx), int(Ly_D * D / dx)
    dom = Domain(nx, ny, nx * dx, ny * dx)
    xc, yc = 8.0 * D, dom.Ly / 2.0
    solid = cylinder_mask(dom, xc, yc, D)

    alpha = nu / Pr
    k = np.full((nx, ny), alpha)
    rho_cp = np.ones((nx, ny))
    # The cylinder is held at T = 1 every step (Solver.T_pin), which is what
    # "isothermal cylinder" means. Its k must still be well above the fluid's,
    # because the harmonic-mean face conductivity then tends to 2*k_fluid and
    # the discrete flux becomes k_f*(T_wall - T_first)/(dx/2) -- the correct
    # half-cell gradient. Leaving k_solid = k_fluid would halve every flux.
    #
    # rho*cp is raised by the same factor so the solid's DIFFUSIVITY matches
    # the fluid's and does not drive the explicit diffusive step limit. Since
    # the pin overwrites the solid temperature every step, its capacity has no
    # physical effect here -- it is purely a numerical stabiliser. Raising k
    # alone cut the time step by 50x and the benchmark simply stopped finishing.
    k[solid] = alpha * 50.0
    rho_cp[solid] = 50.0

    mat = Materials(rho_cp=rho_cp, k=k, solid=solid, nu=nu, rho_f=1.0,
                    beta=0.0, T_ref=0.0, gravity=0.0)
    bc = BC(left="inflow", right="outflow", top="slip", bottom="slip",
            U_in=U, T_in=0.0)
    T0 = np.zeros((nx, ny))
    T0[solid] = 1.0
    s = Solver(dom, mat, bc, T0, cfl=0.4, T_pin=(solid, 1.0))
    s.u[:] = U
    s.u[s.u_solid] = 0.0
    s._apply_velocity_bc()

    probe_i = int((xc + 3.0 * D) / dx)
    probe_j = int((yc + 0.5 * D) / dx)
    t_hist, v_hist, nu_hist = [], [], []
    t_end = t_end_D * D / U
    while s.t < t_end:
        s.step()
        if s.step_count % 10 == 0:
            t_hist.append(s.t)
            v_hist.append(float(s.v[probe_i, probe_j]))
            if s.t > 0.4 * t_end:
                Q, per = s.interface_heat_rate()
                # h = Q / (perimeter * dT); Nu = h D / k_fluid
                nu_hist.append(Q / (np.pi * D * 1.0) * D / alpha)

    t_hist = np.array(t_hist)
    v_hist = np.array(v_hist)

    # Strouhal from the wake probe over the last half of the record.
    half = len(t_hist) // 2
    St = _dominant_frequency(t_hist[half:], v_hist[half:]) * D / U

    # Recirculation length: distance from the rear stagnation point to where
    # the centreline u first turns positive again.
    j = ny // 2
    u_c = 0.5 * (s.u[:-1, j] + s.u[1:, j])
    i0 = int((xc + D / 2) / dx)
    Lr = np.nan
    for i in range(i0, nx - 1):
        if u_c[i] > 0:
            Lr = ((i + 0.5) * dx - (xc + D / 2)) / D
            break

    Nu = float(np.mean(nu_hist[-max(1, len(nu_hist) // 3):])) if nu_hist else np.nan

    import ht
    Nu_cb = ht.Nu_cylinder_Churchill_Bernstein(Re, Pr)
    if verbose:
        print(f"    Re={Re:6.0f} Pr={Pr:4.1f}  {nx}x{ny}  D/dx={D_cells}  "
              f"Lr/D={Lr:5.2f}  St={St:6.4f}  Nu={Nu:7.2f}  "
              f"Nu_CB={Nu_cb:7.2f}  ({(Nu/Nu_cb-1)*100:+6.1f} %)")
    return {"Re": Re, "Pr": Pr, "D_cells": int(D_cells), "nx": nx, "ny": ny,
            "dx_over_D": float(dx / D),
            "Lr_over_D": float(Lr), "St": float(St),
            "Nu": Nu, "Nu_ChurchillBernstein": float(Nu_cb),
            "Nu_error_pct": float((Nu / Nu_cb - 1) * 100),
            "probe_t": t_hist.tolist()[-400:], "probe_v": v_hist.tolist()[-400:]}


def _dominant_frequency(t, y):
    if len(t) < 32:
        return np.nan
    y = y - y.mean()
    if np.abs(y).max() < 1e-8:
        return 0.0
    dt = t[1] - t[0]
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y))))
    f = np.fft.rfftfreq(len(y), dt)
    return float(f[np.argmax(Y[1:]) + 1])


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    t0 = time.time()

    if which in ("all", "cavity"):
        print("1. Lid-driven cavity vs Ghia et al. (1982)")
        out["lid_cavity"] = [lid_driven_cavity(Re, n=128) for Re in (100, 400, 1000)]
        print()

    if which in ("all", "heated"):
        print("2. Differentially heated cavity vs de Vahl Davis (1983)")
        out["heated_cavity"] = [heated_cavity(Ra, n=n) for Ra, n in
                                ((1e3, 48), (1e4, 64), (1e5, 96), (1e6, 160))]
        print()

    if which in ("all", "cyl"):
        print("3. Cylinder cross-flow")
        # Two separate questions, deliberately separated.
        #
        # (a) Is the FLOW right? Wake length at Re=40 and Strouhal number at
        #     Re=100 answer that without touching the energy equation.
        # (b) Is the WALL HEAT FLUX right? At Pr=7 the thermal layer is
        #     thinner than the momentum layer by ~Pr^(1/3), so a first attempt
        #     at D/dx=24 came out 40-70% low and Nu even FELL with Reynolds
        #     number, which is impossible -- a clear signature of an
        #     unresolved thermal layer rather than a wrong model. So the heat
        #     transfer is checked at Pr=0.7, where the two layers are
        #     comparable and the mesh can actually resolve it, plus a
        #     refinement sequence at Pr=7 to show the error shrinking.
        out["cylinder"] = [
            cylinder_flow(40, D_cells=24, Pr=0.7, t_end_D=70, Lx_D=18, Ly_D=9),
            cylinder_flow(100, D_cells=24, Pr=0.7, t_end_D=120, Lx_D=20, Ly_D=10),
            cylinder_flow(200, D_cells=28, Pr=0.7, t_end_D=120, Lx_D=20, Ly_D=10),
        ]
        out["cylinder_refine"] = [
            cylinder_flow(40, D_cells=n, Pr=7.0, t_end_D=70, Lx_D=16, Ly_D=8)
            for n in (16, 24, 36, 48)
        ]
        print()

    print(f"total {time.time()-t0:.0f} s")
    suffix = "" if which == "all" else f"_{which}"
    with open(f"results/data/validation{suffix}.json", "w") as f:
        json.dump(out, f, indent=1)
