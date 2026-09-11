"""
The head-to-head CFD experiment: the same bottle, the same mesh, the same
solver, cooled two ways.

The plane of the simulation is VERTICAL and contains the bottle axis, with
gravity pointing down and (for the creek) the current running left to right.
That choice is deliberate: it is the only plane in which both mechanisms are
in-plane at once, so the box case and the creek case differ *only* in their
boundary conditions -- closed walls versus inflow/outflow. Nothing else about
the two runs is different, which is what makes the comparison a controlled one.

The bottle is modelled conjugately: a glass shell at its own conductivity
around water contents carried at an effective conductivity k_eff = Nu_i * k_w,
with the whole bottle treated as a stationary conducting body.

Leaving the contents as a free fluid was tried first and abandoned. A 3 mm
glass wall is only two or three cells thick, and neither hard masking nor
Brinkman penalisation holds a membrane that thin: hard masking destroys the
divergence-free condition (max|div| went from 1e-14 to O(1)), while
penalisation leaves a residual ~2e-3 m/s inside the wall, enough to carry the
bottle's contents bodily into the bath in twenty seconds. Representing the
interior circulation through k_eff instead -- with the multiplier taken from
the FiPy study in axisym_fipy.py rather than guessed -- keeps the external
flow, the wall conduction and the bath energy balance all exact, which is
where the question actually lives. The cost is that the internal circulation
is modelled rather than resolved, and the sensitivity to that choice is
reported explicitly.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

import numpy as np

from geometry import BOTTLE_500, Bottle
from nsolver import BC, Domain, Materials, Solver
from props import GLASS, Solid, water

G = 9.80665


@dataclass
class CFDCase:
    name: str
    mode: str = "box"                 # 'box' | 'creek'
    U: float = 0.0                    # m/s, creek approach velocity
    bottle: Bottle = field(default_factory=lambda: BOTTLE_500)
    T_bottle0: float = 40.0
    T_bath0: float = 10.0
    dx: float = 1.0e-3                # m
    box_W: float = 0.24               # m, box width  (box mode)
    box_H: float = 0.26               # m, box height (box mode)
    up_D: float = 4.0                 # upstream length in bottle diameters
    down_D: float = 12.0              # downstream length
    half_H: float = 0.13              # half channel height above/below bottle
    t_end: float = 600.0
    isothermal_bottle: bool = False   # fix the bottle surface (for h extraction)
    k_eff_mult: float = 1.0           # effective-conductivity multiplier for
                                      # the water inside the bottle. The whole
                                      # bottle is then a conducting solid and
                                      # internal circulation is represented
                                      # rather than resolved -- see the note in
                                      # the module docstring.
    solid_interior: bool = True
    void_fraction: float = 0.0
    shell: Solid | None = None
    n_proj: int = 2
    eta: float = 1.0e-4          # Brinkman penalisation time, s
    record_every: float = 2.0         # s
    snapshots: tuple = ()

    def build(self):
        B = self.bottle
        if self.shell is not None:
            B = type(B)(**{**B.__dict__, "shell": self.shell})
        dx = self.dx

        if self.mode == "box":
            Lx, Ly = self.box_W, self.box_H
            nx, ny = int(round(Lx / dx)), int(round(Ly / dx))
            xc, yc = Lx / 2.0, Ly / 2.0
            bc = BC(left="wall", right="wall", top="wall", bottom="wall")
        else:
            Lx = (self.up_D + self.down_D) * B.D_outer
            Ly = 2.0 * self.half_H
            nx, ny = int(round(Lx / dx)), int(round(Ly / dx))
            xc, yc = self.up_D * B.D_outer, Ly / 2.0
            bc = BC(left="inflow", right="outflow", top="slip", bottom="slip",
                    U_in=self.U, T_in=self.T_bath0)

        dom = Domain(nx, ny, nx * dx, ny * dx)
        X, Y = dom.centers()

        half_w, half_h = B.D_outer / 2.0, B.H_outer / 2.0
        outer = (np.abs(X - xc) <= half_w) & (np.abs(Y - yc) <= half_h)
        inner = ((np.abs(X - xc) <= half_w - B.wall) &
                 (np.abs(Y - yc) <= half_h - B.wall))
        shell_mask = outer & ~inner

        w_bath = water(self.T_bath0)
        w_hot = water(self.T_bottle0)
        w_film = water(0.5 * (self.T_bath0 + self.T_bottle0))

        rho_cp = np.full((nx, ny), w_bath.rho * w_bath.cp)
        k = np.full((nx, ny), w_bath.k)
        rho_cp[inner] = w_hot.rho * w_hot.cp
        k[inner] = w_hot.k
        rho_cp[shell_mask] = B.shell.rho * B.shell.cp
        k[shell_mask] = B.shell.k

        if self.solid_interior:
            # The whole bottle is treated as a conducting solid: glass shell at
            # its own conductivity, contents at k_eff = mult * k_water.
            k[inner] = w_hot.k * self.k_eff_mult
            solid = outer.copy()
        else:
            solid = shell_mask.copy()
        if self.isothermal_bottle:
            # Replace the whole bottle with a high-conductivity, huge-capacity
            # block so the wetted surface sits at a fixed temperature. Used
            # only for extracting a clean film coefficient.
            solid = outer.copy()
            k[outer] = w_bath.k * 1e4
            rho_cp[outer] = w_bath.rho * w_bath.cp * 1e7

        mat = Materials(rho_cp=rho_cp, k=k, solid=solid,
                        nu=w_film.mu / w_film.rho, rho_f=w_film.rho,
                        beta=w_film.beta, T_ref=self.T_bath0, gravity=G)

        T0 = np.full((nx, ny), self.T_bath0)
        T0[outer] = self.T_bottle0

        return dom, mat, bc, T0, {"xc": xc, "yc": yc, "outer": outer,
                                  "inner": inner, "shell": shell_mask,
                                  "bottle": B}


def run(case: CFDCase, verbose: bool = True, progress_every: float = 60.0):
    dom, mat, bc, T0, geo = case.build()
    B = geo["bottle"]
    s = Solver(dom, mat, bc, T0, cfl=0.35, n_proj=case.n_proj, eta=case.eta)
    if case.mode == "creek":
        s.u[:] = case.U
        s.u[s.u_solid] = 0.0
        s._apply_velocity_bc()

    inner, outer, shell = geo["inner"], geo["outer"], geo["shell"]
    cell = dom.dx * dom.dy
    C_inner = float((mat.rho_cp * inner).sum() * cell)   # J/K per m depth
    C_bottle = float((mat.rho_cp * outer).sum() * cell)
    A_wet = B.A_outer
    # 2-D perimeter -> per-metre-depth area
    per_2d = 2.0 * (B.D_outer + B.H_outer)

    hist = {"t": [], "T_core": [], "T_bottle": [], "T_bath": [],
            "Q": [], "h_eff": [], "T_surf": [], "div": [], "umax": [],
            "T_bath_top": [], "T_bath_bot": [], "E_total": []}
    snaps = {}
    next_rec = 0.0
    next_prog = progress_every
    fluid_out = ~outer
    t0 = time.time()

    while s.t < case.t_end:
        s.step()
        if s.t >= next_rec:
            next_rec += case.record_every
            T = s.T
            T_core = float((T * mat.rho_cp * inner).sum() * cell / C_inner)
            T_bot = float((T * mat.rho_cp * outer).sum() * cell / C_bottle)
            T_bath = float(T[fluid_out].mean())
            Q, per = s.interface_heat_rate(outer)
            # Surface temperature: mean over the outermost shell cells.
            surf = T[outer].max() if case.isothermal_bottle else _surface_T(T, outer)
            dT = surf - T_bath
            h = Q / (per * dT) if abs(dT) > 1e-6 and per > 0 else np.nan
            ny = dom.ny
            hist["t"].append(s.t)
            hist["T_core"].append(T_core)
            hist["T_bottle"].append(T_bot)
            hist["T_bath"].append(T_bath)
            hist["Q"].append(Q)
            hist["h_eff"].append(h)
            hist["T_surf"].append(float(surf))
            hist["div"].append(s.divergence_norm())
            hist["umax"].append(float(max(np.abs(s.u).max(), np.abs(s.v).max())))
            hist["T_bath_top"].append(float(T[:, int(0.9 * ny):][~outer[:, int(0.9*ny):]].mean()))
            hist["T_bath_bot"].append(float(T[:, :int(0.1 * ny)][~outer[:, :int(0.1*ny)]].mean()))
            hist["E_total"].append(float((mat.rho_cp * T).sum() * cell))
        for ts in case.snapshots:
            if ts not in snaps and s.t >= ts:
                snaps[ts] = {"T": s.T.copy(), "u": s.u.copy(), "v": s.v.copy()}
        if verbose and s.t >= next_prog:
            next_prog += progress_every
            print(f"      t={s.t:7.1f}s  steps={s.step_count:7d}  "
                  f"T_core={hist['T_core'][-1]:6.2f}  "
                  f"h={hist['h_eff'][-1]:7.1f}  "
                  f"wall={time.time()-t0:6.0f}s", flush=True)

    out = {k: np.array(v) for k, v in hist.items()}
    out["case"] = case
    out["geo"] = geo
    out["dom"] = dom
    out["mat"] = mat
    out["snaps"] = snaps
    out["final"] = {"T": s.T.copy(), "u": s.u.copy(), "v": s.v.copy()}
    out["C_inner"] = C_inner
    out["A_wet"] = A_wet
    out["per_2d"] = per_2d
    out["wall_time"] = time.time() - t0
    return out


def _surface_T(T, outer):
    """Area-weighted mean temperature of the outermost ring of bottle cells."""
    import scipy.ndimage as ndi
    eroded = ndi.binary_erosion(outer)
    ring = outer & ~eroded
    return float(T[ring].mean())


def summarise(res) -> dict:
    c = res["case"]
    tail = slice(max(1, len(res["t"]) // 2), None)
    return {
        "name": c.name,
        "mode": c.mode,
        "U": c.U,
        "dx_mm": c.dx * 1e3,
        "t_end": float(res["t"][-1]),
        "T_core_final": float(res["T_core"][-1]),
        "T_bath_final": float(res["T_bath"][-1]),
        "T_bath_top_final": float(res["T_bath_top"][-1]),
        "T_bath_bot_final": float(res["T_bath_bot"][-1]),
        "h_mean_late": float(np.nanmean(res["h_eff"][tail])),
        "Q0": float(res["Q"][0]),
        "max_div": float(np.max(res["div"])),
        "energy_drift_pct": float((res["E_total"][-1] / res["E_total"][0] - 1) * 100),
        "wall_time_s": float(res["wall_time"]),
    }
