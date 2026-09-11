"""
Axisymmetric conjugate model of the bottle, in FiPy.

The Navier-Stokes solver in nsolver.py runs in a Cartesian plane, so it gets
the convection right but treats the bottle as an infinite slab in the third
direction. This module fills that gap: a true (r, z) axisymmetric finite-volume
model of water core + glass shell + convective film, on FiPy's CylindricalGrid2D,
where the cell volumes carry the 2*pi*r weighting and the geometry is the real
one.

Internal convection is represented by an effective conductivity in the core,
k_eff = Nu_i * k_water. That is the classical enclosure treatment, and the
multiplier Nu_i is taken from the CFD rather than guessed -- see calibrate().

The external film is represented as a thin extra annulus of conductivity
h * t_film wrapped around the bottle, with its outer face held at the bath
temperature. That is an exact representation of a series film resistance and
avoids hand-rolling a Robin condition, at the cost of two extra cell layers.
"""

from __future__ import annotations

import numpy as np
from fipy import CellVariable, CylindricalGrid2D, DiffusionTerm, TransientTerm
from fipy.tools import numerix

from geometry import BOTTLE_500, Bottle
from props import water


def build(bottle: Bottle = BOTTLE_500, nr: int = 44, nz: int = 100,
          n_film: int = 3, film_thickness: float = 2.0e-3,
          T_bottle0: float = 40.0, T_bath0: float = 10.0,
          h_out: float = 520.0, k_eff_mult: float = 1.0):
    """Assemble the mesh, the property fields and the initial condition."""
    R_out = bottle.D_outer / 2.0
    H_out = bottle.H_outer

    R_dom = R_out + film_thickness
    Z_dom = H_out + 2.0 * film_thickness

    nr_tot = nr + n_film
    nz_tot = nz + 2 * n_film
    dr = R_dom / nr_tot
    dz = Z_dom / nz_tot

    mesh = CylindricalGrid2D(dr=dr, dz=dz, nr=nr_tot, nz=nz_tot)
    r, z = mesh.cellCenters
    z_axis = z - Z_dom / 2.0                      # 0 at mid-height

    in_bottle = (r <= R_out) & (numerix.absolute(z_axis) <= H_out / 2.0)
    in_core = ((r <= R_out - bottle.wall) &
               (numerix.absolute(z_axis) <= H_out / 2.0 - bottle.wall))
    in_shell = in_bottle & ~in_core
    in_film = ~in_bottle

    w_hot = water(T_bottle0)
    w_cold = water(T_bath0)

    # Film layer: conductivity chosen so that t_film / k_film = 1 / h.
    k_film = h_out * film_thickness

    k = CellVariable(mesh=mesh, name="k", value=k_film)
    rc = CellVariable(mesh=mesh, name="rho_cp",
                      value=w_cold.rho * w_cold.cp)
    k.setValue(w_hot.k * k_eff_mult, where=in_core)
    rc.setValue(w_hot.rho * w_hot.cp, where=in_core)
    k.setValue(bottle.shell.k, where=in_shell)
    rc.setValue(bottle.shell.rho * bottle.shell.cp, where=in_shell)

    T = CellVariable(mesh=mesh, name="T", value=T_bath0, hasOld=True)
    T.setValue(T_bottle0, where=in_bottle)

    # The outer face of the film sits at the bath temperature; the axis and
    # nothing else are symmetry planes (FiPy's default zero flux).
    for faces in (mesh.facesRight, mesh.facesTop, mesh.facesBottom):
        T.constrain(T_bath0, faces)

    masks = {"bottle": in_bottle, "core": in_core, "shell": in_shell,
             "film": in_film}
    return mesh, T, k, rc, masks


def simulate(bottle: Bottle = BOTTLE_500, t_end: float = 3600.0,
             dt0: float = 0.5, h_out: float = 520.0, k_eff_mult: float = 1.0,
             T_bottle0: float = 40.0, T_bath0: float = 10.0,
             nr: int = 44, nz: int = 100, verbose: bool = False):
    mesh, T, k, rc, masks = build(bottle, nr=nr, nz=nz, h_out=h_out,
                                  T_bottle0=T_bottle0, T_bath0=T_bath0,
                                  k_eff_mult=k_eff_mult)
    eq = TransientTerm(coeff=rc) == DiffusionTerm(coeff=k.harmonicFaceValue)

    vol = np.asarray(mesh.cellVolumes)
    core = masks["core"]
    bot = masks["bottle"]
    C_core = float((np.asarray(rc.value) * vol * np.asarray(core)).sum())
    C_bot = float((np.asarray(rc.value) * vol * np.asarray(bot)).sum())

    ts, T_core, T_bot, T_centre, T_edge = [], [], [], [], []
    # Centre-of-core and near-wall probes, to expose internal gradients.
    rc_, zc_ = np.asarray(mesh.cellCenters[0]), np.asarray(mesh.cellCenters[1])
    Z = 0.5 * (zc_.max() + zc_.min())
    i_centre = int(np.argmin(rc_**2 + (zc_ - Z) ** 2))
    R_in = bottle.D_outer / 2.0 - bottle.wall
    i_edge = int(np.argmin((rc_ - 0.9 * R_in) ** 2 + (zc_ - Z) ** 2))

    t, dt = 0.0, dt0
    while t < t_end:
        dt = min(dt * 1.03, 20.0, t_end - t)
        T.updateOld()
        for _ in range(2):
            eq.solve(var=T, dt=dt)
        t += dt
        ts.append(t)
        T_core.append(float((np.asarray(T.value) * np.asarray(rc.value) * vol * np.asarray(core)).sum() / C_core))
        T_bot.append(float((np.asarray(T.value) * np.asarray(rc.value) * vol * np.asarray(bot)).sum() / C_bot))
        T_centre.append(float(T.value[i_centre]))
        T_edge.append(float(T.value[i_edge]))
        if verbose and len(ts) % 100 == 0:
            print(f"      t={t:7.1f}  T_core={T_core[-1]:6.2f}")

    return {"t": np.array(ts), "T_core": np.array(T_core),
            "T_bottle": np.array(T_bot), "T_centre": np.array(T_centre),
            "T_edge": np.array(T_edge), "mesh": mesh, "T_final": T.value.copy(),
            "masks": masks, "C_core": C_core, "C_bottle": C_bot}


if __name__ == "__main__":
    import time
    from lumped import Scenario, fmt, simulate as lump_sim, time_to

    print("Axisymmetric FiPy conjugate model vs the lumped model\n")
    print(f"{'k_eff multiplier':>18} {'t->25C':>9} {'t->15C':>9} "
          f"{'core-edge dT @5min':>20}")
    t0 = time.time()
    for mult in (1.0, 5.0, 20.0, 100.0):
        r = simulate(h_out=522.0, k_eff_mult=mult, t_end=3600.0, dt0=0.25)
        i5 = int(np.argmin(np.abs(r["t"] - 300.0)))
        grad = r["T_centre"][i5] - r["T_edge"][i5]

        def tt(target):
            T = r["T_core"]
            if T[-1] > target:
                return np.inf
            i = int(np.argmax(T <= target))
            return float(np.interp(target, [T[i], T[i - 1]], [r["t"][i], r["t"][i - 1]]))

        print(f"{mult:18.0f} {fmt(tt(25)):>9} {fmt(tt(15)):>9} "
              f"{grad:17.2f} K")
    print(f"\n({time.time()-t0:.0f} s)")

    sc = Scenario("lumped, infinite bath at 10 C", U=0.0, bath_volume=None)
    rl = lump_sim(sc, t_end=3600)
    print(f"{'lumped 0-D':>18} {fmt(time_to(rl,25)):>9} {fmt(time_to(rl,15)):>9}")
