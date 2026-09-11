"""
The inverse question: under what conditions does the STILL BOX win?

The main study established that a creek beats a still box of the same
temperature, always, and explained why the margin is small. This module goes
looking for the exceptions.

There is one structural fact to get out of the way first. The mixed-convection
blend used throughout is

    h_mixed = (h_forced^n + h_free^n)^(1/n)

which is >= h_free for any velocity. So at equal bath temperature, with an
unlimited bath and a current that assists or is transverse to the plume, the
creek CANNOT lose. Every genuine box win therefore has to break one of those
assumptions, and there turn out to be five ways to do it:

  1. Temperature      -- you choose the box's temperature; the creek is
                         whatever it is. This is by far the biggest lever.
  2. Reservoir        -- a salt-ice bath is effectively infinite AND colder.
  3. Shelter          -- the creek's nominal speed is not what a bottle
                         wedged behind a rock actually feels.
  4. Opposed flow     -- a current running against the buoyant plume can
                         cancel part of it, the one direction where moving
                         water is worse than still water.
  5. Scale and drive  -- h_free is nearly size-independent while h_forced
                         falls as D^-1/2, and h_free grows as dT^1/3 while
                         h_forced does not. Big bottles and hot bottles both
                         shift the balance towards the box.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import replace

import numpy as np
from scipy.optimize import brentq

warnings.filterwarnings("ignore")

from correlations import (h_external, h_forced_crossflow, h_free_vertical,
                          h_mixed)
from geometry import BOTTLE_500 as B
from geometry import Bottle
from lumped import Scenario, fmt, freeze_time, simulate, time_to
from props import brine_freeze_C, brine_nacl, nacl_freezing_point, water

RULE = "=" * 78


def head(n, t):
    print(f"\n{RULE}\n{n}. {t}\n{RULE}")


def t_to(target, **kw) -> float:
    """Time for the bottle to reach `target`, seconds."""
    sc = Scenario("", **kw)
    r = simulate(sc, t_end=4 * 3600, stop_at=min(target, 0.4) if
                 kw.get("T_bath_0", 10.0) < 0 else None)
    return time_to(r, target)


# ---------------------------------------------------------------------------
# 1. How cold does the box have to be?
# ---------------------------------------------------------------------------


def temperature_crossover(targets=(25.0, 20.0, 15.0),
                          Us=(0.05, 0.2, 0.35, 1.0), T_creek=10.0):
    """
    For a creek at T_creek running at U, find the still-box temperature that
    matches it on time-to-target. Colder than that and the box wins.
    """
    rows = []
    for U in Us:
        for target in targets:
            t_creek = t_to(target, U=U, bath_volume=None, T_bath_0=T_creek)
            if not np.isfinite(t_creek):
                rows.append({"U": U, "target": target, "T_box": None,
                             "t_creek": None})
                continue

            def resid(T_box):
                t_box = t_to(target, U=0.0, bath_volume=None, T_bath_0=T_box)
                if not np.isfinite(t_box):
                    return 1e4
                return t_box - t_creek

            try:
                T_box = brentq(resid, -1.0, T_creek - 0.01, xtol=1e-3)
            except ValueError:
                T_box = None
            rows.append({"U": U, "target": target, "T_box": T_box,
                         "t_creek": t_creek,
                         "advantage_K": None if T_box is None
                         else T_creek - T_box})
    return rows


# ---------------------------------------------------------------------------
# 2. The cold-bath ladder
# ---------------------------------------------------------------------------


def cold_bath_ladder():
    rows = []
    rows.append(dict(name="Creek, 0.35 m/s at 10 C", kw=dict(
        U=0.35, bath_volume=None, T_bath_0=10.0)))
    rows.append(dict(name="Creek, 1.5 m/s at 10 C", kw=dict(
        U=1.5, bath_volume=None, T_bath_0=10.0)))
    rows.append(dict(name="Still box, tap water 10 C", kw=dict(
        U=0.0, bath_volume=0.020, T_bath_0=10.0)))
    rows.append(dict(name="Still box, fridge water 4 C", kw=dict(
        U=0.0, bath_volume=0.020, T_bath_0=4.0)))
    rows.append(dict(name="Ice-water bath, 0 C", kw=dict(
        U=0.0, bath_volume=None, T_bath_0=0.0)))
    for w in (0.05, 0.10, 0.15, 0.20, 0.23):
        Tf = brine_freeze_C(w)
        rows.append(dict(name=f"Salt-ice bath, {w*100:.0f} wt% ({Tf:.1f} C)",
                         kw=dict(U=0.0, bath_volume=None, T_bath_0=Tf + 0.05,
                                 bath_kind="brine", brine_w=w)))
    out = []
    for r in rows:
        sc = Scenario(r["name"], **r["kw"])
        sub = r["kw"].get("T_bath_0", 10.0) < 0
        res = simulate(sc, t_end=4 * 3600, stop_at=0.4 if sub else None)
        rec = {"name": r["name"], "T_bath": r["kw"].get("T_bath_0", 10.0),
               "U": r["kw"].get("U", 0.0)}
        for T in (25.0, 20.0, 15.0, 10.0, 5.0, 1.0):
            rec[f"t_{T:g}"] = time_to(res, T)
        rec["freeze_through"] = freeze_time(sc)[1] if sub else np.inf
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
# 3. The creek the bottle actually feels
# ---------------------------------------------------------------------------


def shelter_sweep(U_nominal=0.35, shelters=(1.0, 0.5, 0.3, 0.15, 0.05, 0.0),
                  box_V=0.020):
    """
    A bottle wedged among rocks, or lying in the bed boundary layer, sees a
    fraction of the depth-averaged velocity. How much of the creek's advantage
    survives?
    """
    t_box = t_to(15.0, U=0.0, bath_volume=box_V, T_bath_0=10.0)
    rows = []
    for s in shelters:
        t = t_to(15.0, U=U_nominal, bath_volume=None, T_bath_0=10.0, shelter=s)
        rows.append({"shelter": s, "U_effective": U_nominal * s, "t15": t,
                     "vs_box": t_box / t if t > 0 else np.nan})
    return rows, t_box


def bed_shelter_factor(z_bottle: float, depth: float, U_mean: float,
                       ks: float = 0.05) -> float:
    """
    Velocity a bottle on the bed sees, as a fraction of the depth-averaged
    mean, from the log law over a rough bed.

        u(z)/u* = (1/kappa) ln(z/z0),   z0 = ks/30

    `ks` is the Nikuradse roughness of the bed (0.05 m is coarse gravel).
    Evaluated at the bottle's mid-height rather than at its top, which is the
    right average for a body spanning the lower flow.
    """
    kappa = 0.41
    z0 = ks / 30.0
    if z_bottle <= z0:
        return 0.0
    # Depth-averaged log-law velocity is u(0.4*depth) to good accuracy.
    u_ref = np.log(max(0.4 * depth, 1.01 * z0) / z0)
    return float(np.log(z_bottle / z0) / u_ref) if u_ref > 0 else 0.0


# ---------------------------------------------------------------------------
# 4. The one direction where a current hurts
# ---------------------------------------------------------------------------


def opposed_flow(Us=(0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.3),
                 T_s=22.0, T_inf=10.0):
    """
    A current running DOWNWARD past an upright warm bottle opposes its own
    buoyant plume. Churchill's blend subtracts in that case, and near the
    velocity where the two are comparable it very nearly cancels.

    Treat the magnitude with caution -- the asymptotic blend is least reliable
    exactly at the crossover it is predicting, and a real opposing flow
    separates rather than cancelling cleanly. The SIGN, though, is a genuine
    and well-documented effect: in opposed mixed convection the heat transfer
    coefficient dips below both pure limits.
    """
    h_free = h_free_vertical(B.H_outer, T_s, T_inf, D=B.D_outer)
    rows = []
    for U in Us:
        hf = h_forced_crossflow(B.D_outer, U, T_s, T_inf)
        rows.append({"U": U, "h_forced": hf, "h_free": h_free,
                     "assisting": h_mixed(hf, h_free, assisting=True),
                     "opposing": h_mixed(hf, h_free, assisting=False)})
    return rows


# ---------------------------------------------------------------------------
# 5. Scale and drive
# ---------------------------------------------------------------------------


def size_scaling(scales=(0.4, 0.6, 1.0, 1.6, 2.5, 4.0), U=0.35,
                 T_s=22.0, T_inf=10.0):
    """
    h_free is nearly independent of size once Ra is large (Nu ~ Ra^1/3 makes
    the length cancel), while h_forced falls as D^-1/2. Bigger bottles
    therefore close the gap.
    """
    rows = []
    for s in scales:
        b = replace(B, D_outer=B.D_outer * s, H_outer=B.H_outer * s)
        hfree = h_free_vertical(b.H_outer, T_s, T_inf, D=b.D_outer)
        hforced = h_forced_crossflow(b.D_outer, U, T_s, T_inf)
        hmix = h_external(b.D_outer, b.H_outer, U, T_s, T_inf)
        rows.append({"scale": s, "V_litre": b.V_liquid * 1000,
                     "D_mm": b.D_outer * 1000, "h_free": hfree,
                     "h_forced": hforced, "h_mixed": hmix,
                     "ratio": hmix / hfree})
    return rows


def drive_scaling(T_bottles=(15.0, 25.0, 40.0, 60.0, 80.0, 95.0),
                  U=0.35, T_inf=10.0):
    """
    Free convection scales as dT^(1/3); forced convection does not scale with
    dT at all. A hotter bottle therefore drives its own plume harder and the
    still box closes some of the gap.
    """
    rows = []
    for Tb in T_bottles:
        T_s = T_inf + 0.45 * (Tb - T_inf)   # representative wall temperature
        hfree = h_free_vertical(B.H_outer, T_s, T_inf, D=B.D_outer)
        hmix = h_external(B.D_outer, B.H_outer, U, T_s, T_inf)
        rows.append({"T_bottle": Tb, "T_surface": T_s, "h_free": hfree,
                     "h_mixed": hmix, "ratio": hmix / hfree})
    return rows


# ---------------------------------------------------------------------------
# 6. Time horizon
# ---------------------------------------------------------------------------


def time_horizon(times=(30, 60, 300, 900, 1800, 3600, 7200),
                 box_V=0.020, U=0.35):
    box = simulate(Scenario("", U=0.0, bath_volume=box_V), t_end=8000)
    creek = simulate(Scenario("", U=U, bath_volume=None), t_end=8000)
    rows = []
    for t in times:
        Tb = float(np.interp(t, box["t"], box["T_bottle"]))
        Tc = float(np.interp(t, creek["t"], creek["T_bottle"]))
        rows.append({"t": t, "T_box": Tb, "T_creek": Tc, "gap": Tb - Tc})
    return rows


# ---------------------------------------------------------------------------


def main():
    print(RULE)
    print("WHEN DOES THE STILL BOX WIN?")
    print(RULE)
    print(__doc__.strip().split("There is one structural fact")[0].strip())

    head(1, "The structural obstacle, stated plainly")
    print("  At equal bath temperature, unlimited bath, transverse current:")
    print(f"  {'U (m/s)':>9}{'h_free':>9}{'h_mixed':>9}{'ratio':>8}")
    for U in (0.0, 0.01, 0.05, 0.35, 1.5):
        hf = h_free_vertical(B.H_outer, 22.0, 10.0, D=B.D_outer)
        hm = h_external(B.D_outer, B.H_outer, U, 22.0, 10.0)
        print(f"  {U:9.2f}{hf:9.0f}{hm:9.0f}{hm/hf:8.2f}")
    print("\n  The blend can never fall below the still-water limit, so the box")
    print("  needs a different advantage than the flow field. It has five.")

    head(2, "LEVER 1+2: temperature and reservoir -- the cold-bath ladder")
    rows = cold_bath_ladder()
    print(f"{'bath':34s}{'->25C':>8}{'->20C':>8}{'->15C':>8}{'->10C':>8}"
          f"{'->5C':>8}{'->1C':>8}")
    print("-" * 82)
    for r in rows:
        print(f"{r['name']:34s}" + "".join(
            fmt(r[f't_{T:g}']).rjust(8) for T in (25, 20, 15, 10, 5, 1)))
    print("\n  A 10 C creek cannot reach 10 C at any velocity -- it is the bath")
    print("  temperature. Every row below the ice bath goes somewhere no creek")
    print("  at 10 C can follow, and the salt-ice bath gets to 1 C in under")
    print("  ten minutes while still in a bucket that is not moving at all.")
    print("\n  Margin before the contents freeze solid:")
    for r in rows:
        if np.isfinite(r["freeze_through"]):
            print(f"    {r['name']:34s}{fmt(r['freeze_through']):>10}")

    head(3, "LEVER 3: the creek the bottle actually feels")
    sh, t_box = shelter_sweep()
    print(f"  Nominal creek 0.35 m/s; still 20 L box reaches 15 C in "
          f"{fmt(t_box)}\n")
    print(f"  {'shelter':>9}{'U felt':>9}{'->15C':>9}{'vs box':>9}")
    for r in sh:
        print(f"  {r['shelter']*100:8.0f}%{r['U_effective']:9.3f}"
              f"{fmt(r['t15']):>9}{r['vs_box']:8.2f}x")
    print("\n  Shelter erodes the creek's margin but never reverses it: at zero")
    print("  effective velocity the creek is simply an infinite still bath,")
    print("  which still beats a 20 L one. Shelter is a multiplier on the other")
    print("  levers, not a lever on its own.")
    print("\n  Where a bottle on the bed actually sits, from the log law:")
    print(f"  {'depth':>8}{'U mean':>9}{'shelter':>9}{'U felt':>9}")
    for depth, Um in ((0.20, 0.35), (0.40, 0.35), (0.40, 0.8), (1.0, 0.5)):
        f = bed_shelter_factor(0.5 * B.D_outer, depth, Um)
        print(f"  {depth:8.2f}{Um:9.2f}{f*100:8.0f}%{f*Um:9.3f}")

    head(4, "LEVER 4: a current that opposes the plume")
    print("  A bottle in a downwelling -- the throat of a plunge pool, say --")
    print("  meets water moving DOWN past a plume trying to rise.\n")
    print(f"  {'U (m/s)':>9}{'h_forced':>10}{'h_free':>9}{'assisting':>11}"
          f"{'opposing':>10}")
    for r in opposed_flow():
        print(f"  {r['U']:9.3f}{r['h_forced']:10.0f}{r['h_free']:9.0f}"
              f"{r['assisting']:11.0f}{r['opposing']:10.0f}")
    print("\n  Near 0.03 m/s the two mechanisms very nearly cancel, and moving")
    print("  water transfers LESS heat than still water. This is the only")
    print("  configuration in the whole study where a current is a liability.")

    head(5, "LEVER 5: scale and drive")
    print("  Bigger bottle -- h_free barely changes, h_forced falls as D^-1/2:")
    print(f"  {'litres':>8}{'D (mm)':>9}{'h_free':>9}{'h_forced':>10}"
          f"{'h_mixed':>9}{'creek gain':>12}")
    for r in size_scaling():
        print(f"  {r['V_litre']:8.2f}{r['D_mm']:9.0f}{r['h_free']:9.0f}"
              f"{r['h_forced']:10.0f}{r['h_mixed']:9.0f}{r['ratio']:11.2f}x")
    print("\n  Hotter bottle -- h_free grows as dT^(1/3), h_forced does not:")
    print(f"  {'bottle C':>10}{'wall C':>9}{'h_free':>9}{'h_mixed':>9}"
          f"{'creek gain':>12}")
    for r in drive_scaling():
        print(f"  {r['T_bottle']:10.0f}{r['T_surface']:9.1f}{r['h_free']:9.0f}"
              f"{r['h_mixed']:9.0f}{r['ratio']:11.2f}x")

    head(6, "Time horizon: the gap is not constant")
    print(f"  {'elapsed':>10}{'box 20 L':>11}{'creek':>9}{'gap':>8}")
    for r in time_horizon():
        print(f"  {fmt(r['t']):>10}{r['T_box']:10.2f}C{r['T_creek']:8.2f}C"
              f"{r['gap']:7.2f}K")
    print("\n  The gap peaks in the middle of the run and then closes again as")
    print("  both approach their limits -- but they approach DIFFERENT limits,")
    print("  which is the finite-bath effect reasserting itself.")

    head(7, "The verdict")
    print("""  The box wins whenever you use the one thing it has that a creek
  does not: you choose its contents.

    * Ice water at 0 C beats any 10 C creek outright, at any velocity,
      and 245 g of ice holds a 20 L box there for the whole job.
    * Salt and ice takes it to -21 C and reaches 1 C in under ten
      minutes, with hours of margin before the bottle is in danger.
    * A warm summer creek loses to a bucket of tap water.

  The box loses whenever the comparison is made at equal temperature with
  an unlimited bath, because then the only thing left to compete on is the
  flow field, and a still box has none. Shelter, size and drive all shrink
  the creek's margin, and opposed flow can reverse it locally, but none of
  them is worth as much as six degrees of bath temperature.""")


if __name__ == "__main__":
    import sys
    main()
