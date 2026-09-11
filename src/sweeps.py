"""
Parameter sweeps: every way the answer can be made to change.

The comparison "box versus creek" hides a dozen choices, and several of them
matter more than the one being asked about. This module varies them one at a
time on the calibrated transient model from lumped.py and records where each
one actually lands.
"""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np

from aeration import h_aerated, induced_liquid_velocity, void_to_superficial
from correlations import _FORCED, h_external, resistance_split
from geometry import BOTTLE_500, Bottle, StillBox
from lumped import Scenario, fmt, half_life, simulate, time_to
from props import ALUMINIUM, GLASS, PET, STAINLESS

T_HOT, T_COLD = 40.0, 10.0
METRICS = (30.0, 25.0, 20.0, 15.0, 12.0)


def metrics(sc: Scenario, t_end: float = 4 * 3600.0) -> dict:
    r = simulate(sc, t_end=t_end)
    out = {"name": sc.name, "T_eq": sc.T_equilibrium()}
    for T in METRICS:
        out[f"t_{T:g}"] = time_to(r, T)
    out["t_half"] = half_life(r)
    out["T_30min"] = float(np.interp(1800, r["t"], r["T_bottle"]))
    out["T_60min"] = float(np.interp(3600, r["t"], r["T_bottle"]))
    out["h_out_0"] = float(r["h_out"][0])
    out["U_0"] = float(r["U"][0])
    out["r"] = r
    return out


def table(rows, cols=("t_25", "t_15", "t_12", "T_30min"), label="scenario",
          width=34):
    hdr = f"{label:<{width}}" + "".join(f"{c:>10}" for c in cols)
    print(hdr)
    print("-" * len(hdr))
    for row in rows:
        line = f"{row['name']:<{width}}"
        for c in cols:
            v = row[c]
            line += f"{fmt(v):>10}" if c.startswith("t_") else f"{v:>9.2f}C"
        print(line)


# ---------------------------------------------------------------------------
# 1. Box volume
# ---------------------------------------------------------------------------


def sweep_box_volume(volumes=(0.002, 0.005, 0.010, 0.020, 0.050, 0.100,
                             0.500, 2.0, 20.0)):
    rows = []
    for V in volumes:
        sc = Scenario(f"box {V*1000:>7.1f} L, still", U=0.0, bath_volume=V)
        rows.append(metrics(sc))
        rows[-1]["V"] = V
    sc = Scenario("box, infinite still reservoir", U=0.0, bath_volume=None)
    rows.append(metrics(sc)); rows[-1]["V"] = np.inf
    return rows


# ---------------------------------------------------------------------------
# 2. Creek velocity
# ---------------------------------------------------------------------------


def sweep_creek_velocity(Us=(0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.6,
                             1.0, 1.5, 2.5)):
    rows = []
    for U in Us:
        sc = Scenario(f"creek {U:4.2f} m/s", U=U, bath_volume=None)
        rows.append(metrics(sc))
        rows[-1]["U_creek"] = U
    return rows


# ---------------------------------------------------------------------------
# 3. Crossover map
# ---------------------------------------------------------------------------


def crossover_map(volumes=(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5, 2.0, np.inf),
                  Us=(0.0, 0.01, 0.03, 0.1, 0.3, 1.0)):
    grid = np.zeros((len(volumes), len(Us)))
    for i, V in enumerate(volumes):
        for j, U in enumerate(Us):
            sc = Scenario("", U=U, bath_volume=None if np.isinf(V) else V)
            grid[i, j] = time_to(simulate(sc, t_end=4 * 3600), 15.0)
    return {"volumes": list(volumes), "Us": list(Us), "t15": grid}


# ---------------------------------------------------------------------------
# 4. Stratification of the box
# ---------------------------------------------------------------------------


def sweep_stratification(V=0.020, fracs=(1.0, 0.6, 0.4, 0.25, 0.15)):
    rows = []
    for f in fracs:
        sc = Scenario(f"box 20 L, {f*100:3.0f}% engaged", U=0.0, bath_volume=V,
                      engaged_fraction=f)
        rows.append(metrics(sc))
        rows[-1]["frac"] = f
    return rows


# ---------------------------------------------------------------------------
# 5. Container material and wall thickness
# ---------------------------------------------------------------------------


def sweep_material():
    rows = []
    for mat in (GLASS, PET, STAINLESS, ALUMINIUM):
        for U, tag in ((0.0, "box 20 L"), (0.35, "creek")):
            b = replace(BOTTLE_500, shell=mat, name=mat.name)
            sc = Scenario(f"{mat.name:<16s} {tag}", bottle=b, U=U,
                          bath_volume=None if U else 0.020)
            rows.append(metrics(sc))
    return rows


def sweep_wall_thickness(walls=(0.001, 0.002, 0.003, 0.005, 0.008)):
    rows = []
    for w in walls:
        for U, tag in ((0.0, "box 20 L"), (0.35, "creek")):
            b = replace(BOTTLE_500, wall=w)
            sc = Scenario(f"wall {w*1000:4.1f} mm, {tag}", bottle=b, U=U,
                          bath_volume=None if U else 0.020)
            rows.append(metrics(sc))
            rows[-1]["wall"] = w
    return rows


# ---------------------------------------------------------------------------
# 6. Bottle size
# ---------------------------------------------------------------------------


def sweep_bottle_size(scales=(0.5, 0.75, 1.0, 1.35, 1.6)):
    rows = []
    for s in scales:
        b = replace(BOTTLE_500, D_outer=BOTTLE_500.D_outer * s,
                    H_outer=BOTTLE_500.H_outer * s)
        for U, tag in ((0.0, "box 20 L"), (0.35, "creek")):
            sc = Scenario(f"{b.V_liquid*1e6:6.0f} mL, {tag}", bottle=b, U=U,
                          bath_volume=None if U else 0.020)
            rows.append(metrics(sc))
            rows[-1]["V_liquid"] = b.V_liquid
    return rows


# ---------------------------------------------------------------------------
# 7. Orientation
# ---------------------------------------------------------------------------


def sweep_orientation():
    rows = []
    for orient in ("vertical", "horizontal"):
        for U, tag in ((0.0, "box 20 L"), (0.35, "creek")):
            sc = Scenario(f"{orient:<10s} {tag}", U=U, orientation=orient,
                          bath_volume=None if U else 0.020)
            rows.append(metrics(sc))
    return rows


# ---------------------------------------------------------------------------
# 8. Aeration
# ---------------------------------------------------------------------------


def sweep_aeration(voids=(0.0, 0.002, 0.01, 0.05, 0.15)):
    """
    Aeration enters the transient model through an effective velocity (the
    agitation term) and through the bath film properties (the penalty term).
    """
    B = BOTTLE_500
    rows = []
    for U, tag, V in ((0.0, "box 20 L", 0.020), (0.35, "creek", None)):
        for a in voids:
            j_g = void_to_superficial(a)
            U_ind = induced_liquid_velocity(j_g, B.H_outer)
            U_eff = float(np.hypot(U, U_ind))
            sc = Scenario(f"{tag}, void {a*100:4.1f}%", U=U_eff,
                          bath_volume=V, void_fraction=a)
            rows.append(metrics(sc))
            rows[-1].update({"void": a, "U_eff": U_eff, "base_U": U})
    return rows


# ---------------------------------------------------------------------------
# 9. Bath temperature difference
# ---------------------------------------------------------------------------


def sweep_bath_temperature():
    rows = []
    for Tb in (4.0, 7.0, 10.0, 13.0, 16.0):
        for U, tag, V in ((0.0, "box", 0.020), (0.35, "creek", None)):
            sc = Scenario(f"{tag} at {Tb:4.1f} C", U=U, bath_volume=V,
                          T_bath_0=Tb)
            rows.append(metrics(sc))
            rows[-1]["T_bath"] = Tb
    return rows


# ---------------------------------------------------------------------------
# 10. Correlation uncertainty
# ---------------------------------------------------------------------------


def correlation_spread(U=0.35):
    rows = []
    for m in _FORCED:
        sc = Scenario(f"creek, {m}", U=U, bath_volume=None, forced_method=m)
        rows.append(metrics(sc))
    return rows


# ---------------------------------------------------------------------------
# 11. A box that is not insulated, and a box with ice
# ---------------------------------------------------------------------------


def sweep_box_realism():
    rows = []
    rows.append(metrics(Scenario("box 20 L, insulated", U=0.0, bath_volume=0.020)))
    rows.append(metrics(Scenario("box 20 L, open to a 20 C room", U=0.0,
                                 bath_volume=0.020, T_room=20.0, UA_room=2.0)))
    rows.append(metrics(Scenario("box 20 L, in a 2 C fridge", U=0.0,
                                 bath_volume=0.020, T_room=2.0, UA_room=2.0)))
    rows.append(metrics(Scenario("box 20 L, gently stirred (0.05 m/s)", U=0.05,
                                 bath_volume=0.020)))
    rows.append(metrics(Scenario("box 20 L, ice-water at 0 C", U=0.0,
                                 bath_volume=0.020, T_bath_0=0.0)))
    # Melting ice pins the bath at 0 C until the ice is gone, which makes a
    # small box behave like an infinite reservoir. The latent heat needed is
    # tiny: cooling a 500 mL bottle from 40 C to 5 C releases ~73 kJ, which
    # melts only 220 g of ice.
    rows.append(metrics(Scenario("box + melting ice, pinned 0 C", U=0.0,
                                 bath_volume=None, T_bath_0=0.0)))
    rows.append(metrics(Scenario("creek 0.35 m/s at 10 C", U=0.35,
                                 bath_volume=None)))
    return rows


if __name__ == "__main__":
    import sys
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    store = {}

    def emit(key, rows, cols=("t_25", "t_15", "t_12", "T_30min"), width=34):
        store[key] = [{k: (None if isinstance(v, float) and not np.isfinite(v)
                           else v)
                       for k, v in r.items() if k != "r"} for r in rows]
        table(rows, cols=cols, width=width)
        print()

    if what in ("all", "1"):
        print("=" * 78); print("1. BOX VOLUME  (still water, no current)"); print("=" * 78)
        emit("box_volume", sweep_box_volume())

    if what in ("all", "2"):
        print("=" * 78); print("2. CREEK VELOCITY  (infinite reservoir at 10 C)"); print("=" * 78)
        emit("creek_velocity", sweep_creek_velocity())

    if what in ("all", "4"):
        print("=" * 78); print("4. STRATIFICATION OF THE BOX"); print("=" * 78)
        emit("stratification", sweep_stratification())

    if what in ("all", "5"):
        print("=" * 78); print("5. CONTAINER MATERIAL"); print("=" * 78)
        emit("material", sweep_material())
        print("=" * 78); print("5b. WALL THICKNESS"); print("=" * 78)
        emit("wall", sweep_wall_thickness())

    if what in ("all", "6"):
        print("=" * 78); print("6. BOTTLE SIZE"); print("=" * 78)
        emit("size", sweep_bottle_size())

    if what in ("all", "7"):
        print("=" * 78); print("7. ORIENTATION"); print("=" * 78)
        emit("orientation", sweep_orientation())

    if what in ("all", "8"):
        print("=" * 78); print("8. AERATION (entrained air)"); print("=" * 78)
        emit("aeration", sweep_aeration())

    if what in ("all", "9"):
        print("=" * 78); print("9. BATH TEMPERATURE"); print("=" * 78)
        emit("bath_T", sweep_bath_temperature())

    if what in ("all", "10"):
        print("=" * 78); print("10. CORRELATION SPREAD (same case, six correlations)"); print("=" * 78)
        emit("corr", correlation_spread())

    if what in ("all", "11"):
        print("=" * 78); print("11. MORE REALISTIC BOXES"); print("=" * 78)
        emit("realism", sweep_box_realism())

    if what in ("all", "3"):
        print("=" * 78); print("3. CROSSOVER MAP: time to 15 C"); print("=" * 78)
        cm = crossover_map()
        print(f"{'box volume':>14}" + "".join(f"{u:>10.2f}" for u in cm["Us"]))
        print(f"{'':>14}" + "".join(f"{'m/s':>10}" for _ in cm["Us"]))
        for i, V in enumerate(cm["volumes"]):
            lbl = "infinite" if np.isinf(V) else f"{V*1000:.0f} L"
            print(f"{lbl:>14}" + "".join(fmt(t).rjust(10) for t in cm["t15"][i]))
        store["crossover"] = {"volumes": [None if np.isinf(v) else v for v in cm["volumes"]],
                              "Us": cm["Us"],
                              "t15": [[None if not np.isfinite(x) else x for x in row]
                                      for row in cm["t15"]]}
        print()

    with open(f"results/data/sweeps_{what}.json", "w") as f:
        json.dump(store, f, indent=1, default=float)
    print(f"written results/data/sweeps_{what}.json")
