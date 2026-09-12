"""Emit a single JSON blob of every series the report page plots."""

from __future__ import annotations

import json
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from aeration import h_aerated
from correlations import (_FORCED, h_external, h_forced_crossflow,
                          h_free_vertical, h_internal, resistance_split)
from geometry import BOTTLE_500 as B
from lumped import Scenario, simulate, time_to
from props import bubbly_mixture, do_saturation, water, water_with_dissolved_o2

out = {}

# --- 1. resistance chain -----------------------------------------------
hi = h_internal(B.D_inner, B.H_inner, 30.0, 22.0)
chain = []
for nm, U in (("Still box", 0.0), ("Creek 0.05 m/s", 0.05),
              ("Creek 0.35 m/s", 0.35), ("Creek 1.5 m/s", 1.5)):
    ho = h_external(B.D_outer, B.H_outer, U, 22.0, 10.0)
    sp = resistance_split(B, ho, hi)
    chain.append({"name": nm, "U": U, "h_out": ho,
                  "inside": sp["inside"] * 100, "wall": sp["wall"] * 100,
                  "outside": sp["outside"] * 100,
                  "U_overall": 1.0 / (sp["R_total"] * B.A_outer)})
out["chain"] = chain
out["wall_equiv_h"] = B.U_wall

# --- 2. cooling curves --------------------------------------------------
curves = []
for nm, kw in (("Box 5 L", dict(U=0.0, bath_volume=0.005)),
               ("Box 20 L", dict(U=0.0, bath_volume=0.020)),
               ("Box 200 L", dict(U=0.0, bath_volume=0.200)),
               ("Creek 0.35 m/s", dict(U=0.35, bath_volume=None))):
    r = simulate(Scenario(nm, **kw), t_end=3600, n_out=400)
    idx = np.linspace(0, len(r["t"]) - 1, 90).astype(int)
    curves.append({"name": nm,
                   "t": [round(float(x), 1) for x in r["t"][idx]],
                   "T": [round(float(x), 3) for x in r["T_bottle"][idx]],
                   "Tbath": [round(float(x), 3) for x in r["T_bath"][idx]]})
out["curves"] = curves

# --- 3. h versus velocity ----------------------------------------------
U = np.concatenate([[0.0], np.logspace(-2.3, 0.48, 46)])
hv = {"U": [float(u) for u in U],
      "mixed": [h_external(B.D_outer, B.H_outer, u, 22.0, 10.0) for u in U],
      "lo": [], "hi": []}
# The band is the spread of the SIX correlations after each is blended with
# buoyancy, i.e. the uncertainty in the quantity actually plotted. Showing the
# raw forced-convection spread instead would send the band diving away from the
# curve at low velocity, where every pure-forced correlation tends to zero but
# the real transfer is floored by natural convection.
from correlations import h_mixed
h_free_ref = h_free_vertical(B.H_outer, 22.0, 10.0, D=B.D_outer)
for u in U:
    vals = [h_mixed(h_forced_crossflow(B.D_outer, u, 22.0, 10.0, m), h_free_ref)
            for m in _FORCED]
    hv["lo"].append(min(vals)); hv["hi"].append(max(vals))
hv["free"] = h_free_ref
out["h_vs_U"] = hv

# --- 4. aeration --------------------------------------------------------
voids = np.linspace(0, 0.30, 31)
out["aeration"] = {
    "void": [float(v) for v in voids],
    "series": [{"name": nm, "net": [h_aerated(B.D_outer, B.H_outer, u, a,
                                              22.0, 10.0)["net_effect_pct"]
                                    for a in voids]}
               for u, nm in ((0.0, "Still box"), (0.05, "Creek 0.05 m/s"),
                             (0.35, "Creek 0.35 m/s"), (1.0, "Creek 1.0 m/s"))],
    "split": {k: [h_aerated(B.D_outer, B.H_outer, 0.35, a, 22.0, 10.0)[k]
                  for a in voids]
              for k in ("property_effect_pct", "agitation_effect_pct",
                        "net_effect_pct")},
}

# --- 5. dissolved oxygen / bubbles -------------------------------------
base = water(10.0)
o2 = water_with_dissolved_o2(10.0, saturation=1.0)
b1 = bubbly_mixture(10.0, 0.01)
out["oxygen"] = {
    "saturation_table": [{"T": T, "mgL": do_saturation(T)}
                         for T in (0, 5, 10, 15, 20, 25)],
    "props": [{"name": n,
               "do": abs(fn(o2) / fn(base) - 1) * 100,
               "bub": abs(fn(b1) / fn(base) - 1) * 100}
              for n, fn in (
                  ("density  \u03c1", lambda f: f.rho),
                  ("volumetric capacity  \u03c1c\u209a",
                   lambda f: f.rho * f.cp),
                  ("conductivity  k", lambda f: f.k),
                  ("viscosity  \u03bc", lambda f: f.mu))],
}

# --- 6. crossover grid --------------------------------------------------
vols = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5, None]
Us = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
grid = []
for V in vols:
    row = []
    for u in Us:
        t = time_to(simulate(Scenario("", U=u, bath_volume=V), t_end=3 * 3600), 15.0)
        row.append(None if not np.isfinite(t) else round(t / 60.0, 1))
    grid.append(row)
out["crossover"] = {"volumes": vols, "Us": Us, "t15_min": grid}

# --- 7. density anomaly -------------------------------------------------
out["anomaly"] = [{"Tb": Tb, "Ts": Ts, "rho_b": water(Tb).rho,
                   "rho_s": water(Ts).rho,
                   "h": h_free_vertical(B.H_outer, Ts, Tb, D=B.D_outer)}
                  for Tb, Ts in ((10, 25), (10, 15), (4, 10), (0, 8), (2, 6))]

# --- 8. headline table --------------------------------------------------
head = []
for nm, kw in (("Still box, 20 L at 10 °C", dict(U=0.0, bath_volume=0.020)),
               ("Creek, 0.35 m/s at 10 °C", dict(U=0.35, bath_volume=None)),
               ("Still box, 2 L at 10 °C", dict(U=0.0, bath_volume=0.002)),
               ("Ice-water box, 0 °C", dict(U=0.0, bath_volume=None, T_bath_0=0.0))):
    sc = Scenario(nm, **kw)
    r = simulate(sc, t_end=4 * 3600)
    head.append({"name": nm,
                 "t25": time_to(r, 25), "t15": time_to(r, 15),
                 "t12": time_to(r, 12),
                 "T30": float(np.interp(1800, r["t"], r["T_bottle"])),
                 "Teq": sc.T_equilibrium()})
out["headline"] = [{k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                    for k, v in h.items()} for h in head]

# --- 9. levers ----------------------------------------------------------
from dataclasses import replace
from props import ALUMINIUM, PET
base_t = time_to(simulate(Scenario("", U=0.0, bath_volume=0.020), t_end=4*3600), 15)
levers = [("Baseline: still 20 L box", dict(U=0.0, bath_volume=0.020))]
levers += [
    ("Swap glass for aluminium", dict(U=0.0, bath_volume=0.020,
                                      bottle=replace(B, shell=ALUMINIUM))),
    ("Use 4 °C water, not 10 °C", dict(U=0.0, bath_volume=0.020, T_bath_0=4.0)),
    ("Move it to a 0.35 m/s creek", dict(U=0.35, bath_volume=None)),
    ("Bubble air through the box", dict(U=0.215, bath_volume=0.020,
                                        void_fraction=0.01)),
    ("Use 250 mL, not 500 mL", dict(U=0.0, bath_volume=0.020,
                                    bottle=replace(B, D_outer=B.D_outer*0.794,
                                                   H_outer=B.H_outer*0.794))),
    ("Stir the box by hand", dict(U=0.1, bath_volume=0.020)),
    ("Halve the wall to 1.5 mm", dict(U=0.0, bath_volume=0.020,
                                      bottle=replace(B, wall=0.0015))),
    ("Swap glass for PET", dict(U=0.0, bath_volume=0.020,
                                bottle=replace(B, shell=PET))),
]
out["levers"] = []
for nm, kw in levers:
    t = time_to(simulate(Scenario("", **kw), t_end=4 * 3600), 15)
    out["levers"].append({"name": nm, "t15_min": round(t / 60, 1),
                          "speedup": round(base_t / t, 2)})

# --- 10. bottle facts ---------------------------------------------------
out["bottle"] = {"D": B.D_outer, "H": B.H_outer, "wall": B.wall,
                 "V": B.V_liquid, "A": B.A_outer,
                 "C_liquid": B.C_liquid(), "C_shell": B.C_shell(),
                 "R_wall": B.R_wall, "U_wall": B.U_wall}

# --- 11. CFD results, if present ---------------------------------------
for tag in ("grid", "boxsize", "creek", "h2h"):
    try:
        path = f"results/data/cfd_{tag}.json"
        if not os.path.exists(path):
            path = f"results/data/cfd_{tag}_partial.json"
        d = json.load(open(path))
        for r in d["records"]:
            if tag != "h2h":
                r.pop("series", None)
        out[f"cfd_{tag}"] = d
    except FileNotFoundError:
        pass
for tag in ("cavity", "heated", "cyl"):
    path = f"results/data/validation_{tag}.json"
    if not os.path.exists(path):
        path = f"results/data/validation_{tag}_partial.json"
    if os.path.exists(path):
        out[f"val_{tag}"] = json.load(open(path))

# --- 12. The inverse study: when does the box win? ---------------------
import box_wins as bw
import other_factors as of
from props import brine_freeze_C

out["inverse"] = {
    "crossover": [
        {k: (None if isinstance(v, float) and not np.isfinite(v) else v)
         for k, v in r.items()} for r in bw.temperature_crossover()],
    "ladder": [
        {k: (None if isinstance(v, float) and not np.isfinite(v) else v)
         for k, v in r.items()} for r in bw.cold_bath_ladder()],
    "shelter": bw.shelter_sweep()[0],
    "shelter_box_t15": bw.shelter_sweep()[1],
    "opposed": bw.opposed_flow(),
    "size": bw.size_scaling(),
    "drive": bw.drive_scaling(),
    "horizon": bw.time_horizon(),
}
out["other"] = {
    "mixing": of.internal_mixing(),
    "fill": of.fill_level(),
    "contents": of.contents_sweep(),
    "evaporative": of.evaporative(),
    "wet_bulb": of.wet_bulb(),
}
for tag in ("placement", "aspect"):
    path = f"results/data/cfd_{tag}.json"
    if os.path.exists(path):
        d = json.load(open(path))
        for r in d["records"]:
            r.pop("series", None)
        out[f"cfd_{tag}"] = d

json.dump(out, open("results/data/web.json", "w"), default=float)
print("written results/data/web.json",
      f"({len(json.dumps(out, default=float))/1024:.0f} kB)")
print("keys:", ", ".join(out.keys()))
