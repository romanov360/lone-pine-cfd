"""
The CFD campaign. Each entry point writes a .npz to results/data/.

Split into short, independently runnable studies so they can be launched in
parallel and so a failure in one does not cost the others.
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np

from case_bottle import CFDCase, run, summarise


def _save(tag, records, extra=None):
    path = f"results/data/cfd_{tag}.json"
    with open(path, "w") as f:
        json.dump({"records": records, "extra": extra or {}}, f, indent=1,
                  default=float)
    print(f"  -> {path}")


def _series(res, keys=("t", "T_core", "T_bath", "h_eff", "Q", "T_surf",
                       "T_bath_top", "T_bath_bot")):
    return {k: np.asarray(res[k]).tolist() for k in keys}


# ---------------------------------------------------------------------------
# A. Grid convergence, box, isothermal bottle
# ---------------------------------------------------------------------------


def grid_study(dxs=(2.0e-3, 1.5e-3, 1.0e-3, 0.75e-3), t_end=60.0):
    print("A. Grid convergence: bottle at 25 C in a still 10 C box")
    recs = []
    for dx in dxs:
        c = CFDCase(f"grid dx={dx*1e3:.2f}mm", mode="box", dx=dx, t_end=t_end,
                    isothermal_bottle=True, T_bottle0=25.0, T_bath0=10.0,
                    record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.6 * t_end
        h = float(np.nanmean(np.asarray(r["h_eff"])[tail]))
        s = summarise(r)
        s["h_quasi_steady"] = h
        s["series"] = _series(r)
        recs.append(s)
        print(f"   dx={dx*1e3:5.2f} mm  {s['dx_mm']:.2f}  h={h:7.1f} W/m2K  "
              f"div={s['max_div']:.1e}  dE={s['energy_drift_pct']:+.3f}%  "
              f"({s['wall_time_s']:.0f} s)", flush=True)

    # Richardson extrapolation on the three finest grids.
    if len(recs) >= 3:
        h3, h2, h1 = [recs[i]["h_quasi_steady"] for i in (-3, -2, -1)]
        r_ratio = recs[-3]["dx_mm"] / recs[-2]["dx_mm"]
        try:
            p = np.log(abs((h3 - h2) / (h2 - h1))) / np.log(r_ratio)
            h_exact = h1 + (h1 - h2) / (r_ratio**p - 1)
            gci = 1.25 * abs((h1 - h2) / h1) / (r_ratio**p - 1) * 100
        except Exception:
            p, h_exact, gci = np.nan, np.nan, np.nan
        print(f"   observed order p = {p:.2f}, h(dx->0) = {h_exact:.1f} W/m2K, "
              f"GCI = {gci:.1f} %")
        _save("grid", recs, {"p": p, "h_extrapolated": h_exact, "GCI_pct": gci})
    else:
        _save("grid", recs)
    return recs


# ---------------------------------------------------------------------------
# B. Box size, isothermal bottle
# ---------------------------------------------------------------------------


def box_size_study(sizes=((0.12, 0.20), (0.18, 0.26), (0.26, 0.34),
                          (0.40, 0.46)), dx=1.25e-3, t_end=180.0):
    print("B. Box size: how much does confinement matter?")
    recs = []
    for W, H in sizes:
        c = CFDCase(f"box {W*100:.0f}x{H*100:.0f} cm", mode="box", dx=dx,
                    box_W=W, box_H=H, t_end=t_end, isothermal_bottle=True,
                    T_bottle0=25.0, T_bath0=10.0, record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.5 * t_end
        s = summarise(r)
        h_tail = np.asarray(r["h_eff"])[tail]
        s["h_quasi_steady"] = float(np.nanmean(h_tail))
        s["h_std"] = float(np.nanstd(h_tail))
        s["W"], s["H"] = W, H
        s["gap_D"] = (W - 0.070) / 2 / 0.070
        s["series"] = _series(r)
        recs.append(s)
        print(f"   {W*100:5.1f} x {H*100:5.1f} cm  side gap = {s['gap_D']:4.2f} D  "
              f"h={s['h_quasi_steady']:6.1f}+/-{s['h_std']:4.1f}  bath +{s['T_bath_final']-10:.3f} K  "
              f"top-bottom split {s['T_bath_top_final']-s['T_bath_bot_final']:+.3f} K  "
              f"({s['wall_time_s']:.0f} s)", flush=True)
    _save("boxsize", recs)
    return recs


# ---------------------------------------------------------------------------
# C. Creek velocity, isothermal bottle
# ---------------------------------------------------------------------------


def creek_study(Us=(0.0, 0.005, 0.01, 0.02, 0.05, 0.10), dx=1.2e-3, t_end=45.0):
    print("C. Creek velocity: h versus current")
    import correlations as co
    recs = []
    for U in Us:
        c = CFDCase(f"creek U={U}", mode="creek", U=U, dx=dx, t_end=t_end,
                    isothermal_bottle=True, T_bottle0=25.0, T_bath0=10.0,
                    up_D=2.5, down_D=6.5, half_H=0.13,
                    record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.55 * t_end
        s = summarise(r)
        h_tail = np.asarray(r["h_eff"])[tail]
        s["h_quasi_steady"] = float(np.nanmean(h_tail))
        s["h_std"] = float(np.nanstd(h_tail))
        s["h_correlation"] = co.h_external(0.070, 0.161, U, 25.0, 10.0)
        s["series"] = _series(r)
        recs.append(s)
        print(f"   U={U:6.3f} m/s  h_CFD={s['h_quasi_steady']:7.1f}  "
              f"h_corr={s['h_correlation']:7.1f}  "
              f"ratio={s['h_quasi_steady']/s['h_correlation']:5.2f}  "
              f"({s['wall_time_s']:.0f} s)", flush=True)
    _save("creek", recs)
    return recs


# ---------------------------------------------------------------------------
# D. The head-to-head transient
# ---------------------------------------------------------------------------


def head_to_head(dx=1.5e-3, t_end=360.0, k_eff_mult=5.0):
    print("D. Transient head-to-head, conjugate bottle, identical mesh")
    cases = [
        CFDCase("still box", mode="box", U=0.0, dx=dx, box_W=0.24, box_H=0.30,
                t_end=t_end, k_eff_mult=k_eff_mult, record_every=2.0,
                snapshots=(20.0, 120.0, 300.0)),
        CFDCase("creek 0.02 m/s", mode="creek", U=0.02, dx=dx, t_end=t_end,
                up_D=3.0, down_D=8.0, half_H=0.15, k_eff_mult=k_eff_mult,
                record_every=2.0, snapshots=(20.0, 120.0, 300.0)),
    ]
    recs = []
    for c in cases:
        print(f"   {c.name} ...", flush=True)
        r = run(c, verbose=True, progress_every=120.0)
        s = summarise(r)
        s["series"] = _series(r)
        s["k_eff_mult"] = k_eff_mult
        recs.append(s)
        np.savez_compressed(
            f"results/data/field_{c.name.replace(' ','_').replace('.','p')}.npz",
            T=r["final"]["T"], u=r["final"]["u"], v=r["final"]["v"],
            outer=r["geo"]["outer"], dx=c.dx, nx=r["dom"].nx, ny=r["dom"].ny,
            **{f"snapT_{int(k)}": v["T"] for k, v in r["snaps"].items()},
            **{f"snapu_{int(k)}": v["u"] for k, v in r["snaps"].items()},
            **{f"snapv_{int(k)}": v["v"] for k, v in r["snaps"].items()})
        print(f"   {c.name:16s} T_core {40.0} -> {s['T_core_final']:.2f} C, "
              f"bath {s['T_bath_final']:.2f} C, "
              f"stratification {s['T_bath_top_final']-s['T_bath_bot_final']:+.2f} K, "
              f"dE={s['energy_drift_pct']:+.3f}%  ({s['wall_time_s']/60:.1f} min)",
              flush=True)
    _save("head2head", recs)
    return recs


if __name__ == "__main__":
    what = sys.argv[1]
    t0 = time.time()
    {"grid": grid_study, "boxsize": box_size_study, "creek": creek_study,
     "h2h": head_to_head}[what]()
    print(f"total {time.time()-t0:.0f} s")
