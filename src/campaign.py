"""
The CFD campaign. Each entry point writes a .npz to results/data/.

Split into short, independently runnable studies so they can be launched in
parallel and so a failure in one does not cost the others.
"""

from __future__ import annotations

import json
import os
import re
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


def creek_study(Us=(0.0, 0.01, 0.02, 0.05, 0.10, 0.20), dx=1.2e-3,
                up_D=2.0, down_D=4.0):
    """
    Film coefficient against current, bottle held isothermal.

    Run length is set from the FLOW-THROUGH time, not fixed. A first attempt
    used a constant 45 s for every velocity and produced nonsense at the slow
    end -- at U = 0.005 m/s the domain needs 126 s just to exchange its
    contents once, so the reported h was measured on a flow that had not yet
    arrived, and came out at half the still-water value. A creek that slow is
    below what this domain can resolve in reasonable time and is dropped; the
    rest are given at least three and a half flow-throughs before averaging.
    """
    import correlations as co
    D = 0.070
    Lx = (up_D + down_D) * D
    recs = []
    for U in Us:
        t_end = 60.0 if U == 0 else max(60.0, 3.5 * Lx / U)
        c = CFDCase(f"creek U={U}", mode="creek", U=U, dx=dx, t_end=t_end,
                    isothermal_bottle=True, T_bottle0=25.0, T_bath0=10.0,
                    up_D=up_D, down_D=down_D, half_H=0.13,
                    record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.65 * t_end
        s = summarise(r)
        h_tail = np.asarray(r["h_eff"])[tail]
        s["h_quasi_steady"] = float(np.nanmean(h_tail))
        s["h_std"] = float(np.nanstd(h_tail))
        s["flow_throughs"] = t_end * U / Lx if U else 0.0
        s["h_correlation"] = co.h_external(0.070, 0.161, U, 25.0, 10.0)
        s["series"] = _series(r)
        recs.append(s)
        print(f"   U={U:6.3f} m/s  t_end={t_end:6.0f}s "
              f"({s['flow_throughs']:4.1f} flow-throughs)  "
              f"h_CFD={s['h_quasi_steady']:6.1f}+/-{s['h_std']:4.1f}  "
              f"h_corr={s['h_correlation']:7.1f}  "
              f"ratio={s['h_quasi_steady']/s['h_correlation']:5.2f}  "
              f"({s['wall_time_s']:.0f} s)", flush=True)
    _save("creek", recs)
    return recs


# ---------------------------------------------------------------------------
# B2. Where to put the bottle, and what shape the box should be
# ---------------------------------------------------------------------------


def placement_study(offsets=(-0.30, -0.15, 0.0, 0.15, 0.30), dx=1.25e-3,
                    t_end=180.0, W=0.24, H=0.34):
    """
    A still box stratifies: the bottle's own plume pools warm water at the top.
    So height within the box should matter, and this measures how much.
    """
    print("B2. Bottle height within the box")
    recs = []
    for off in offsets:
        c = CFDCase(f"offset {off:+.2f}", mode="box", dx=dx, box_W=W, box_H=H,
                    y_offset=off, t_end=t_end, isothermal_bottle=True,
                    T_bottle0=25.0, T_bath0=10.0, record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.5 * t_end
        s = summarise(r)
        h_tail = np.asarray(r["h_eff"])[tail]
        s["h_quasi_steady"] = float(np.nanmean(h_tail))
        s["h_std"] = float(np.nanstd(h_tail))
        s["y_offset"] = off
        s["series"] = _series(r)
        recs.append(s)
        print(f"   offset {off:+5.2f} ({'low' if off < 0 else 'high' if off > 0 else 'centre':>6})  "
              f"h={s['h_quasi_steady']:6.1f}+/-{s['h_std']:4.1f}  "
              f"bath +{s['T_bath_final']-10:.2f} K  "
              f"top-bottom {s['T_bath_top_final']-s['T_bath_bot_final']:+.2f} K  "
              f"({s['wall_time_s']:.0f} s)", flush=True)
    _save("placement", recs)
    return recs


def aspect_study(shapes=((0.14, 0.514), (0.19, 0.379), (0.24, 0.300),
                         (0.31, 0.232), (0.40, 0.180)),
                 dx=1.25e-3, t_end=180.0):
    """
    Box shape at (very nearly) constant volume. Tall and narrow gives the plume
    a long chimney but squeezes the sides; short and wide does the reverse.
    """
    print("B3. Box shape at constant volume")
    recs = []
    for W, H in shapes:
        c = CFDCase(f"box {W*100:.0f}x{H*100:.0f}", mode="box", dx=dx,
                    box_W=W, box_H=H, t_end=t_end, isothermal_bottle=True,
                    T_bottle0=25.0, T_bath0=10.0, record_every=1.0, n_proj=1)
        r = run(c, verbose=False)
        tail = np.asarray(r["t"]) > 0.5 * t_end
        s = summarise(r)
        h_tail = np.asarray(r["h_eff"])[tail]
        s["h_quasi_steady"] = float(np.nanmean(h_tail))
        s["h_std"] = float(np.nanstd(h_tail))
        s["W"], s["H"] = W, H
        s["area_m2"] = W * H
        s["side_gap_D"] = (W - 0.070) / 2 / 0.070
        s["head_room_H"] = (H - 0.16142) / 2 / 0.16142
        s["series"] = _series(r)
        recs.append(s)
        print(f"   {W*100:4.0f} x {H*100:5.1f} cm (area {W*H*1e4:5.0f} cm2)  "
              f"side gap {s['side_gap_D']:4.2f}D  head room {s['head_room_H']:4.2f}H  "
              f"h={s['h_quasi_steady']:6.1f}+/-{s['h_std']:4.1f}  "
              f"({s['wall_time_s']:.0f} s)", flush=True)
    _save("aspect", recs)
    return recs


# ---------------------------------------------------------------------------
# D. The head-to-head transient
# ---------------------------------------------------------------------------


def head_to_head(dx=1.5e-3, t_end=360.0, k_eff_mult=5.0, only=None):
    """
    Same bottle, same mesh, same solver; only the boundary conditions differ.

    The creek domain is shorter than it looks like it should be, on purpose.
    Its top is a free surface (free-slip), so the bottle's warm plume spreads
    along it as a gravity current with no wall friction to slow it down --
    reaching ~0.19 m/s against the 0.03 m/s peak in the closed box, and
    dragging the CFL-limited time step down with it. That is real physics, not
    a numerical artefact, but it makes a long channel very expensive for what
    it adds.
    """
    print("D. Transient head-to-head, conjugate bottle, identical mesh")
    cases = [
        CFDCase("still box", mode="box", U=0.0, dx=dx, box_W=0.24, box_H=0.30,
                t_end=t_end, k_eff_mult=k_eff_mult, record_every=2.0,
                snapshots=(20.0, 120.0, 300.0)),
        CFDCase("creek 0.02 m/s", mode="creek", U=0.02, dx=dx, t_end=t_end,
                up_D=2.0, down_D=4.5, half_H=0.13, k_eff_mult=k_eff_mult,
                record_every=2.0, snapshots=(20.0, 120.0, 300.0)),
    ]
    if only:
        cases = [c for c in cases if only in c.name]
        if not cases:
            raise SystemExit(f"no head-to-head case matching {only!r}")
    recs = []
    for c in cases:
        print(f"   {c.name} ...", flush=True)
        r = run(c, verbose=True, progress_every=120.0)
        s = summarise(r)
        s["series"] = _series(r)
        s["k_eff_mult"] = k_eff_mult
        recs.append(s)
        tag = re.sub(r"[^A-Za-z0-9]+", "_", c.name).strip("_")
        np.savez_compressed(
            f"results/data/field_{tag}.npz",
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
    # Merge rather than overwrite: the box and creek cases are expensive
    # enough that they are often run separately, and the second run should not
    # discard the first.
    path = "results/data/cfd_h2h.json"
    existing = []
    if os.path.exists(path):
        existing = json.load(open(path)).get("records", [])
    names = {r["name"] for r in recs}
    merged = [r for r in existing if r["name"] not in names] + recs
    _save("h2h", merged)
    return recs


if __name__ == "__main__":
    what = sys.argv[1]
    t0 = time.time()
    fns = {"grid": grid_study, "boxsize": box_size_study, "creek": creek_study,
           "h2h": head_to_head, "placement": placement_study,
           "aspect": aspect_study}
    if what == "h2h" and len(sys.argv) > 2:
        # e.g.  python3 src/campaign.py h2h creek
        head_to_head(only=sys.argv[2])
    else:
        fns[what]()
    print(f"total {time.time()-t0:.0f} s")
