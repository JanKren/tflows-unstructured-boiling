#!/usr/bin/env python3
"""
Scriven campaign analysis: current vs previous d_f, and the extended stencil.

Every run directory keeps its previous-formulation trace as
`bench-data-olddf.dat` alongside the current `bench-data.dat`, so each level and
each stencil configuration carries a paired comparison in which the mesh, the
time step and everything else are identical.

Three questions, all answered from the radius traces:

  1. Does the interfacial distance change the structured results at all?
     On a Cartesian mesh the cell-centre connection is axis-aligned, so the
     normal projection coincides with simpler constructions EXCEPT where the
     direction cosine falls below the floor c_min = 0.05 -- a band of roughly
     2.9 degrees about each great circle, so a few per cent of interface faces.
     The floor raises d_f there, which weakens the coupling, so the expected
     sign is a slightly SMALLER radius.  Measured, not assumed.

  2. What does the extended stencil buy, coupled rather than frozen?
     The frozen-state estimate was -24 to -38% on the shape distortion.  A
     coupled run is the honest test.

  3. Do the observed convergence orders change?

Radius convention: bench-data column 4 is the volume-equivalent radius
R_vol = (3V/4pi)^(1/3).  R_front is larger (3.3% structured, 4.9% polyhedral)
and is the physically correct one, but the figures in the manuscript were built
from column 4, so orders quoted here use it for comparability.

TRAP: check for splicing first.  `Scriven-Poly-100/bench-data.dat` holds three
concatenated runs and `np.interp` over it returns nonsense silently.  This
script de-splices by selecting the segment whose step size matches `--dt`,
falling back to the last (most recent) segment, and it always says so.  Taking
the LONGEST segment would be wrong for exactly that file: its longest pass ran
at dt = 1 us while the run its control describes is the shorter 2 us one.

Usage:  python3 analyze_dfix_campaign.py <dir-with-run-subdirs> [--t 1.0e-3]
"""

import argparse
import os
import glob

import numpy as np

BETA = 4.06022
R0 = 50.0e-6
ALPHA_L = 1.6753e-7
T0 = R0 ** 2 / (4.0 * BETA ** 2 * ALPHA_L)     # virtual origin, 226.3 us


def scriven_radius(t):
    return 2.0 * BETA * np.sqrt(ALPHA_L * (t + T0))


def load_clean(path, want_dt=None, label=""):
    """Load a bench-data file, returning one monotone-in-time segment.

    Picking the LONGEST segment is wrong: Scriven-Poly-100 holds three runs and
    its longest segment is a dt = 1 us pass, while the segment its own control
    describes is the shorter dt = 2 us one.  Prefer the segment whose step size
    matches `want_dt`; failing that take the LAST, which is the most recent
    pass.  Splicing is always announced, never silently resolved.
    """
    if not os.path.exists(path):
        return None
    b = np.loadtxt(path)
    if b.ndim == 1 or len(b) < 3:
        return None
    resets = np.where(np.diff(b[:, 0]) <= 0)[0]
    if len(resets) == 0:
        return b
    segs = np.split(b, resets + 1)
    segs = [s for s in segs if len(s) > 2]
    dts = [float(np.median(np.diff(s[:, 0]))) for s in segs]
    pick = len(segs) - 1
    if want_dt is not None:
        for i, dt in enumerate(dts):
            if abs(dt - want_dt) < 1e-12:
                pick = i
    print(f"    ! {label}: SPLICED, {len(segs)} segments "
          f"(dt {', '.join(f'{d:.1e}' for d in dts)}); using segment {pick+1} "
          f"({len(segs[pick])} rows, dt={dts[pick]:.1e})")
    return segs[pick]


def radius_at(b, t):
    if b is None or t > b[-1, 0] or t < b[0, 0]:
        return np.nan
    return float(np.interp(t, b[:, 0], b[:, 3]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="directory containing the run subdirectories")
    ap.add_argument("--t", type=float, default=1.0e-3,
                    help="evaluation time [s] (default 1 ms, the analysis time)")
    ap.add_argument("--dt", type=float, default=2.0e-6,
                    help="expected time step; selects the right segment of a "
                         "spliced bench-data file (default 2 us)")
    a = ap.parse_args()

    runs = sorted(d for d in glob.glob(os.path.join(a.root, "*"))
                  if os.path.isdir(d)
                  and os.path.exists(os.path.join(d, "bench-data.dat")))
    if not runs:
        raise SystemExit(f"no run directories with bench-data.dat under {a.root}")

    Ra = scriven_radius(a.t)
    print(f"\n{'='*94}")
    print(f"  Scriven campaign at t = {a.t*1e3:.3f} ms   "
          f"R_Scriven = {Ra*1e6:.3f} um")
    print(f"{'='*94}")
    print(f"  {'run':<22s} {'R prev':>9s} {'R curr':>9s} {'change':>9s} "
          f"{'err prev':>9s} {'err curr':>9s} {'gain':>8s}")
    print("  " + "-" * 90)

    table = {}
    for d in runs:
        name = os.path.basename(d)
        cur = load_clean(os.path.join(d, "bench-data.dat"), a.dt, name+"/curr")
        old = load_clean(os.path.join(d, "bench-data-olddf.dat"), a.dt, name+"/prev")
        Rc, Ro = radius_at(cur, a.t), radius_at(old, a.t)
        if np.isnan(Rc):
            print(f"  {name:<22s} {'(not yet at this time)':>50s}")
            continue
        ec = Rc / Ra - 1.0
        table[name] = (Rc, ec)
        if np.isnan(Ro):
            print(f"  {name:<22s} {'--':>9s} {Rc*1e6:9.3f} {'--':>9s} "
                  f"{'--':>9s} {ec*100:+8.2f}% {'--':>8s}")
        else:
            eo = Ro / Ra - 1.0
            print(f"  {name:<22s} {Ro*1e6:9.3f} {Rc*1e6:9.3f} "
                  f"{100*(Rc/Ro-1):+8.3f}% {eo*100:+8.2f}% {ec*100:+8.2f}% "
                  f"{100*(abs(eo)-abs(ec)):+7.2f}pp")

    # ---- convergence order per stencil configuration -----------------------
    H = {"075": 4.0e-6, "100": 3.0e-6, "125": 2.4e-6, "150": 2.0e-6}
    cfgs = sorted({n.split("-")[-1] for n in table if n.count("-") >= 2})
    if cfgs:
        print(f"\n  observed order of the radius error (fitted over the levels present)")
        print("  " + "-" * 62)
        for c in cfgs:
            hs, es, lv = [], [], []
            for L in ("075", "100", "125", "150"):
                for n, (R, e) in table.items():
                    if n.endswith("-" + c) and L in n:
                        hs.append(H[L]); es.append(abs(e)); lv.append(L)
            if len(hs) >= 3:
                p = np.polyfit(np.log(hs), np.log(es), 1)[0]
                print(f"    {c:<10s} levels {','.join(lv):<18s} p = {p:+.3f}   "
                      + "  ".join(f"{e*100:.2f}%" for e in es))

    # ---- what the extended stencil buys, at matched level ------------------
    print(f"\n  extended stencil against its own baseline, same level")
    print("  " + "-" * 62)
    for L in ("075", "100", "125", "150"):
        base = next((v for n, v in table.items() if L in n and n.endswith("-off")), None)
        if base is None:
            continue
        for c in ("extonly", "extw1"):
            v = next((v for n, v in table.items() if L in n and n.endswith("-" + c)), None)
            if v is None:
                continue
            print(f"    level {L}  {c:<9s} err {abs(base[1])*100:6.3f}% -> "
                  f"{abs(v[1])*100:6.3f}%   ({(abs(v[1])/abs(base[1])-1)*100:+.1f}% relative)")
    print()


if __name__ == "__main__":
    main()
