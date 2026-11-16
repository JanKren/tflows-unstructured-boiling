#!/usr/bin/env python3
"""
Screen every run directory for spliced/mixed output before it is used.

Three datasets in this project have turned out to contain output from more than
one run interleaved in the same directory (Scriven-Poly-75, Scriven-Poly-100 and
the Sucking validation case).  Read as time series they give wrong answers, and
in one case a published number moved by two orders of magnitude once the
contamination was removed.  This screen is cheap and should precede any
time-series analysis.

Two independent tests, both necessary:

  physical : the vapour volume of an evaporating case must increase
             monotonically.  A decrease is proof of splicing (or of a diverged
             solver) and does not depend on file metadata at all.

  metadata : snapshots written by one run have monotonically increasing mtimes
             and a single write date.  Several write dates, or an mtime that
             goes backwards with increasing timestep, indicate splicing.
             Weaker than the physical test -- rsync and cp do not preserve
             order -- so it is reported but never used alone.

Usage:  python3 screen_run_integrity.py [--sample 10]
"""

import argparse
import glob
import os
import time
from collections import Counter

import numpy as np
import pyvista as pv

DIRS = [
    ("Stefan", "/home/user/runs/tflows-vof/Validation/Stefan", "*-ts*"),
    ("Stefan_budget", "Data/Stefan_budget", "stefan-ts*"),
    ("Sucking", "/home/user/runs/tflows-vof/Validation/Sucking", "*-ts*"),
    ("Sucking_clean", "Data/Sucking_clean", "*-ts*"),
    ("Sucking_Fine", "/home/user/runs/tflows-vof/Validation/Sucking_Fine", "*-ts*"),
]
for lvl in (75, 100, 125, 150):
    DIRS.append((f"Tests Poly-{lvl}",
                 f"/home/user/runs/tflows-vof/Validation/Scriven/Scriven-Poly-{lvl}",
                 f"Scriven-Poly-{lvl}_dual-ts*"))
    DIRS.append((f"Merlin Poly-{lvl}",
                 f"/home/user/archive/Scriven-Merlin/Scriven-Poly-{lvl}", "*-ts*"))
    DIRS.append((f"aniso-{lvl:03d}-off", f"Data/Aniso/aniso-{lvl:03d}-off", "*-ts*"))
# the 075 run re-recorded cold to 600 steps, so Figure 8(b) reaches 1.2 ms on
# measured data rather than on an interpolation between the 1 and 10 us runs
DIRS.append(("aniso-075-off-600", "Data/Aniso/aniso-075-off-600", "*-ts*"))
DIRS.append(("Merlin Struct-125",
             "/home/user/archive/Scriven-Merlin/Scriven-Struct-125-smallDt", "bubble-ts*"))


def files_of(d, pat):
    fs = [f for f in glob.glob(os.path.join(d, pat + ".pvtu"))
          + glob.glob(os.path.join(d, pat + ".vtu"))
          if "bnd" not in os.path.basename(f)
          and "front" not in os.path.basename(f)
          and "plane" not in os.path.basename(f)]
    out = []
    for f in fs:
        try:
            out.append((int(f.split("-ts")[1].split(".")[0]), f))
        except (IndexError, ValueError):
            pass
    return sorted(out)


def screen_bench_data(d, r0_expect=50.0e-6, rtol=0.05):
    """First-row physical check on bench-data.dat.

    A run built against the wrong case's User_Mod initialises NO bubble and
    then runs to completion: interfacial area is exactly 0.0 at every step and
    the vapour volume sits at the full domain volume, constant.  The
    monotonicity test below does not catch that -- a frozen record never
    decreases -- and it cost 110 core-hours on 30 July 2026 before anyone
    looked at row one (results/dfix-rerun-failure.md).

    Costs nothing and is decisive, so it runs first:
      - interfacial area at t=0 must be non-zero (there must BE an interface)
      - R at t=0 must match the initialisation radius
      - R must actually change over the record (not frozen)

    Columns: time, interfacial area, vapour volume, radius.
    """
    fn = os.path.join(d, "bench-data.dat")
    if not os.path.exists(fn):
        return None
    try:
        b = np.loadtxt(fn)
    except Exception:
        return "unreadable"
    if b.ndim == 1:
        b = b[None, :]
    if b.shape[0] < 2 or b.shape[1] < 4:
        return "too short"
    area0, r0, r_last = b[0, 1], b[0, 3], b[-1, 3]
    if area0 == 0.0:
        return f"DEAD: area=0 at t=0, R frozen at {r_last*1e6:.2f} um"
    if abs(r0 / r0_expect - 1.0) > rtol:
        return f"DEAD: R(t=0) = {r0*1e6:.2f} um, expected {r0_expect*1e6:.0f}"
    if r_last == r0:
        return f"DEAD: R never changed from {r0*1e6:.2f} um"

    # Splicing shows up in bench-data.dat before it shows up anywhere else:
    # concatenated segments leave the time column non-monotone, duplicated, or
    # carrying more than one step size.  Scriven-Poly-100 held three runs
    # (dt = 1 us twice, then 2 us) and passed both the first-row and the
    # vapour-volume tests, because row one is fine and the volume only dips at
    # the two joins.  np.interp on such a file returns nonsense silently.
    tcol = b[:, 0]
    d = np.diff(tcol)
    n_reset = int((d <= 0).sum())
    n_dup = len(tcol) - len(np.unique(np.round(tcol, 12)))
    steps = np.unique(np.round(d[d > 0], 12))
    if n_reset:
        return (f"SPLICED: {n_reset} time reset(s), {n_dup} duplicate times, "
                f"{len(steps)} distinct dt -> de-splice before use")
    if len(steps) > 1:
        return (f"SUSPECT: {len(steps)} distinct time steps "
                f"({', '.join(f'{v:.1e}' for v in steps[:4])})")
    return f"ok: R {r0*1e6:.2f} -> {r_last*1e6:.2f} um, dt={steps[0]:.1e}, {len(b)} rows"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=10)
    ap.add_argument("--bench-only", action="store_true",
                    help="run only the cheap bench-data.dat first-row screen")
    a = ap.parse_args()

    print(f"\n{'='*100}\n  bench-data.dat first-row screen "
          f"(catches a run that initialised no bubble)\n{'='*100}")
    # the 50 um initial radius is Scriven-specific; Stefan and Sucking are
    # planar and write a bench-data.dat with an unrelated meaning
    for label, d, _ in DIRS:
        if not any(k in label for k in ("Poly", "Struct", "Scriven", "aniso")):
            continue
        v = screen_bench_data(d)
        if v is not None:
            print(f"  {label:<20s} {v}")
    if a.bench_only:
        return

    print(f"\n{'='*100}")
    print("  run-integrity screen: physical monotonicity + write-metadata")
    print(f"{'='*100}")
    print("  %-20s %6s %16s %7s %8s  %s" %
          ("dataset", "snaps", "ts range", "dates", "mtime-br", "vapour volume"))
    print("  " + "-" * 98)

    for label, d, pat in DIRS:
        fs = files_of(d, pat)
        if not fs:
            print("  %-20s %6s  (no snapshots)" % (label, "-"))
            continue
        ts = [t for t, _ in fs]
        dates = Counter(time.strftime("%m-%d", time.localtime(os.path.getmtime(f)))
                        for _, f in fs)
        mt = [os.path.getmtime(f) for _, f in fs]
        br = sum(1 for i in range(1, len(mt)) if mt[i] < mt[i - 1])

        sel = np.linspace(0, len(fs) - 1, min(a.sample, len(fs))).astype(int)
        vols, bad = [], 0
        for i in sel:
            try:
                m = pv.read(fs[i][1])
                al = np.asarray(m.cell_data["Vof Sharp [1]"], float)
                vv = np.asarray(m.cell_data["Grid Cell Volume [m^3]"], float)
                vols.append(float(((1.0 - al) * vv).sum()))
            except Exception:
                bad += 1
        drops = sum(1 for i in range(1, len(vols)) if vols[i] < vols[i - 1] * 0.999)
        verdict = ("SPLICED" if drops else
                   "suspect" if (len(dates) > 1 or br) else "clean")
        print("  %-20s %6d %16s %7d %8d  %s" %
              (label, len(fs), f"{ts[0]}..{ts[-1]}", len(dates), br,
               f"{'OK' if not drops else str(drops)+' DECREASES'}  -> {verdict}"
               + (f"  [{bad} unreadable]" if bad else "")))

    print("\n  SPLICED  = vapour volume decreases -> two or more runs interleaved,")
    print("             or a diverged solver. Do not use as a time series.")
    print("  suspect  = metadata inconsistent but physically monotonic; usually")
    print("             an rsync/cp artefact. Safe for single-snapshot use.")


if __name__ == "__main__":
    main()
