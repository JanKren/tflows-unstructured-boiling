#!/usr/bin/env python3
"""Achieved linear residuals and iteration counts, read from T-Flows run logs.

    python3 analyze_linear_residuals.py LABEL=LOG[,LOG...] [...]
            [--tol U=1e-8,T=1e-8,PP=1e-8,VOF=1e-8] [--cap U=5,T=5,PP=200,VOF=6]

Every outer iteration of a time step prints one line

    # U : n res | V : n res | W : n res | PP : n res | VOL : res | T : n res #

and the first iteration a second one for the volume fraction,  # VOF : n res.
VOL is not a solver residual: it is the largest volume imbalance of a cell
after the velocity correction, |sum of face fluxes - source| dt / V_cell
(Correct_Velocity.f90), the fraction of a cell volume left unbalanced per step.
The outer loop stops when max(VOL, res_U, res_V, res_W, res_T) falls below
TOLERANCE_FOR_SIMPLE_ALGORITHM (Iter_Mod/Max_Fields_Residual.f90), but not
before MIN_SIMPLE_ITERATIONS.
`n` is the number of Krylov iterations the solver took and `res` the residual
it stopped at.  The residual is the one of Native_Mod (Bicg.f90, Cg.f90):

    res = sqrt( sum_i ( r_i / a_ii )^2 ),      r = b - A x,

the Euclidean norm over all cells of the residual divided by the diagonal, so
it carries the unit of the unknown (K, m/s, Pa, -) and is an absolute measure
of the error left in it.  The solver stops when res < TOLERANCE_FOR_*_SOLVER
or after MAX_ITERATIONS_FOR_*_SOLVER iterations.

For each log the script reports, for the LAST outer iteration of every time
step and for ALL outer iterations: median, 95th percentile and maximum of the
residual and of the iteration count, the fraction of solves that stopped on
the cap with the residual still above the tolerance (when --tol and --cap are
given), and the number of outer iterations per step.
"""
import re
import sys

import numpy as np

PAIR = re.compile(r"(U|V|W|PP|T|VOF)\s*:\s*(\d+)\s+([0-9.+\-Ee]+)")
VOL = re.compile(r"VOL\s*:\s*([0-9.+\-Ee]+)")
STEP = re.compile(r"Time step\s*:\s*(\d+)")
ITER = re.compile(r"Iteration:\s*(\d+)")


def parse(paths):
    """Return {var: list over steps of list over iterations of (n, res)}."""
    steps = []
    cur = None
    for p in paths:
        with open(p, errors="replace") as f:
            for ln in f:
                if "Time step" in ln:
                    m = STEP.search(ln)
                    if m:
                        cur = {"U": [], "V": [], "W": [], "PP": [], "T": [],
                               "VOF": [], "VOL": []}
                        steps.append(cur)
                    continue
                if cur is None or "#" not in ln or ":" not in ln:
                    continue
                if " U   :" in ln or "# U  " in ln or "PP " in ln:
                    for var, n, r in PAIR.findall(ln):
                        cur[var].append((int(n), float(r)))
                    m = VOL.search(ln)
                    if m:
                        cur["VOL"].append((0, float(m.group(1))))
                elif "VOF :" in ln:
                    for var, n, r in PAIR.findall(ln):
                        cur[var].append((int(n), float(r)))
    return [s for s in steps if s["PP"]]


def stats(x):
    x = np.asarray(x, dtype=float)
    return np.median(x), np.percentile(x, 95), x.max()


def report(label, steps, tol, cap):
    n_it = np.array([len(s["PP"]) for s in steps])
    print(f"\n=== {label}: {len(steps)} time steps, outer iterations per step "
          f"median {np.median(n_it):.0f}, min {n_it.min()}, max {n_it.max()}")
    print(f"  {'eq.':<5}{'which':<7}{'res median':>12}{'res p95':>12}"
          f"{'res max':>12}{'it med':>8}{'it p95':>8}{'it max':>8}"
          f"{'cap & res>tol':>16}")
    for var in ("U", "V", "W", "PP", "T", "VOF"):
        for which in ("last", "all"):
            if var == "VOF" and which == "all":
                continue
            if which == "last":
                rows = [s[var][-1] for s in steps if s[var]]
            else:
                rows = [x for s in steps for x in s[var]]
            if not rows:
                continue
            n = np.array([r[0] for r in rows])
            r = np.array([r[1] for r in rows])
            if not np.any(n) and not np.any(r):
                continue
            rm, rp, rx = stats(r)
            im, ip, ix = stats(n)
            hit = ""
            if var in tol and var in cap:
                bad = np.mean((n >= cap[var]) & (r > tol[var]))
                hit = f"{100 * bad:9.2f} %"
            print(f"  {var:<5}{which:<7}{rm:12.3e}{rp:12.3e}{rx:12.3e}"
                  f"{im:8.0f}{ip:8.0f}{ix:8.0f}{hit:>16}")
    vol = np.array([s["VOL"][-1][1] for s in steps if s["VOL"]])
    if len(vol):
        print(f"  VOL  last   {np.median(vol):12.3e}{np.percentile(vol, 95):12.3e}"
              f"{vol.max():12.3e}   (largest cell volume imbalance after "
              f"the last outer iteration, |sum of face fluxes - source| dt / V)")


def kv(s, cast):
    return {k: cast(v) for k, v in (x.split("=") for x in s.split(","))}


def main():
    args = sys.argv[1:]
    tol, cap, jobs = {}, {}, []
    i = 0
    while i < len(args):
        if args[i] == "--tol":
            tol = kv(args[i + 1], float); i += 2
        elif args[i] == "--cap":
            cap = kv(args[i + 1], int); i += 2
        else:
            jobs.append(args[i]); i += 1
    if not jobs:
        sys.exit(__doc__)
    for j in jobs:
        label, paths = j.split("=", 1)
        report(label, parse(paths.split(",")), tol, cap)


if __name__ == "__main__":
    main()
