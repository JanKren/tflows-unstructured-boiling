#!/usr/bin/env python3
"""Kinematic radial-expansion test: does the transport alone deform the bubble?

The volume fraction of a sphere is advected by the production CICSAM in the
exact flux field of a point source at the bubble centre, on the Scriven cube
and meshes, with no energy equation, no mass transfer and no momentum or
pressure solve.  The source strength follows the Scriven law, so the exact
sphere has R_S(t) = 2 beta sqrt(D_l (t + t0)) and the interface Courant number
has the production history.  Run directories and the solver hooks are under
/home/user/runs/caf-review-2026-11/kinematic (results/caf-kinematic-advection.md).

Run names are <mesh>-dt<step>-<model>, the model being where the point-source
flux is switched off:

  sol    nowhere: the liquid velocity continued through the interface band,
         which is solenoidal;
  step   on faces whose two cells have a mean volume fraction below 1/2: the
         velocity outside the computed bubble and none inside it;
  mask   on faces inside the exact sphere R_S(t);
  mix    scaled by the mean volume fraction of the two cells (a ramp).

For every saved front this prints the measures of the manuscript, computed by
the manuscript's own functions:

  * body - axis: 100 (R_body / R_axis - 1) from the area-weighted regression of
    the front radius on the cubic invariant s = nx^4 + ny^4 + nz^4
    (plot_noadvection.front / fit, the metric of Figure fig:noadvection and of
    Table tab:hyperc_ladder);
  * the azimuthal spectrum of the front radius, area-weighted least squares on
    every front element (analyze_hyperc_shape_m4.spectrum), with m = 2 and
    m = 6 as the cubic-symmetry check.

Screens, printed before any shape number: the time column is monotone, the
sphere volume follows the exact increment, the volume fraction stays in [0, 1],
the discrete divergence of the imposed flux, and the corrector sweep counts of
the log.

Usage:  python3 analyze_kinematic_expansion.py [--csv results/caf-kinematic-expansion.csv]
"""

import argparse
import csv
import glob
import json
import os
import re

import numpy as np

from plot_noadvection import front, fit
from analyze_hyperc_shape_m4 import spectrum

ROOT = "/home/user/runs/caf-review-2026-11/kinematic"
BETA, D_L, R0 = 4.0647, 1.680e-7, 50.0e-6
T0 = R0**2 / (4.0 * BETA**2 * D_L)
L = 300.0e-6

# run -> (family, level, nominal structured h [m], dt [s], note)
H = {75: L / 74, 100: L / 99, 125: L / 124, 150: L / 149}
VARIANT = {
    "sol":  "solenoidal: the liquid velocity (point source) everywhere",
    "mix":  "one-fluid velocity vof*u_l, face vof = mean of the two cells",
    "step": "one-fluid velocity, step at vof = 1/2",
    "mask": "no flux inside the exact sphere R_S(t)",
}


def find_runs():
    """(name, family, level, dt, note) of every run directory present."""
    out = []
    for d in sorted(glob.glob(os.path.join(ROOT, "[sp][to][0-9]*-dt*"))):
        name = os.path.basename(d)
        m = re.match(r"(st|po)(\d+)-dt(\w+?)-(\w+?)(-samp)?$", name)
        if not m:
            continue
        ctl = open(os.path.join(d, "control")).read()
        dt = float(re.search(r"^\s*TIME_STEP\s+(\S+)", ctl, re.M).group(1))
        note = VARIANT.get(m.group(4), m.group(4))
        if m.group(5):
            note += "; sub-sampled initial field"
        out.append((name, "struct" if m.group(1) == "st" else "poly",
                    int(m.group(2)), dt, note))
    fam = {"struct": 0, "poly": 1}
    return sorted(out, key=lambda r: (r[0].split("-", 2)[2], fam[r[1]],
                                      r[3] != 2e-6, r[2], r[3]))


RUNS = find_runs()


def r_scriven(t):
    return 2.0 * BETA * np.sqrt(D_L * (t + T0))


def screen(d):
    """The integrity screens of one run; returns a dict or None."""
    bench = os.path.join(d, "bench-data.dat")
    init = os.path.join(d, "kin-init.dat")
    if not (os.path.exists(bench) and os.path.exists(init)):
        return None
    b = np.loadtxt(bench)
    v0 = np.loadtxt(init)[2]
    t, vol = b[:, 0], b[:, 2]
    out = dict(n=len(t), t_end=t[-1])
    out["monotone"] = bool((np.diff(t) > 0).all())
    exact = 4.0 * np.pi / 3.0 * (r_scriven(t)**3 - R0**3)
    out["vol_defect"] = float(np.abs((vol - v0) - exact).max() / exact[-1])
    out["vol_defect_end"] = float(((vol[-1] - v0) - exact[-1]) / exact[-1])
    out["rvol0_um"] = (3.0 * v0 / (4.0 * np.pi))**(1.0 / 3.0) * 1e6
    out["fmin"] = float(b[:, 5].min())
    out["fmax1"] = float(b[:, 6].max())
    out["rfront_end_um"] = float(b[-1, 4] * 1e6)
    out["rvol_end_um"] = float(b[-1, 3] * 1e6)
    flux = os.path.join(d, "kin-flux.dat")
    if os.path.exists(flux):
        f = np.loadtxt(flux)
        out["div_defect"] = float(f[:, 5].max())
        out["src_sum"] = float(f[0, 7])
    # corrector sweeps and unbounded cells, from the log itself
    log = open(os.path.join(d, "run.log"), errors="replace").read()
    sw = re.findall(r"VOF_CORR: scheme beta\s+sweeps=(\d+) of max (\d+)\s+"
                    r"unbounded lo=(\d+) hi=(\d+)", log)
    out["nsub"] = len(sw)
    out["first"] = sum(1 for x in sw if int(x[0]) == 1)
    out["exhausted"] = sum(1 for x in sw if int(x[0]) > int(x[1]))
    out["unb_max"] = max([int(x[2]) + int(x[3]) for x in sw] or [0])
    out["subcycled"] = os.path.exists(os.path.join(d, "alert-dt-vof.dat"))
    out["finished"] = os.path.exists(os.path.join(d, "run.exit")) and \
        open(os.path.join(d, "run.exit")).read().strip() == "exit 0"
    return out


def shapes(d, dt, h):
    """Shape measures of every saved front of one run.  The per-front numbers
    are cached beside the run (shape_cache.json), keyed by file and its
    modification time, because reading a few hundred fronts takes minutes."""
    cpath = os.path.join(d, "shape_cache.json")
    cache = json.load(open(cpath)) if os.path.exists(cpath) else {}
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "*front-ts*.pvtu")),
                    key=lambda p: int(re.search(r"-ts(\d+)", p).group(1))):
        ts = int(re.search(r"-ts(\d+)", f).group(1))
        key = f"{os.path.basename(f)}:{os.path.getmtime(f):.0f}"
        if key not in cache:
            got = front(f)
            if got is None:
                continue
            a, r, s = got
            c0, c1, mean = fit(a, r, s)
            axis, body = c0 + c1, c0 + c1 / 3.0
            r0, sp, n = spectrum(f)
            cache[key] = dict(R=float(mean),
                              dev=float(100.0 * (body / axis - 1.0)),
                              dRds=float(c1 / mean), n=int(n),
                              m4=float(sp[4][0]), phi4=float(sp[4][1]),
                              m8=float(sp[8][0]), m2=float(sp[2][0]),
                              m6=float(sp[6][0]))
        c = cache[key]
        rows.append(dict(ts=ts, t=ts * dt, Rh=c["R"] / h,
                         Rex=r_scriven(ts * dt), **c))
    json.dump(cache, open(cpath, "w"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/caf-kinematic-expansion.csv")
    args = ap.parse_args()

    print("=" * 100)
    print("SCREENS")
    print("=" * 100)
    print(f"{'run':<17}{'rows':>6}{'t_end ms':>10}{'mono':>6}{'R_vol(0)':>10}"
          f"{'dV/exact max':>14}{'div defect':>12}{'min a':>10}{'max a-1':>10}"
          f"{'done':>6}")
    scr = {}
    for name, fam, lev, dt, note in RUNS:
        d = os.path.join(ROOT, name)
        s = screen(d)
        if s is None:
            continue
        scr[name] = s
        print(f"{name:<17}{s['n']:>6}{s['t_end']*1e3:>10.4f}"
              f"{str(s['monotone']):>6}{s['rvol0_um']:>10.4f}"
              f"{s['vol_defect']:>14.2e}{s.get('div_defect', np.nan):>12.2e}"
              f"{s['fmin']:>10.1e}{s['fmax1']:>10.1e}"
              f"{str(s['finished']):>6}")
        if s["nsub"]:
            print(f"{'':<17}   corrector: bounded after the first sweep on "
                  f"{100.0*s['first']/s['nsub']:.1f} % of {s['nsub']} VOF "
                  f"steps, exhausted on {100.0*s['exhausted']/s['nsub']:.1f} %"
                  f" (at most {s['unb_max']} cells then clipped)")
        if s["subcycled"]:
            print(f"{'':<17}   NOTE: alert-dt-vof.dat exists (sub-cycling)")

    print()
    print("=" * 100)
    print("SHAPE  (body - axis from the cubic-invariant regression; azimuthal"
          " spectrum of the front radius)")
    print("=" * 100)
    allrows = []
    for name, fam, lev, dt, note in RUNS:
        d = os.path.join(ROOT, name)
        if name not in scr:
            continue
        rows = shapes(d, dt, H[lev])
        if not rows:
            continue
        print(f"\n{name}   {fam} {lev}^3, dt = {dt*1e6:g} us   {note}")
        print(f"{'ts':>7}{'t ms':>7}{'R um':>9}{'R_ex um':>9}{'R/h':>7}"
              f"{'body-axis %':>13}{'dR/ds':>9}{'m4 %':>8}{'phi4':>7}"
              f"{'m8 %':>7}{'m2 %':>7}{'m6 %':>7}{'n':>7}")
        for r in rows:
            print(f"{r['ts']:>7}{r['t']*1e3:>7.3f}{r['R']*1e6:>9.2f}"
                  f"{r['Rex']*1e6:>9.2f}{r['Rh']:>7.1f}{r['dev']:>+13.3f}"
                  f"{r['dRds']:>+9.4f}{r['m4']:>8.3f}{r['phi4']:>7.2f}"
                  f"{r['m8']:>7.3f}{r['m2']:>7.3f}{r['m6']:>7.3f}{r['n']:>7}")
            allrows.append(dict(run=name, family=fam, level=lev, dt=dt,
                                note=note, **r))

    if allrows:
        os.makedirs(os.path.dirname(args.csv) or ".", exist_ok=True)
        keys = ["run", "family", "level", "dt", "note", "ts", "t", "R", "Rex",
                "Rh", "dev", "dRds", "m4", "phi4", "m8", "m2", "m6", "n"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in allrows:
                w.writerow({k: r[k] for k in keys})
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
