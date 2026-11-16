#!/usr/bin/env python3
"""Standard advection tests with the production CICSAM, beside the literature.

Run directories: /home/user/runs/caf-review-2026-11/kinematic/adv (set up by
setup_adv.sh there; results/caf-kinematic-advection.md).  The face fluxes are
imposed and discretely solenoidal, and momentum, pressure and energy are not
solved, so these are tests of the volume-fraction transport alone.

  disc   two-dimensional disc, R = 0.25, in the uniform flow U = (1, 0.5) for
         t = 4 on square cells (one cell deep): the test of Roenby, Bredmose
         and Jasak, R. Soc. Open Sci. 3 (2016) 160405, Section 4.1, whose
         Table 1 lists isoAdvector, MULES, HRIC and their own CICSAM.
  sph    sphere, R = 0.25, translated along the body diagonal of a Cartesian
         mesh, U = (1, 1, 1), over sqrt(3) (one way) or there and back.
  def    sphere, R = 0.15, in the deformation field of LeVeque with period
         T = 6, returning at t = 3 (their Section 4.4, Table 8, isoAdvector).

Error norms as defined there: E1 = sum V |a - a_exact| / sum V a_exact (a the
volume fraction of the advected body), dV_rel the relative change of its
volume, min(a) and 1 - max(a), and dW_rel, the volume of the cells with
0.01 < a < 0.99 relative to that of the exact solution, minus one.  (Roenby et
al. measure the volume between the two isosurfaces; cells are counted here.)

The Courant number is Co = dt * sum_f max(F_f, 0) / V, i.e. (|u|+|v|+|w|) dt/h
in a uniform flow.  The solver sub-cycles the volume-fraction step above
MAX_COURANT_VOF (0.25 by default); runs named *-sub keep that default, the
others raise it to 1 so that the step is taken at the stated Co.

Usage:  python3 analyze_advection_tests.py
"""

import glob
import os
import re

import numpy as np

ROOT = "/home/user/runs/caf-review-2026-11/kinematic/adv"

# Roenby et al. (2016), Table 1(a)-(e): (Nx, Co) -> isoAdvector, MULES, HRIC, CICSAM
ROENBY_T1 = {
    "E1":   {(10, 0.5): (0.11, 0.76, 0.86, 1.0),
             (20, 0.5): (0.035, 0.43, 0.56, 0.73),
             (40, 0.5): (0.021, 0.32, 0.38, 0.5),
             (40, 0.2): (0.014, 0.2, 0.22, 0.27),
             (40, 0.1): (0.017, 0.18, 0.19, 0.15)},
    "dV":   {(10, 0.5): (0, -0.0019, -0.026, 0.0024),
             (20, 0.5): (-6.7e-14, -0.00018, -0.001, 0.0042),
             (40, 0.5): (-3.6e-10, -1.1e-05, -1e-05, 0.0012),
             (40, 0.2): (-4.7e-13, -5.2e-10, 1.5e-05, 0.00044),
             (40, 0.1): (-2.1e-13, -2.7e-13, -8.8e-05, 0.00019)},
    "min":  {(10, 0.5): (-3.9e-15, 0, 0, -0.084),
             (20, 0.5): (0, 0, 0, -0.24),
             (40, 0.5): (0, 0, -5.9e-13, -0.23),
             (40, 0.2): (0, 0, 0, -0.091),
             (40, 0.1): (0, 0, 0, -0.047)},
    "dW":   {(10, 0.5): (-0.0003, 0.34, 1.4, 0.21),
             (20, 0.5): (-0.0051, 0.44, 1.9, 0.2),
             (40, 0.5): (0.026, 0.51, 3, 0.27),
             (40, 0.2): (0.012, 0.32, 1.7, 0.14),
             (40, 0.1): (0.0065, 0.35, 1.2, 0.036)},
}
# Roenby et al. (2016), Table 8: Nx -> isoAdvector E1 at Co = 0.5
ROENBY_T8 = {64: 0.22, 128: 0.047, 256: 0.012}


def read(run):
    d = os.path.join(ROOT, run)
    fin = os.path.join(d, "kin-final.dat")
    if not os.path.exists(fin):
        return None
    v = np.loadtxt(fin)
    out = dict(t=v[0], e1_init=v[1], e1_exact=v[2], dv=v[3], amin=v[4],
               amax1=v[5], dw=v[6])
    # the advected body is 1 - vof: its minimum is 1 - max(vof)
    out["body_min"] = -v[5]
    out["body_1max"] = v[4]
    log = open(os.path.join(d, "run.log")).read()
    sw = re.findall(r"VOF_CORR: scheme beta\s+sweeps=(\d+) of max (\d+)\s+"
                    r"unbounded lo=(\d+) hi=(\d+)", log)
    out["nsub"] = len(sw)
    out["exhausted"] = sum(1 for s in sw if int(s[0]) > int(s[1]))
    out["unb"] = sum(int(s[2]) + int(s[3]) for s in sw)
    ctl = open(os.path.join(d, "control")).read()
    out["nsteps"] = int(re.search(r"NUMBER_OF_TIME_STEPS\s+(\d+)", ctl).group(1))
    out["dt"] = float(re.search(r"^\s*TIME_STEP\s+(\S+)", ctl, re.M).group(1))
    out["comax"] = float(re.search(r"^\s*MAX_COURANT_VOF\s+(\S+)", ctl, re.M).group(1))
    fl = np.loadtxt(os.path.join(d, "kin-flux.dat"))
    out["defect"] = float(fl[:, 5].max())
    return out


def co_of(run):
    m = re.search(r"-co(\d+)", run)
    s = m.group(1)
    return float("0." + s[1:]) if s.startswith("0") else float(s)


def main():
    runs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "*-co*")))

    print("=" * 108)
    print("DISC IN UNIFORM FLOW U = (1, 0.5), t = 4  --  Roenby et al. (2016)"
          " Section 4.1, Table 1")
    print("=" * 108)
    print(f"{'run':<18}{'(Nx,Co)':>10}{'substeps':>9}{'E1':>9}{'dV_rel':>11}"
          f"{'min(a)':>10}{'1-max(a)':>10}{'dW_rel':>8}{'unb':>6}  |"
          f"{'  Roenby E1: isoAdv  MULES   HRIC  CICSAM':>40}")
    for r in runs:
        if not r.startswith("disc"):
            continue
        o = read(r)
        if o is None:
            continue
        nx = int(re.match(r"disc(\d+)", r).group(1))
        co = co_of(r)
        ref = ROENBY_T1["E1"].get((nx, co))
        refs = ("" if ref is None else
                f"{ref[0]:>12}{ref[1]:>7}{ref[2]:>7}{ref[3]:>8}")
        print(f"{r:<18}{'(%d,%g)' % (nx, co):>10}{o['nsub']/o['nsteps']:>9.2f}"
              f"{o['e1_exact']:>9.3f}{o['dv']:>+11.1e}{o['body_min']:>10.1e}"
              f"{o['body_1max']:>10.1e}{o['dw']:>8.2f}{o['unb']:>6}  |{refs}")
    print("\n  Roenby et al. Table 1 for their CICSAM, same rows: "
          "dV_rel 2.4e-3, 4.2e-3, 1.2e-3, 4.4e-4, 1.9e-4;  "
          "min(a) -0.084, -0.24, -0.23, -0.091, -0.047;\n"
          "  dW_rel 0.21, 0.2, 0.27, 0.14, 0.036.")

    for key, title in (("sph", "SPHERE ALONG THE BODY DIAGONAL, U = (1,1,1), "
                                "R = 0.25, R/h = 8"),
                       ("def", "SPHERE IN THE DEFORMATION FIELD, T = 6, t = 3, "
                               "R = 0.15  --  Roenby et al. Section 4.4")):
        print()
        print("=" * 108)
        print(title)
        print("=" * 108)
        print(f"{'run':<18}{'mesh':>6}{'Co':>6}{'substeps':>9}{'E1 exact':>10}"
              f"{'E1 init':>9}{'dV_rel':>11}{'min(a)':>10}{'1-max(a)':>10}"
              f"{'dW_rel':>8}{'unb':>7}{'defect':>10}")
        for r in runs:
            if not re.match(r"p?" + key, r):
                continue
            o = read(r)
            if o is None:
                continue
            mesh = "poly" if r.startswith("p") else "hex"
            e1x = f"{o['e1_exact']:>10.3f}" if o["e1_exact"] >= 0 else f"{'--':>10}"
            print(f"{r:<18}{mesh:>6}{co_of(r):>6g}{o['nsub']/o['nsteps']:>9.2f}"
                  f"{e1x}{o['e1_init']:>9.3f}{o['dv']:>+11.1e}"
                  f"{o['body_min']:>10.1e}{o['body_1max']:>10.1e}"
                  f"{o['dw']:>8.2f}{o['unb']:>7}{o['defect']:>10.1e}")
        if key == "def":
            print("\n  Roenby et al. Table 8, isoAdvector at Co = 0.5: "
                  "E1 = 0.22 (Nx = 64), 0.047 (128), 0.012 (256).")
        else:
            print("\n  E1 exact: against the sphere at its final position "
                  "(one-way runs).  E1 init: against the initial field\n"
                  "  (there-and-back runs, *-rev, where the two coincide).")


if __name__ == "__main__":
    main()
