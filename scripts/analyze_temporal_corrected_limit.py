#!/usr/bin/env python3
"""Structured spatial limit after removing the temporal residual level by level.

Section 3.2.1 quotes a doubly extrapolated structured mean limit.  That number
had been computed ad hoc and never committed, so it could not be reproduced
after the renormalisation of 2026-08.  This script is that calculation, written
down.

Method.  Table 3 gives the mean radius error at dt = 2 us on the four levels.
Section 3.2.1 gives the shift in R(1 ms) on halving the step, level by level, as
a fraction of the dt = 2 us radius.  Those shifts are ratios of two measured
radii, so they are invariant under the choice of R_Scriven and did not move with
the renormalisation.  With the first order in dt that the splitting argument
fixes, the deficit at dt = 2 us is twice the 2 -> 1 us shift.  The 75^3 level is
taken from its own three-step extrapolation in Table 2 instead, that level
having been run at 0.1 us as well.

Eq. (richardson) is then refitted to the corrected sequence.  Both the exponent
and the limit are reported, with a drop-one jackknife over the four exactly
determined three-point fits, as in plot_scriven_extrapolation.py.
"""
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit

H = np.array([4.05, 3.03, 2.42, 2.01])      # nominal cell size [um]
HT = H / H[-1]

# Table 3, structured mean error at dt = 2 us, after the 2026-08 renormalisation
E_DT2 = np.array([1.96, 7.11, 9.57, 10.75])   # Table 2, R_Scriven = 116.642 um

# Section 3.2.1, shift in R(1 ms) on halving the step, as a fraction of the
# dt = 2 us radius.  Invariant under R_Scriven.
SHIFT = np.array([1.63, 2.83, 3.78, 4.59]) / 100.0

# Table 2, structured dt -> 0 mean error on 75^3 from its own three steps
E_LIMIT_75 = 5.48

def model(x, a, p, b):
    return a * x**p + b

def fit(y):
    q, _ = curve_fit(model, HT, y, p0=[-0.05, 2.0, 1.1], maxfev=400000)
    j = [curve_fit(model, HT[list(i)], y[list(i)], p0=q, maxfev=400000)[0]
         for i in combinations(range(4), 3)]
    return q, (min(x[1] for x in j), max(x[1] for x in j)), \
              (min(x[2] for x in j), max(x[2] for x in j))

def main():
    deficit = 2.0 * SHIFT * (1.0 + E_DT2 / 100.0) * 100.0
    corr = E_DT2 + deficit
    corr[0] = E_LIMIT_75
    print("level        dt=2us      deficit    corrected")
    for h, e, d, c in zip(H, E_DT2, deficit, corr):
        print(f"  h={h:4.2f}   {e:+7.2f}    {d:+7.2f}    {c:+7.2f}")
    def fixed(x, a, b):            # the order used in the manuscript
        return a * x**2.0 + b
    print()
    for lab, y in (("uncorrected", 1 + E_DT2 / 100.0), ("dt-corrected", 1 + corr / 100.0)):
        q, _ = curve_fit(fixed, HT, y, p0=[-0.05, 1.1], maxfev=400000)
        jb = [curve_fit(fixed, HT[list(i)], y[list(i)], p0=q, maxfev=400000)[0][1]
              for i in __import__("itertools").combinations(range(4), 3)]
        print(f"  p = 2 fixed   {lab:14s} b = {q[1]:.4f} ({(q[1]-1)*100:+.1f}%)"
              f"  range {min(jb):.3f} to {max(jb):.3f}")
    print()
    for lab, y in (("uncorrected", 1 + E_DT2 / 100.0), ("dt-corrected", 1 + corr / 100.0)):
        q, pr, br = fit(y)
        print(f"\n{lab:>12}:  p = {q[1]:.2f} [{pr[0]:.2f}, {pr[1]:.2f}]   "
              f"b = {q[2]:.4f} [{br[0]:.3f}, {br[1]:.3f}]")

if __name__ == "__main__":
    main()
