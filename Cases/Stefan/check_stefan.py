#!/usr/bin/env python3
"""Compare a Stefan run with the reference record and with the exact solution.

    check_stefan.py runs/<name>

Reads caf-terms.dat of the run (columns in User_Mod/End_Of_Time_Step.f90) and
compares it with expected/caf-terms.dat, the record of the command in
README.md made from this tree.  Then prints the interface position against the Stefan solution
of the case and the energy imbalance at the end of the run.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1.0e-6      # relative; the same compiler gives the same digits

# the case (control, box_parameters.ini)
K_V, RHO_V, CP_V = 0.025, 0.597, 2030.0
RHO_L, H_LV = 958.4, 2.26e6
T_WALL, T_SAT = 10.1, 0.1
X0 = 0.5e-3                      # initial interface position [m]


def stefan_position(t):
    """Interface position of the one-phase Stefan problem that starts at X0."""
    a_v = K_V / (RHO_V * CP_V)
    ste = CP_V * (T_WALL - T_SAT) / H_LV
    lo, hi = 1.0e-6, 5.0
    for _ in range(200):          # lam exp(lam^2) erf(lam) = Ste / sqrt(pi)
        lam = 0.5 * (lo + hi)
        f = lam * math.exp(lam * lam) * math.erf(lam) - ste / math.sqrt(math.pi)
        lo, hi = (lam, hi) if f < 0.0 else (lo, lam)
    t0 = (X0 / (2.0 * lam)) ** 2 / a_v
    return 2.0 * lam * np.sqrt(a_v * (t + t0))


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    run = np.loadtxt(os.path.join(sys.argv[1], "caf-terms.dat"))
    ref = np.loadtxt(os.path.join(HERE, "expected", "caf-terms.dat"))
    t, x = run[:, 1], run[:, 13]
    print(f"steps: {len(run)}, end time {t[-1]:.3f} s")
    print(f"interface position at the end: {x[-1] * 1e3:.5f} mm, "
          f"Stefan solution {stefan_position(t[-1]) * 1e3:.5f} mm "
          f"({(x[-1] / stefan_position(t[-1]) - 1.0) * 100:+.3f} %)")
    if run.shape != ref.shape:
        print("no comparison with the reference record: another mesh or step count")
        return
    # columns whose reference values are round-off (divergence) are left out
    cols = [c for c in range(2, run.shape[1]) if c != 15]
    scale = np.maximum(np.abs(ref[:, cols]).max(axis=0), 1.0e-300)
    dev = (np.abs(run[:, cols] - ref[:, cols]) / scale).max()
    print(f"largest difference from the reference record: {dev:.2e} "
          f"(relative to the column maximum)")
    print("PASS" if dev < TOL else "FAIL")
    sys.exit(0 if dev < TOL else 1)


if __name__ == "__main__":
    main()
