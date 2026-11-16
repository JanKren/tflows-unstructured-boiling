#!/usr/bin/env python3
"""Compare a radial-expansion run with the reference record and the exact radius.

    check_expansion.py <run folder> [expected/st075-dt2-step]

Reads bench-data.dat and kin-flux.dat of the run (columns in
User_Mod/End_Of_Time_Step.f90 and Beginning_Of_Compute_Vof.f90), prints the
radius of the vapour volume and of the front against the exact radius at the
end, and, when a reference record is given or the default one fits, the
largest difference from it.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1.0e-9      # relative; the same compiler and decomposition give the same digits


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    run = sys.argv[1]
    ref_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, "expected", "st075-dt2-step")
    bench = np.loadtxt(os.path.join(run, "bench-data.dat"))
    flux = np.loadtxt(os.path.join(run, "kin-flux.dat"))
    r_exact = flux[-1, 3]
    print(f"steps: {len(bench)}, end time {bench[-1, 0]:.4e} s")
    print(f"exact radius            {r_exact * 1e6:9.3f} um")
    print(f"radius of the volume    {bench[-1, 3] * 1e6:9.3f} um "
          f"({(bench[-1, 3] / r_exact - 1.0) * 100:+.2f} %)")
    print(f"radius of the front     {bench[-1, 4] * 1e6:9.3f} um "
          f"({(bench[-1, 4] / r_exact - 1.0) * 100:+.2f} %)")
    print(f"largest flux defect     {flux[:, 5].max():.2e} of a cell volume per step")
    ref = np.loadtxt(os.path.join(ref_dir, "bench-data.dat"))
    if ref.shape != bench.shape:
        print("no comparison with the reference record: another step count")
        return
    cols = [0, 1, 2, 3, 4, 7, 8]          # 5 and 6 are round-off of the bounds
    scale = np.abs(ref[:, cols]).max(axis=0)
    dev = (np.abs(bench[:, cols] - ref[:, cols]) / scale).max()
    print(f"largest difference from {os.path.basename(ref_dir)}: {dev:.2e} "
          f"(relative to the column maximum)")
    print("PASS" if dev < TOL else "DIFFERENT")
    sys.exit(0 if dev < TOL else 1)


if __name__ == "__main__":
    main()
