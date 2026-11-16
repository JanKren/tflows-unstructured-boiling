#!/usr/bin/env python3
"""
Static-droplet benchmark table from the solver's own spurious.dat  (R2 #7).

Supersedes the .vtu route in analyze_static_droplet.py for the summary table.
Both are valid, but they measure the pressure jump differently:

  spurious.dat  p_in, p_out are volume-weighted averages over PURE-phase
                cells (alpha == max / min), written every time step by
                User_Mod/End_Of_Time_Step.f90
  .vtu route    p_in, p_out are split GEOMETRICALLY by radius

Mixing the two in one table would confound mesh topology with the
measurement definition.  The polyhedral runs abort in the finalisation
save path (a defect in the polyhedral branch of Form_Cells_Comm, unrelated
to the physics) and so may leave no .vtu at all, while spurious.dat is
complete for every run -- hence spurious.dat is the common source.

Column layout, from End_Of_Time_Step.f90:
    time   a_vof   u_rms   u_max   p_max-p_min   p_in-p_out

Reported quantities:
    Ca      = mu * max|u| / sigma          spurious current, should -> 0
    dp err  = (p_in - p_out) / (sigma/R) - 1
The droplet is a CYLINDER with its axis along y, so the exact jump is the
cylindrical value sigma/R, not the spherical 2*sigma/R.

Usage:
  python3 analyze_droplet_spurious.py <dir-with-*.dat> [--sigma S] [--mu M]
                                      [--radius R] [--tail-frac F]
"""

import argparse
import glob
import os

import numpy as np

# run name -> (label, representative cell size h)
ORDER = [
    ("cavity32",       "hex  32",   2.0 / 32),
    ("cavity32_corr",  "hex  32 +RC", 2.0 / 32),
    ("cavity64",       "hex  64",   2.0 / 64),
    ("cavity64_corr",  "hex  64 +RC", 2.0 / 64),
    ("cavity128",      "hex 128",   2.0 / 128),
    ("cavity128_corr", "hex 128 +RC", 2.0 / 128),
    ("poly32_off",     "poly 32",   2.0 / 32),
    ("poly32_corr",    "poly 32 +RC", 2.0 / 32),
    ("poly64_off",     "poly 64",   2.0 / 64),
    ("poly64_corr",    "poly 64 +RC", 2.0 / 64),
    # curvature-smoothing series: recomputed curvature (+RC) with N smoothing
    # cycles.  sm0 isolates the correction from the smoothing; the Scriven
    # production runs use 12 cycles on both mesh families.
    ("poly32_rc_sm0",  "poly 32 sm0", 2.0 / 32),
    ("poly32_rc_sm2",  "poly 32 sm2", 2.0 / 32),
    ("poly32_rc_sm4",  "poly 32 sm4", 2.0 / 32),
    ("poly32_rc_sm8",  "poly 32 sm8", 2.0 / 32),
    ("poly32_nocorr_sm12", "poly 32 sm12", 2.0 / 32),
    ("poly64_rc_sm4",  "poly 64 sm4", 2.0 / 64),
    ("hex32_rc_sm2",   "hex  32 sm2", 2.0 / 32),
    ("hex32_nocorr_sm12", "hex  32 sm12", 2.0 / 32),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--sigma", type=float, default=15.0e3)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--radius", type=float, default=0.4)
    ap.add_argument("--tail-frac", type=float, default=0.1,
                    help="average over this trailing fraction of the record, "
                         "so the value is not one noisy final step")
    a = ap.parse_args()

    dp_exact = a.sigma / a.radius          # cylindrical Laplace
    print(f"\n{'='*78}")
    print("  STATIC DROPLET  --  spurious currents and Laplace jump")
    print(f"{'='*78}")
    print(f"  sigma = {a.sigma:g}, mu = {a.mu:g}, R = {a.radius:g}"
          f"   ->  exact dp = sigma/R = {dp_exact:,.1f} Pa")
    print(f"  values averaged over the final {a.tail_frac*100:.0f}% "
          f"of each record\n")
    print(f"  {'case':<14s} {'steps':>6s} {'a_vof':>10s} {'Ca':>11s}"
          f" {'dp [Pa]':>11s} {'dp err':>9s}")
    print("  " + "-" * 66)

    results = {}
    for key, label, h in ORDER:
        fn = os.path.join(a.dir, f"{key}.dat")
        if not os.path.exists(fn):
            print(f"  {label:<14s} {'--':>6s}   (not available yet)")
            continue
        d = np.loadtxt(fn)
        if d.ndim == 1:
            d = d[None, :]
        n = len(d)
        tail = max(1, int(a.tail_frac * n))
        u_max = d[-tail:, 3].mean()
        dp = d[-tail:, 5].mean()
        a_vof = d[-1, 1]
        ca = a.mu * u_max / a.sigma
        err = dp / dp_exact - 1.0
        results[key] = (ca, dp, err, h)
        print(f"  {label:<14s} {n:6d} {a_vof:10.6f} {ca:11.3e}"
              f" {dp:11.1f} {err*100:8.2f}%")

    print("\n  effect of the Rhie--Chow corrections (Ca reduction factor)")
    for base in ("cavity32", "cavity64", "cavity128",
                 "poly32_off", "poly64_off"):
        corr = base.replace("_off", "") + "_corr" if base.endswith("_off") \
            else base + "_corr"
        if base in results and corr in results:
            f = results[base][0] / results[corr][0]
            print(f"    {base:<14s} -> {corr:<14s} {f:8.1f}x")
    print()


if __name__ == "__main__":
    main()
