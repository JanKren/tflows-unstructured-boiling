#!/usr/bin/env python3
"""Is the review-branch binary the published scheme?  Compares the first rows
of the pilots with the published traces (same mesh, control and time step) and
with the October T-Flows-jan runs at four cycles.

    python3 compare_pilot.py
"""
import os
import numpy as np

P = "/home/user/runs/caf-review-2026-11/static-drop/traces_pilot"
PAIRS = [("m_hex_s12_p32", "/home/user/runs/Droplet3D/hex64_result/spurious.dat",
          "published hex64, 12 cycles (paper branch)"),
         ("m_poly_s4_p32", "/home/user/runs/Droplet3D/poly64_result/spurious.dat",
          "published poly64, 4 cycles (paper branch)"),
         ("m_hex_s4_p32", "/home/user/runs/FlowBoiling/static_sphere_merlin_20261008_v1/"
          "full/cases/hex64_cicsam_csf/spurious.dat",
          "October hex64, 4 cycles (T-Flows-jan)"),
         ("m_poly_s4_p32", "/home/user/runs/FlowBoiling/static_sphere_merlin_20261008_v1/"
          "full/cases/poly64_cicsam_csf/spurious.dat",
          "October poly64, 4 cycles (T-Flows-jan)")]
names = ["time", "volume", "u_rms", "u_max", "p_max-p_min", "p_in-p_out"]
for pilot, ref, lab in PAIRS:
    f = os.path.join(P, pilot, "spurious.dat")
    if not os.path.exists(f):
        print(pilot, "missing")
        continue
    a = np.loadtxt(f)
    b = np.loadtxt(ref, max_rows=len(a))
    n = min(len(a), len(b))
    rel = np.abs(a[:n] / b[:n] - 1.0)
    print(f"{pilot} against {lab}, {n} rows")
    for j, nm in enumerate(names):
        print(f"    {nm:<12} max rel diff {rel[:, j].max():.3e}   at row "
              f"{n}: {a[n-1, j]:.6e} vs {b[n-1, j]:.6e}")
