#!/usr/bin/env python3
"""Body - axis of the kinematic runs interpolated to common radii, beside the
coupled runs (numbers of the coupled runs computed with the same functions
from Data/unif-transport and Data/Aniso/aniso-*-off).  Reads the CSV written
by analyze_kinematic_expansion.py."""
import csv, sys, numpy as np
CSV = "/home/user/papers/UnstructuredBoiling/results/caf-kinematic-expansion.csv"
runs = {}
for r in csv.DictReader(open(CSV)):
    runs.setdefault(r["run"], []).append((float(r["R"]) * 1e6, float(r["dev"]), float(r["m4"]), float(r["phi4"]), float(r["Rh"])))
RT = [60, 70, 80, 90, 100, 110, 118]
print("body - axis [%] at fixed radius R [um]")
print(f"{'run':<22}" + "".join(f"{x:>8}" for x in RT))
for name in sorted(runs):
    d = sorted(runs[name]); R = np.array([x[0] for x in d]); v = np.array([x[1] for x in d])
    print(f"{name:<22}" + "".join((f"{np.interp(x, R, v):>+8.2f}" if R.min() <= x <= R.max() else f"{'--':>8}") for x in RT))
print()
CU = [(66.39, 4.235), (81.31, 4.909), (95.05, 5.152), (107.87, 5.857), (119.94, 7.126)]
print("against the coupled run at a uniform flux (75^3, dt = 1 us):")
for name in ("st075-dt2-step", "st075-dt1-step", "st075-dt01-step"):
    if name not in runs: continue
    d = sorted(runs[name]); R = np.array([x[0] for x in d]); v = np.array([x[1] for x in d])
    print(f"  {name:<18}" + "  ".join(f"R={x:6.2f}: {np.interp(x, R, v):+.2f} vs {c:+.2f} ({100*np.interp(x, R, v)/c:.0f} %)" for x, c in CU if R.min() <= x <= R.max()))
