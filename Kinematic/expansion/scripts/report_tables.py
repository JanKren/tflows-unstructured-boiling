#!/usr/bin/env python3
"""Compact tables for the note and the report from the CSV written by
analyze_kinematic_expansion.py: per run, at t = 0.2, 0.4, 0.6, 0.8, 1.0 ms,
front radius, R/h, body - axis, shape m = 4, phi_4, m = 2, m = 6."""
import csv, numpy as np
CSV = "/home/user/papers/UnstructuredBoiling/results/caf-kinematic-expansion.csv"
runs = {}
for r in csv.DictReader(open(CSV)):
    runs.setdefault(r["run"], []).append(r)
order = ["st075-dt2-step", "st100-dt2-step", "st125-dt2-step", "st150-dt2-step", "st075-dt1-step",
         "st075-dt2-step-samp", "st075-dt2-sol", "st100-dt2-sol", "st125-dt2-sol", "st075-dt2-mask",
         "st075-dt2-mix", "po075-dt2-step", "po100-dt2-step", "po125-dt2-step", "po075-dt2-sol", "po100-dt2-sol"]
T = [0.2e-3, 0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3]
print("| run | " + " | ".join(f"t = {t*1e3:.1f} ms" for t in T) + " | phi4 range | max m2, m6 [%] |")
print("|---|" + "---|" * (len(T) + 2))
for name in order:
    if name not in runs: continue
    d = runs[name]
    cells = []
    for t in T:
        m = [r for r in d if abs(float(r["t"]) - t) < 1e-9]
        if not m: cells.append("--"); continue
        r = m[0]
        cells.append(f"R {float(r['R'])*1e6:.1f}, R/h {float(r['Rh']):.1f}: {float(r['dev']):+.2f} / {float(r['m4']):.3f}")
    ph = [float(r["phi4"]) for r in d if float(r["m4"]) > 0.1]
    phs = f"{min(ph):.2f}-{max(ph):.2f}" if ph else "scattered"
    mm = max(max(float(r["m2"]), float(r["m6"])) for r in d)
    print(f"| `{name}` | " + " | ".join(cells) + f" | {phs} | {mm:.3f} |")
print("\n(each cell: front radius [um], R/h: body - axis [%] / shape m = 4 [%];", "last saved front of unfinished runs:",
      ", ".join(f"{n}: ts {max(int(r['ts']) for r in runs[n])}" for n in order if n in runs and max(float(r['t']) for r in runs[n]) < 0.999e-3), ")")
