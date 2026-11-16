#!/usr/bin/env python3
"""Radial profile of the volume fraction by direction, for one saved volume
snapshot of a kinematic expansion run: is liquid left behind the alpha = 1/2
level, and where?   Usage: profile_diag.py RUN_DIR TS"""
import sys, glob, numpy as np, pyvista as pv
d, ts = sys.argv[1], int(sys.argv[2])
H_ARG = float(sys.argv[3]) if len(sys.argv) > 3 else 300e-6 / 74
f = glob.glob(f"{d}/*-ts{ts:06d}.pvtu")
f = [x for x in f if "front" not in x][0]
m = pv.read(f)
names = list(m.cell_data.keys())
key = [n for n in names if "vof" in n.lower() or "Vof" in n][0]
a = np.asarray(m.cell_data[key], float)
c = np.asarray(m.cell_centers().points, float)
vol = np.asarray(m.compute_cell_sizes(length=False, area=False, volume=True).cell_data["Volume"], float)
r = np.linalg.norm(c, axis=1)
ok = r > 0
s = np.where(ok, ((c / np.where(ok, r, 1)[:, None])**4).sum(1), 1.0)
h = H_ARG
ff = glob.glob(f"{d}/*front-ts{ts:06d}.pvtu")[0]
fm = pv.read(ff)
fa = np.asarray(fm.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"], float)
fc = np.asarray(fm.cell_centers().points, float)
fr = np.linalg.norm(fc, axis=1)
fs = ((fc / fr[:, None])**4).sum(1)
rfront = (fa * fr).sum() / fa.sum()
vvap = (vol * (1 - a)).sum()
rvol = (3 * vvap / (4 * np.pi))**(1 / 3)
print(f"{d.split('/')[-1]}  ts {ts}: field '{key}', R_front = {rfront*1e6:.2f} um, R_vol = {rvol*1e6:.2f} um, ratio {rfront/rvol:.4f}")
classes = (("axis      (s > 0.90)", s > 0.90, fs > 0.90),
           ("face diag (0.45-0.55)", (s > 0.45) & (s < 0.55), (fs > 0.45) & (fs < 0.55)),
           ("body diag (s < 0.36)", s < 0.36, fs < 0.36))
print(f"{'direction':<24}{'R_front(dir)':>13}{'liquid behind front':>21}{'thickness 0.01-0.99':>21}")
for lab, msk, fmsk in classes:
    rf = (fa[fmsk] * fr[fmsk]).sum() / fa[fmsk].sum()
    # liquid volume behind the local front, as an equivalent layer thickness
    inside = msk & (r < rf)
    omega_frac = vol[msk & (r < 0.8 * rf)].sum() / (4 / 3 * np.pi * (0.8 * rf)**3)
    liq = (vol[inside] * a[inside]).sum()
    area = 4 * np.pi * rf**2 * omega_frac
    band = msk & (a > 0.01) & (a < 0.99)
    thick = vol[band].sum() / area
    print(f"{lab:<24}{rf*1e6:>10.2f} um{liq/area/h:>16.3f} h{thick/h:>19.2f} h")
print("\nmean alpha against (r - R_front(dir)) / h:")
edges = np.arange(-8.0, 3.01, 1.0)
print(f"{'bin':>12}" + "".join(f"{lab.split('(')[0].strip():>12}" for lab, _, _ in classes))
rows = []
for lab, msk, fmsk in classes:
    rf = (fa[fmsk] * fr[fmsk]).sum() / fa[fmsk].sum()
    x = (r[msk] - rf) / h
    idx = np.digitize(x, edges) - 1
    row = []
    for i in range(len(edges) - 1):
        k = idx == i
        row.append((vol[msk][k] * a[msk][k]).sum() / max(vol[msk][k].sum(), 1e-300))
    rows.append(row)
for i in range(len(edges) - 1):
    print(f"{edges[i]:>5.0f}..{edges[i+1]:<5.0f}" + "".join(f"{rows[j][i]:>12.4f}" for j in range(3)))
