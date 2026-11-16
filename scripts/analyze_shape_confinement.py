#!/usr/bin/env python3
"""dR/ds at matched confinement -- the source of Table tab:shape_anisotropy.

USE THIS, NOT analyze_confinement.py, for anything that feeds that table.
analyze_confinement.py answers a different question (is the deformation caused
by the domain boundary?) and its shape number is an UNWEIGHTED fit over cells
with 0.35 < alpha < 0.65.  That is not the paper's shape metric, and it comes
out 39 to 86 per cent low.  Filling two blank columns with it would have put
numbers in the same rows as values it disagrees with by a factor of three.

The metric here is the paper's own, as defined in check_shape_metric_zero.py:
the area-weighted regression of front radius on the cubic invariant
s = nx^4+ny^4+nz^4, normalised by the area-weighted mean radius.  Its validation
is that it reproduces the published STRUCTURED column -- level 125 returns
-0.1648 against the tabulated -0.1647.

Confinement is R_front/(L/2), the FRONT radius, not the volume-equivalent one.
The two differ by about 5 per cent on these runs, which is the difference
between level 100 covering the R/(L/2) = 0.765 row and appearing to fall short
of it.

Data sources are not interchangeable between levels.  Poly-75's copy under
runs/ is spliced (radius drops 105.6 -> 65.1 um at ts130) and the archived copy
is clean; Poly-100 is the other way round, and its archived volume pieces were
deleted in the July 2026 cleanup, leaving readable .pvtu headers with no data
behind them.  Provenance attaches to each product, not to the directory.

Usage:  python3 analyze_shape_confinement.py
"""
import glob, numpy as np, pyvista as pv

TARGETS = (0.50, 0.60, 0.70, 0.765)
HALF = 150e-6
CASES = [
    ("75",  "/home/user/archive/Scriven-Merlin/Scriven-Poly-75"),
    ("100", "/home/user/runs/tflows-vof/Validation/Scriven/Scriven-Poly-100"),
    ("125", "/home/user/archive/Scriven-Merlin/Scriven-Poly-125"),
    ("150", "/home/user/archive/Scriven-Merlin/Scriven-Poly-150"),
]
PUB = {"75": (None,)*4, "100": (None,)*4,
       "125": (-0.0146, -0.0247, -0.0357, -0.0452),
       "150": (-0.0170, -0.0271, -0.0372, -0.0450)}

def shape(f):
    # some archived fronts have lost their piece files; skip rather than abort
    try:
        m = pv.read(f)
    except Exception:
        return None
    if m.n_cells == 0: return None
    z = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(z.cell_data["Area"], float)
    c = np.asarray(z.cell_centers().points, float)
    c = c - np.average(c, axis=0, weights=a)
    r = np.linalg.norm(c, axis=1); s = ((c/r[:, None])**4).sum(axis=1)
    W = a/a.sum(); sm, rm = (W*s).sum(), (W*r).sum()
    return rm, (W*(s-sm)*(r-rm)).sum()/(W*(s-sm)**2).sum()/rm

out = {}
for lv, d in CASES:
    fs = sorted(glob.glob(f"{d}/*front-ts*.pvtu"),
                key=lambda f: int(f.split("-ts")[1].split(".")[0]))
    rows, bad = [], 0
    for f in fs[::max(1, len(fs)//30)]:   # ~30 points is plenty to interpolate
        o = shape(f)
        if o is None:
            bad += 1; continue
        rows.append((o[0]/HALF, o[1]))
    if not rows:
        print(f"  Poly-{lv}: all {len(fs)} fronts unreadable", flush=True)
        out[lv] = [None]*4; continue
    rows.sort()
    c = np.array([r[0] for r in rows]); s = np.array([r[1] for r in rows])
    vals = []
    for t in TARGETS:
        vals.append(float(np.interp(t, c, s)) if c.min() <= t <= c.max() else None)
    out[lv] = vals
    print(f"  Poly-{lv:<4} {len(rows):>3} fronts ({bad} unreadable), confinement {c.min():.3f}-{c.max():.3f}",
          flush=True)
    for t, v, p in zip(TARGETS, vals, PUB[lv]):
        note = "" if p is None else f"   published {p:+.4f}  ratio {v/p:.2f}" if v else ""
        print(f"      {t:.3f}  " + (f"{v:+.4f}" if v is not None else "  ---  ") + note,
              flush=True)
print("\n  === recomputed polyhedral block ===", flush=True)
for i, t in enumerate(TARGETS):
    cells = []
    for lv in ("75", "100", "125", "150"):
        v = out[lv][i]
        cells.append(f"{v:+.4f}" if v is not None else "  ---  ")
    r = (out["150"][i]/out["75"][i]) if (out["150"][i] and out["75"][i]) else None
    print(f"  {t:<6} " + " ".join(f"{c:>9}" for c in cells)
          + (f"   150/75 = {r:.2f}x" if r else ""), flush=True)
