#!/usr/bin/env python3
"""Where the spurious velocity sits and how noisy the curvature is, from one
saved field of a static-droplet run.

    python3 field_stats.py <label> <file.pvtu> [...]
"""
import sys
import numpy as np
import pyvista as pv

R0, KEX = 0.4, 2.0 / 0.4
for lab, f in zip(sys.argv[1::2], sys.argv[2::2]):
    m = pv.read(f)
    cc = m.cell_centers().points
    u = np.linalg.norm(m.cell_data["Velocity [m/s]"], axis=1)
    a = m.cell_data["Vof Sharp [1]"]
    k = m.cell_data["Vof Curvature [1/m]"]
    vol = m.cell_data["Grid Cell Volume [m^3]"]
    r = np.linalg.norm(cc - 1.0, axis=1)
    dw = np.minimum(cc, 2.0 - cc).min(axis=1)
    i = int(np.argmax(u))
    print(f"{lab}: {m.n_cells} cells")
    print(f"  global max |u| {u[i]:.3f} at r = {r[i]:.3f} from the centre, "
          f"{dw[i]:.3f} from a wall, fun = {a[i]:.3f}")
    for name, mask in (("interface cells", (a > 1e-3) & (a < 1 - 1e-3)),
                       ("band |r-R| < R/4", np.abs(r - R0) < 0.25 * R0),
                       ("inside r < 0.3", r < 0.3),
                       (">0.2 from walls", dw > 0.2),
                       ("<0.1 from walls", dw < 0.1)):
        print(f"    {name:<18} n = {mask.sum():6d}  max {u[mask].max():8.3f}  "
              f"p99 {np.percentile(u[mask], 99):8.3f}  median "
              f"{np.median(u[mask]):8.4f}")
    ic = (a > 0.05) & (a < 0.95) & (np.abs(k) > 1e-12)
    kk = np.abs(k[ic])
    print(f"  curvature in cells 0.05 < fun < 0.95 (n = {ic.sum()}): mean "
          f"{kk.mean():.4f} ({100 * (kk.mean() / KEX - 1):+.2f} % of 2/R), "
          f"std {100 * kk.std() / KEX:.2f} %, range "
          f"{100 * (kk.min() / KEX - 1):+.1f} .. {100 * (kk.max() / KEX - 1):+.1f} %")
    cen = (cc * (a * vol)[:, None]).sum(0) / (a * vol).sum()
    print(f"  centroid moved {np.linalg.norm(cen - 1.0) / R0:.4f} R; cell "
          f"volume min/mean {vol.min() / vol.mean():.3f}")
