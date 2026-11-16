#!/usr/bin/env python3
"""The cells that carry the global velocity maximum of the published polyhedral
run: how many, where, how small."""
import numpy as np, pyvista as pv
m = pv.read("/home/user/runs/Droplet3D/poly64_result/drop3d_poly64_dual-ts033255.pvtu")
cc = m.cell_centers().points
u = np.linalg.norm(m.cell_data["Velocity [m/s]"], axis=1)
a = m.cell_data["Vof Sharp [1]"]
vol = m.cell_data["Grid Cell Volume [m^3]"]
p = m.cell_data["Pressure [Pa]"]
r = np.linalg.norm(cc - 1.0, axis=1)
dw = np.minimum(cc, 2.0 - cc).min(axis=1)
far = r > 0.8
print("cells further than 0.8 from the centre: %d, of which |u| > 5: %d, > 1: %d"
      % (far.sum(), (far & (u > 5)).sum(), (far & (u > 1)).sum()))
for i in np.argsort(np.where(far, u, 0))[::-1][:12]:
    corner = np.minimum(cc[i], 2 - cc[i])
    print("  |u| %6.2f  pos %s  dist to the three nearest walls %s  vol/mean %.3f  fun %.3f  p %.0f"
          % (u[i], np.round(cc[i], 3), np.round(np.sort(corner), 3), vol[i] / vol.mean(), a[i], p[i]))
print("pressure outside (fun < 1e-6, r > 0.6): mean %.1f  std %.1f" % (p[(a < 1e-6) & (r > 0.6)].mean(), p[(a < 1e-6) & (r > 0.6)].std()))
print("max |u| over all cells excluding those within 0.05 of a wall: %.3f" % u[dw > 0.05].max())
