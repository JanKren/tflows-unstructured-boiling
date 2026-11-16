#!/usr/bin/env python3
"""Distribution of the CICSAM blend weight gamma_f on the structured Scriven run.

    python3 analyze_gamma_distribution.py

The question this answers: is gamma_f close to one, on the grounds that the
interface normal and the flow are parallel for a radially expanding bubble?

It is not, and the reason is what gamma_f is actually measured against.
Predict_Beta.f90:61-108 computes

    dotprod  = grad_alpha(donor) . (dx(s), dy(s), dz(s))
    ang      = dotprod / (|grad_alpha(donor)| * d(s))
    gamma_f  = min(ang**2, 1)

so the angle is between the interface normal and the CELL-TO-CELL CONNECTION
VECTOR of the face -- a property of the mesh -- not between the normal and the
velocity. Were it the latter, a radially expanding bubble with a radial normal
would indeed give gamma_f = 1 everywhere and the blend would be inert.

On a Cartesian mesh the connection vector is axis-aligned, so for a face whose
normal direction is i,

    gamma_f = n_i^2,      i in {x, y, z}

and the three values available at any interface cell sum to one. gamma_f can
therefore only be near one where the interface normal happens to align with a
mesh axis, and is 1/3 on the body diagonals where n_x^2 = n_y^2 = n_z^2. The
sphere-averaged value is <n_i^2> = 1/3 exactly.

That is the whole mechanism behind the four-fold anisotropy: the blend samples
the interface orientation against six fixed directions, so its imprint carries
the cubic symmetry of the mesh however accurate the normal is.

The gradient is recomputed here from the saved `Vof Sharp` field by central
differences, which is what Predict_Beta differentiates; the saved
`Vof SurfaceNormals` is the smoothed normal used for curvature and is a
different field.

Run directory is spliced -- multiples of 100 are the coherent run B (see the
project notes) -- so only those timesteps are read.
"""

import glob
import os
import re
import sys

import numpy as np
import pyvista as pv

CASE = "/home/user/archive/Scriven-Merlin/Scriven-Struct-125-smallDt"
N_SNAP = 6               # snapshots sampled across the record
A_LO, A_HI = 0.01, 0.99  # interface cells


def load_grid(path):
    """Volume fraction on the structured index grid, plus the spacing."""
    m = pv.read(path)
    c = np.asarray(m.cell_centers().points, float)
    a = np.asarray(m.cell_data["Vof Sharp [1]"], float)
    ux, iy = np.unique(np.round(c[:, 0], 10), return_inverse=True)
    uy, jy = np.unique(np.round(c[:, 1], 10), return_inverse=True)
    uz, kz = np.unique(np.round(c[:, 2], 10), return_inverse=True)
    n = (len(ux), len(uy), len(uz))
    g = np.full(n, np.nan)
    g[iy, jy, kz] = a
    return g, float(ux[1] - ux[0])


def gammas(g, h):
    """gamma_f = n_i^2 at interface cells, for each of the three face directions."""
    gx, gy, gz = (np.gradient(g, h, axis=0),
                  np.gradient(g, h, axis=1),
                  np.gradient(g, h, axis=2))
    mag = np.sqrt(gx**2 + gy**2 + gz**2)
    sel = (g > A_LO) & (g < A_HI) & (mag > 0)
    nx, ny, nz = gx[sel] / mag[sel], gy[sel] / mag[sel], gz[sel] / mag[sel]
    return nx**2, ny**2, nz**2


def main():
    files = sorted((f for f in glob.glob(os.path.join(CASE, "bubble-ts*.pvtu"))
                    if "bnd" not in f and "front" not in f),
                   key=lambda f: int(re.search(r"ts0*(\d+)", f).group(1)))
    files = [f for f in files
             if int(re.search(r"ts0*(\d+)", f).group(1)) % 100 == 0]
    if not files:
        sys.exit(f"no coherent snapshots in {CASE}")
    files = files[:: max(1, len(files) // N_SNAP)]

    print("gamma_f = cos^2(angle between interface normal and the cell-connection")
    print("vector). On a Cartesian mesh that is n_i^2 for a face along axis i.\n")
    print(f"{'ts':>6} {'n_iface':>9} {'mean':>7} {'median':>8} "
          f"{'>0.9':>7} {'>0.5':>7} {'<0.1':>7} {'sum':>6}")
    for f in files:
        ts = int(re.search(r"ts0*(\d+)", f).group(1))
        g, h = load_grid(f)
        gx, gy, gz = gammas(g, h)
        allg = np.concatenate([gx, gy, gz])
        print(f"{ts:>6} {len(gx):>9} {allg.mean():>7.4f} "
              f"{np.median(allg):>8.4f} {100*(allg>0.9).mean():>6.1f}% "
              f"{100*(allg>0.5).mean():>6.1f}% {100*(allg<0.1).mean():>6.1f}% "
              f"{(gx+gy+gz).mean():>6.3f}")

    print("\n  'sum' is the mean of n_x^2+n_y^2+n_z^2 over interface cells; it is 1 by")
    print("  construction and is printed only as a check on the gradient.")
    print("  A sphere gives <n_i^2> = 1/3 exactly, so a mean near 0.333 means the")
    print("  blend is nowhere near unity: it is compressive only on the faces whose")
    print("  axis the interface normal happens to align with.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
