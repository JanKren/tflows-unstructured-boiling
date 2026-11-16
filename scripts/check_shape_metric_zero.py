#!/usr/bin/env python3
"""Zero check for the shape metric: does it read anisotropy off a perfect sphere?

The metric of Section 4 regresses the front radius on the cubic invariant
s = nx^4+ny^4+nz^4. If the extraction of a spherical interface on a Cartesian
mesh biased that regression, every shape percentage in the paper would inherit
the bias, so the metric has to be shown to return zero on a sphere it did not
help create.

The check cannot be done from solver output: no front is written before the
first step, and by the first structured output at t = 0.2 ms the interface has
already deformed (dR/ds = -0.0975). It is therefore reproduced here from an
analytically initialised sphere, extracted the same way -- an alpha = 0.5
isosurface on a uniform Cartesian mesh -- and measured with the same
area-weighted regression the paper uses.

This reproduces the extraction step, not T-Flows' isoap polygon: the isosurface
is triangulated here rather than assembled as a single polygon per cell. What
it tests is whether an alpha = 0.5 surface through a Cartesian mesh carries a
four-fold radial signature, which is the question the check exists to answer.

Run at the four Scriven resolutions with the production radius, R0 = 50 um in a
300 um cube.
"""

import sys

import numpy as np
import pyvista as pv

L = 300.0e-6        # cube side [m]
R0 = 50.0e-6        # initial bubble radius [m]
SUB = 6             # sub-samples per direction when computing cell volume fractions


def volume_fraction(n_cells, sub=SUB):
    """Exact-ish liquid volume fraction per cell for a sphere of radius R0.

    alpha = 1 outside the bubble, 0 inside, matching the paper's convention that
    alpha is the liquid fraction. Sub-sampling each cell on a sub^3 lattice keeps
    the initialisation free of any directional preference of its own.
    """
    h = L / n_cells
    c1 = (np.arange(n_cells) + 0.5) * h - L / 2          # cell-centre coordinates
    off = (np.arange(sub) + 0.5) * h / sub - h / 2       # sub-sample offsets
    ox, oy, oz = np.meshgrid(off, off, off, indexing="ij")
    o2 = (ox**2 + oy**2 + oz**2).ravel()
    ox, oy, oz = ox.ravel(), oy.ravel(), oz.ravel()

    X, Y, Z = np.meshgrid(c1, c1, c1, indexing="ij")
    alpha = np.empty(X.shape, dtype=np.float64)
    flat = [a.ravel() for a in (X, Y, Z)]
    n = flat[0].size
    out = np.empty(n)
    chunk = 200_000
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        d2 = ((flat[0][a:b, None] + ox) ** 2
              + (flat[1][a:b, None] + oy) ** 2
              + (flat[2][a:b, None] + oz) ** 2)
        out[a:b] = (d2 > R0 * R0).mean(axis=1)           # liquid fraction
    alpha[:] = out.reshape(X.shape)
    return alpha, h


def measure(n_cells):
    alpha, h = volume_fraction(n_cells)
    grid = pv.ImageData(dimensions=(n_cells,) * 3, spacing=(h,) * 3,
                        origin=(-L / 2 + h / 2,) * 3)
    grid.point_data["alpha"] = alpha.ravel(order="F")
    surf = grid.contour([0.5], scalars="alpha")
    if surf.n_cells == 0:
        return None
    sized = surf.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(sized.cell_data["Area"], float)
    c = np.asarray(sized.cell_centers().points, float)
    c = c - np.average(c, axis=0, weights=a)             # centre on the surface itself
    r = np.linalg.norm(c, axis=1)
    s = ((c / r[:, None]) ** 4).sum(axis=1)
    W = a / a.sum()
    sm, rm = (W * s).sum(), (W * r).sum()
    b1 = (W * (s - sm) * (r - rm)).sum() / (W * (s - sm) ** 2).sum()
    return rm, b1 / rm, surf.n_cells


def main():
    print(f"{'mesh':>10} {'h [um]':>8} {'elements':>10} {'R [um]':>9} {'dR/ds':>10}")
    print("-" * 52)
    worst = 0.0
    for n in (74, 99, 124, 149):
        out = measure(n)
        if out is None:
            print(f"{n:>10} extraction failed")
            continue
        rm, d, ne = out
        worst = max(worst, abs(d))
        print(f"{str(n)+'^3':>10} {L/n*1e6:>8.2f} {ne:>10} {rm*1e6:>9.3f} {d:>+10.5f}")
    print(f"\n  largest |dR/ds| on a perfect sphere: {worst:.5f}")
    print("  for comparison, the structured runs reach -0.10 to -0.18 (Table 12)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
