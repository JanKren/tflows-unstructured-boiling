#!/usr/bin/env python3
"""Is the reconstructed interfacial AREA orientation-biased on a Cartesian mesh?

The question behind this is whether a *uniform* interfacial mass flux would move
a sphere uniformly.  A cell receives mass in proportion to the interfacial area
the solver reconstructs in it, and the interface then advances by however far
that mass carries it.  If the reconstructed area is the area actually swept, the
two cancel and v_n = mdot"/rho everywhere.  If the reconstruction is biased by
an amount depending on orientation, the cancellation fails and the interface
speed picks up a cubic modulation -- with no advection scheme involved.

Measuring this on saved run output does not work: the bubble is already deformed
by 6.8% at the earliest snapshot, and area per steradian scales with the local
r^2, so the deformation swamps the quantity of interest.  There is no
front output at t = 0.  So the test is built here from scratch, on an EXACT
sphere, where there is nothing to correct for.

Method, per near-interface cell of a uniform Cartesian lattice:

  A_true  the exact area of the sphere inside the cell.  Obtained by sampling
          points uniformly on the sphere and binning them by cell, so each cell
          gets 4 pi R^2 * (its share of the samples).

  A_plic  the area the PLIC reconstruction assigns.  The interface is one plane
          per cell, oriented along the exact radial normal and positioned to
          reproduce the cell's exact volume fraction.  For a plane cutting a
          cube, the cut area is exactly dV/dd -- moving the plane by dd sweeps
          A dd -- so both the volume and the area follow in closed form from the
          same standard expression, with no polygon clipping needed.

The ratio A_plic/A_true is then regressed on the cubic invariant
s = nx^4+ny^4+nz^4 (1 on a coordinate axis, 1/3 on a body diagonal).  A flat
ratio means a uniform flux moves a sphere uniformly and the area is exonerated;
a slope means it does not.

Usage:  python3 analyze_plic_area_bias.py [--n 74] [--R 110e-6]
"""

import argparse

import numpy as np


def cube_volume_below(d, n):
    """Volume of the unit cube below the plane n.x = d, n unit and n_i > 0."""
    n1, n2, n3 = n
    s = lambda x: np.maximum(x, 0.0) ** 3
    num = (s(d) - s(d - n1) - s(d - n2) - s(d - n3)
           + s(d - n1 - n2) + s(d - n1 - n3) + s(d - n2 - n3)
           - s(d - n1 - n2 - n3))
    return num / (6.0 * n1 * n2 * n3)


def cube_cut_area(d, n):
    """Area of the plane n.x = d inside the unit cube.  Equals dV/dd exactly,
    because translating the plane by dd sweeps exactly (area) x dd."""
    n1, n2, n3 = n
    s = lambda x: np.maximum(x, 0.0) ** 2
    num = (s(d) - s(d - n1) - s(d - n2) - s(d - n3)
           + s(d - n1 - n2) + s(d - n1 - n3) + s(d - n2 - n3)
           - s(d - n1 - n2 - n3))
    return num / (2.0 * n1 * n2 * n3)


def plane_offset(alpha, n, iters=60):
    """Offset d with cube_volume_below(d, n) = alpha.  Bisection: V is monotone
    in d, and a vectorised bisection is both simpler and safer here than a
    piecewise-analytic inversion, which has to branch on the ordering of n."""
    lo = np.zeros_like(alpha)
    hi = np.full_like(alpha, float(n[0] + n[1] + n[2]))
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        v = cube_volume_below(mid, n)
        too_small = v < alpha
        lo = np.where(too_small, mid, lo)
        hi = np.where(too_small, hi, mid)
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=74, help="cells per side")
    ap.add_argument("--box", type=float, default=300e-6)
    ap.add_argument("--R", type=float, default=110e-6)
    ap.add_argument("--samples", type=int, default=20_000_000)
    a = ap.parse_args()

    h = a.box / a.n
    R = a.R
    print(f"\n  Cartesian lattice {a.n}^3, h = {h*1e6:.3f} um, exact sphere "
          f"R = {R*1e6:.1f} um  (R/h = {R/h:.1f})")

    # ---- A_true: sample the sphere surface, bin by cell -------------------
    rng = np.random.default_rng(0)
    v = rng.normal(size=(a.samples, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    pts = v * R
    idx = np.floor(pts / h).astype(np.int64)          # cell index, origin-centred
    key = (idx[:, 0] + 1000) * 4_000_000 + (idx[:, 1] + 1000) * 2000 + (idx[:, 2] + 1000)
    uniq, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    A_true = 4.0 * np.pi * R ** 2 * cnt / a.samples

    # cell centre and the exact radial normal there
    ctr = np.zeros((len(uniq), 3))
    for k in range(3):
        ctr[:, k] = np.bincount(inv, weights=pts[:, k], minlength=len(uniq)) / cnt
    cell_lo = np.floor(ctr / h) * h
    cell_ctr = cell_lo + 0.5 * h
    rc = np.linalg.norm(cell_ctr, axis=1)
    nrm = cell_ctr / rc[:, None]

    # ---- alpha: exact volume fraction of the sphere in each of those cells --
    sub = 12
    g = (np.arange(sub) + 0.5) / sub
    gx, gy, gz = np.meshgrid(g, g, g, indexing="ij")
    off = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1) * h
    alpha = np.empty(len(uniq))
    CH = 20000
    for i0 in range(0, len(uniq), CH):
        i1 = min(i0 + CH, len(uniq))
        p = cell_lo[i0:i1, None, :] + off[None, :, :]
        alpha[i0:i1] = (np.linalg.norm(p, axis=2) < R).mean(axis=1)

    band = (alpha > 1e-6) & (alpha < 1 - 1e-6) & (A_true > 0)
    alpha, nrm, A_true, cell_ctr = alpha[band], nrm[band], A_true[band], cell_ctr[band]
    print(f"  interface cells: {band.sum():,}   (sphere sampled with "
          f"{a.samples/1e6:.0f}M points, {a.samples/band.sum():.0f} per cell)")

    # ---- A_plic, cell by cell (the normal differs per cell, so loop) -------
    A_plic = np.empty(len(alpha))
    na = np.abs(nrm)
    for i in range(len(alpha)):
        nn = np.sort(na[i])
        nn = np.maximum(nn, 1e-9)
        nn = nn / np.linalg.norm(nn)
        d = plane_offset(np.array([alpha[i]]), nn)
        A_plic[i] = cube_cut_area(d, nn)[0] * h * h

    ratio = A_plic / A_true
    s = (nrm ** 4).sum(-1)

    # area-weighted regression of the ratio on the cubic invariant
    W = np.sqrt(A_true)
    M = np.vstack([W, W * s]).T
    c0, c1 = np.linalg.lstsq(M, W * ratio, rcond=None)[0]
    axis_v, face_v, body_v = c0 + c1, c0 + c1 / 2, c0 + c1 / 3

    print(f"\n  total  A_plic / A_true      = {A_plic.sum()/A_true.sum():.5f}")
    print(f"  ratio at s=1   (axis)       = {axis_v:.5f}")
    print(f"  ratio at s=1/2 (face diag)  = {face_v:.5f}")
    print(f"  ratio at s=1/3 (body diag)  = {body_v:.5f}")
    print(f"  body - axis                 = {100*(body_v/axis_v-1):+.3f} %")
    print(f"  scatter (std of ratio)      = {100*ratio.std():.2f} %")

    # binned, to show the shape rather than trust one linear fit
    print(f"\n  {'s range':>14}{'cells':>9}{'<A_plic/A_true>':>18}")
    edges = np.linspace(1 / 3, 1, 9)
    for i in range(len(edges) - 1):
        m = (s >= edges[i]) & (s < edges[i + 1])
        if m.sum() > 20:
            w = A_true[m]
            print(f"  {edges[i]:.3f}-{edges[i+1]:.3f}{m.sum():9d}"
                  f"{(ratio[m]*w).sum()/w.sum():18.5f}")
    print("\n  A flat column means a uniform flux moves a sphere uniformly:")
    print("  the area cancels out of v_n and cannot produce a four-fold mode.\n")


if __name__ == "__main__":
    main()
