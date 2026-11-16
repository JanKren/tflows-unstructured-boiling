#!/usr/bin/env python3
"""Reproduce the Table `tab:cicsam_beta` evaluation on an EXACT sphere and test
how the four-fold content of beta_f depends on how faces are aggregated.

Lattice: cube of side L = 300 um, N^3 cells (N = 124 for the '125^3' mesh), so
h = 2.42 um.  Sphere radius R = 110 um (R/h = 45).  alpha = liquid fraction
(0 inside, 1 outside), exact by sub-sampling in the interface band.  Velocity =
prescribed radial expansion, so the donor of every face is the inner cell.
beta_f transcribed from analyze_cicsam_beta.py (itself from Predict_Beta.f90).

Aggregates compared:
  (A) paper's: per-cell unweighted mean of beta over 'valid' faces, cells in
      the band |r - R| < 1.5 h, LSQ Fourier fit in azimuth over all such cells
  (B) same but weighted by the face's share of the cell's outward flux |n_i|
  (C) restricted to front cells 0.3 < alpha < 0.7, flux-weighted
  (D) front cells, unweighted
"""
import numpy as np, time, sys

N = 124; L = 300e-6; h = L / N; R = 110e-6; CO = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
SUB = 10

t0 = time.time()
x1 = (np.arange(N) + 0.5) * h - L / 2
X, Y, Z = np.meshgrid(x1, x1, x1, indexing="ij")
rr = np.sqrt(X**2 + Y**2 + Z**2)
alpha = (rr > R).astype(float)                      # liquid outside
near = np.abs(rr - R) < 1.0 * h                     # cells that may be cut
idx = np.argwhere(near)
print(f"cells: {N**3:,}, cut candidates: {len(idx):,}")
# sub-sample
off = (np.arange(SUB) + 0.5) / SUB - 0.5
ox, oy, oz = np.meshgrid(off, off, off, indexing="ij")
ox, oy, oz = ox.ravel() * h, oy.ravel() * h, oz.ravel() * h
chunk = 4000
for c0 in range(0, len(idx), chunk):
    ii = idx[c0:c0 + chunk]
    cx = X[ii[:, 0], ii[:, 1], ii[:, 2]][:, None] + ox[None, :]
    cy = Y[ii[:, 0], ii[:, 1], ii[:, 2]][:, None] + oy[None, :]
    cz = Z[ii[:, 0], ii[:, 1], ii[:, 2]][:, None] + oz[None, :]
    frac = (np.sqrt(cx**2 + cy**2 + cz**2) > R).mean(axis=1)
    alpha[ii[:, 0], ii[:, 1], ii[:, 2]] = frac
print(f"alpha built in {time.time() - t0:.1f}s; vapour volume error "
      f"{((1 - alpha).sum() * h**3) / (4 / 3 * np.pi * R**3) - 1:+.2e}")

def grad_central(f):
    return np.stack([np.gradient(f, h, axis=k) for k in range(3)], -1)

g = grad_central(alpha)

def beta_faces(axis, gamma_fix=None):
    lo = [slice(None)] * 3; hi = [slice(None)] * 3
    lo[axis] = slice(0, -1); hi[axis] = slice(1, None); lo, hi = tuple(lo), tuple(hi)
    fc = [0.5 * (c[lo] + c[hi]) for c in (X, Y, Z)]
    fr = np.sqrt(sum(f**2 for f in fc))
    vflux = fc[axis] / np.where(fr > 0, fr, 1.0)     # u.n_face for radial expansion
    pos = vflux > 0
    pick = lambda a, b: np.where(pos, a, b)
    ad = pick(alpha[lo], alpha[hi]); aa = pick(alpha[hi], alpha[lo])
    signo = np.where(pos, 1.0, -1.0)
    gd = np.stack([pick(g[lo][..., k], g[hi][..., k]) for k in range(3)], -1)
    gmag = np.linalg.norm(gd, axis=-1)
    dot = signo * gd[..., axis] * h
    au = np.clip(aa - 2.0 * dot, 0.0, 1.0)
    den = aa - au
    good = np.abs(den) > 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        adt = np.where(good, (ad - au) / np.where(good, den, 1.0), 0.0)
        inr = (adt >= 0) & (adt <= 1)
        cbc = np.where(inr, np.minimum(1.0, adt / CO), adt)
        uq = np.where(inr, np.minimum(CO * adt + (1 - CO) * (6 * adt + 3) / 8, cbc), adt)
        ang = np.where(gmag * h > 1e-30, dot / np.where(gmag * h > 1e-30, gmag * h, 1.0), 0.0)
        gam = np.minimum(ang**2, 1.0) if gamma_fix is None else np.full_like(ang, gamma_fix)
        aft = gam * cbc + (1 - gam) * uq
        db = 1 - adt
        beta = np.where(np.abs(db) > 1e-12, (aft - adt) / np.where(np.abs(db) > 1e-12, db, 1.0), 0.0)
    beta = np.clip(np.nan_to_num(beta), 0, 1)
    valid = good & inr & (gmag > 1e-30)
    w = np.abs(vflux)                                 # |n_i| of the face: flux share
    return beta, valid, w, lo, hi

def fourier_fit(theta, y, M=8):
    A = [np.ones_like(theta)]
    for m in range(1, M + 1):
        A += [np.cos(m * theta), np.sin(m * theta)]
    A = np.stack(A, 1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    a0 = c[0]; out = {}
    for m in range(1, M + 1):
        a, b = c[2 * m - 1], c[2 * m]
        out[m] = (np.hypot(a, b) / a0, (np.rad2deg(np.arctan2(b, a) / m)) % (360 / m))
    return a0, out

theta = np.arctan2(Y, X)
band = np.abs(rr - R) < 1.5 * h
front = (alpha > 0.3) & (alpha < 0.7)
front_liq = (alpha > 0.5) & (alpha < 1.0)          # a liquid-side band, crude

for label, gfix in (("blend cos^2 (as implemented)", None), ("gamma=1 (Hyper-C)", 1.0),
                    ("gamma=0.5", 0.5), ("gamma=0 (UQ)", 0.0)):
    acc = np.zeros(alpha.shape); cnt = np.zeros(alpha.shape)
    accw = np.zeros(alpha.shape); cntw = np.zeros(alpha.shape)
    for axis in range(3):
        b, v, w, lo, hi = beta_faces(axis, gfix)
        for tgt in (lo, hi):
            acc[tgt] += np.where(v, b, 0.0); cnt[tgt] += v
            accw[tgt] += np.where(v, w * b, 0.0); cntw[tgt] += np.where(v, w, 0.0)
    meanb = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    meanw = np.where(cntw > 0, accw / np.where(cntw > 0, cntw, 1), np.nan)
    print(f"\n### {label}   (Co = {CO})")
    for agg, mask, arr in (("(A) band, unweighted  [paper]", band & (cnt > 0), meanb),
                           ("(B) band, flux-weighted", band & (cntw > 0), meanw),
                           ("(C) front 0.3<a<0.7, flux-weighted", front & (cntw > 0), meanw),
                           ("(D) front 0.3<a<0.7, unweighted", front & (cnt > 0), meanb)):
        a0, fo = fourier_fit(theta[mask], arr[mask])
        print(f"  {agg:38s} n={mask.sum():6d}  <beta>={a0:.4f}  "
              f"m4={100 * fo[4][0]:6.2f}% @ {fo[4][1]:5.1f}   m2={100 * fo[2][0]:5.2f}%  "
              f"m6={100 * fo[6][0]:5.2f}%  m8={100 * fo[8][0]:5.2f}%")

# ---------------------------------------------------------------- belt check --
print("\n=== equatorial belt |lat| < 15 deg, binned azimuthal profile of the per-cell beta ===")
lat = np.arcsin(np.clip(Z / np.maximum(rr, 1e-30), -1, 1))
belt = np.abs(lat) < np.deg2rad(15)
edges = np.deg2rad(np.arange(0, 91, 7.5))
def profile(mask, arr):
    th = np.mod(theta[mask], np.pi / 2)              # fold onto [0, 90) by cubic symmetry
    v = arr[mask]
    return [np.nanmean(v[(th >= edges[k]) & (th < edges[k + 1])]) for k in range(len(edges) - 1)]
for label, gfix in (("blend cos^2", None), ("gamma=1 (Hyper-C)", 1.0), ("gamma=0 (UQ)", 0.0)):
    acc = np.zeros(alpha.shape); cnt = np.zeros(alpha.shape)
    accw = np.zeros(alpha.shape); cntw = np.zeros(alpha.shape)
    for axis in range(3):
        b, v, w, lo, hi = beta_faces(axis, gfix)
        for tgt in (lo, hi):
            acc[tgt] += np.where(v, b, 0.0); cnt[tgt] += v
            accw[tgt] += np.where(v, w * b, 0.0); cntw[tgt] += np.where(v, w, 0.0)
    meanb = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
    meanw = np.where(cntw > 0, accw / np.where(cntw > 0, cntw, 1), np.nan)
    print(f"\n{label}: azimuth bins 0-7.5, ..., 82.5-90 deg (folded); 0 = axis, 45 = face diagonal")
    for agg, mask, arr in (("band, unweighted [paper]", band & belt & (cnt > 0), meanb),
                           ("band, flux-weighted", band & belt & (cntw > 0), meanw),
                           ("front, unweighted", front & belt & (cnt > 0), meanb),
                           ("front, flux-weighted", front & belt & (cntw > 0), meanw)):
        p = profile(mask, arr)
        print(f"  {agg:26s} " + " ".join(f"{x:5.3f}" for x in p))
