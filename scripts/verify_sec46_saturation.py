#!/usr/bin/env python3
"""Fraction of interfacial faces on which Hyper-C is saturated (alpha_tilde_D > Co),
rebuilt on the exact-sphere lattice of verify_sec46_beta_aggregation.py, following
plot_ramp_sampling.py::saturated_fraction_report (all faces with a gradient and
0 <= alpha_tilde_D <= 1)."""
import numpy as np, sys
sys.argv = [sys.argv[0], "0.2"]
exec(open("verify_sec46_beta_aggregation.py").read().split("def fourier_fit")[0])   # reuse lattice + beta_faces
vals = []
for axis in range(3):
    lo = [slice(None)] * 3; hi = [slice(None)] * 3
    lo[axis] = slice(0, -1); hi[axis] = slice(1, None); lo, hi = tuple(lo), tuple(hi)
    fc = [0.5 * (c[lo] + c[hi]) for c in (X, Y, Z)]
    fr = np.sqrt(sum(f**2 for f in fc)); pos = fc[axis] / np.where(fr > 0, fr, 1) > 0
    pick = lambda a, b: np.where(pos, a, b)
    ad = pick(alpha[lo], alpha[hi]); aa = pick(alpha[hi], alpha[lo]); sg = np.where(pos, 1.0, -1.0)
    gd = np.stack([pick(g[lo][..., k], g[hi][..., k]) for k in range(3)], -1)
    dot = sg * gd[..., axis] * h
    au = np.clip(aa - 2 * dot, 0, 1); den = aa - au; ok = np.abs(den) > 1e-12
    adt = np.where(ok, (ad - au) / np.where(ok, den, 1), 0.0)
    sel = ok & (adt >= 0) & (adt <= 1) & (np.linalg.norm(gd, axis=-1) > 1e-30)
    vals.append(adt[sel])
v = np.concatenate(vals)
print(f"interfacial faces: {v.size:,}   median alpha_tilde_D = {np.median(v):.3f}")
for co in (0.014, 0.020, 0.023, 0.026, 0.032, 0.047, 0.2, 0.3, 0.44, 0.6):
    print(f"  Co = {co:5.3f}   saturated fraction (alpha_tilde_D > Co) = {100 * (v > co).mean():5.1f}%")
