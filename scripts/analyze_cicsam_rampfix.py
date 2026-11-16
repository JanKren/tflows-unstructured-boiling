#!/usr/bin/env python3
"""Can the four-fold mode be removed by fixing the sampling rate, not the normal?

    python3 analyze_cicsam_rampfix.py

Section 4 locates the anisotropy in the normalised donor value: the volume
fraction of a cubic cell crosses 0 to 1 over a normal distance W(n) = h||n||_1,
so three cells sampled at spacing h read that ramp at an orientation-dependent
rate. Writing that rate out,

    per-step fraction of the ramp   Delta = |n_i| / ||n||_1

which is 1 on an axis, 1/2 on a face diagonal, 1/3 on the body diagonal --- the
factor of three the section quotes.

The point of this script is that Delta is not an abstraction: it is exactly the
quantity `dotprod` in Predict_Beta. There

    dotprod = dalpha/dx_i * h        and, in the transition, |grad alpha| = 1/W,

so dotprod = (n_i / (h||n||_1)) * h = |n_i| / ||n||_1 = Delta, and the upwind
value is reconstructed as alfa_u = alfa_a - 2*dotprod. The orientation
dependence therefore enters through a single scalar, one line above the place
the limiter is usually blamed.

That suggests a fix the mechanism dictates rather than one guessed at: hold the
sampling rate constant. Variants tested here, all on the same exact spherical
volume fraction used for Table 8, so the numbers are directly comparable:

  baseline           as implemented                     (expect m=4 ~ 6.4%)
  gamma = 0.5        blend weight frozen                (expect m=4 ~ 5.7%)
  dotprod = const    sampling rate made orientation-free
  |grad| * h         magnitude instead of axis component
  Co * ||n||_1       ramp-normalised Courant in the CBC bound
  Co / ||n||_1       the same, opposite sense

A fix that works must drive m=4 towards zero while leaving m=2 and m=6 small;
cubic symmetry is the check that the evaluation has not been broken instead.
"""

import sys

import numpy as np

from analyze_anisotropy_attribution import (
    load_lattice, fit_bubble, fourier_lsq_phase, DEFAULT_CASE,
)
from analyze_extraction_error import exact_sphere_vof
from analyze_cicsam_beta import grad_central

CO = 0.2               # matches Table 8


def beta_variant(vof, grad_a, h, co, axis, ctr, X, Y, Z,
                 gamma_fix=None, dot_fix=None, dot_mode=None, co_scale=None):
    """Predict_Beta transcription with one knob turned at a time."""
    lo = [slice(None)] * 3;  hi = [slice(None)] * 3
    lo[axis] = slice(0, -1); hi[axis] = slice(1, None)
    lo, hi = tuple(lo), tuple(hi)

    fc = [0.5 * (c[lo] + c[hi]) for c in (X, Y, Z)]
    fr = np.sqrt(sum((fc[k] - ctr[k]) ** 2 for k in range(3)))
    v_flux = (fc[axis] - ctr[axis]) / np.where(fr > 0, fr, 1.0)
    pos = v_flux > 0.0
    pick = lambda a, b: np.where(pos, a, b)

    alfa_d = pick(vof[lo], vof[hi])
    alfa_a = pick(vof[hi], vof[lo])
    signo  = np.where(pos, 1.0, -1.0)

    g_d = np.stack([pick(grad_a[lo][..., k], grad_a[hi][..., k])
                    for k in range(3)], axis=-1)
    gmag = np.linalg.norm(g_d, axis=-1)

    dotprod = signo * g_d[..., axis] * h          # == Delta for a plane
    if dot_mode == "mag":                          # |grad| h  (drops the n_i)
        dotprod = signo * gmag * h
    if dot_fix is not None:                        # orientation-free sampling
        dotprod = signo * dot_fix
    prodmag = gmag * h

    alfa_u = np.clip(alfa_a - 2.0 * dotprod, 0.0, 1.0)
    den = alfa_a - alfa_u
    good = np.abs(den) > 1e-12
    with np.errstate(divide="ignore", invalid="ignore"):
        alfa_d_til = np.where(good, (alfa_d - alfa_u)
                              / np.where(good, den, 1.0), 0.0)
        l1 = np.abs(g_d).sum(-1) / np.where(gmag > 0, gmag, 1.0)   # ||n||_1
        cod = np.full_like(alfa_d_til, min(1.0, co))
        if co_scale == "mul":
            cod = np.minimum(1.0, co * l1)
        elif co_scale == "div":
            cod = np.minimum(1.0, co / np.where(l1 > 0, l1, 1.0))
        inrange = (alfa_d_til >= 0.0) & (alfa_d_til <= 1.0)
        alfa_cbc = np.where(inrange,
                            np.minimum(1.0, alfa_d_til / np.maximum(cod, 1e-30)),
                            alfa_d_til)
        alfa_uq = np.where(inrange,
                           np.minimum(cod * alfa_d_til + (1.0 - cod)
                                      * (6.0 * alfa_d_til + 3.0) / 8.0,
                                      alfa_cbc),
                           alfa_d_til)
        ang = np.where(prodmag > 1e-30,
                       signo * g_d[..., axis] * h
                       / np.where(prodmag > 1e-30, prodmag, 1.0), 0.0)
        gamma_f = (np.minimum(ang ** 2, 1.0) if gamma_fix is None
                   else np.full_like(ang, gamma_fix))
        alfa_f_til = gamma_f * alfa_cbc + (1.0 - gamma_f) * alfa_uq
        db = 1.0 - alfa_d_til
        beta = np.where(np.abs(db) > 1e-12,
                        (alfa_f_til - alfa_d_til)
                        / np.where(np.abs(db) > 1e-12, db, 1.0), 0.0)
    beta = np.clip(np.nan_to_num(beta), 0.0, 1.0)
    return beta, good & inrange & (gmag > 1e-30)


def evaluate(vof, grad_a, h, ctr, X, Y, Z, R, shape, **kw):
    """Per-cell aggregate of beta_f over all six faces, then fit m = 2,4,6.

    Transcribed from analyze_cicsam_beta.py.  beta_f on a single axis is not
    cubically symmetric by construction, so m2 and m6 blow up if the six faces
    are not aggregated onto the cell first -- which is the sanity check that
    caught an earlier wrong version of this script.
    """
    acc = np.zeros(shape); cnt = np.zeros(shape)
    for axis in range(3):
        b, valid = beta_variant(vof, grad_a, h, CO, axis, ctr, X, Y, Z, **kw)
        lo = [slice(None)] * 3; hi = [slice(None)] * 3
        lo[axis] = slice(0, -1); hi[axis] = slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        for tgt in (lo, hi):                 # both cells share the face
            acc[tgt] += np.where(valid, b, 0.0)
            cnt[tgt] += valid.astype(float)
    mean_b = np.where(cnt > 0, acc / np.where(cnt > 0, cnt, 1.0), np.nan)

    dx, dy = X - ctr[0], Y - ctr[1]
    rr = np.sqrt(dx**2 + dy**2 + (Z - ctr[2])**2)
    ok = (np.abs(rr - R) < 1.5 * h) & (cnt > 0)
    if ok.sum() < 100:
        return None
    amps, phas, a0 = fourier_lsq_phase(np.arctan2(dy, dx)[ok], mean_b[ok])
    return a0, amps, phas, int(ok.sum())


def main():
    d = load_lattice(DEFAULT_CASE)
    fit = fit_bubble(d)
    ctr = np.asarray(fit["ctr"], float)
    R, h = fit["R"], d["h"]
    X, Y, Z = d["X"], d["Y"], d["Z"]
    vof = exact_sphere_vof(d, ctr, R, 12, band=3.0)
    grad_a = grad_central(vof, h)

    print(f"exact sphere, R/h = {R/h:.1f}, Co = {CO}\n")
    print(f"  {'variant':<34}{'<beta>':>8}{'m=4 [%]':>10}{'phi_4':>8}"
          f"{'m=2':>7}{'m=6':>7}")
    cases = [
        ("baseline (as implemented)",        {}),
        ("gamma_f frozen at 0.5",            dict(gamma_fix=0.5)),
        ("dotprod fixed = 1.0",              dict(dot_fix=1.0)),
        ("dotprod fixed = 0.5",              dict(dot_fix=0.5)),
        ("dotprod fixed = 0.25",             dict(dot_fix=0.25)),
        ("dotprod = |grad a| h",             dict(dot_mode="mag")),
        ("Co * ||n||_1 in the CBC bound",    dict(co_scale="mul")),
        ("Co / ||n||_1 in the CBC bound",    dict(co_scale="div")),
        ("dotprod=0.5 AND gamma_f=0.5",      dict(dot_fix=0.5, gamma_fix=0.5)),
    ]
    for name, kw in cases:
        out = evaluate(vof, grad_a, h, ctr, X, Y, Z, R, d['shape'], **kw)
        if out is None:
            print(f"  {name:<34}{'no faces':>8}"); continue
        a0, amps, phas, n = out
        m4, m2, m6 = 100.0*amps[4], 100.0*amps[2], 100.0*amps[6]
        print(f"  {name:<34}{a0:>8.4f}{m4:>10.2f}{phas[4]%90:>8.2f}"
              f"{m2:>7.2f}{m6:>7.2f}")
    print("\n  Delta = |n_i|/||n||_1 is the per-step fraction of the ramp, and is")
    print("  exactly `dotprod` in Predict_Beta: 1 on an axis, 1/2 on a face")
    print("  diagonal, 1/3 on the body diagonal.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
