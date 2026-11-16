#!/usr/bin/env python3
"""
Does an extended gradient stencil reduce the anisotropy on the ACTUAL
simulation state, or only on an idealised one?

Figure 13(c) and stencil_parameter_scan.py both show that adding edge and
corner neighbours suppresses the stencil's own anisotropy on an analytic
field.  analyze_anisotropy_attribution.py shows the stencil contributes
only ~0.4% of the ~7% measured on the real state.  Those two facts are
compatible with either conclusion:

  - the extended stencil is irrelevant (it fixes a 0.4% effect), or
  - the extended stencil helps anyway, because a wider support also
    low-pass filters the state-induced anisotropy it is handed.

This script settles it by running the extended stencil on the real
simulation state, using the same bit-exact reimplementation of T-Flows'
gradient that was validated in the attribution study.

Stencil families (unique connection directions, each giving a +/- pair):
  face6                 3 directions -> 6 neighbours   (what T-Flows does)
  face6_edge12          9 directions -> 18 neighbours
  face6_edge12_corner8 13 directions -> 26 neighbours

Weighting w = |r|^-p, with p = 0 reproducing T-Flows' unweighted LSQ.

Usage:  python3 analyze_extended_stencil.py [case.pvtu]
"""

import sys

import numpy as np
from scipy.special import erf

from analyze_anisotropy_attribution import (
    DEFAULT_CASE, load_lattice, fit_bubble, fourier_lsq, fourier_lsq_phase,
)


# ======================================================================
#  Stencil families
# ======================================================================
def pair_directions(family):
    """Unique connection directions; o and -o are the same face pair."""
    face = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    edge = [(1, 1, 0), (1, -1, 0), (1, 0, 1),
            (1, 0, -1), (0, 1, 1), (0, 1, -1)]
    corner = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]
    out = list(face)
    if family in ("face6_edge12", "face6_edge12_corner8"):
        out += edge
    if family == "face6_edge12_corner8":
        out += corner
    return [np.array(o, dtype=float) for o in out]


def _slices(o):
    """Slice pair (lo, hi) such that hi is lo shifted by offset o."""
    lo, hi = [], []
    for k in range(3):
        if o[k] > 0:
            lo.append(slice(0, -1)); hi.append(slice(1, None))
        elif o[k] < 0:
            lo.append(slice(1, None)); hi.append(slice(0, -1))
        else:
            lo.append(slice(None)); hi.append(slice(None))
    return tuple(lo), tuple(hi)


# ======================================================================
#  Generalised interface-modified LSQ gradient
# ======================================================================
def gradient(d, fit, T, family="face6", p_weight=0.0,
             use_ideal_interface=False):
    """Interface-modified LSQ gradient with an arbitrary stencil family.

    With family='face6' and p_weight=0 this reduces exactly to
    Calculate_Grad_Matrix_With_Front + Grad_Component_No_Refresh_With_Front.
    """
    h, shape = d["h"], d["shape"]
    vof, r, R, T_sat = d["vof"], fit["r"], fit["R"], fit["T_sat"]

    G = [np.zeros(shape) for _ in range(6)]
    G0 = [np.zeros(shape) for _ in range(6)]
    b = [np.zeros(shape) for _ in range(3)]
    IJ = [(0, 0, 0), (1, 1, 1), (2, 2, 2), (3, 0, 1), (4, 0, 2), (5, 1, 2)]

    for o in pair_directions(family):
        lo, hi = _slices(o)
        conn = o * h
        w = np.linalg.norm(conn) ** (-p_weight) if p_weight else 1.0

        a1, a2 = vof[lo], vof[hi]
        T1, T2 = T[lo], T[hi]

        if use_ideal_interface:
            s1, s2 = r[lo] - R, r[hi] - R
            crosses = (s1 * s2) < 0.0
            den = s1 - s2
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(crosses & (np.abs(den) > 1e-30),
                                s1 / np.where(np.abs(den) > 1e-30, den, 1.0),
                                0.0)
        else:
            crosses = ((a1 - 0.5) * (a2 - 0.5)) < 0.0
            den = np.abs(a1 - a2)
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(crosses & (den > 1e-15),
                                np.abs(a1 - 0.5) / np.where(den > 1e-15,
                                                            den, 1.0), 0.0)
        frac = np.clip(np.nan_to_num(frac), 0.0, 1.0)

        f = frac[..., None]
        r1 = np.where(crosses[..., None], f * conn, conn)
        r2 = np.where(crosses[..., None], (1.0 - f) * conn, conn)
        d1 = np.where(crosses, T_sat - T1, T2 - T1)
        d2 = np.where(crosses, T2 - T_sat, T2 - T1)

        # per-face weights: shortened vectors get their own |r|^-p
        if p_weight:
            n1 = np.linalg.norm(r1, axis=-1)
            n2 = np.linalg.norm(r2, axis=-1)
            w1 = np.where(n1 > 1e-30, n1, 1.0) ** (-p_weight)
            w2 = np.where(n2 > 1e-30, n2, 1.0) ** (-p_weight)
            w1 = np.where(n1 > 1e-30, w1, 0.0)
            w2 = np.where(n2 > 1e-30, w2, 0.0)
        else:
            w1 = w2 = 1.0

        for k, i, j in IJ:
            G[k][lo] += w1 * r1[..., i] * r1[..., j]
            G[k][hi] += w2 * r2[..., i] * r2[..., j]
            G0[k][lo] += w * conn[i] * conn[j]
            G0[k][hi] += w * conn[i] * conn[j]
        for i in range(3):
            b[i][lo] += w1 * d1 * r1[..., i]
            b[i][hi] += w2 * d2 * r2[..., i]

    def jac_of(g):
        return (g[0] * g[1] * g[2] - g[0] * g[5] ** 2 - g[3] ** 2 * g[2]
                + 2.0 * g[3] * g[4] * g[5] - g[4] ** 2 * g[1])

    jac = jac_of(G)
    bad = jac <= 0.0
    n_fb = int(bad.sum())
    for k in range(6):
        G[k] = np.where(bad, G0[k], G[k])
    jac = jac_of(G)
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = [
            (G[1] * G[2] - G[5] ** 2) / jac,
            (G[0] * G[2] - G[4] ** 2) / jac,
            (G[0] * G[1] - G[3] ** 2) / jac,
            -(G[3] * G[2] - G[4] * G[5]) / jac,
            (G[3] * G[5] - G[4] * G[1]) / jac,
            -(G[0] * G[5] - G[3] * G[4]) / jac,
        ]
    MAP = [[0, 3, 4], [3, 1, 5], [4, 5, 2]]
    grad = np.stack(
        [sum(inv[MAP[i][j]] * b[j] for j in range(3)) for i in range(3)],
        axis=-1)
    return np.nan_to_num(grad), n_fb


def liquid_band(d, fit, use_ideal):
    is_vap = (fit["r"] < fit["R"]) if use_ideal else (d["vof"] < 0.5)
    nb = np.zeros(d["shape"], dtype=bool)
    for axis in range(3):
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis], hi[axis] = slice(0, -1), slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        nb[lo] |= is_vap[hi]
        nb[hi] |= is_vap[lo]
    return (~is_vap) & nb


def measure(d, fit, grad, mask, g_ex):
    q = (fit["P"] - fit["ctr"])[mask]
    rr = np.linalg.norm(q, axis=-1)
    nrm = q / rr[:, None]
    g_n = np.abs(np.sum(grad[mask] * nrm, axis=-1))
    amps, phas, a0 = fourier_lsq_phase(np.arctan2(q[:, 1], q[:, 0]), g_n)
    return amps, a0, a0 / g_ex, phas


def main(fname):
    print("=" * 80)
    print("  EXTENDED STENCIL, TESTED ON THE ACTUAL SIMULATION STATE")
    print("=" * 80)
    print(f"  case: {fname}\n")

    d = load_lattice(fname)
    fit = fit_bubble(d)
    g_ex = 2.0 * fit["dT"] / (np.sqrt(np.pi) * fit["delta_T"])
    print(f"  R = {fit['R']*1e6:.2f} um, h = {d['h']*1e6:.4f} um, "
          f"delta_T = {fit['delta_T']/d['h']:.2f} h")
    print(f"  exact interface gradient = {g_ex:,.0f} K/m\n")

    x = fit["r"] - fit["R"]
    T_ideal = np.where(x > 0.0,
                       fit["T_sat"] + fit["dT"] * erf(x / fit["delta_T"]),
                       fit["T_sat"])

    # ---- sanity: face6/p=0 must still reproduce the VTU -----------------
    band = liquid_band(d, fit, False)
    g0, _ = gradient(d, fit, d["T"], "face6", 0.0, False)
    a = np.linalg.norm(g0[band], axis=-1)
    bm = np.linalg.norm(d["grad_saved"][band], axis=-1)
    ok = bm > 0
    err = float(np.median(np.abs(a[ok] - bm[ok]) / bm[ok]))
    print(f"  sanity check, face6/p=0 vs saved VTU: "
          f"median rel. error {100*err:.4f}%  "
          f"{'OK' if err < 1e-3 else 'MISMATCH'}\n")
    if err >= 1e-3:
        sys.exit("generalised stencil does not reduce to T-Flows; aborting")

    fams = ["face6", "face6_edge12", "face6_edge12_corner8"]
    ps = [0.0, 0.5, 1.0]

    for ideal, title in [(False, "ACTUAL simulation state "
                                 "(actual T, actual alpha interface)"),
                         (True, "IDEALISED state "
                                "(erf field, ideal sphere interface)")]:
        mask = liquid_band(d, fit, ideal)
        T = T_ideal if ideal else d["T"]
        print("=" * 80)
        print(f"  {title}")
        print("=" * 80)
        print(f"  {'stencil':<24s} {'p':>4s} {'a4 [%]':>8s} {'phi4':>8s} "
              f"{'max':<5s} {'m2':>6s} {'m6':>6s} {'m8':>6s} {'/exact':>7s} "
              f"{'fb':>6s}")
        print("  " + "-" * 86)
        for fam in fams:
            for p in ps:
                g, nfb = gradient(d, fit, T, fam, p, ideal)
                amps, a0, ratio, phas = measure(d, fit, g, mask, g_ex)
                tag = " <- T-Flows" if (fam == "face6" and p == 0.0
                                        and not ideal) else ""
                print(f"  {fam:<24s} {p:>4.1f} {100*amps[4]:8.4f} "
                      f"{phas[4]:8.2f} "
                      f"{'diag' if 30<=phas[4]<=60 else 'axes' if (phas[4]<15 or phas[4]>75) else 'mix':<5s} "
                      f"{100*amps[2]:6.3f} {100*amps[6]:6.3f} "
                      f"{100*amps[8]:6.2f} {ratio:7.3f} "
                      f"{nfb:6d}{tag}")
        print()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CASE)
