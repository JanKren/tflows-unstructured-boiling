#!/usr/bin/env python3
"""
Figure: the substitution test as four polar panels (brief section 6).

The 2x2 attribution is currently a table.  As four polar plots on shared axes the
argument is immediate: case 4, fed self-consistent analytical inputs, is a circle;
case 1, the real state, carries a visible four-fold mode; and cases 2 and 3 -- one
input idealised each -- are worse than either, and are four-fold in ANTIPHASE with
one another, which is what makes their cancellation in case 1 visible rather than
merely asserted.

Panels show |grad T . n| on the liquid-side interface band, normalised by the mean
of case 1 so the four panels are directly comparable in radius, against azimuth.
Grey radial guides mark the coordinate axes; the fitted m=4 mode is overlaid.

Usage:  python3 plot_substitution_polar.py [case.pvtu] [--out PATH] [--bins 72]
"""

import argparse
import os

import numpy as np
from scipy.special import erf

from analyze_anisotropy_attribution import (
    DEFAULT_CASE, load_lattice, fit_bubble, modified_gradient,
    fourier_lsq_phase, M_MAX,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="?", default=DEFAULT_CASE)
    ap.add_argument("--out", default="results/figs/T0.2_substitution_polar.png")
    ap.add_argument("--bins", type=int, default=72)
    a = ap.parse_args()

    d = load_lattice(a.case)
    fit = fit_bubble(d)
    vof = d["vof"]
    x = fit["r"] - fit["R"]
    T_ideal = np.where(x > 0.0,
                       fit["T_sat"] + fit["dT"] * erf(x / fit["delta_T"]),
                       fit["T_sat"])

    cases = [
        ("(a) case 1", "simulation $T$, $\\alpha$ interface", d["T"], False),
        ("(b) case 2", "analytical $T$, $\\alpha$ interface", T_ideal, False),
        ("(c) case 3", "simulation $T$, exact sphere", d["T"], True),
        ("(d) case 4", "analytical $T$, exact sphere", T_ideal, True),
    ]

    edges = np.linspace(-np.pi, np.pi, a.bins + 1)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    out = []
    ref = None
    for tag, sub, Tf, ideal in cases:
        is_vap = (fit["r"] < fit["R"]) if ideal else (vof < 0.5)
        nb = np.zeros_like(is_vap)
        for ax in range(3):
            nb |= np.roll(is_vap, 1, ax) | np.roll(is_vap, -1, ax)
        mask = (~is_vap) & nb

        grad, _ = modified_gradient(d, fit, Tf, ideal)
        q = (fit["P"] - fit["ctr"])[mask]
        rr = np.linalg.norm(q, axis=-1)
        gn = np.abs(np.sum(grad[mask] * (q / rr[:, None]), axis=-1))
        th = np.arctan2(q[:, 1], q[:, 0])

        amps, phas, a0 = fourier_lsq_phase(th, gn, M_MAX)
        if ref is None:
            ref = a0
        idx = np.digitize(th, edges) - 1
        prof = np.array([gn[idx == b].mean() if (idx == b).any() else np.nan
                         for b in range(a.bins)]) / ref
        fit4 = (a0 / ref) * (1.0 + amps[4] * np.cos(4 * (ctr - np.radians(phas[4]))))
        out.append((tag, sub, prof, fit4, 100 * amps[4], phas[4]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'lines.linewidth': 1.3})

    rmax = max(np.nanmax(p) for _, _, p, _, _, _ in out) * 1.08
    fig, axes = plt.subplots(1, 4, figsize=(12.4, 3.5),
                             subplot_kw={'projection': 'polar'})
    for ax, (tag, sub, prof, fit4, a4, p4) in zip(axes, out):
        t = np.concatenate([ctr, ctr[:1]])
        ax.plot(t, np.concatenate([prof, prof[:1]]), color='0.35', lw=1.0)
        ax.plot(t, np.concatenate([fit4, fit4[:1]]), color='C3', lw=1.6)
        ax.set_ylim(0, rmax)
        ax.set_yticklabels([])
        ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
        ax.set_xticklabels(['0', '', '90', '', '180', '', '270', ''],
                           fontsize=7)
        ax.grid(alpha=0.3)
        ax.set_title(f"{tag}\n{sub}\n$a_4$={a4:.2f}%, $\\varphi_4$={p4:.0f}$^\\circ$",
                     fontsize=8, pad=12)
    fig.suptitle(r"Interfacial gradient $|\nabla T\cdot\mathbf{n}|$ by azimuth, "
                 r"normalised by the case-1 mean "
                 r"(grey: binned data, red: fitted $m=4$)", fontsize=9, y=1.06)
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    print(f"\n  figure -> {a.out}")
    for tag, sub, _, _, a4, p4 in out:
        print(f"    {tag} {sub:38s} a4 = {a4:7.4f}%   phi4 = {p4:6.2f} deg")


if __name__ == "__main__":
    main()
