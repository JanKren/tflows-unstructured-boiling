#!/usr/bin/env python3
"""The kinematic radial-expansion test (Reviewer 2, comment 2.8).

The volume fraction of a sphere is advected by the production CICSAM in a
prescribed radial flow on the Scriven meshes: no energy equation, no mass
transfer, no momentum or pressure solve.  What is prescribed is the face flux,
the exact flux of a point source at the bubble centre, in four velocity models
(results/caf-kinematic-advection.md):

  solenoidal     the liquid velocity Q/(4 pi r^2) continued through the
                 interface band;
  jump           that velocity outside the computed bubble and none inside it:
                 the flux of a face is kept where the mean volume fraction of
                 its two cells exceeds 1/2;
  jump, exact    the same jump placed on the exact sphere R_S(t) instead of on
                 the computed interface;
  ramp           the one-fluid velocity vof * u_l, with the face value of vof
                 the mean of the two cells.

Panel (a): body - axis against radius on the structured 75^3 mesh, with the
coupled run at a uniform interfacial flux (Data/unif-transport) for
comparison.  Panel (b): the jump model on every mesh computed.

The numbers come from analyze_kinematic_expansion.py, which writes
results/caf-kinematic-expansion.csv with the manuscript's own shape functions.

Usage:  python3 plot_kinematic_expansion.py [--out pics/kinematic_expansion.pdf]
"""

import argparse
import csv
import os

import numpy as np

from plot_noadvection import series, UNIF_ON

CSV = "results/caf-kinematic-expansion.csv"


def load():
    runs = {}
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            runs.setdefault(r["run"], []).append(
                {k: (float(v) if k not in ("run", "family", "note") else v)
                 for k, v in r.items()})
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/kinematic_expansion.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True,
                         'grid.alpha': 0.3})

    runs = load()
    coupled = series(UNIF_ON)

    fig = plt.figure(figsize=(7.4, 3.3))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    # ---- (a) the three velocity models on the structured 75^3 -------------
    def xy(name):
        d = runs.get(name, [])
        return [r["R"] * 1e6 for r in d], [r["dev"] for r in d]

    ax1.plot([d["R"] * 1e6 for d in coupled], [d["dev"] for d in coupled],
             '^--', color='0.65', ms=4, lw=1.1,
             label=r"coupled solver, uniform $\dot m$")
    for name, col, mark, fill, lab in (
            ("st075-dt2-step", 'tab:grey', 'o', 'tab:grey',
             "velocity jump at the computed interface"),
            ("st075-dt2-mask", 'tab:grey', 'v', 'none',
             "velocity jump at the exact sphere"),
            ("st075-dt2-mix", 'tab:grey', 'D', 'none',
             "velocity ramp across the band"),
            ("st075-dt2-sol", 'tab:red', 's', 'tab:red',
             "solenoidal (liquid velocity throughout)")):
        x, y = xy(name)
        if x:
            ax1.plot(x, y, '-', marker=mark, color=col, mfc=fill, ms=3.5,
                     lw=1.3, label=lab)
    ax1.axhline(0.0, color='k', lw=0.6, ls=':')
    ax1.set_xlabel(r"interface radius $R$ [$\mu$m]")
    ax1.set_ylabel(r"body diagonal $-$ axis [%]")
    ax1.legend(fontsize=6.4, loc='upper left', framealpha=0.92,
               handlelength=1.8, borderpad=0.4)
    lo, hi = ax1.get_ylim()
    ax1.set_ylim(lo, hi + 0.55 * (hi - lo))
    ax1.set_title(r"(a) structured $75^3$, kinematic unless stated",
                  fontsize=9, pad=6)

    # ---- (b) the jump model on every mesh --------------------------------
    col_s = {75: '0.75', 100: '0.55', 125: '0.35', 150: '0.0'}
    for lev in (75, 100, 125, 150):
        d = runs.get(f"st{lev:03d}-dt2-step", [])
        if d:
            ax2.plot([r["R"] * 1e6 for r in d], [r["dev"] for r in d], 'o-',
                     color=col_s[lev], ms=3, lw=1.2,
                     label=rf"structured ${lev}^3$")
    col_p = {75: '#f4a582', 100: '#d6604d', 125: '#b2182b', 150: '#67001f'}
    for lev in (75, 100, 125, 150):
        d = runs.get(f"po{lev:03d}-dt2-step", [])
        if d:
            ax2.plot([r["R"] * 1e6 for r in d], [r["dev"] for r in d], 's-',
                     color=col_p[lev], ms=3, lw=1.2,
                     label=rf"polyhedral ${lev}^3$")
    ax2.axhline(0.0, color='k', lw=0.6, ls=':')
    ax2.set_xlabel(r"interface radius $R$ [$\mu$m]")
    ax2.set_ylabel(r"body diagonal $-$ axis [%]")
    ax2.legend(fontsize=6.4, loc='upper left', framealpha=0.92, ncol=2,
               handlelength=1.6, borderpad=0.4, columnspacing=1.0)
    lo, hi = ax2.get_ylim()
    ax2.set_ylim(lo, hi + 0.30 * (hi - lo))
    ax2.set_title(r"(b) velocity jump at the computed interface",
                  fontsize=9, pad=6)

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    fig.savefig(a.out.replace('.pdf', '.png'), bbox_inches='tight')
    print(f"wrote {a.out} and {a.out.replace('.pdf', '.png')}")


if __name__ == "__main__":
    main()
