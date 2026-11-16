#!/usr/bin/env python3
"""Standard advection tests of the production CICSAM against the literature.

Panel (a): the disc of Roenby, Bredmose and Jasak (2016), Section 4.1, after
translation by sixteen radii in U = (1, 0.5) on 40 cells per unit length: the
alpha = 1/2 contour at three Courant numbers against the exact circle.
Panel (b): the shape error E1 of that test on their five (Nx, Co) rows, the
present solver beside the four schemes of their Table 1.

Data: /home/user/runs/caf-review-2026-11/kinematic/adv (analyze_advection_tests.py
prints the full table, with the other norms and the three-dimensional tests).

Usage:  python3 plot_advection_tests.py [--out pics/advection_tests.pdf]
"""

import argparse
import glob
import os

import numpy as np
import pyvista as pv

from analyze_advection_tests import ROOT, ROENBY_T1, read

ROWS = [(10, 0.5, "disc10-co05"), (20, 0.5, "disc20-co05"),
        (40, 0.5, "disc40-co05"), (40, 0.2, "disc40-co02"),
        (40, 0.1, "disc40-co01")]


def field(run):
    """Cell centres and volume fraction of the last saved snapshot."""
    # a serial run writes one .vtu, a parallel one a .pvtu and its pieces
    f = sorted(glob.glob(os.path.join(ROOT, run, "*-ts*.pvtu"))) or \
        sorted(glob.glob(os.path.join(ROOT, run, "*-ts*.vtu")))
    m = pv.read(f[-1])
    key = [n for n in m.cell_data.keys() if "vof" in n.lower()][0]
    c = np.asarray(m.cell_centers().points, float)
    a = np.asarray(m.cell_data[key], float)
    return c[:, 0], c[:, 1], a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/advection_tests.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True,
                         'grid.alpha': 0.3})

    fig = plt.figure(figsize=(7.4, 3.3))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    # ---- (a) final contours on Nx = 40 ------------------------------------
    th = np.linspace(0, 2 * np.pi, 400)
    ax1.plot(4.5 + 0.25 * np.cos(th), 2.5 + 0.25 * np.sin(th), '-',
             color='k', lw=0.8, label="exact")
    for run, col, ls, lab in (("disc40-co05", '0.65', '--', r"$\mathrm{Co}=0.5$"),
                              ("disc40-co02", 'tab:grey', '-', r"$\mathrm{Co}=0.2$"),
                              ("disc40-co01", 'tab:red', '-', r"$\mathrm{Co}=0.1$")):
        x, y, v = field(run)
        xs, ys = np.unique(np.round(x, 9)), np.unique(np.round(y, 9))
        grid = np.full((len(ys), len(xs)), np.nan)
        ix = np.searchsorted(xs, np.round(x, 9))
        iy = np.searchsorted(ys, np.round(y, 9))
        grid[iy, ix] = v
        cs = ax1.contour(xs, ys, grid, levels=[0.5], colors=[col],
                         linewidths=1.2, linestyles=ls)
        ax1.plot([], [], ls, color=col, lw=1.2, label=lab)
        del cs
    ax1.set_xlim(4.5 - 0.42, 4.5 + 0.42)
    ax1.set_ylim(2.5 - 0.42, 2.5 + 0.42)
    ax1.set_aspect('equal')
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$y$")
    ax1.legend(fontsize=6.6, loc='upper right', framealpha=0.92,
               handlelength=1.8, borderpad=0.4)
    ax1.set_title(r"(a) disc at $t = 4$, $N_x = 40$", fontsize=9, pad=6)

    # ---- (b) E1 beside the four schemes of Roenby et al. ------------------
    lab = [rf"$({nx},{co:g})$" for nx, co, _ in ROWS]
    x = np.arange(len(ROWS))
    ref = np.array([ROENBY_T1["E1"][(nx, co)] for nx, co, _ in ROWS])
    own = np.array([read(run)["e1_exact"] for _, _, run in ROWS])
    w = 0.16
    for k, (name, col, hatch) in enumerate((("isoAdvector", '0.85', ''),
                                            ("MULES", '0.65', ''),
                                            ("HRIC", '0.45', ''),
                                            ("CICSAM, Roenby et al.", '0.2', ''))):
        ax2.bar(x + (k - 2) * w, ref[:, k], w, color=col, label=name)
    ax2.bar(x + 2 * w, own, w, color='tab:red', label="CICSAM, present solver")
    ax2.set_yscale('log')
    ax2.set_xticks(x)
    ax2.set_xticklabels(lab, fontsize=7.5)
    ax2.set_xlabel(r"$(N_x,\ \mathrm{Co})$")
    ax2.set_ylabel(r"shape error $E_1$")
    ax2.set_ylim(0.008, 12.0)
    ax2.legend(fontsize=6.2, loc='upper center', framealpha=0.92, ncol=2,
               handlelength=1.4, borderpad=0.4, columnspacing=1.0)
    ax2.grid(axis='x', visible=False)
    ax2.set_title("(b) disc in uniform flow", fontsize=9, pad=6)

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    fig.savefig(a.out.replace('.pdf', '.png'), bbox_inches='tight')
    print(f"wrote {a.out} and {a.out.replace('.pdf', '.png')}")


if __name__ == "__main__":
    main()
