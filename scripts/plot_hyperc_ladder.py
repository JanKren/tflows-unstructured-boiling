#!/usr/bin/env python3
"""Structured refinement ladder: pure Hyper-C against the CICSAM blend.

Panel (a) is the result.  R_vol at t = 1 ms is plotted against 1/n, so the
continuum limit sits at the left axis and a convergent scheme walks toward it.
The blend walks the other way.  Both sequences fit `R = R_inf - C n^-p` to
better than 0.03 um, so the extrapolated intercepts (open markers on the axis)
are meaningful, and they land on opposite sides of Scriven: -4.3% for Hyper-C,
+19.7% for the blend.

Two conventions matter here and both are easy to get wrong.

  * R_vol, not R_front.  Scriven predicts a volume-equivalent radius, while
    R_front/R_vol is scheme-dependent (1.007 Hyper-C, 1.035 blend, 1.077 UQ)
    because the blend weight sets the effective interface thickness.  Plotting
    R_front against Scriven would flatter Hyper-C for a reason that has nothing
    to do with its accuracy.  Column 4 of bench-data.dat is R_vol.
  * dt is fixed at 1 us across the ladder.  That is deliberate -- Hyper-C is
    Courant-dependent through alfa_cbc = min(1, alfa_d_til/Co), so dt must not
    move along a refinement study of it.  The consequence is that any temporal
    error is constant along the ladder and lands entirely in the intercept,
    which is the leading explanation on offer for Hyper-C's -4.3% residual.

Panel (b) is the counterweight, and the reason this figure has two panels
instead of one.  The four-fold deformation *rises* with refinement for both
schemes, and Hyper-C's advantage at matched radius narrows from 2.6x to 1.9x.
Showing (a) alone would sell a convergence fix as an isotropy fix.

Usage:  python3 plot_hyperc_ladder.py [--out pics/hyperc_ladder.pdf]
"""

import argparse
import glob
import os

import numpy as np
from scipy.optimize import curve_fit

BETA, AL, R0 = 4.06022, 1.6753e-7, 50.0e-6
T0 = R0 ** 2 / (4 * BETA ** 2 * AL)
R_SCRIVEN = 2 * BETA * np.sqrt(AL * (1e-3 + T0)) * 1e6      # um, volume-equivalent

# block N in a .dom yields N-1 cells per direction
LEVELS = [75, 100, 125, 150]
NCELL = np.array([l - 1 for l in LEVELS], float)

BENCH = {"hyperc": ["Data/hyperc-bench/pure-075.dat",
                    "Data/hyperc-bench/pure-100.dat",
                    "Data/hyperc-bench/pure-125.dat",
                    "Data/hyperc-bench/pure-150.dat"],
         "blend":  ["Data/blend-bench/st75-dt1000.dat",
                    "Data/blend-bench/st100-dt1000.dat",
                    "Data/blend-bench/st125-dt1000.dat",
                    "Data/hyperc-bench/base-150-dt1.dat"]}

FRONTS = {"hyperc": ["Data/pure075", "Data/pure100", "Data/pure125", "Data/pure150"],
          "blend":  ["Data/gamma-dt1000", "Data/gamma-st100",
                     "Data/gamma-st125", "Data/gamma-st150"]}

STYLE = {"hyperc": (r"pure Hyper-C, $\gamma_f=1$", "tab:red", "o", "-"),
         "blend":  (r"blend, $\gamma_f=\cos^2\theta$", "tab:grey", "s", "-")}


def r_at_1ms(path):
    """R_vol at t = 1 ms, with the integrity screen this project insists on."""
    a = np.loadtxt(path)
    dt = np.diff(a[:, 0])
    assert (dt > 0).all(), f"{path}: time column not monotone (spliced?)"
    assert abs(a[0, 3] * 1e6 - 50.0) < 0.1, f"{path}: R(0) is not 50 um"
    assert (np.diff(a[:, 2]) > 0).all(), f"{path}: vapour volume not monotone"
    return a[999, 3] * 1e6


def deviation(directory, ts):
    """Body-diagonal minus axis radius, area-weighted, as a percentage."""
    import pyvista as pv
    f = glob.glob(os.path.join(directory, f"*front-ts{ts:06d}.pvtu"))
    if not f:
        return None
    m = pv.read(f[0])
    s = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(s.cell_data["Area"], float)
    c = np.asarray(s.cell_centers().points, float)
    ctr = (a[:, None] * c).sum(0) / a.sum()
    q = c - ctr
    r = np.linalg.norm(q, axis=1)
    ok = r > 0
    a, q, r = a[ok], q[ok], r[ok]
    sv = ((q / r[:, None]) ** 4).sum(-1)          # cubic invariant, 1 axis .. 1/3 body
    w = np.sqrt(a)
    c0, c1 = np.linalg.lstsq(np.vstack([w, w * sv]).T, w * r, rcond=None)[0]
    return float((a * r).sum() / a.sum() * 1e6), 100 * ((c0 + c1 / 3) / (c0 + c1) - 1)


def fit(n, y, p0):
    f = lambda x, ri, c, p: ri - c * x ** -p
    po, pc = curve_fit(f, n, y, p0=p0, maxfev=400000)
    return po, np.sqrt(np.diag(pc)), np.abs(y - f(n, *po)).max()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/hyperc_ladder.pdf")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True, 'grid.alpha': 0.3})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.3))
    print(f"\n  Scriven R_vol(1 ms) = {R_SCRIVEN:.2f} um\n")

    x = 1.0 / NCELL
    for key, p0 in (("blend", [139., -5e4, 1.3]), ("hyperc", [111., 5e4, 2.3])):
        lab, col, mk, ls = STYLE[key]
        R = np.array([r_at_1ms(p) for p in BENCH[key]])
        (ri, c, p), se, resid = fit(NCELL, R, p0)
        xs = np.linspace(0, x.max() * 1.06, 400)
        with np.errstate(divide='ignore'):
            ys = np.where(xs > 0, ri - c * (1 / np.maximum(xs, 1e-12)) ** -p, ri)
        ax1.plot(xs, ys, ls, color=col, lw=1.0, alpha=0.55)
        ax1.plot(x, R, mk, color=col, ms=5.5, label=lab)
        ax1.plot([0], [ri], mk, color=col, ms=7, mfc='none', mew=1.4)
        ax1.annotate(rf"$R_\infty={ri:.1f}$" "\n" rf"$({100*(ri/R_SCRIVEN-1):+.1f}\%)$,"
                     "\n" rf"$p={p:.2f}$",
                     xy=(0, ri), xytext=(0.055 * x.max(), ri),
                     fontsize=6.4, color=col, va='center')
        print(f"  {key:<7} p={p:5.2f}+/-{se[2]:.2f}  R_inf={ri:7.2f}+/-{se[0]:.2f} um"
              f"  ({100*(ri/R_SCRIVEN-1):+6.2f}%)  max resid {resid:.4f} um")
        print("          " + "  ".join(f"{l}:{v:7.3f}" for l, v in zip(LEVELS, R)))

    ax1.axhline(R_SCRIVEN, ls='--', color='k', lw=1.3, zorder=1)
    ax1.text(x.max() * 0.98, R_SCRIVEN + 0.9, "Scriven", fontsize=7,
             ha='right', va='bottom')
    ax1.set_xlim(-0.0006, x.max() * 1.06)
    # headroom for the blend's intercept annotation, which sits near the top
    ax1.set_ylim(104, 145)
    ax1.set_xlabel(r"$1/n$   (continuum limit at 0)")
    ax1.set_ylabel(r"$R_\mathrm{vol}$ at $t=1$ ms  [$\mu$m]")
    ax1.set_title("(a) trajectory converges only for Hyper-C", fontsize=9)
    # upper right is the empty corner: the blend curve is already below 128 there
    ax1.legend(fontsize=7, loc='upper right')
    sec = ax1.secondary_xaxis('top', functions=(lambda v: v, lambda v: v))
    sec.set_xticks(x)
    sec.set_xticklabels([str(l) for l in LEVELS], fontsize=6.5)
    sec.set_xlabel("mesh level", fontsize=7.5)

    print()
    for key in ("hyperc", "blend"):
        lab, col, mk, ls = STYLE[key]
        dev = []
        for d in FRONTS[key]:
            g = deviation(d, 1000)
            dev.append(np.nan if g is None else g[1])
        ax2.plot(LEVELS, dev, ls, color=col, marker=mk, ms=5.5, lw=1.4, label=lab)
        print(f"  {key:<7} m=4 at 1 ms  " +
              "  ".join(f"{l}:{v:+6.2f}%" for l, v in zip(LEVELS, dev)))
    ax2.set_xticks(LEVELS)
    ax2.set_xlabel("mesh level  (cells per direction $+1$)")
    ax2.set_ylabel(r"body diagonal $-$ axis  [\%]".replace("\\%", "%"))
    ax2.set_title("(b) but the anisotropy still grows", fontsize=9)
    ax2.set_ylim(0, None)
    # lower right is the only empty corner: both curves rise to the right
    ax2.legend(fontsize=7, loc='lower right')

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, bbox_inches='tight')
    fig.savefig(args.out.replace('.pdf', '.png'), bbox_inches='tight')
    print(f"\n  wrote {args.out} and {args.out.replace('.pdf', '.png')}\n")


if __name__ == "__main__":
    main()
