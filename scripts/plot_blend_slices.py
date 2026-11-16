#!/usr/bin/env python3
"""Equatorial cross-section of the three blend settings at t = 1 ms.

The 3-D render (`plot_blend_bubbles.py`) shows the deformation as colour; this
shows it as shape, and separates the two errors onto two panels:

  (a) the z = 0 outlines at TRUE scale.  The size error dominates and is the
      point of the panel -- pure Hyper-C undergrown, pure UQ overgrown, the
      blend between them.  All three still look like circles, which is why the
      four-fold mode has to be extracted by regression rather than seen.
  (b) each outline with the departure from its OWN mean radius magnified ten
      times.  Size is divided out, so only shape remains: the blend and UQ
      strongly four-fold, Hyper-C markedly rounder.

The analytical Scriven solution is drawn in both panels, but it needs care.
Scriven predicts a VOLUME-EQUIVALENT radius, 116.39 um at 1 ms, while the solid
outlines are FRONT radii from the reconstructed interface.  The two diverge by
the effective interface thickness, which is what the blend weight controls and
which differs per scheme: 1.007 for Hyper-C, 1.035 for the blend, 1.077 for UQ.
Laying one volume-equivalent circle against three front outlines would compare
three different conventions against one and would flatter Hyper-C, whose front
nearly coincides with its own R_vol, while penalising UQ twice over.

Panel (a) therefore draws both for each scheme: the solid front outline, and a
thin dashed circle at that scheme's own R_vol.  The black dashed circle is
Scriven.  Comparing dashed against black is the like-for-like accuracy
statement; the solid outline shows where the reconstruction puts the interface.

Panel (b) needs no such care -- the analytical shape is a sphere, so the
reference is a circle at 1.0 whatever the radius.

Usage:  python3 plot_blend_slices.py [--out pics/blend_slices.pdf]
"""

import argparse
import glob
import os

import numpy as np
import pyvista as pv

CASES = [("Data/gamma-dt1000", "st75-dt1000", r"blend, $\gamma_f=\cos^2\theta$", "tab:grey"),
         ("Data/gamma-hyperc", "st75-hyperc", r"pure Hyper-C, $\gamma_f=1$", "tab:red"),
         ("Data/gamma-uq", "st75-uq", r"pure UQ, $\gamma_f=0$", "tab:blue")]
BETA, AL, R0 = 4.06022, 1.6753e-7, 50.0e-6
T0 = R0 ** 2 / (4 * BETA ** 2 * AL)
R_SCRIVEN = 2 * BETA * np.sqrt(AL * (1e-3 + T0)) * 1e6      # um, volume-equivalent
TS = 1000
# these deviations are ~6%, eight times the no-transport case the x10 in
# plot_bubble_slice.py was sized for; x10 here drives the curve past the radial
# limits and matplotlib wraps it through the origin
AMP = 4.0


def outline(directory, ts=TS, slab=0.06):
    f = glob.glob(os.path.join(directory, f"*front-ts{ts:06d}.pvtu"))
    if not f:
        return None
    m = pv.read(f[0])
    sized = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(sized.cell_data["Area"], float)
    c = np.asarray(sized.cell_centers().points, float)
    ctr = (a[:, None] * c).sum(0) / a.sum()
    q = c - ctr
    R = float((a * np.linalg.norm(q, axis=1)).sum() / a.sum())
    sel = np.abs(q[:, 2]) < slab * R
    q = q[sel]
    th = np.arctan2(q[:, 1], q[:, 0])
    r = np.hypot(q[:, 0], q[:, 1])
    o = np.argsort(th)
    return th[o], r[o], R


def binned(th, r, nbin=180):
    edges = np.linspace(-np.pi, np.pi, nbin + 1)
    idx = np.clip(np.digitize(th, edges) - 1, 0, nbin - 1)
    num = np.bincount(idx, weights=r, minlength=nbin)
    cnt = np.bincount(idx, minlength=nbin)
    ok = cnt > 0
    mid = 0.5 * (edges[:-1] + edges[1:])
    t, v = mid[ok], num[ok] / cnt[ok]
    return np.append(t, t[0] + 2 * np.pi), np.append(v, v[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/blend_slices.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True,
                         'grid.alpha': 0.3})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.9),
                                   subplot_kw={'projection': 'polar'})
    print(f"\n  equatorial slice at t = 1 ms (ts{TS})")
    print(f"  analytical Scriven R_vol = {R_SCRIVEN:.2f} um")
    span = []
    # volume-equivalent radii, read from the solver's own per-step record
    rvol = {}
    for _, bench, _, _ in CASES:
        arr = np.loadtxt(f"Data/blend-bench/{bench}.dat")
        rvol[bench] = arr[TS - 1, 3] * 1e6

    circ = np.linspace(-np.pi, np.pi, 361)
    for d, bench, lab, col in CASES:
        ls = "-"
        g = outline(d)
        if g is None:
            print(f"  {lab}: no ts{TS} front")
            continue
        th, r, R = g
        t, v = binned(th, r)
        ax1.plot(t, v * 1e6, ls, color=col, lw=1.6,
                 label=rf"{lab}: $R_\mathrm{{front}}={R*1e6:.0f}$, "
                       rf"$R_\mathrm{{vol}}={rvol[bench]:.0f}\,\mu$m")
        # that scheme's own volume-equivalent radius: this is the quantity
        # Scriven predicts, so dashed-against-black is the fair comparison
        ax1.plot(circ, np.full_like(circ, rvol[bench]), '--', color=col, lw=0.9,
                 alpha=0.85)
        vm = v.mean()
        mag = (vm + AMP * (v - vm)) / vm
        span.append((mag.min(), mag.max()))
        ax2.plot(t, mag, ls, color=col, lw=1.6)
        d45 = np.abs(((np.degrees(t) % 90) - 45)) < 12
        d00 = (np.degrees(t) % 90) < 12
        print(f"    {lab:<34} R_front={R*1e6:7.2f}  R_vol={rvol[bench]:7.2f} "
              f"({100*(rvol[bench]/R_SCRIVEN-1):+6.2f}%)   "
              f"face-diag minus axis {100*(v[d45].mean()/v[d00].mean()-1):+6.2f}%")

    ax1.plot(circ, np.full_like(circ, R_SCRIVEN), '--', color='k', lw=1.6,
             label=rf"Scriven, $R_\mathrm{{vol}}={R_SCRIVEN:.0f}\,\mu$m",
             zorder=5)
    ax1.set_title("(a) true scale", fontsize=9, pad=14)
    ax1.set_rlabel_position(112)
    ax1.tick_params(labelsize=7)
    ax1.legend(fontsize=6.3, loc='lower center', bbox_to_anchor=(0.5, -0.34),
               framealpha=0.92, handlelength=1.5)

    ax2.set_title(rf"(b) departure from own $\langle R\rangle$, $\times{AMP:.0f}$",
                  fontsize=9, pad=14)
    # size the radial axis from the data, with a margin, so nothing is clipped
    lo = min(x[0] for x in span); hi = max(x[1] for x in span)
    pad = 0.08 * (hi - lo)
    ax2.set_rlim(lo - pad, hi + pad)
    # the inner tick lands where the four lobes cross; drop it and put the
    # labels in the gap between lobes rather than on one
    ax2.set_rticks([1.0, round(hi, 2)])
    ax2.set_rlabel_position(100)
    ax2.tick_params(labelsize=7)
    # the analytical shape is a sphere, so the reference is a circle at 1.0
    ax2.plot(circ, np.ones_like(circ), '--', color='k', lw=1.6, zorder=5)
    ax2.text(np.radians(20), 1.0, "Scriven", fontsize=6.5, color='k',
             ha='left', va='bottom')
    for ang in (45, 135, 225, 315):
        ax2.plot([np.radians(ang)] * 2, ax2.get_ylim(), ':', color='0.6', lw=0.8)
    ax2.text(0.5, -0.13, r"dotted: face diagonals ($s=1/2$)",
             transform=ax2.transAxes, ha='center', va='top',
             fontsize=7, color='0.35')

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    fig.savefig(a.out.replace('.pdf', '.png'), bbox_inches='tight')
    print(f"\n  wrote {a.out} and {a.out.replace('.pdf', '.png')}\n")


if __name__ == "__main__":
    main()
