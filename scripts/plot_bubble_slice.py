#!/usr/bin/env python3
"""Equatorial cross-section of the interface, with and without the transport.

The rendered surfaces show the deformation as colour; this shows it as shape.
A cross-section is the more honest of the two, because it also shows how small
the effect is in absolute terms -- at 8% the unmodified bubble still looks
round to the eye, which is exactly why the mode went unnoticed for so long.

Two panels:

  (a) the z = 0 outline at true scale.  Both curves are visibly circular.  This
      panel exists to stop the second one being read as a literal picture of
      the interface.
  (b) the same outlines with the radial departure from the mean amplified ten
      times about a circle.  The four-fold shape appears: the unmodified bubble
      is drawn out at 45, 135, 225, 315 degrees and pulled in on the axes,
      while the transport-removed one stays round.

CAREFUL with what those 45-degree directions are.  In the z = 0 plane they are
FACE diagonals, s = nx^4+ny^4+nz^4 = 1/2, not body diagonals -- the body
diagonals do not lie in this plane at all.  Since R(s) is linear in s to good
approximation, the contrast visible here is at most (1 - 1/2)/(1 - 1/3) = 3/4
of the axis-to-body-diagonal figure quoted everywhere else.  Measured over
+-12 degree bins it comes out lower still -- +5.64% against +0.68% -- because
those bins average s over a range and compress the contrast.  Quote the 3-D
regression for numbers; use this figure to show the shape.

Points are taken from the front within a thin slab about z = 0 rather than by
interpolating a cut plane, so the curve is the reconstructed interface itself
and carries its own scatter rather than a smoothed version of it.

Usage:  python3 plot_bubble_slice.py [--out pics/bubble_slice.pdf]
"""

import argparse
import glob
import os

import numpy as np
import pyvista as pv

NOADV = "Data/inflate-unif"
FULL = "Data/unif-transport"         # uniform mdot, transport ACTIVE
AMP = 10.0


def outline(path, slab=0.06):
    """(theta, r) of the front within a thin equatorial slab, sorted by angle."""
    m = pv.read(path)
    sized = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(sized.cell_data["Area"], float)
    c = np.asarray(sized.cell_centers().points, float)
    ctr = (a[:, None] * c).sum(0) / a.sum()
    q = c - ctr
    r3 = np.linalg.norm(q, axis=1)
    R = float((a * r3).sum() / a.sum())
    sel = np.abs(q[:, 2]) < slab * R          # slab thickness scales with radius
    q = q[sel]
    th = np.arctan2(q[:, 1], q[:, 0])
    r = np.hypot(q[:, 0], q[:, 1])
    o = np.argsort(th)
    return th[o], r[o], R


def smooth_polar(th, r, nbin=180):
    """Mean radius per angular bin, wrapped, so the curve closes."""
    edges = np.linspace(-np.pi, np.pi, nbin + 1)
    idx = np.clip(np.digitize(th, edges) - 1, 0, nbin - 1)
    num = np.bincount(idx, weights=r, minlength=nbin)
    cnt = np.bincount(idx, minlength=nbin)
    ok = cnt > 0
    mid = 0.5 * (edges[:-1] + edges[1:])
    t, v = mid[ok], num[ok] / cnt[ok]
    return np.append(t, t[0] + 2 * np.pi), np.append(v, v[0])


def pick(directory, target):
    best = None
    for f in sorted(glob.glob(os.path.join(directory, "*front-ts*.pvtu"))):
        try:
            th, r, R = outline(f)
        except Exception:
            continue
        if len(r) < 50:
            continue
        if best is None or abs(R - target) < abs(best[3] - target):
            best = (f, th, r, R)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/bubble_slice.pdf")
    ap.add_argument("--target", type=float, default=119.0e-6)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True,
                         'grid.alpha': 0.3})

    fn, thn, rn, Rn = pick(NOADV, a.target)
    ff, thf, rf, Rf = pick(FULL, a.target)
    print(f"\n  transport active  : {os.path.basename(ff)}  R = {Rf*1e6:.2f} um, "
          f"{len(rf)} slab points")
    print(f"  transport removed : {os.path.basename(fn)}  R = {Rn*1e6:.2f} um, "
          f"{len(rn)} slab points")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7),
                                   subplot_kw={'projection': 'polar'})

    series = ((thf, rf, Rf, 'tab:grey', 'transport active'),
              (thn, rn, Rn, 'tab:red', 'transport removed'))

    # ---- (a) true scale ---------------------------------------------------
    for th, r, R, col, lab in series:
        t, v = smooth_polar(th, r)
        ax1.plot(t, v * 1e6, '-', color=col, lw=1.6, label=lab)
    ax1.set_title("(a) true scale", fontsize=9, pad=14)
    ax1.set_rlabel_position(112)
    ax1.tick_params(labelsize=7)
    ax1.legend(fontsize=7, loc='lower center', bbox_to_anchor=(0.5, -0.30),
               ncol=1, framealpha=0.9)

    # ---- (b) departure from the mean, amplified ---------------------------
    for th, r, R, col, lab in series:
        t, v = smooth_polar(th, r)
        vm = v.mean()
        ax2.plot(t, (vm + AMP * (v - vm)) * 1e6, '-', color=col, lw=1.6)
    ax2.set_title(rf"(b) departure from $\langle R\rangle$, $\times{AMP:.0f}$",
                  fontsize=9, pad=14)
    ax2.set_rlabel_position(112)
    ax2.tick_params(labelsize=7)
    for ang in (45, 135, 225, 315):
        ax2.plot([np.radians(ang)] * 2, ax2.get_ylim(), ':', color='0.6', lw=0.8)
    ax2.annotate("face diagonals ($s=1/2$)", xy=(np.radians(45), ax2.get_ylim()[1]),
                 xytext=(0.5, -0.30), textcoords='axes fraction',
                 fontsize=7, ha='center', color='0.35')

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    fig.savefig(a.out.replace('.pdf', '.png'), bbox_inches='tight')
    for th, r, R, col, lab in series:
        t, v = smooth_polar(th, r)
        ax_r = v[np.argmin(np.abs(((t + np.pi/4) % (np.pi/2)) - np.pi/4))]
        # equatorial contrast, measured: mean at 45-ish vs mean at 0/90-ish
        d45 = np.abs(((np.degrees(t) % 90) - 45)) < 12
        d00 = (np.degrees(t) % 90) < 12
        if d45.any() and d00.any():
            print(f"  {lab:<18} face-diagonal minus axis = "
                  f"{100*(v[d45].mean()/v[d00].mean() - 1):+.3f}%  "
                  f"(binned; below the 3-D fit, which is the number to quote)")
    print(f"  wrote {a.out} and {a.out.replace('.pdf', '.png')}")


if __name__ == "__main__":
    main()
