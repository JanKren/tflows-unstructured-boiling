#!/usr/bin/env python3
"""The three bubbles at t = 1 ms: blend, pure Hyper-C, pure ULTIMATE-QUICKEST.

Renders the reconstructed interface itself, coloured by the local departure
from that bubble's own mean radius, so both errors are visible at once:

  * the SHAPE, as colour -- the four-fold pattern, depressed on the coordinate
    axes and raised towards the cube corners;
  * the SIZE, as the actual extent of each sphere -- the three are drawn at
    true relative scale on a common camera, so Hyper-C is visibly the smallest
    and UQ the largest.

That second point is the reason not to normalise each bubble to unit radius,
which would make the shapes easier to compare and would hide the trade the
figure exists to show: freezing the blend weight improves the shape and wrecks
the growth rate, in opposite directions for the two limits.

Colour scale is common to all three panels.  Per-panel scales would make
Hyper-C look as deformed as the blend.

Usage:  python3 plot_blend_bubbles.py [--out pics/blend_bubbles.png]
"""

import argparse
import glob
import os

import numpy as np
import pyvista as pv

BETA, AL, R0 = 4.06022, 1.6753e-7, 50.0e-6
T0 = R0 ** 2 / (4 * BETA ** 2 * AL)
R_SCRIVEN = 2 * BETA * np.sqrt(AL * (1e-3 + T0))          # at t = 1 ms

CASES = [("Data/gamma-dt1000", r"blend, $\gamma_f=\cos^2\theta$"),
         ("Data/gamma-hyperc", r"pure Hyper-C, $\gamma_f=1$"),
         ("Data/gamma-uq", r"pure UQ, $\gamma_f=0$")]
TS = 1000                                                  # dt = 1 us -> 1 ms


def load(directory, ts=TS):
    f = glob.glob(os.path.join(directory, f"*front-ts{ts:06d}.pvtu"))
    if not f:
        return None
    m = pv.read(f[0]).clean()
    sized = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(sized.cell_data["Area"], float)
    cc = np.asarray(sized.cell_centers().points, float)
    ctr = (a[:, None] * cc).sum(0) / a.sum()
    rc = np.linalg.norm(cc - ctr, axis=1)
    R = float((a * rc).sum() / a.sum())
    n = (cc - ctr) / rc[:, None]
    s = (n ** 4).sum(-1)
    w = np.sqrt(a)
    c0, c1 = np.linalg.lstsq(np.vstack([w, w * s]).T, w * rc, rcond=None)[0]
    dev = 100 * ((c0 + c1 / 3) / (c0 + c1) - 1)
    p = np.asarray(m.points, float)
    m.point_data["dev"] = 100.0 * (np.linalg.norm(p - ctr, axis=1) - R) / R
    m.points = p - ctr                       # centred, NOT scaled
    return m, R, dev


def render(mesh, clim, half_extent, out):
    """Common camera and a fixed view extent, so relative sizes are preserved."""
    pl = pv.Plotter(off_screen=True, window_size=(900, 900))
    pl.add_mesh(mesh, scalars="dev", cmap="coolwarm", clim=clim,
                show_scalar_bar=False, smooth_shading=True,
                specular=0.25, specular_power=18)
    pl.set_background("white")
    pl.view_vector((1.0, 0.85, 0.75), viewup=(0, 0, 1))
    pl.camera.SetParallelProjection(True)    # parallel: size on screen == size
    pl.camera.SetParallelScale(half_extent)
    pl.screenshot(out, transparent_background=False)
    pl.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/blend_bubbles.png")
    a = ap.parse_args()

    pv.start_xvfb()
    pv.OFF_SCREEN = True

    got = []
    for d, lab in CASES:
        g = load(d)
        if g is None:
            raise SystemExit(f"no ts{TS} front under {d}")
        got.append((lab,) + g)

    Rmax = max(g[2] for g in got)
    lim = max(abs(np.percentile(g[1].point_data["dev"], [1, 99])).max() for g in got)
    lim = float(np.ceil(lim))
    print(f"\n  R_Scriven(1 ms) = {R_SCRIVEN*1e6:.2f} um")
    print(f"  common colour scale: +-{lim:.0f}%   view half-extent {Rmax*1.08*1e6:.0f} um")

    tmp = "/tmp/_blendbub"
    os.makedirs(tmp, exist_ok=True)
    for i, (lab, m, R, dev) in enumerate(got):
        render(m, (-lim, lim), Rmax * 1.08, f"{tmp}/{i}.png")
        print(f"    {lab:<34} R={R*1e6:7.2f} um ({100*(R/R_SCRIVEN-1):+6.2f}%)"
              f"   body-axis {dev:+6.2f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300})

    fig = plt.figure(figsize=(7.4, 3.0))
    for i, (lab, m, R, dev) in enumerate(got):
        ax = fig.add_subplot(1, 3, i + 1)
        ax.imshow(plt.imread(f"{tmp}/{i}.png"))
        ax.set_title(lab, fontsize=8.5, pad=3)
        # two lines: three single-line captions collide at this width.  R is
        # the FRONT radius, which is not R_vol -- the two differ by the
        # effective interface thickness, and that differs per scheme (1.007
        # for Hyper-C, 1.035 blend, 1.077 UQ), so the convention must be stated
        ax.text(0.5, 0.02, rf"$R_\mathrm{{front}}={R*1e6:.0f}\,\mu$m",
                transform=ax.transAxes, ha='center', va='top', fontsize=7.5)
        ax.text(0.5, -0.06, rf"body$-$axis ${dev:+.1f}\%$",
                transform=ax.transAxes, ha='center', va='top', fontsize=7.5)
        ax.axis("off")

    cax = fig.add_axes([0.32, 0.02, 0.36, 0.035])
    cb = fig.colorbar(cm.ScalarMappable(colors.Normalize(-lim, lim), "coolwarm"),
                      cax=cax, orientation="horizontal")
    cb.set_ticks([-lim, 0, lim])
    cb.ax.tick_params(labelsize=7, length=2, pad=1)
    cb.set_label(r"$(r-\langle R\rangle)/\langle R\rangle$ [%]",
                 fontsize=7.5, labelpad=2)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight")
    fig.savefig(a.out.replace(".png", ".pdf"), bbox_inches="tight")
    print(f"\n  wrote {a.out} and {a.out.replace('.png', '.pdf')}\n")


if __name__ == "__main__":
    main()
