#!/usr/bin/env python3
"""Anisotropy at t = 1 ms for the three settings of the CICSAM blend weight.

gamma_f = cos^2(angle between grad(alpha) and the cell connection) is the only
explicitly orientation-dependent term in the limiter.  Freezing it removes that
dependence: 1 gives pure Hyper-C, 0 pure ULTIMATE-QUICKEST.  Both endpoints are
markedly less anisotropic than the blend that interpolates between them, which
is the signature of the *variation* of gamma_f being the problem rather than
either scheme.

Panel (a) is the paper's shape metric -- interface radius against the cubic
invariant s = nx^4+ny^4+nz^4, unity along a coordinate axis and 1/3 along a body
diagonal, so a negative slope is a bubble bulging towards the corners.
Panel (b) is the same data as normalised radius against azimuth, where the
four-fold mode is directly visible.

CAREFUL: at fixed time the three runs are at different radii (R_vol = 116.7,
107.0 and 132.3 um), because freezing the blend also changes the effective
interface thickness and hence the mass transfer rate.  The four-fold mode
collapses on R/h, so this figure compares at matched TIME and not at matched
R/h; the radius of each is printed in the legend for that reason, and the
matched-radius comparison is the one quoted in the text.

Usage:  python3 plot_blend_limits.py [--out pics/blend_limits.pdf]
"""

import argparse
import os

import numpy as np
import pyvista as pv

CASES = [("Data/gamma-dt1000", r"blend, $\gamma_f=\cos^2\theta$", "tab:grey", "o"),
         ("Data/gamma-hyperc", r"pure Hyper-C, $\gamma_f=1$", "tab:red", "s"),
         ("Data/gamma-uq", r"pure UQ, $\gamma_f=0$", "tab:blue", "^")]
TS = 1000          # dt = 1 us, so this is exactly t = 1 ms


def front(directory, ts=TS):
    import glob
    f = glob.glob(os.path.join(directory, f"*front-ts{ts:06d}.pvtu"))
    if not f:
        return None
    m = pv.read(f[0])
    if m.n_cells == 0:
        return None
    sized = m.compute_cell_sizes(length=False, area=True, volume=False)
    a = np.asarray(sized.cell_data["Area"], float)
    c = np.asarray(sized.cell_centers().points, float)
    ctr = (a[:, None] * c).sum(0) / a.sum()
    q = c - ctr
    r = np.linalg.norm(q, axis=1)
    ok = r > 0
    a, q, r = a[ok], q[ok], r[ok]
    n = q / r[:, None]
    s = (n ** 4).sum(-1)
    theta = np.arctan2(q[:, 1], q[:, 0])
    # equatorial band for the polar panel, so it is one shell not a projection
    band = np.abs(n[:, 2]) < 0.12
    w = np.sqrt(a)
    c0, c1 = np.linalg.lstsq(np.vstack([w, w * s]).T, w * r, rcond=None)[0]
    mean = float((a * r).sum() / a.sum())
    dev = 100 * ((c0 + c1 / 3) / (c0 + c1) - 1)
    return dict(a=a, r=r, s=s, theta=theta[band], rb=r[band], ab=a[band],
                c0=c0, c1=c1, mean=mean, dev=dev)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="pics/blend_limits.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'font.family': 'serif',
                         'savefig.dpi': 300, 'axes.grid': True,
                         'grid.alpha': 0.3})

    fig = plt.figure(figsize=(7.4, 3.4))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection='polar')

    edges = np.linspace(1 / 3, 1, 13)
    mid = 0.5 * (edges[:-1] + edges[1:])
    tb = np.linspace(-np.pi, np.pi, 49)
    tmid = 0.5 * (tb[:-1] + tb[1:])

    print(f"\n  body-axis at t = 1 ms (dt = 1 us, ts{TS})")
    for d, lab, col, mark in CASES:
        g = front(d)
        if g is None:
            print(f"  {lab}: no ts{TS} front")
            continue
        norm = g["mean"]
        # ---- (a) cubic invariant ------------------------------------------
        idx = np.clip(np.digitize(g["s"], edges) - 1, 0, len(mid) - 1)
        num = np.bincount(idx, weights=g["a"] * g["r"], minlength=len(mid))
        den = np.bincount(idx, weights=g["a"], minlength=len(mid))
        m = den > 0
        ax1.plot(mid[m], (num[m] / den[m]) / norm, mark, color=col, ms=3.5, lw=0,
                 label=rf"{lab}, $R={norm*1e6:.0f}\,\mu$m")
        ss = np.linspace(1 / 3, 1, 50)
        ax1.plot(ss, (g["c0"] + g["c1"] * ss) / norm, '-', color=col, lw=1.2)
        # ---- (b) polar ----------------------------------------------------
        j = np.clip(np.digitize(g["theta"], tb) - 1, 0, len(tmid) - 1)
        nu = np.bincount(j, weights=g["ab"] * g["rb"], minlength=len(tmid))
        de = np.bincount(j, weights=g["ab"], minlength=len(tmid))
        k = de > 0
        t, v = tmid[k], (nu[k] / de[k]) / (nu[k].sum() / de[k].sum() if False else norm)
        ax2.plot(np.append(t, t[0] + 2 * np.pi), np.append(v, v[0]), '-',
                 color=col, lw=1.5)
        print(f"    {lab:<34} R={norm*1e6:7.2f} um   body-axis {g['dev']:+7.3f}%")

    ax1.axhline(1.0, color='k', lw=0.6, ls=':')
    ax1.set_xlabel(r"cubic invariant $s = n_x^4+n_y^4+n_z^4$")
    ax1.set_ylabel(r"$R(s)\,/\,\langle R\rangle$")
    ax1.set_title("(a) shape at $t = 1$ ms", fontsize=9, pad=6)
    lo, hi = ax1.get_ylim()
    ax1.set_ylim(lo - 0.10 * (hi - lo), hi)
    ax1.text(0.02, 0.03, "body diagonal", transform=ax1.transAxes,
             fontsize=6.5, ha='left', va='bottom', color='0.35')
    ax1.text(0.98, 0.03, "axis", transform=ax1.transAxes,
             fontsize=6.5, ha='right', va='bottom', color='0.35')
    ax1.legend(fontsize=6.4, loc='upper right', framealpha=0.92,
               handlelength=1.4, borderpad=0.35)

    ax2.set_title(r"(b) $R(\theta)/\langle R\rangle$, equatorial band",
                  fontsize=9, pad=14)
    # the deviations are a few per cent, so the radial axis has to be zoomed
    # onto them; drawn from 0 they are invisible and the panel says nothing
    ax2.set_rlim(0.94, 1.05)
    ax2.set_rticks([0.96, 1.00, 1.04])
    ax2.set_rlabel_position(112)
    ax2.tick_params(labelsize=7)
    ax2.grid(alpha=0.35)

    fig.tight_layout()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, bbox_inches='tight')
    fig.savefig(a.out.replace('.pdf', '.png'), bbox_inches='tight')
    print(f"\n  wrote {a.out} and {a.out.replace('.pdf', '.png')}\n")


if __name__ == "__main__":
    main()
