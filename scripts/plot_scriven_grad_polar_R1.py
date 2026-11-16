#!/usr/bin/env python3
"""Figure 14 (pics/scriven_grad_polar_R1.{png,pdf}), revision R1.

Drawn from the SAME runs and by the SAME procedure as Table 4
(tab:m4_convergence), so that the figure and the table cannot disagree:

  runs       structured 75^3 and 125^3 at dt = 2 us (Data/Aniso/aniso-*-off),
             polyhedral 125^3 at dt = 2 us, all at t = 1.0 ms (ts = 500)
  sample     liquid-side interface cells over the whole bubble: alpha > 0.5
             with at least one face neighbour at alpha < 0.5
  quantity   |grad T . n|, n the radial unit vector from the vapour centroid
  modes      least-squares fit of m = 0..8 to the scattered (theta, g) points,
             theta the azimuth about the z axis; amplitudes relative to m = 0

The submitted figure used the dt = 10 us structured runs, a one-cell equatorial
slice within two cells of the ANALYTICAL radius (both phases), 5 degree bins
and an FFT of the binned means.  That is the procedure the measurement
appendix calls unreliable, and it gave 15.4 / 11.1 % where the table has
10.18 / 14.39 % (Reviewer 1, comment 5).

On the polyhedral mesh a face neighbour is a cell sharing at least three
vertices, found on a shell about the interface only.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_anisotropy_attribution import load_lattice, fit_bubble

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'savefig.dpi': 300,
    'axes.grid': False,
    'lines.linewidth': 1.5,
})

HERE = os.path.dirname(os.path.abspath(__file__))
ANISO = os.path.join(HERE, "Data", "Aniso")
POLY = ("/home/user/runs/tflows-vof/Validation/Scriven/Scriven-Poly-125/"
        "Scriven-Poly-125_dual-ts000500.pvtu")
OUT = os.path.join(HERE, "pics", "scriven_grad_polar_R1.png")
M_MAX = 8


def fourier_lsq(theta, g, m_max=8):
    """Same fit as measure_sigma_sweep.fourier_lsq (the Table 4 procedure);
    repeated here because importing that script runs its sweep."""
    cols = [np.ones_like(theta)]
    for m in range(1, m_max + 1):
        cols += [np.cos(m * theta), np.sin(m * theta)]
    coef, *_ = np.linalg.lstsq(np.stack(cols, axis=1), g, rcond=None)
    a0 = coef[0]
    amp, pha = {}, {}
    for m in range(1, m_max + 1):
        c, sn = coef[2 * m - 1], coef[2 * m]
        amp[m] = float(np.hypot(c, sn) / abs(a0))
        pha[m] = float(np.degrees(np.arctan2(sn, c)) / m) % (360.0 / m)
    return amp, pha


def fourier_curve(theta, g, th_plot):
    """The fitted series itself, for drawing."""
    cols = [np.ones_like(theta)]
    for m in range(1, M_MAX + 1):
        cols += [np.cos(m * theta), np.sin(m * theta)]
    coef, *_ = np.linalg.lstsq(np.stack(cols, axis=1), g, rcond=None)
    out = np.full_like(th_plot, coef[0])
    for m in range(1, M_MAX + 1):
        out += coef[2 * m - 1] * np.cos(m * th_plot)
        out += coef[2 * m] * np.sin(m * th_plot)
    return out


def measure_structured(path):
    d = load_lattice(path)
    fit = fit_bubble(d)
    is_vap = d["vof"] < 0.5
    nb = np.zeros_like(is_vap)
    for ax in range(3):
        nb |= np.roll(is_vap, 1, ax) | np.roll(is_vap, -1, ax)
    mask = (~is_vap) & nb
    q = (fit["P"] - fit["ctr"])[mask]
    g = d["grad_saved"][mask]
    return q, g


def measure_polyhedral(path):
    mesh = pv.read(path)
    if mesh.n_cells == 0:
        sys.exit(f"{path} opens but yields no cells")
    vof = np.asarray(mesh.cell_data["Vof Sharp [1]"], dtype=np.float64)
    vol = np.asarray(mesh.cell_data["Grid Cell Volume [m^3]"], dtype=np.float64)
    grad = np.asarray(mesh.cell_data["Temperature Gradients [K/m]"],
                      dtype=np.float64)
    cc = mesh.cell_centers().points

    vap = vof < 0.5
    ctr = np.average(cc[vap], axis=0, weights=vol[vap])
    r = np.linalg.norm(cc - ctr, axis=1)
    R = (3.0 * (vol * (1.0 - vof)).sum() / (4.0 * np.pi)) ** (1.0 / 3.0)
    h = float(np.cbrt(np.median(vol)))

    # shell about the interface, wide enough to hold every cut cell and its
    # neighbours whatever the deformation
    shell = np.where(np.abs(r - R) < 8.0 * h)[0]
    sub = mesh.extract_cells(shell)
    ids = np.asarray(sub.cell_data["vtkOriginalCellIds"])

    # cell-to-vertex incidence of the shell
    rows, cols = [], []
    for i in range(sub.n_cells):
        p = sub.get_cell(i).point_ids
        rows.extend([i] * len(p))
        cols.extend(p)
    C = sparse.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)),
                          shape=(sub.n_cells, sub.n_points))
    shared = (C @ C.T).tocoo()
    face = (shared.data >= 3) & (shared.row != shared.col)
    i, j = shared.row[face], shared.col[face]

    vap_s = vof[ids] < 0.5
    has_vap_nb = np.zeros(sub.n_cells, dtype=bool)
    np.logical_or.at(has_vap_nb, i, vap_s[j])
    band = ids[(~vap_s) & has_vap_nb]

    return cc[band] - ctr, grad[band]


def reduce(q, g):
    r = np.linalg.norm(q, axis=1)
    n = q / r[:, None]
    theta = np.arctan2(q[:, 1], q[:, 0])
    gn = np.abs((g * n).sum(-1))
    amp, pha = fourier_lsq(theta, gn, m_max=M_MAX)
    return dict(theta=theta, gn=gn, amp=amp, pha=pha, n=len(r), R=r.mean(),
                mean=gn.mean())


def main():
    cases = [
        (r"Structured $75^3$", "s",
         os.path.join(ANISO, "aniso-075-off",
                      "Scriven-Struct-075-ts000500.pvtu")),
        (r"Structured $125^3$", "s",
         os.path.join(ANISO, "aniso-125-off",
                      "Scriven-Struct-125-ts000500.pvtu")),
        (r"Polyhedral $125^3$", "p", POLY),
    ]
    res = []
    for label, kind, path in cases:
        q, g = (measure_structured if kind == "s" else measure_polyhedral)(path)
        d = reduce(q, g)
        d["label"] = label
        res.append(d)
        print(f"{label:22s} n = {d['n']:6d}  <R> = {d['R']*1e6:7.2f} um  "
              f"<|gT.n|> = {d['mean']:.3e} K/m  "
              + "  ".join(f"m{m}={100*d['amp'][m]:.2f}%"
                          for m in range(1, M_MAX + 1))
              + f"  phi4 = {d['pha'][4]:.2f} deg")

    fig = plt.figure(figsize=(16, 8.6))
    fig.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.07,
                        hspace=0.42, wspace=0.28)
    th = np.linspace(-np.pi, np.pi, 721)
    rng = np.random.default_rng(0)

    for k, d in enumerate(res):
        ax = fig.add_subplot(2, 3, k + 1, projection='polar')
        # every band cell enters the fit; a fixed random third is drawn
        show = rng.random(d["n"]) < 1.0 / 3.0
        ax.scatter(d["theta"][show], d["gn"][show], s=2, alpha=0.25,
                   color='0.45', label='cells', rasterized=True)
        ax.plot(th, fourier_curve(d["theta"], d["gn"], th), color='red',
                lw=2, label='fitted modes $m \\leq 8$')
        rmax = 1.02 * float(np.percentile(d["gn"], 99.5))
        ax.set_rlim(0.0, rmax)
        for step in (2.0e4, 2.5e4, 4.0e4, 5.0e4, 1.0e5):
            ticks = np.arange(step, 0.9 * rmax, step)
            if len(ticks) <= 4:
                break
        ax.set_rticks(ticks)
        ax.set_yticklabels([])
        for t in ticks:
            ax.text(np.radians(22.5), t, f"{t/1e4:g}", fontsize=8,
                    ha='center', va='center', zorder=5,
                    bbox=dict(facecolor='white', edgecolor='none',
                              alpha=0.85, pad=0.15))
        ax.set_title(f'{d["label"]}\n'
                     f'mean $|\\nabla T\\cdot\\mathbf{{n}}|$ = '
                     f'{d["mean"]:.2e} K/m', pad=18, fontsize=10)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=2,
                  fontsize=8, frameon=False)
        ax.text(0.5, -0.19,
                r'rings: $|\nabla T\cdot\mathbf{n}|$ in $10^{4}$ K/m',
                transform=ax.transAxes, ha='center', va='top', fontsize=8,
                color='0.25')

    for k, d in enumerate(res):
        ax = fig.add_subplot(2, 3, k + 4)
        a = np.array([100 * d["amp"][m] for m in range(1, M_MAX + 1)])
        ax.bar(range(1, M_MAX + 1), a,
               color=['C3' if m == 4 else 'C0' for m in range(1, M_MAX + 1)],
               edgecolor='black', linewidth=0.5)
        ax.set_xlabel('Azimuthal mode $m$')
        ax.set_ylabel('Relative amplitude [%]')
        ax.set_title(d["label"], fontsize=10)
        ax.set_xticks(range(1, M_MAX + 1))
        ax.set_ylim(0, a.max() * 1.3 + 0.3)
        ax.annotate(f'm=4: {a[3]:.2f}%', xy=(4, a[3]),
                    xytext=(5.6, a[3] * 1.12 + 0.1), fontsize=9, ha='center',
                    arrowprops=dict(arrowstyle='->', color='C3', lw=1))

    fig.savefig(OUT, bbox_inches='tight')
    fig.savefig(OUT.replace('.png', '.pdf'), bbox_inches='tight')
    print("wrote", OUT)


if __name__ == "__main__":
    main()
