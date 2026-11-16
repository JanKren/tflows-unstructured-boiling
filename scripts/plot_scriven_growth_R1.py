#!/usr/bin/env python3
"""Figure 8 (pics/scriven_growth_R1.{png,pdf}), revision R1.

Bubble radius against time on the 75^3 and 125^3 meshes of BOTH families at
the SAME time step, dt = 2 us, against the analytical Scriven solution and the
geometric-VOF data of Bures & Sato (2021).

The submitted figure drew the structured curves at dt = 10 us and the
polyhedral ones at 2 us (Reviewer 2), and it mixed radius definitions: the
structured curves were the volume-equivalent radius of bench-data.dat, the
polyhedral ones the front radius.  Here every T-Flows curve is R_front, the
radius the paper uses throughout, reconstructed exactly as in
plot_scriven_sensitivity.py (Figures 9 and 10), from the same records:

  structured  Data/blend-dt2/aniso-{075,125}-off-600.dat + their fronts
  polyhedral  Data/polydfix/run-{75,125}.dat + their fronts

The reference curves are the volume-equivalent radius of an octant domain,
shifted to their own virtual origin as in the submitted figure.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import scriven_reference as ref
from plot_scriven_sensitivity import series_from_ratio, ROOT

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'savefig.dpi': 300,
    'axes.linewidth': 0.6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'lines.linewidth': 1.2,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
})

BLUE, ORANGE, GREEN = '#2171b5', '#d94801', '#238b45'


def bures_sato(fname):
    """Volume-equivalent radius of the octant-domain reference runs."""
    t, r = [], []
    with open(os.path.join(ROOT, "Data", fname)) as f:
        for line in f:
            parts = line.strip().split('=')
            if len(parts) == 2:
                v = parts[1].split()
                if len(v) == 3:
                    t.append(float(v[0]))
                    r.append(2.0 * (3.0 * float(v[2]) / (4.0 * np.pi)) ** (1 / 3))
    t, r = np.array(t), np.array(r)
    t0 = r[0] ** 2 / (4.0 * ref.BETA ** 2 * ref.D_L)
    return t + t0, r


def main():
    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.7))

    ta = np.linspace(0.0, 1.5e-3, 600)
    ax.plot(ta * 1e3, 2.0 * ref.BETA * np.sqrt(ref.D_L * ta) * 1e6, 'k-',
            lw=1.5, label='Analytical')

    rows = []
    for lv, ls in (("075", '--'), ("125", '-')):
        t, r = series_from_ratio(
            None, 2e-6,
            bench=f"{ROOT}/Data/blend-dt2/aniso-{lv}-off-600.dat",
            fronts=f"{ROOT}/Data/Aniso/aniso-{lv}-off-600/*front*.pvtu")
        ax.plot((t + ref.T0) * 1e3, r, color=BLUE, ls=ls,
                label=f"Struct. ${int(lv)}^3$")
        rows.append((f"struct-{lv}", t, r))
    for lv, ls in (("75", '--'), ("125", '-')):
        t, r = series_from_ratio(f"run-{lv}", 2e-6)
        ax.plot((t + ref.T0) * 1e3, r, color=ORANGE, ls=ls,
                label=f"Poly. ${lv}^3$")
        rows.append((f"poly-{lv}", t, r))
    for fname, lab, ls in (("vol-cart2.out", r"Bure\v{s} & Sato, $h = 3.91~\mu$m", ':'),
                           ("vol-cart4.out", r"Bure\v{s} & Sato, $h = 1.95~\mu$m", '-.')):
        t, r = bures_sato(fname)
        ax.plot(t * 1e3, r * 1e6, color=GREEN, ls=ls, lw=1.0,
                label=lab.replace(r"\v{s}", "š"))
        rows.append((fname, t - ref.T0, r * 1e6))

    ax.set_xlabel('Time [ms]')
    ax.set_ylabel(r'Bubble radius [$\mu$m]')
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(0.0, 150.0)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.legend(fontsize=7, loc='lower right')
    fig.tight_layout()
    out = os.path.join(ROOT, "pics", "scriven_growth_R1")
    fig.savefig(out + ".png", bbox_inches='tight')
    fig.savefig(out + ".pdf", bbox_inches='tight')
    print("wrote", out + ".{png,pdf}")

    # the numbers the text quotes: error against the analytical radius at
    # fixed simulation time (t = 0 at the initialised 50 um bubble)
    print(f"\n{'curve':16s}" + "".join(f"{'t=%.1f ms' % x:>12s}"
                                       for x in (0.4, 0.6, 0.8, 1.0)))
    for name, t, r in rows:
        line = f"{name:16s}"
        for x in (0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3):
            if t[0] <= x <= t[-1]:
                e = np.interp(x, t, r) / ref.r_scriven_um(x) - 1.0
                line += f"{100*e:+11.2f}%"
            else:
                line += f"{'--':>12s}"
        print(line)


if __name__ == "__main__":
    main()
