#!/usr/bin/env python3
"""Stefan energy closure against mesh size and time step.

    python3 plot_stefan_energy_convergence.py [--runs DIR]

Reads the per-step records of the series under
`~/runs/caf-review-2026-11/stefan/runs/` through
`analyze_stefan_energy_terms.analyse`, prints the table of
`results/caf-energy-mass-solver.md` and writes
`pics/stefan_energy_convergence.pdf`.

Two configurations of the solver of the submitted manuscript:

    as submitted               its control file
    energy equation converged  no under-relaxation of the energy equation and
                               upwind convection in it, which removes the two
                               defects named in the note through the control
                               file alone

Left panel, |closure| of the whole record against the cell size, at
dt = 5 ms and with dt proportional to h (5 ms at 20 um); right panel, against
the time step on the 100-cell mesh.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import analyze_stefan_energy_terms as st

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = "/home/user/runs/caf-review-2026-11/stefan/runs"

plt.rcParams.update({
    "font.size": 11, "font.family": "serif", "savefig.dpi": 300,
    "axes.grid": True, "grid.alpha": 0.3,
})

# configuration -> (label, {N: run} at dt = 5 ms, {dt[ms]: run} at N = 100,
#                   {N: run} with dt proportional to h)
CONF = {
    "P": ("as submitted",
          {40: "P_N040", 100: "base_paper", 200: "P_N200", 400: "P_N400"},
          {10.0: "P_dt10", 5.0: "base_paper", 2.5: "paper_dt2"},
          {40: "P_s040", 100: "base_paper", 200: "P_s200"}),
    "Pfix": ("energy equation converged",
             {40: "Pfix_N040", 100: "Pfix_N100", 200: "Pfix_N200",
              400: "Pfix_N400"},
             {10.0: "Pfix_dt10", 5.0: "Pfix_N100", 2.5: "Pfix_dt2.5"},
             {40: "Pfix_s040", 100: "Pfix_N100", 200: "Pfix_s200"}),
}
STYLE = {"P": ("k", "o"), "Pfix": ("C0", "s")}
KEYS = ["closure", "(1) vapour cells", "(2) Eq.11 - Eq.7, vapour",
        "(3) estimate into liquid", "(4) iteration lag",
        "(5) booked - produced", "liquid sink"]


def get(run, window, n_expected=None):
    d = os.path.join(RUNS, run)
    f = os.path.join(d, "caf-terms.dat")
    if not os.path.exists(f):
        return None
    r = st.analyse(d, window)
    return r


def order(x, y):
    """Observed order between consecutive levels, on |y|."""
    out = []
    for i in range(len(x) - 1):
        if y[i] == 0 or y[i + 1] == 0:
            out.append(np.nan)
        else:
            out.append(np.log(abs(y[i]) / abs(y[i + 1]))
                       / np.log(x[i] / x[i + 1]))
    return out


def main():
    global RUNS
    if "--runs" in sys.argv:
        RUNS = sys.argv[sys.argv.index("--runs") + 1]
    window = (4.5, 9.0)

    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3))
    for c, (lab, mesh, step, scaled) in CONF.items():
        col, mk = STYLE[c]
        for axis, series, xname, xof, ls, fill, tag in (
                (ax[0], mesh, "N", lambda n: 2.0e3 / n, "-", col,
                 r", $\Delta t = 5$ ms"),
                (ax[0], scaled, "N (dt ~ h)", lambda n: 2.0e3 / n, ":",
                 "none", r", $\Delta t \propto h$"),
                (ax[1], step, "dt [ms]", lambda d: d, "-", col, "")):
            xs, rows = [], []
            for k in sorted(series, reverse=(xname == "dt [ms]")):
                r = get(series[k], window)
                if r is None:
                    continue
                if abs(r["n"] * r["dt"] - 9.0) > 1e-6:
                    print(f"  ({series[k]} incomplete: t_end = "
                          f"{r['n'] * r['dt']:.3f} s, skipped)")
                    continue
                xs.append(k)
                rows.append(r)
            if not rows:
                continue
            print(f"\n{lab}: series in {xname}")
            hdr = f"  {xname:>10}" + "".join(f"{k[:16]:>18}" for k in KEYS) \
                + f"{'window':>9}{'x_front':>9}{'max':>7}{'x_vol':>8}" \
                + f"{'mass(lag)':>11}"
            print(hdr)
            for k, r in zip(xs, rows):
                w = r["whole"]
                print(f"  {k:>10g}" + "".join(f"{100 * w[q]:>18.3f}"
                                               for q in KEYS)
                      + f"{100 * r['window']['closure']:>9.3f}"
                      + f"{100 * r['x']['front, final']:>9.3f}"
                      + f"{100 * r['x']['front, max |e|']:>7.3f}"
                      + f"{100 * r['x']['volume, final']:>8.3f}"
                      + f"{r['mass']['lagged']:>11.2e}")
            xv = np.array([xof(k) for k in xs])
            yv = np.array([100 * r["whole"]["closure"] for r in rows])
            o = order(xv, yv)
            print("  observed order of |closure| between levels: "
                  + ", ".join(f"{v:.2f}" for v in o))
            axis.loglog(xv, np.abs(yv), color=col, marker=mk, ls=ls, lw=1.3,
                        ms=5, mfc=fill, label=lab + tag)
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for a, ticks in ((ax[0], [5, 10, 20, 50]), (ax[1], [2.5, 5, 10])):
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_major_formatter(plain)
        a.xaxis.set_minor_formatter(NullFormatter())
        a.yaxis.set_major_formatter(plain)
        a.set_ylim(0.08, 40)
    ax[0].set_xlabel(r"cell size $h$ [$\mu$m]")
    ax[0].set_ylabel("energy closure error [%]")
    ax[0].set_title("(a)", loc="left", fontsize=11)
    ax[1].set_xlabel(r"time step $\Delta t$ [ms]")
    ax[1].set_title(r"(b) $h = 20~\mu$m", loc="left", fontsize=11)
    h, l = ax[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.16, 1, 1])
    fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.0))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "pics",
                                 f"stefan_energy_convergence.{ext}"),
                    bbox_inches="tight")
    print("\nwritten: pics/stefan_energy_convergence.{pdf,png}")


if __name__ == "__main__":
    main()
