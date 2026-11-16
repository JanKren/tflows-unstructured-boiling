#!/usr/bin/env python3
"""Figure 13 (pics/mass_history_R1.{pdf,png}), revision R1.

Same records and same drawing as plot_mass_history.py, with the outflow
re-integrated by the quadrature of the SCHEME for all three cases.

The hooks accumulated the egress with the trapezoidal rule.  The solver
advances the volume fraction of step n with the face fluxes of step n-1, so the
consistent rule is the lagged one, E_n = sum_{k<=n} out_{k-1} dt.  The choice
is not cosmetic: the submitted statement that halving the time step doubles the
conservation error compared one run integrated by the trapezoidal rule with one
integrated at the end of the step (results/caf-energy-mass-solver.md, C).

The outflow rate is recovered from the recorded cumulative egress
E_n = M_s - M_d:  out_1 = E_1/dt (the hook seeds out_prev with out_1) and
out_n = 2 (E_n - E_{n-1})/dt - out_{n-1}.

With the lagged rule the error is 2 rho_v/rho_l times the expelled mass
fraction on every record, which is the signature of the sign of the liquid
source term -m_dot/rho_l in the volume-fraction equation.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import plot_mass_history as base       # CASES, ROOT, rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_V, RHO_L = 0.597, 958.4


def lagged(a):
    """t, M_d/M_0 and the conservation error with the lagged quadrature."""
    t, m_d, m_s = a[:, 0], a[:, 1], a[:, 2]
    dt = t[1] - t[0]
    if not np.allclose(np.diff(t), dt, rtol=1e-6):
        sys.exit("time column is not uniform: the record is spliced")
    e = m_s - m_d
    out = np.empty_like(e)
    out[0] = e[0] / dt
    for n in range(1, len(e)):
        out[n] = 2.0 * (e[n] - e[n - 1]) / dt - out[n - 1]
    prev = np.concatenate(([out[0]], out[:-1]))
    m0 = m_d[0] + prev[0] * dt
    err = (m_d + np.cumsum(prev) * dt) / m0 - 1.0
    return t, m_d / m0, err


def main():
    have = [c for c in base.CASES
            if os.path.exists(os.path.join(base.ROOT, c[1], "mass-balance.dat"))]
    if len(have) != len(base.CASES):
        sys.exit("a mass-balance.dat record is missing")

    fig, axes = plt.subplots(1, len(have), figsize=(4.7 * len(have), 3.8))
    print(f"  {'case':<34}{'expelled':>10}{'final error':>13}"
          f"{'error/expelled':>16}   (2 rho_v/rho_l = {2*RHO_V/RHO_L:.4e})")
    for ax, (lab, d, tol, unit, scale) in zip(axes, have):
        a = np.loadtxt(os.path.join(base.ROOT, d, "mass-balance.dat"))
        t, md, err = lagged(a)
        t = t * scale

        ax.axhline(0.0, color="k", lw=0.9)
        ax.plot(t, err, color="C0", lw=1.9)
        lim = 1.15 * np.abs(err).max()
        ax.set_ylim(-lim, lim)
        ax.set_xlim(0, t.max())
        ax.set_xlabel(f"time [{unit}]")
        ax.set_ylabel(r"$M_s/M_0 - 1$", color="C0")
        ax.tick_params(axis="y", labelcolor="C0")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        tw = ax.twinx()
        tw.plot(t, md, color="C3", lw=1.5, ls="--")
        tw.set_ylim(min(md) - 0.02, 1.02)
        tw.set_ylabel(r"$M_d/M_0$", color="C3")
        tw.tick_params(axis="y", labelcolor="C3")
        tw.grid(False)
        ax.set_title(lab, fontsize=16, loc="left")

        frac = 1.0 - md[-1]
        print(f"  {lab:<34}{100*frac:>9.1f}%{err[-1]:>13.4e}{err[-1]/frac:>16.4e}")

    key = [Line2D([], [], color="C0", lw=1.9, label=r"$M_s/M_0 - 1$"),
           Line2D([], [], color="C3", lw=1.5, ls="--", label=r"$M_d/M_0$"),
           Line2D([], [], color="k", lw=0.9, label="exact conservation")]
    fig.tight_layout(rect=[0, 0.15, 1, 1])
    fig.legend(handles=key, loc="lower center", ncol=3, frameon=False,
               fontsize=17, handlelength=1.6, columnspacing=1.6,
               bbox_to_anchor=(0.5, 0.0))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "pics", f"mass_history_R1.{ext}"),
                    dpi=200, bbox_inches="tight")
    print("written: pics/mass_history_R1.{pdf,png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
