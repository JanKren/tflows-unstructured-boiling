#!/usr/bin/env python3
"""
T2.3 -- evaporation rate and interfacial heat flux vs the Scriven solution
        (Reviewer 2, point 4).

The residual radius error mixes two distinct failures, and the referee's point is
that they should be separated.  This script does that by comparing the evaporation
rate against the analytical rate evaluated at TWO radii:

  ratio_local = mdot_sim(t) / mdot_Scriven(R_sim(t))
      the rate the analytical solution would produce for a bubble of the size the
      simulation actually has.  Deviation from 1 is a failure of the interfacial
      CLOSURE -- the gradient, the extrapolation, the boundary-layer resolution --
      independent of how the bubble got to that size.

  ratio_traj  = R_sim(t) / R_Scriven(t)
      the accumulated trajectory error.  A run can have a perfect closure and still
      sit off the analytical curve because of its history.

A run whose closure is exact but whose trajectory is not has an advection or
history problem; a run whose closure is off has a gradient or resolution problem.

Analytical relations (virtual origin t0 = R0^2/(4 beta^2 alpha_l)):

    R(t)      = 2 beta sqrt(alpha_l (t + t0))
    dm/dt     = rho_v 4 pi R^2 dR/dt = 2 pi rho_v R^3 / (t + t0)
    Q_int     = h_lv dm/dt                      (energy jump condition)

The solver's `Vof MassTransfer [kg/m^3/s]` is MISLABELLED: Save_Vtu_Fields writes
Vof % m_dot raw, and Mass_Transfer_Estimate sets it to (q_1-q_0)/h_lv with q in
watts, so it is kg/s PER CELL, already integrated.  The domain total is a plain
sum; weighting by cell volume under-predicts by ~16 orders of magnitude.

Usage:
  python3 analyze_scriven_evap_rate.py [--fig results/figs/T2.3_evap_rate.png]
"""

import argparse
import glob
import os

import numpy as np
import pyvista as pv

RHO_V, RHO_L, CP_L, K_L = 0.597, 958.4, 4216.0, 0.679
H_LV = 2.26e6
BETA = 4.06022
R0 = 5.0e-5
ALPHA_L = K_L / (RHO_L * CP_L)
T0 = R0 ** 2 / (4.0 * BETA ** 2 * ALPHA_L)


def r_scriven(t):
    return 2.0 * BETA * np.sqrt(ALPHA_L * (t + T0))


def mdot_scriven_at_R(R, t):
    """Analytical evaporation rate for a bubble of radius R at time t."""
    return 2.0 * np.pi * RHO_V * R ** 3 / (t + T0)


def scan(pattern, dt, label, stride=1):
    files = sorted((f for f in glob.glob(pattern)
                    if "bnd" not in os.path.basename(f)
                    and "front" not in os.path.basename(f)),
                   key=lambda f: int(f.split("-ts")[1].split(".")[0]))
    files = files[::stride]
    rows, bad = [], 0
    for fn in files:
        ts = int(fn.split("-ts")[1].split(".")[0])
        if ts == 0:
            continue
        try:
            m = pv.read(fn)
            _ = m.cell_data["Vof Sharp [1]"]      # truncated files fail here
        except Exception:
            bad += 1
            continue
        cd = m.cell_data
        vof = np.asarray(cd["Vof Sharp [1]"], float)
        vol = np.asarray(cd["Grid Cell Volume [m^3]"], float)
        md = np.asarray(cd["Vof MassTransfer [kg/m^3/s]"], float)
        v_vap = ((1.0 - vof) * vol).sum()
        R = (3.0 * v_vap / (4.0 * np.pi)) ** (1.0 / 3.0)
        rows.append((ts * dt, R, md.sum(), v_vap))
    if bad:
        print(f"  [{label}: skipped {bad} unreadable/truncated snapshot(s)]")
    return label, np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", default=None)
    ap.add_argument("--stride", type=int, default=6,
                    help="subsample the dense polyhedral output")
    a = ap.parse_args()

    S = "/home/user/runs/tflows-vof/Validation/Scriven"
    cases = []
    for lvl in (75, 125):
        p = f"Data/Aniso/aniso-{lvl:03d}-off/*-ts*.pvtu"
        if glob.glob(p):
            cases.append(scan(p, 2.0e-6, f"structured {lvl}"))
    # Scriven-Poly-75 in the Tests tree mixes two runs (radius drops
    # 112.45 -> 66.98 um at ts300, mtime break at ts190) and is excluded;
    # Poly-100 is substituted. Poly-100/125 verified monotonic in R.
    for lvl in (100, 125):
        p = f"{S}/Scriven-Poly-{lvl}/{'Scriven-Poly-%d_dual' % lvl}-ts*.pvtu"
        if glob.glob(p):
            cases.append(scan(p, 2.0e-6, f"polyhedral {lvl}", stride=a.stride))

    print(f"\n{'='*94}")
    print("  T2.3  evaporation rate vs the Scriven solution")
    print(f"{'='*94}")
    print(f"  t0 = {T0*1e6:.2f} us, alpha_l = {ALPHA_L:.6e} m2/s, "
          f"rho_v = {RHO_V} kg/m3, h_lv = {H_LV:.3e} J/kg\n")
    print("  %-16s %8s %10s %10s %9s %13s %13s %9s %9s" %
          ("case", "t [ms]", "R_sim[um]", "R_Scr[um]", "R ratio",
           "mdot_sim[kg/s]", "mdot_Scr(R)", "LOCAL", "Q_int[W]"))
    print("  " + "-" * 108)
    store = {}
    for label, arr in cases:
        if arr.size == 0:
            continue
        t, R, md, _ = arr.T
        rs = r_scriven(t)
        ms = mdot_scriven_at_R(R, t)
        store[label] = (t, R / rs, md / ms)
        for i in range(len(t)):
            print("  %-16s %8.4f %10.4f %10.4f %9.5f %13.6e %13.6e "
                  "%9.5f %9.4e" %
                  (label if i == 0 else "", t[i] * 1e3, R[i] * 1e6,
                   rs[i] * 1e6, R[i] / rs[i], md[i], ms[i], md[i] / ms[i],
                   H_LV * md[i]))
        print()

    print("  LOCAL = mdot_sim / mdot_Scriven(R_sim): the closure error at the")
    print("  bubble's actual size. R ratio: the accumulated trajectory error.\n")
    print("  %-16s %26s %26s" % ("case", "LOCAL closure (mean, last)",
                                 "trajectory R/R_Scr (last)"))
    print("  " + "-" * 72)
    for label, (t, rr, lr) in store.items():
        print("  %-16s %13.5f %12.5f %20.5f" %
              (label, lr.mean(), lr[-1], rr[-1]))

    if a.fig and store:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10, 'font.family': 'serif',
                             'savefig.dpi': 300, 'axes.grid': False,
                             'lines.linewidth': 1.4})
        fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.4))
        mk = {'structured 75': ('C0', 'o-'), 'structured 125': ('C1', 's-'),
              'polyhedral 75': ('C2', '^--'), 'polyhedral 125': ('C3', 'v--')}
        for label, (t, rr, lr) in store.items():
            c, s = mk.get(label, ('k', '-'))
            ax[0].plot(t * 1e3, lr, s, color=c, label=label, ms=4)
            ax[1].plot(t * 1e3, rr, s, color=c, label=label, ms=4)
        for k, ttl in ((0, r"local closure  $\dot m_{sim}/\dot m_{Scriven}(R_{sim})$"),
                       (1, r"trajectory  $R_{sim}/R_{Scriven}(t)$")):
            ax[k].axhline(1.0, color='k', lw=0.8, ls=':')
            ax[k].set_xlabel("t [ms]"); ax[k].set_title(ttl, fontsize=9)
        ax[0].legend(fontsize=7.5)
        fig.tight_layout()
        os.makedirs(os.path.dirname(a.fig), exist_ok=True)
        fig.savefig(a.fig, bbox_inches='tight')
        print(f"\n  figure -> {a.fig}")


if __name__ == "__main__":
    main()
