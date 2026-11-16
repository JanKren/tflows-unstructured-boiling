#!/usr/bin/env python3
"""Figure 3 (pics/stefan_sucking_interface_R1.{png,pdf}), revision R1.

The Stefan panels are drawn from the runs with the energy equation converged:
no under-relaxation of the energy equation and first-order upwind convection
for it, the two control settings that remove the two causes of the -3.5 %
energy imbalance of the submitted run (results/caf-energy-mass-solver.md, B).
Two time steps are shown, 5 ms (the step of the paper) and 2.5 ms, because the
remaining position error is first order in the interface Courant number.

The Sucking panels are unchanged: the record of the submitted run.

Inputs
  Stefan     ~/runs/caf-review-2026-11/stefan/runs/Pfix_N100, Pfix_dt2.5
  Sucking    Unstructured-Boiling-Paper-Repo/Sucking/interface_position.{num,exa}
  analytical Unstructured-Boiling-Paper-Repo/Stefan/interface_position.dat
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_stefan_sucking as base

REPO = "/home/user/src/Unstructured-Boiling-Paper-Repo"
RUNS = "/home/user/runs/caf-review-2026-11/stefan/runs"
OUT = "/home/user/papers/UnstructuredBoiling/pics/stefan_sucking_interface_R1"

base.SUCKING_DIR = os.path.join(REPO, "Sucking")


def stefan(run):
    s = np.loadtxt(os.path.join(RUNS, run, "stefans_solution.dat"))
    return s[:, 0], s[:, 1]


def main():
    ana = np.loadtxt(os.path.join(REPO, "Stefan", "interface_position.dat"),
                     comments='#')
    t_suck, x_suck, ta_suck, xa_suck = base.load_sucking()

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(ana[:, 0], ana[:, 1] * 1e3, 'k-', linewidth=2, label='Analytical')
    for run, lab, sty in (("Pfix_N100", r'$\Delta t = 5$ ms', 'C0--'),
                          ("Pfix_dt2.5", r'$\Delta t = 2.5$ ms', 'C1:')):
        t, x = stefan(run)
        ax.plot(t, x * 1e3, sty, linewidth=1.5, label=lab)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Interface position [mm]')
    ax.set_title('(a) Stefan problem')
    ax.legend(loc='lower right')

    ax = axes[0, 1]
    for run, lab, sty in (("Pfix_N100", r'$\Delta t = 5$ ms', 'C0-'),
                          ("Pfix_dt2.5", r'$\Delta t = 2.5$ ms', 'C1-')):
        t, x = stefan(run)
        xa = np.interp(t, ana[:, 0], ana[:, 1])
        e = 100.0 * (x - xa) / xa
        ax.plot(t, e, sty, linewidth=1.5, label=lab)
        print(f"Stefan {run:11s}: final {e[-1]:+.2f} %, mean |e| "
              f"{np.mean(np.abs(e)):.2f} %, max |e| {np.max(np.abs(e)):.2f} %")
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Relative error [%]')
    ax.set_title('(b) Stefan position error')
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.legend(loc='lower right')

    ax = axes[1, 0]
    ax.plot(ta_suck, xa_suck * 1e3, 'k-', linewidth=2, label='Analytical')
    step = max(1, len(t_suck) // 500)
    ax.plot(t_suck[::step], x_suck[::step] * 1e3, 'C0--', linewidth=1.5,
            label='T-Flows')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Interface position [mm]')
    ax.set_title('(c) Sucking problem')
    ax.legend(loc='lower right')

    ax = axes[1, 1]
    xa = np.interp(t_suck, ta_suck, xa_suck)
    e = 100.0 * (x_suck - xa) / xa
    ax.plot(t_suck[::step], e[::step], 'C0-', linewidth=1.5)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Relative error [%]')
    ax.set_title('(d) Sucking position error')
    ax.axhline(y=0, color='k', linewidth=0.5)
    print(f"Sucking: mean error {np.mean(e):+.2f} %, mean |e| "
          f"{np.mean(np.abs(e)):.2f} %, max |e| {np.max(np.abs(e)):.2f} %")

    plt.tight_layout()
    for ext in ('.png', '.pdf'):
        plt.savefig(OUT + ext, bbox_inches='tight')
    print("wrote", OUT + ".{png,pdf}")


if __name__ == "__main__":
    main()
