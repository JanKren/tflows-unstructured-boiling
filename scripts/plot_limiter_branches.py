#!/usr/bin/env python3
"""The two branches of the CICSAM blend against the normalised donor value.

    python3 plot_limiter_branches.py

Hyper-C takes the CBC bound min(1, alpha_D_tilde/Co), which saturates once
alpha_D_tilde > Co: on those faces the face value is the acceptor's whatever the
stencil reads, so the orientation dependence of alpha_D_tilde does not reach
beta_f.  ULTIMATE-QUICKEST is nearly linear in alpha_D_tilde and passes it on.
That is the mechanism behind the ordering of Table `tab:cicsam_beta`, 3.21% of
four-fold content at gamma_f = 1 rising to 9.71% at gamma_f = 0.
Self-contained (numpy + matplotlib); Co = 0.2 as in that table.
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CO = 0.2

plt.rcParams.update({"font.size": 10, "font.family": "serif",
                     "axes.grid": True, "grid.alpha": 0.3})


def main():
    fig, ax = plt.subplots(figsize=(5.0, 3.9))
    ad = np.linspace(0, 1, 400)
    hyperc = np.minimum(1.0, ad / CO)
    uq = np.minimum(CO * ad + (1 - CO) * (6 * ad + 3) / 8, hyperc)
    ax.plot(ad, hyperc, color="C0", lw=2.0,
            label=r"Hyper-C: $\min(1,\ \tilde\alpha_D/\mathrm{Co})$")
    ax.plot(ad, uq, color="C3", lw=2.0, label="ULTIMATE-QUICKEST")
    ax.axvspan(CO, 1.0, color="C0", alpha=0.10)
    ax.axvline(CO, color="k", ls="--", lw=1.0, alpha=0.7)
    ax.annotate(r"$\mathrm{Co}$", xy=(CO, 0.05), xytext=(CO + 0.03, 0.05), fontsize=9)
    ax.annotate("Hyper-C saturated:\nface value = acceptor value",
                xy=(0.60, 0.30), ha="center", va="center", fontsize=9, color="C0")
    ax.set_xlabel(r"$\tilde\alpha_D$")
    ax.set_ylabel(r"$\tilde\alpha_f$")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.08)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "pics", f"limiter_branches.{ext}"),
                    dpi=200, bbox_inches="tight")
    print("written: pics/limiter_branches.{pdf,png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
