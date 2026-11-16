"""
poster_figure_style.py
-----------------------
Re-export the poster's embedded matplotlib figures at A0-legible sizes.

Why: the figures currently dropped into ICBCHT12_poster_*.pptx were saved at
paper font sizes. On a printed A0 viewed from ~1 m, the in-figure text (tick
labels, legends, axis labels, and especially the polar-plot titles such as
"min/max ratio = 0.597") is unreadable. The body text and equations on the
poster are fine — only the rasterised plots need re-exporting.

How to use: import this in your existing analyze_*.py scripts BEFORE creating
the figure, call apply_poster_style(), and size each figure to the box it
occupies on the poster (see SIZES below). Re-export, then swap the PNGs in the
deck. Source data (VTU fields, *.out) is NOT in this folder, so this has to be
run where the simulation output lives.

Rule of thumb: save the figure at the SAME physical size it occupies on the
poster (in inches) and set fonts to ~22-26 pt. Then there is no up-scaling, and
what you see in the figure is exactly what prints.
"""

import matplotlib as mpl


def apply_poster_style(base=24):
    """Large-font defaults tuned for A0 viewing distance. `base` is the body
    font size in points at final printed size."""
    mpl.rcParams.update({
        "font.size":        base,
        "axes.titlesize":   base + 2,
        "axes.labelsize":   base,
        "xtick.labelsize":  base - 4,
        "ytick.labelsize":  base - 4,
        "legend.fontsize":  base - 4,
        "figure.titlesize": base + 4,
        # thicker strokes so curves read from a distance
        "lines.linewidth":  3.0,
        "lines.markersize": 10,
        "axes.linewidth":   1.6,
        "xtick.major.width": 1.6,
        "ytick.major.width": 1.6,
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "legend.frameon":   False,
        "savefig.dpi":      300,
        "savefig.bbox":     "tight",
        "savefig.pad_inches": 0.05,
        "font.family":      "sans-serif",
        "mathtext.fontset": "dejavusans",
    })


# Target figure sizes (width, height) in INCHES — match the poster box so the
# PNG is printed 1:1 with no up-scaling. Measured from v5 (A0 portrait, 33x47").
SIZES = {
    "stefan_sucking":   (9.0, 6.0),   # col 1, 1-D benchmarks (2x2 panel)
    "scriven_sensitivity": (9.0, 5.0),  # col 2, 4-panel convergence/overshoot
    "scriven_convergence": (9.0, 4.0),
    "annular_quantitative": (9.0, 4.0),  # col 2, (a) MTR vs film, (b) by phase
    "grad_polar":       (9.0, 5.0),   # col 3 + header, polar |grad T| + Fourier
    "headline_polar":   (12.5, 3.2),  # top-right "Headline Result" trio
}

# Per-figure notes:
#   headline_polar / grad_polar:
#     - Drop the tiny per-axes titles ("Structured 75^3 / mean |grad T| = ... /
#       min/max ratio = ...") from inside the plot. Put those numbers in the
#       poster caption instead, or as ONE large annotation per panel (>=22 pt).
#     - Increase polar tick label size and the radial-grid label size.
#     - Make the red binned-mean line width >= 4; thin grey scatter is fine.
#   stefan_sucking / scriven_*:
#     - subplot labels (a,b,c,d): >=26 pt, bold, placed inside each axes corner.
#     - legends: one per panel, >=20 pt; shorten entries ("poly 125^3").
#   annular_quantitative:
#     - the "r = -0.62" / bar value labels ("4.1x") should be >=24 pt.

# Minimal example -------------------------------------------------------------
if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    apply_poster_style(base=24)
    fig, ax = plt.subplots(figsize=SIZES["scriven_convergence"])
    h = np.array([4.0, 3.0, 2.4, 2.0])
    err = np.array([7.1, 5.6, 4.6, 3.9])
    ax.loglog(h, err, "o-", color="#6E14DC", label="polyhedral")
    ax.loglog(h, 1.2 * h**1.3, "--", color="0.5", label=r"$p=1.3$")
    ax.set_xlabel(r"cell size $h$ [$\mu$m]")
    ax.set_ylabel("radius error [%]")
    ax.legend()
    fig.savefig("scriven_convergence_poster.png")
    print("wrote scriven_convergence_poster.png at", SIZES["scriven_convergence"], "in")
