#!/usr/bin/env python3
"""Static droplet over several viscous times, and at a high density ratio.

    python3 analyze_static_drop_long.py            # table + figure
    python3 analyze_static_drop_long.py --no-plot

Answers the referee's static-drop comment (CAF-D-26-01303, Reviewer 2): the
published record of 0.12 s is 0.75 of one viscous time rho R^2 / mu = 0.16 s,
and both phases carry the same density.  The runs read here extend the R/h =
12.8 level to 0.5 s (3.1 viscous times, 242 capillary times sqrt(rho R^3 /
sigma)) on both mesh families, and repeat it for a vapour-like sphere in a
liquid, rho_l/rho_v = 1000 and mu_l/mu_v = 100.

Traces (written every step by User_Mod/End_Of_Time_Step.f90):

  spurious.dat      time | volume of fun | u_rms | u_max | p_max-p_min | p_in-p_out
                    the six columns of the published runs
  spurious_ext.dat  time | u_max | u_max interface cells | u_max band |r-R| < R/4
                    | u_max > 0.2 from every wall | u_max sphere phase
                    | u_max outer phase | x y z of the global maximum
                    | kinetic energy | u_rms (volume weighted) | centroid x y z
                    | min(fun) | max(fun) | p_in-p_out

Ca = mu_l max|u| / sigma throughout, mu_l = 1 being the viscosity of the outer
(liquid) phase in every case.  The published definition takes the maximum over
the whole domain.  On the polyhedral dual that maximum sits in a box corner,
0.01 from a wall and 1.7 from the droplet, so the interface-cell maximum is
reported beside it wherever the extended trace exists.
"""

import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = "/home/user/runs/caf-review-2026-11/static-drop/traces"
OLD = "/home/user/runs/Droplet3D"
ARC = "/home/user/archive/Merlin7-scratch-20260930/Droplet3D"
OCT = ("/home/user/runs/FlowBoiling/static_sphere_merlin_20261008_v1/"
       "full/cases")

SIGMA, MU_L, R0 = 15.0e3, 1.0, 0.4
EXACT_DP = 2.0 * SIGMA / R0              # 75 000 Pa
T_MU = 1.0 * R0**2 / MU_L                # rho_l R^2 / mu_l = 0.16 s
T_SIG = np.sqrt(1.0 * R0**3 / SIGMA)     # sqrt(rho_l R^3 / sigma) = 2.07e-3 s

# label, trace directory, mesh, smoothing cycles, property set, solver
CASES = [
    ("published hex 12",   f"{OLD}/hex64_result",       "hex",  12, "matched", "paper"),
    ("published poly 4",   f"{OLD}/poly64_result",      "poly",  4, "matched", "paper"),
    ("Oct hex 4",          f"{OCT}/hex64_cicsam_csf",   "hex",   4, "matched", "T-Flows-jan"),
    ("Oct poly 4",         f"{OCT}/poly64_cicsam_csf",  "poly",  4, "matched", "T-Flows-jan"),
    ("m_hex_s12",          f"{NEW}/m_hex_s12",          "hex",  12, "matched", "review"),
    ("m_hex_s4",           f"{NEW}/m_hex_s4",           "hex",   4, "matched", "review"),
    ("m_poly_s4",          f"{NEW}/m_poly_s4",          "poly",  4, "matched", "review"),
    ("m_poly_s12",         f"{NEW}/m_poly_s12",         "poly", 12, "matched", "review"),
    ("h_hex_s12",          f"{NEW}/h_hex_s12",          "hex",  12, "high",    "review"),
    ("h_hex_s4",           f"{NEW}/h_hex_s4",           "hex",   4, "high",    "review"),
    ("h_poly_s4",          f"{NEW}/h_poly_s4",          "poly",  4, "high",    "review"),
    # coarse preview, R/h = 6.4, not a level the manuscript tabulates
    ("c32m_hex_s4",        f"{NEW}/c32m_hex_s4",        "hex",   4, "matched", "review"),
    ("c32m_poly_s4",       f"{NEW}/c32m_poly_s4",       "poly",  4, "matched", "review"),
    ("c32h_hex_s4",        f"{NEW}/c32h_hex_s4",        "hex",   4, "high",    "review"),
    ("c32h_poly_s4",       f"{NEW}/c32h_poly_s4",       "poly",  4, "high",    "review"),
]

# exact-curvature runs (EXACT_CURVATURE_TEST 1: kappa = 2/R in every cell) and
# the run with the solver's own curvature each is compared with
EXACT = [("ex1_m_hex_s4",  "m_hex_s4",  "hex 64^3, matched"),
         ("ex1_m_poly_s4", "m_poly_s4", "poly 64, matched"),
         ("ex1_h_hex_s4",  "h_hex_s4",  "hex 64^3, 1000:100"),
         ("ex1_h_poly_s4", "h_poly_s4", "poly 64, 1000:100")]
PILOT = "/home/user/runs/caf-review-2026-11/static-drop/traces_pilot"

# windows, in viscous times: [centre - half, centre + half]
MARKS = [0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
HALF = 0.05


def load(path):
    """Both traces of one run; None where a file is missing or too short."""
    out = {}
    for key, name in (("base", "spurious.dat"), ("ext", "spurious_ext.dat")):
        f = os.path.join(path, name)
        if not os.path.exists(f):
            out[key] = None
            continue
        try:
            d = np.loadtxt(f)
        except ValueError:                       # a row cut by a live writer
            d = np.genfromtxt(f, invalid_raise=False)
        out[key] = d if d.ndim == 2 and len(d) > 100 else None
    return out


def screen(t):
    """A spliced or restarted trace has a non-monotone time column."""
    dt = np.diff(t)
    return int((dt <= 0).sum()), float(np.median(dt))


def window(t, y, centre, half=HALF):
    m = (t >= (centre - half) * T_MU) & (t < (centre + half) * T_MU)
    return (y[m].mean(), y[m].std()) if m.sum() > 10 else (np.nan, np.nan)


def trend(t, y, t0):
    """e-folding time of y for t > t0, from a straight-line fit to log y on
    window means of 0.1 t_mu.  Positive: growth.  inf: no measurable change."""
    edges = np.arange(t0, t[-1] + 1e-12, 0.1 * T_MU)
    if len(edges) < 4:
        return np.nan
    tc, yc = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b)
        if m.sum() > 10:
            tc.append(0.5 * (a + b))
            yc.append(y[m].mean())
    if len(tc) < 3:
        return np.nan
    slope = np.polyfit(tc, np.log(yc), 1)[0]
    return np.inf if abs(slope) < 1e-12 else 1.0 / slope / T_MU


def report():
    rows = []
    print(f"exact jump {EXACT_DP:.0f} Pa, viscous time {T_MU:.3f} s, capillary "
          f"time {T_SIG:.3e} s, published record 0.12 s = {0.12 / T_MU:.2f} "
          f"viscous = {0.12 / T_SIG:.0f} capillary times\n")
    for lab, path, mesh, cyc, props, solver in CASES:
        d = load(path)
        b, e = d["base"], d["ext"]
        if b is None:
            print(f"{lab:<18} trace missing ({path})")
            continue
        t = b[:, 0]
        nonmono, dt = screen(t)
        ca = MU_L * b[:, 3] / SIGMA
        dp = 100.0 * (b[:, 5] - EXACT_DP) / EXACT_DP
        ca_if = MU_L * e[:, 2] / SIGMA if e is not None else None
        te = e[:, 0] if e is not None else None
        rows.append((lab, mesh, cyc, props, t, ca, te, ca_if, e))

        print(f"{lab:<18} {mesh:<4} {cyc:>2} cycles, {props}, {solver}: "
              f"{len(t)} rows, dt {dt:.4e}, t_end {t[-1]:.4f} s = "
              f"{t[-1] / T_MU:.2f} t_mu = {t[-1] / T_SIG:.0f} t_sigma"
              + (f"  ** {nonmono} non-monotone time rows **" if nonmono else ""))

        # the published statistic: second half of the first 0.12 s
        m = (t >= 0.06) & (t <= 0.12)
        print(f"    0.06-0.12 s (the published window): Ca {ca[m].mean():.3e}"
              f"   dp {dp[m].mean():+.2f} +- {dp[m].std():.2f} %"
              + (f"   Ca_interface {ca_if[(te >= 0.06) & (te <= 0.12)].mean():.3e}"
                 if ca_if is not None else ""))

        line = "    t/t_mu      " + "".join(f"{c:>10.2f}" for c in MARKS)
        print(line)
        print("    Ca global   " + "".join(
            f"{window(t, ca, c)[0]:>10.2e}" for c in MARKS))
        if ca_if is not None:
            print("    Ca interface" + "".join(
                f"{window(te, ca_if, c)[0]:>10.2e}" for c in MARKS))
            print("    Ca core     " + "".join(
                f"{window(te, MU_L * e[:, 4] / SIGMA, c)[0]:>10.2e}"
                for c in MARKS))
            print("    u_max sphere" + "".join(
                f"{window(te, e[:, 5], c)[0]:>10.2e}" for c in MARKS))
            print("    u_max outer " + "".join(
                f"{window(te, e[:, 6], c)[0]:>10.2e}" for c in MARKS))
            print("    E_kin       " + "".join(
                f"{window(te, e[:, 10], c)[0]:>10.2e}" for c in MARKS))
        print("    u_rms       " + "".join(
            f"{window(t, b[:, 2], c)[0]:>10.2e}" for c in MARKS))
        print("    dp err [%]  " + "".join(
            f"{window(t, dp, c)[0]:>+10.2f}" for c in MARKS))
        print("    dp std [%]  " + "".join(
            f"{window(t, dp, c)[1]:>10.2f}" for c in MARKS))

        drift = (b[-1, 1] - b[0, 1]) / b[0, 1]
        extra = ""
        if e is not None:
            cen = e[-1, 12:15] - 1.0
            xm = e[-1, 7:10]
            dwall = min(xm.min(), (2.0 - xm).min())
            extra = (f", centroid moved {np.linalg.norm(cen) / R0:.4f} R, "
                     f"fun in [{e[:, 15].min():.1e}, 1{e[:, 16].max() - 1:+.1e}], "
                     f"global maximum at the end {np.linalg.norm(xm - 1.0):.2f} "
                     f"from the centre and {dwall:.3f} from a wall")
        print(f"    volume drift {drift:+.2e}{extra}")
        if t[-1] > 1.3 * T_MU:
            tau = trend(t, ca, 1.0 * T_MU)
            print(f"    Ca global, t > t_mu: e-folding time {tau:+.1f} t_mu "
                  "(negative = decay)")
            if ca_if is not None:
                tau = trend(te, ca_if, 1.0 * T_MU)
                print(f"    Ca interface, t > t_mu: e-folding time {tau:+.1f} "
                      "t_mu")
                tau = trend(te, e[:, 10], 1.0 * T_MU)
                print(f"    E_kin, t > t_mu: e-folding time {tau:+.1f} t_mu")
        print()
    return rows


def report_exact():
    """Francois-type test: with the exact curvature in every cell the surface
    force is a pure gradient, and a balanced-force discretisation leaves no
    flow.  Compared at equal step counts with the solver's own curvature."""
    marks = [1, 10, 100, 400, 1000, 2000, 3000]
    print("exact curvature (kappa = 2/R everywhere) against the solver's "
          "curvature, u_max [m/s] in interface cells / anywhere, by step")
    print("    step" + "".join(f"{m:>22d}" for m in marks))
    for ex, base, lab in EXACT:
        for tag, path in (("exact", f"{NEW}/{ex}"),
                          ("solver", f"{NEW}/{base}"),
                          ("solver (pilot)", f"{PILOT}/{base}_p32")):
            e = load(path)["ext"]
            if e is None:
                continue
            row = "".join(
                f"{e[m - 1, 2]:>11.3e}/{e[m - 1, 1]:<10.3e}" if len(e) >= m
                else f"{'-':>22}" for m in marks)
            jump = 100.0 * (e[-1, 17] - EXACT_DP) / EXACT_DP
            print(f"  {lab:<20} {tag:<15}{row}   dp err at the last row "
                  f"{jump:+.4f} %")
            if tag == "solver":
                break
    print()


def smooth(t, y, width):
    """Mean over consecutive windows of the given width in time."""
    edges = np.arange(t[0], t[-1] + width, width)
    idx = np.searchsorted(t, edges)
    tc, yc = [], []
    for a, b in zip(idx[:-1], idx[1:]):
        if b > a:
            tc.append(t[a:b].mean())
            yc.append(y[a:b].mean())
    return np.array(tc), np.array(yc)


def plot(rows, prefixes=("m_", "h_"), name="static_drop_ca_history"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif",
                         "axes.grid": True, "grid.alpha": 0.3,
                         "savefig.dpi": 300})
    style = {("hex", 12): ("C0", "hexahedral, 12 cycles"),
             ("hex", 4): ("C2", "hexahedral, 4 cycles"),
             ("poly", 4): ("C3", "polyhedral, 4 cycles"),
             ("poly", 12): ("C1", "polyhedral, 12 cycles")}
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True)
    titles = {"matched": r"(a) $\rho_l/\rho_v = 1$, $\mu_l/\mu_v = 1$",
              "high": r"(b) $\rho_l/\rho_v = 1000$, $\mu_l/\mu_v = 100$"}
    for ax, props in zip(axes, ("matched", "high")):
        for lab, mesh, cyc, p, t, ca, te, ca_if, e in rows:
            if p != props or not lab.startswith(prefixes):
                continue
            col, leg = style[(mesh, cyc)]
            if ca_if is not None:
                x, y = smooth(te, ca_if, 0.02 * T_MU)
                ax.semilogy(x / T_MU, y, color=col, lw=1.3, label=leg)
                x, y = smooth(t, ca, 0.02 * T_MU)
                ax.semilogy(x / T_MU, y, color=col, lw=0.8, ls=":")
            else:
                x, y = smooth(t, ca, 0.02 * T_MU)
                ax.semilogy(x / T_MU, y, color=col, lw=1.3, label=leg)
        ax.axvline(0.12 / T_MU, color="0.4", lw=0.8, ls="--")
        ax.set_xlabel(r"$t\,\mu_l/(\rho_l R^2)$")
        ax.set_title(titles[props], fontsize=10, loc="left")
    axes[0].set_ylabel(r"$\mathrm{Ca} = \mu_l \max|\mathbf{u}|/\sigma$")
    axes[0].legend(fontsize=8, loc="best")
    fig.tight_layout()
    out = os.path.join(HERE, "pics", name)
    fig.savefig(out + ".pdf")
    fig.savefig(out + ".png")
    print("wrote", out + ".pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    r = report()
    report_exact()
    if not a.no_plot:
        # the figure of the letter: R/h = 12.8, only once those runs exist
        if any(x[0].startswith(("m_", "h_")) for x in r):
            plot(r)
        if any(x[0].startswith("c32") for x in r):
            plot(r, ("c32m_", "c32h_"), "static_drop_ca_history_coarse")
