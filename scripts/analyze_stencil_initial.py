#!/usr/bin/env python3
"""The interface-modified stencil's own anisotropy, measured on the initial condition.

The question is how much four-fold structure the gradient operator manufactures
by itself, as opposed to how much it reads off a state that already carries it.
That needs the operator applied to inputs free of mesh-aligned error, and the
initial condition supplies them without any substitution: at t = 0 the interface
IS the initialised sphere and the temperature IS the analytical profile, so the
state is ideal by construction rather than by replacement.

Two columns are reported per mesh and the comparison between them is the point.

  as initialised   the real alpha and the real T that the solver wrote at ts 0
  ideal / ideal    the exact sphere of the same volume, and an erf profile of
                   the fitted layer thickness about it

They agree to a few tenths of a percentage point, which is what says the
initialisation is clean: substituting analytical inputs for the solver's own
changes almost nothing, because at t = 0 there is nothing to substitute away.
The mid-run version of this test had to replace both inputs on a deformed state,
where the two are co-adapted and cannot be separated without manufacturing a
four-fold mismatch of the substitution's own making.

WHAT THIS DOES NOT SHOW.  The layer spans delta_T/h = 1.7 to 3.3 cells here,
against 5.5 to 6.4 at t = 1 ms, so the operator is measured under a harder
resolution than the one the reported results are produced at.  The amplitudes
below are therefore an upper bound on its contribution at the production state,
not an estimate of it.

Usage:  python3 analyze_stencil_initial.py
"""

import os
import sys

import numpy as np
from scipy.special import erf

from analyze_anisotropy_attribution import (load_lattice, fit_bubble,
                                            modified_gradient, metrics)

ROOT = "/home/user/runs/tflows-vof/Validation/Scriven"
LEVELS = ("75", "100", "125", "150")
VALIDATE_AT = 10        # first saved step; ts 0 carries no computed gradient


def liquid_band(d, fit, use_ideal):
    """Liquid-side interface band, built from the interface the gradient used.

    Mixing the two -- an alpha-derived mask against an ideal-sphere gradient --
    selects mismatched cells and manufactures artefacts, so each column gets
    its own band, exactly as analyze_anisotropy_attribution does.
    """
    is_vap = (fit["r"] < fit["R"]) if use_ideal else (d["vof"] < 0.5)
    nb = np.zeros(d["shape"], dtype=bool)
    for axis in range(3):
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis], hi[axis] = slice(0, -1), slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        nb[lo] |= is_vap[hi]
        nb[hi] |= is_vap[lo]
    return (~is_vap) & nb


def ideal_field(d, fit):
    """Spherically symmetric erf profile about the fitted sphere.

    Any spherically symmetric T(r) has exactly zero azimuthal content, so the
    radial shape does not enter the m = 4 result; it is fitted to the layer
    only so that the magnitude column reads against the same state.
    """
    x = fit["r"] - fit["R"]
    return np.where(x > 0.0,
                    fit["T_sat"] + fit["dT"] * erf(x / fit["delta_T"]),
                    fit["T_sat"])


def main():
    rows = []
    print(f"  {'mesh':>6} {'cells':>11} {'h [um]':>8} {'R [um]':>8} "
          f"{'dT [um]':>8} {'dT/h':>6} {'m4 init':>9} {'m4 ideal':>9} "
          f"{'|gT.n| init':>9} {'ideal':>9} {'corr':>7}")
    for lv in LEVELS:
        path = os.path.join(ROOT, f"Scriven-Struct-{lv}-smallDt",
                            "bubble-ts000000.pvtu")
        # A missing input stops the run rather than dropping a mesh silently.
        assert os.path.exists(path), f"missing initial condition: {path}"
        d = load_lattice(path)
        fit = fit_bubble(d)
        h = d["h"]
        g_ex = 2.0 * fit["dT"] / (np.sqrt(np.pi) * fit["delta_T"])

        out = {}
        for tag, T, ideal in (("init", d["T"], False),
                              ("ideal", ideal_field(d, fit), True)):
            grad, _ = modified_gradient(d, fit, T, ideal)
            amps, a0, _ = metrics(d, fit, grad, liquid_band(d, fit, ideal))
            out[tag] = (100.0 * amps[4], a0 / g_ex, grad)

        # The reimplementation is checked against the gradient T-Flows wrote --
        # but not on this snapshot.  The temperature gradient is formed inside
        # the energy solve, and ts 0 is written before the time loop runs, so
        # its "Temperature Gradients" array holds initialisation residue: it is
        # non-zero in about a tenth of the domain, peaks in a far corner and
        # sits a factor six below the interfacial value.  The check is made on
        # the first saved step instead, where the field is complete and the
        # bubble has grown by one and a half per cent, and it is a check on the
        # OPERATOR, which does not depend on the instant it is applied at.
        dv = load_lattice(os.path.join(ROOT, f"Scriven-Struct-{lv}-smallDt",
                                       f"bubble-ts{VALIDATE_AT:06d}.pvtu"))
        fv = fit_bubble(dv)
        gv, _ = modified_gradient(dv, fv, dv["T"], False)
        m = liquid_band(dv, fv, False)
        a = np.linalg.norm(gv[m], axis=-1)
        b = np.linalg.norm(dv["grad_saved"][m], axis=-1)
        ok = b > 0
        corr = float(np.corrcoef(a[ok], b[ok])[0, 1])

        rows.append((lv, d["vof"].size, h, fit["R"], fit["delta_T"],
                     fit["delta_T"] / h, out["init"][0], out["ideal"][0],
                     out["init"][1], out["ideal"][1], corr))
        print(f"  {lv:>6} {d['vof'].size:>11,} {h*1e6:>8.3f} "
              f"{fit['R']*1e6:>8.2f} {fit['delta_T']*1e6:>8.3f} "
              f"{fit['delta_T']/h:>6.2f} {out['init'][0]:>8.2f}% "
              f"{out['ideal'][0]:>8.2f}% {out['init'][1]:>9.3f} "
              f"{out['ideal'][1]:>9.3f} {corr:>7.4f}")

    bad = [r[0] for r in rows if r[10] < 0.99]
    if bad:
        print(f"\n  SUSPECT: reimplementation did not validate on {bad}; "
              f"treat the amplitudes as unverified.")
        return 1

    m4i = [r[7] for r in rows]
    print(f"\n  The operator's own m = 4 falls {m4i[0]:.2f}% to {m4i[-1]:.2f}% "
          f"as delta_T/h rises {rows[0][5]:.1f} to {rows[-1][5]:.1f}.")
    print(f"  As-initialised and ideal agree to "
          f"{max(abs(r[6]-r[7]) for r in rows):.2f} percentage points, so the "
          f"substitution changes almost nothing at t = 0.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
