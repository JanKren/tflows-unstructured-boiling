#!/usr/bin/env python3
"""Liquid left behind the alpha = 1/2 level, by direction.

The kinematic expansion test (results/caf-kinematic-advection.md) reproduces
the corner bulge only when the transporting velocity drops to zero behind the
computed interface.  This script measures what that leaves behind: the mean
volume fraction inside the bubble, more than two cells behind the local front,
resolved by direction (the cubic invariant s of the position vector: 1 on a
mesh axis, 1/2 on a face diagonal, 1/3 on a body diagonal).

It is applied to the kinematic runs and, unchanged, to field snapshots of the
coupled production runs that were already on disk, so the comparison needs no
new coupled run.

The last column is a volume-conservation estimate, not a fit.  If the vapour
produced per unit solid angle is the same in every direction, then a direction
which keeps a liquid fraction a_d behind the front must carry its front at

    R_d^3 - R_0^3 = (R_axis^3 - R_0^3) (1 - a_axis) / (1 - a_d),

and "explained" is (R_d - R_axis) from that relation over the measured one.

Usage:  python3 analyze_interior_residue.py
"""

import glob
import os

import numpy as np
import pyvista as pv

L = 300.0e-6
R0 = 50.0e-6
KIN = "/home/user/runs/caf-review-2026-11/kinematic"
CASES = [
    ("kinematic, velocity jump, struct 75^3, t = 0.5 ms",
     f"{KIN}/diag075-step-keep", 250, L / 74),
    ("kinematic, solenoidal, struct 75^3, t = 0.5 ms",
     f"{KIN}/diag075-sol-keep", 250, L / 74),
    ("coupled production (computed mdot), struct 75^3, t = 0.6 ms",
     "Data/Aniso/aniso-075-off", 300, L / 74),
    ("coupled production (computed mdot), struct 100^3, t = 0.6 ms",
     "Data/Aniso/aniso-100-off", 300, L / 99),
    ("coupled production (computed mdot), struct 125^3, t = 0.6 ms",
     "Data/Aniso/aniso-125-off", 300, L / 124),
    ("coupled production (computed mdot), struct 150^3, t = 0.6 ms",
     "Data/Aniso/aniso-150-off", 300, L / 149),
    ("kinematic, velocity jump, poly 75^3, t = 0.5 ms",
     f"{KIN}/diagpo075-step-keep", 250, L / 74),
    ("coupled production (computed mdot), poly 75^3, t = 0.6 ms",
     "/home/user/archive/Scriven-Merlin/Scriven-Poly-75", 300, L / 74),
]
CLASSES = (("axis", lambda s: s > 0.90),
           ("face diagonal", lambda s: (s > 0.45) & (s < 0.55)),
           ("body diagonal", lambda s: s < 0.36))


def one(d, ts, h):
    vf = [x for x in glob.glob(os.path.join(d, f"*-ts{ts:06d}.pvtu"))
          if "front" not in x]
    ff = glob.glob(os.path.join(d, f"*front-ts{ts:06d}.pvtu"))
    if not vf or not ff:
        return None
    m = pv.read(vf[0])
    if m.n_cells == 0:
        return None
    key = [n for n in m.cell_data.keys() if "vof" in n.lower()][0]
    a = np.asarray(m.cell_data[key], float)
    c = np.asarray(m.cell_centers().points, float)
    vol = np.asarray(m.compute_cell_sizes(length=False, area=False,
                                          volume=True).cell_data["Volume"], float)
    r = np.linalg.norm(c, axis=1)
    r = np.where(r > 0, r, 1e-30)
    s = ((c / r[:, None])**4).sum(1)

    fm = pv.read(ff[0])
    fa = np.asarray(fm.compute_cell_sizes(length=False, area=True,
                                          volume=False).cell_data["Area"], float)
    fc = np.asarray(fm.cell_centers().points, float)
    fr = np.linalg.norm(fc, axis=1)
    fs = ((fc / fr[:, None])**4).sum(1)

    out = dict(rfront=(fa * fr).sum() / fa.sum(),
               rvol=(3 * (vol * (1 - a)).sum() / (4 * np.pi))**(1 / 3))
    # The solver's own record of the vapour volume, where the run has one
    # (all snapshots here are at dt = 2 us).  Preferred over the sum above,
    # because VTK's volume of a polyhedral cell is not reliable.
    bench = os.path.join(d, "bench-data.dat")
    if os.path.exists(bench):
        b = np.loadtxt(bench)
        i = int(np.argmin(np.abs(b[:, 0] - ts * 2.0e-6)))
        if abs(b[i, 0] - ts * 2.0e-6) < 1.0e-6:
            out["rvol"] = float(b[i, 3])
    for lab, sel in CLASSES:
        rf = (fa[sel(fs)] * fr[sel(fs)]).sum() / fa[sel(fs)].sum()
        # more than two cells behind the local front, outside the initial sphere
        k = sel(s) & (r < rf - 2 * h) & (r > R0 + 2 * h)
        out[lab] = (rf, (vol[k] * a[k]).sum() / vol[k].sum())
    return out


def main():
    print("| case | R_front/R_vol | axis: R [um] | residue | face diag.: R | residue "
          "| body diag.: R | residue | explained |")
    print("|---|---|---|---|---|---|---|---|---|")
    for lab, d, ts, h in CASES:
        o = one(d, ts, h)
        if o is None:
            print(f"| {lab} | (no readable volume snapshot) |")
            continue
        ra, aa = o["axis"]
        rb, ab = o["body diagonal"]
        pred = ((ra**3 - R0**3) * (1 - aa) / (1 - ab) + R0**3)**(1 / 3)
        expl = (f"{100 * (pred - ra) / (rb - ra):.0f} %"
                if abs(rb - ra) > 0.2e-6 else "--")
        print(f"| {lab} | {o['rfront']/o['rvol']:.4f} | " +
              " | ".join(f"{o[c][0]*1e6:.2f} | {o[c][1]:.4f}" for c, _ in CLASSES) +
              f" | {expl} |", flush=True)


if __name__ == "__main__":
    main()
