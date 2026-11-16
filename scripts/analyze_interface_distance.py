#!/usr/bin/env python3
"""How close does the interface come to a cell centre, and how large does the
interface coefficient a_f = k A_int,f / d_f become relative to an ordinary
conduction coefficient?

    python3 analyze_interface_distance.py <structured Scriven .pvtu> [...]

On every face whose two cells straddle alpha = 1/2 the solver cuts the
cell-centre connection at

    w = |alpha_P - 1/2| / |alpha_P - alpha_Q|          (Eq. front_interpolation)

so cell P sees the interface at |r| = w h and cell Q at (1 - w) h, and the
distance in the coefficient is d_f = |r| max(|r.n|/|r|, c_min), c_min = 0.05.
On a Cartesian mesh the connection is a mesh axis e_i, so |r.n|/|r| = |n_i|,
and the area projected on the face is A |n_i|.  Relative to the conduction
coefficient k A / h of an uncut face of the same cell,

    a_f / (k A / h) = |n_i| / ( w max(|n_i|, c_min) ),

which is 1 / w wherever |n_i| > c_min.  The script reports the distribution
of w (both sides of every cut face), of d_f / h and of that ratio, and the
share of the cell's matrix diagonal that a_f takes at the given time step.
"""
import sys

import numpy as np
import pyvista as pv

K_L, K_V = 0.679, 0.025
RCP_L, RCP_V = 958.4 * 4216.0, 0.597 * 2030.0
C_MIN = 0.05


def main():
    for fn in sys.argv[1:]:
        m = pv.read(fn)
        if m.n_cells == 0:
            sys.exit(f"{fn}: no cells behind the header")
        xc = m.cell_centers().points
        al = np.asarray(m.cell_data["Vof Sharp [1]"], dtype=np.float64)
        # structured block: sort into a 3-D array
        ux, uy, uz = (np.unique(np.round(xc[:, i], 12)) for i in range(3))
        n = len(ux)
        assert n * len(uy) * len(uz) == m.n_cells, "not a structured block"
        h = ux[1] - ux[0]
        ix = np.searchsorted(ux, np.round(xc[:, 0], 12))
        iy = np.searchsorted(uy, np.round(xc[:, 1], 12))
        iz = np.searchsorted(uz, np.round(xc[:, 2], 12))
        a = np.empty((n, len(uy), len(uz)))
        a[ix, iy, iz] = al
        g = np.stack(np.gradient(a, h), axis=-1)
        nrm = g / (np.linalg.norm(g, axis=-1, keepdims=True) + 1e-300)

        w_all, d_all, ratio_all, kc = [], [], [], []
        for ax in range(3):
            sl1 = [slice(None)] * 3
            sl2 = [slice(None)] * 3
            sl1[ax] = slice(0, -1)
            sl2[ax] = slice(1, None)
            a1, a2 = a[tuple(sl1)], a[tuple(sl2)]
            cut = (a1 - 0.5) * (a2 - 0.5) < 0
            w = np.abs(a1[cut] - 0.5) / np.abs(a1[cut] - a2[cut])
            nf = 0.5 * (nrm[tuple(sl1)][cut] + nrm[tuple(sl2)][cut])
            nf /= np.linalg.norm(nf, axis=-1, keepdims=True)
            ni = np.abs(nf[:, ax])
            for ww, aa in ((w, a1[cut]), (1.0 - w, a2[cut])):
                d = ww * h * np.maximum(ni, C_MIN)
                w_all.append(ww)
                d_all.append(d / h)
                ratio_all.append(ni * h / np.maximum(d, 1e-300))
                kc.append(np.where(aa > 0.5, 1, 0))
        w = np.concatenate(w_all)
        d = np.concatenate(d_all)
        r = np.concatenate(ratio_all)
        liq = np.concatenate(kc).astype(bool)

        print(f"\n{fn}\n  {m.n_cells} cells, h = {h * 1e6:.3f} um, "
              f"{len(w) // 2} cut faces, {len(w)} one-sided coefficients")
        pct = [0, 0.1, 1, 5, 25, 50]
        print("  percentile        " + "".join(f"{p:>10g}" for p in pct))
        print("  w = |r|/h         " + "".join(
            f"{np.percentile(w, p):>10.2e}" for p in pct))
        print("  d_f/h             " + "".join(
            f"{np.percentile(d, p):>10.2e}" for p in pct))
        pct2 = [50, 75, 95, 99, 99.9, 100]
        print("  percentile        " + "".join(f"{p:>10g}" for p in pct2))
        print("  a_f / (k A / h)   " + "".join(
            f"{np.percentile(r, p):>10.2e}" for p in pct2))
        for lim in (10, 100, 1000):
            print(f"  share of coefficients with a_f > {lim:4d} k A/h: "
                  f"{100 * np.mean(r > lim):.3f} %")
        # share of the cell diagonal: rho cp V / dt + 6 k A / h + a_f
        for dt in (2e-6, 1e-6):
            for name, sel, k, rcp in (("liquid", liq, K_L, RCP_L),
                                      ("vapour", ~liq, K_V, RCP_V)):
                cond = k * h                      # k A / h with A = h^2
                inert = rcp * h ** 3 / dt
                share = r[sel] * cond / (inert + 5 * cond + r[sel] * cond)
                print(f"  dt = {dt * 1e6:.0f} us, {name} side: a_f is "
                      f"{100 * np.median(share):.1f} % (median), "
                      f"{100 * np.percentile(share, 99):.1f} % (99th pct) "
                      f"of the diagonal; inertia/conduction = "
                      f"{inert / cond:.2f}")


if __name__ == "__main__":
    main()
