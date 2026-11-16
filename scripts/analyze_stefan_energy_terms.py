#!/usr/bin/env python3
"""Term-by-term energy and mass budget of the Stefan problem.

    python3 analyze_stefan_energy_terms.py LABEL=RUN_DIR [LABEL=RUN_DIR ...]
            [--window 4.5 9.0] [--csv out.csv]

Reads `caf-terms.dat` and `caf-vof.dat`, written once per time step by the
hooks in `instrumentation/caf_stefan_terms/` (End_Of_Time_Step.f90 documents
the columns).  Rates are recorded, so the integrals below use the quadrature
of the scheme itself: the energy equation is implicit Euler, so every heat
rate is the end-of-step value times dt, and the volume-fraction step n+1 is
advanced with the mass transfer estimated from the temperature of step n.

ENERGY.  With E(t) the time integral of a rate from 0 to t,

    supplied   E_wall                                  heat through the wall
    stored     dH_v + dH_l + E_lat                     sensible + latent
               E_lat = h_lv rho_v [V_vap(t) - V_vap(0)]

and the closure error of the paper is (stored - supplied)/supplied.  The gap
supplied - stored is split exactly into

    (1) E_wall - dH_v - E_int_v     vapour cells: heat in, minus stored,
                                    minus removed through the interface
                                    coupling a_f (T_P - T_sat)       [Eq. 11]
    (2) E_int_v - E_qv              the same vapour-side heat as the
                                    mass-transfer estimate sees it,
                                    k_v grad T . S of the extrapolated
                                    least-squares gradient            [Eq. 7]
    (3) E_ql                        heat the estimate sends into the liquid,
                                    k_l grad T . S on the liquid side; it
                                    lowers m_dot and is zero in the exact
                                    problem, the liquid being at T_sat
    (4) E_mdot(end) - E_mdot(vof)   the estimate of the last outer iteration
                                    against the one the next volume-fraction
                                    step is advanced with
    (5) E_mdot(vof) - E_lat         latent heat booked by m_dot against the
                                    vapour actually produced
    (6) -dH_l                       sensible heat of the liquid cells

with h_lv m_dot = E_qv - E_ql.  (1)+...+(6) = supplied - stored identically.
`liquid sink` is the heat the coupling feeds into the liquid cells,
-E_int_l - dH_l, which the liquid energy equation loses.

MASS.  M_d is the mass in the domain and out_rate the mass flux through the
outflow.  The volume-fraction step n is advanced with the face fluxes of step
n-1, so the egress consistent with the scheme is out_rate(n-1) dt, and

    M_s(n) = M_d(n) + sum_{k<=n} out_rate(k-1) dt

is the quantity that should stay at M_0.  The error is reported with that
quadrature, with the trapezoidal rule and with the end-of-step rate, and
divided by the expelled mass fraction, for comparison with 2 rho_v/rho_l.
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf

RHO_V, RHO_L = 0.597, 958.4
CP_V, CP_L = 2030.0, 4216.0
K_V = 0.025
H_LV = 2.26e6
DT_WALL = 10.0
X0 = 0.5e-3           # initial interface position [m]
AREA = 1.0e-4         # cross-section [m^2]; caf_area.txt in the run overrides it

D_V = K_V / (RHO_V * CP_V)
STE = CP_V * DT_WALL / H_LV
XI = brentq(lambda x: x * np.exp(x * x) * erf(x) - STE / np.sqrt(np.pi),
            1e-6, 1.0)
T0 = X0 ** 2 / (4.0 * XI ** 2 * D_V)


def x_exact(t):
    return 2.0 * XI * np.sqrt(D_V * (t + T0))


def load(d):
    a = np.loadtxt(os.path.join(d, "caf-terms.dat"))
    v = np.loadtxt(os.path.join(d, "caf-vof.dat"))
    # one VOF record per time step: keep the last of each step (sub-cycling)
    last = {}
    for row in v:
        last[int(row[0])] = row
    steps = a[:, 0].astype(int)
    mv = np.array([last[s][2] for s in steps])
    return a, mv


def analyse(d, window):
    a, mdot_vof = load(d)
    area = AREA
    if os.path.exists(os.path.join(d, "caf_area.txt")):
        area = float(open(os.path.join(d, "caf_area.txt")).read())
    t = a[:, 1]
    dt = t[1] - t[0]
    m_d, v_vap, out, q_out, mdot_end = (a[:, i] for i in range(2, 7))
    q_wall, q_iv, q_il, h_v, h_l, x_f = (a[:, i] for i in range(8, 14))
    q0, q1 = a[:, 18], a[:, 19]

    v0 = X0 * area
    hv0 = RHO_V * CP_V * 0.5 * DT_WALL * v0       # linear initial profile

    def E(rate):
        return np.cumsum(rate) * dt

    e_wall, e_iv, e_il = E(q_wall), E(q_iv), E(q_il)
    e_qv, e_ql = E(-q0), E(-q1)
    e_m_end = E(mdot_end) * H_LV
    e_m_vof = E(mdot_vof) * H_LV
    e_lat = H_LV * RHO_V * (v_vap - v0)
    dh_v, dh_l = h_v - hv0, h_l

    def at(arr, tt):
        return arr[np.argmin(np.abs(t - tt))]

    def span(arr, t1, t2):
        lo = at(arr, t1) if t1 > 0 else 0.0
        return at(arr, t2) - lo

    res = {}
    for name, (t1, t2) in (("whole", (0.0, t[-1])), ("window", window)):
        s = lambda arr: span(arr, t1, t2)
        sup = s(e_wall)
        terms = {
            "closure": (s(dh_v) + s(dh_l) + s(e_lat) - sup) / sup,
            "(1) vapour cells": (s(e_wall) - s(dh_v) - s(e_iv)) / sup,
            "(2) Eq.11 - Eq.7, vapour": (s(e_iv) - s(e_qv)) / sup,
            "(3) estimate into liquid": s(e_ql) / sup,
            "(4) iteration lag": (s(e_m_end) - s(e_m_vof)) / sup,
            "(5) booked - produced": (s(e_m_vof) - s(e_lat)) / sup,
            "(6) -dH_l": -s(dh_l) / sup,
            "liquid sink": (-s(e_il) - s(dh_l)) / sup,
            "latent share of stored": s(e_lat) / (s(dh_v) + s(dh_l) + s(e_lat)),
            "supplied [J]": sup,
        }
        res[name] = terms

    # interface position, the front intersection and the volume-equivalent one
    xe = x_exact(t)
    ex_front = (x_f - xe) / xe
    ex_vol = (v_vap / area - xe) / xe
    res["x"] = {"front, final": ex_front[-1],
                "front, mean |e|": np.mean(np.abs(ex_front)),
                "front, max |e|": np.max(np.abs(ex_front)),
                "volume, final": ex_vol[-1],
                "volume, mean |e|": np.mean(np.abs(ex_vol)),
                "volume, max |e|": np.max(np.abs(ex_vol))}

    # mass: three quadratures of the egress
    prev = np.concatenate(([out[0]], out[:-1]))
    m0 = m_d[0] + prev[0] * dt
    frac = 1.0 - m_d / m0
    quad = {"lagged": prev, "trapezoid": 0.5 * (out + prev), "end of step": out}
    res["mass"] = {"expelled fraction": frac[-1]}
    for k, r in quad.items():
        err = (m_d + np.cumsum(r) * dt) / m0 - 1.0
        res["mass"][k] = err[-1]
        res["mass"][k + " / expelled"] = err[-1] / frac[-1]
    res["mass"]["2 rho_v/rho_l"] = 2.0 * RHO_V / RHO_L
    # global continuity of the last outer iteration: Q_out against the
    # volume produced by the mass transfer of that iteration
    res["mass"]["Q_out / [mdot (1/rho_v - 1/rho_l)] - 1, median"] = np.median(
        q_out[5:] / (mdot_end[5:] * (1 / RHO_V - 1 / RHO_L)) - 1.0)
    res["dt"] = dt
    res["n"] = len(t)
    res["T_min - T_sat"] = a[:, 17].min() - 0.1
    return res


def main():
    args = sys.argv[1:]
    window = (4.5, 9.0)
    csv = None
    jobs = []
    i = 0
    while i < len(args):
        if args[i] == "--window":
            window = (float(args[i + 1]), float(args[i + 2])); i += 3
        elif args[i] == "--csv":
            csv = args[i + 1]; i += 2
        else:
            jobs.append(args[i].split("=", 1)); i += 1
    if not jobs:
        sys.exit(__doc__)

    print(f"Stefan similarity solution: Ste = {STE:.6f}, xi = {XI:.6f}, "
          f"D_v = {D_V:.4e} m2/s, t0 = {T0:.4f} s")
    out = {lab: analyse(d, window) for lab, d in jobs}
    labs = [lab for lab, _ in jobs]
    w = max(12, max(len(x) for x in labs) + 1)

    def table(title, key, fmt, scale=1.0, skip=()):
        print(f"\n{title}")
        print(f"  {'':<30}" + "".join(f"{x:>{w}}" for x in labs))
        for name in out[labs[0]][key]:
            if name in skip:
                continue
            print(f"  {name:<30}" + "".join(
                f"{fmt.format(scale * out[x][key][name]):>{w}}" for x in labs))

    for key, title in (("whole", "ENERGY, whole record, % of the heat supplied"),
                       ("window", f"ENERGY, t in [{window[0]:g}, {window[1]:g}] s,"
                                  " % of the heat supplied")):
        table(title, key, "{:+.3f}", 100.0,
              skip=("supplied [J]", "latent share of stored"))
        print(f"  {'latent share of stored [%]':<30}" + "".join(
            f"{100 * out[x][key]['latent share of stored']:>{w}.2f}" for x in labs))
        print(f"  {'supplied [J]':<30}" + "".join(
            f"{out[x][key]['supplied [J]']:>{w}.5e}" for x in labs))
    table("INTERFACE POSITION, relative error against the similarity solution [%]",
          "x", "{:+.3f}", 100.0)
    table("MASS, M_s/M_0 - 1 at the end of the record", "mass", "{:+.4e}")
    print(f"\n  {'dt [s]':<30}" + "".join(f"{out[x]['dt']:>{w}.3e}" for x in labs))
    print(f"  {'steps':<30}" + "".join(f"{out[x]['n']:>{w}d}" for x in labs))
    print(f"  {'min T - T_sat [K]':<30}" + "".join(
        f"{out[x]['T_min - T_sat']:>{w}.3e}" for x in labs))

    if csv:
        with open(csv, "w") as f:
            f.write("label,section,term,value\n")
            for lab in labs:
                for sec in ("whole", "window", "x", "mass"):
                    for k, v in out[lab][sec].items():
                        f.write(f"{lab},{sec},{k},{v:.10e}\n")


if __name__ == "__main__":
    main()
