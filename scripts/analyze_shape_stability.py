#!/usr/bin/env python3
"""Linear shape stability of the Scriven bubble at the paper's parameters.

Reviewer 2 (CAF-D-26-01303): "A vapour bubble growing in uniformly superheated
liquid without Gibbs-Thomson correction is susceptible to diffusive
(Mullins-Sekerka-type) shape instability ... The authors ... do not estimate
the growth rate of the l = 4 mode at the present Ja, R and sigma."

Two mechanisms are estimated, on the analytical trajectory
R(t) = 2 beta sqrt(D_l (t + t0)) of the benchmark.

1. The hydrodynamic problem.  The amplitude a_n of the shape mode n of a
   bubble whose interface moves with the liquid obeys (Plesset 1954;
   Prosperetti 1977; Prosperetti & Seminara 1978)

       a'' + B(t) a' - A(t) a = 0
       A = (n-1) R''/R - (n-1)(n+1)(n+2) sigma/(rho_l R^3)
           - (n-1)(n+2) 2 nu R'/R^3
       B = 3 R'/R + (n+2)(2n+1) 2 nu / R^2

   (thin-viscous-layer form; the boundary-layer correction of Prosperetti is
   switched on with BL = True and changes the damping by about a tenth).
   It is integrated from a seed a(0) = a0, a'(0) = 0.

2. The diffusive (Mullins-Sekerka) mechanism.  A protrusion of the interface
   into the thermal layer steepens the local gradient and raises the local
   mass flux.  For a solidification front the interface velocity IS that flux
   divided by the latent heat, which gives the growth rate (n-1) R'/R.  For a
   vapour bubble the interface moves with the liquid, and the local flux
   displaces it relative to the liquid only at  m''/rho_l = (rho_v/rho_l) R'.
   The diffusive growth rate is therefore at most

       s_MS = (rho_v/rho_l) g R'/R ,   g <= R/delta_T

   with g the sensitivity of the flux to a unit relative displacement
   (g = n-1 in the quasi-static limit, R/delta_T if the isotherms do not
   follow the interface at all).  It is compared with the damping rate B/2.
   The vapour recoil pressure rho_v R'^2 enters the normal-stress balance and
   is compared with the capillary pressure of the same perturbation.
"""

import numpy as np
from scipy.integrate import solve_ivp

import scriven_reference as ref

RHO_L, RHO_V = ref.RHO_L, ref.RHO_V
MU_L = 2.82e-4
SIGMA = 0.059
NU = MU_L / RHO_L
BETA, D_L, T0 = ref.BETA, ref.D_L, ref.T0
T_END = 1.0e-3


def R(t):
    return 2.0 * BETA * np.sqrt(D_L * (t + T0))


def Rd(t):
    return BETA * np.sqrt(D_L / (t + T0))


def Rdd(t):
    return -0.5 * BETA * np.sqrt(D_L) * (t + T0) ** -1.5


def omega(n, r):
    return np.sqrt((n - 1) * (n + 1) * (n + 2) * SIGMA / (RHO_L * r ** 3))


def coeffs(n, t, bl):
    r, rd, rdd = R(t), Rd(t), Rdd(t)
    d = min(np.sqrt(NU / omega(n, r)), r / (2 * n)) / r if bl else 0.0
    A = ((n - 1) * rdd / r - omega(n, r) ** 2
         - ((n - 1) * (n + 2) + 2 * n * (n + 2) * (n - 1) * d)
         * 2 * NU * rd / r ** 3)
    B = (3 * rd / r
         + ((n + 2) * (2 * n + 1) - 2 * n * (n + 2) ** 2 * d) * 2 * NU / r ** 2)
    return A, B


def integrate(n, bl):
    def rhs(t, y):
        A, B = coeffs(n, t, bl)
        return [y[1], A * y[0] - B * y[1]]
    sol = solve_ivp(rhs, (0.0, T_END), [1.0, 0.0], method="DOP853",
                    rtol=1e-10, atol=1e-14, dense_output=True)
    t = np.linspace(0.0, T_END, 200001)
    a = sol.sol(t)[0]
    return t, a


def envelope(t, a, t_lo, t_hi):
    m = (t >= t_lo) & (t <= t_hi)
    return np.abs(a[m]).max()


def main():
    print(f"beta = {BETA:.4f}  D_l = {D_L:.4e} m2/s  t0 = {T0*1e6:.1f} us  "
          f"nu_l = {NU:.3e} m2/s")
    print(f"R: {R(0)*1e6:.2f} -> {R(T_END)*1e6:.2f} um over {T_END*1e3:.1f} ms\n")

    ja = RHO_L * ref.CP_L * ref.D_T / (RHO_V * ref.H_LV)
    for t in (0.0, T_END):
        r, rd = R(t), Rd(t)
        delta = ja * D_L / rd            # k_l dT / (rho_v h_lv R')
        print(f"t = {t*1e3:.1f} ms: R'/R = {rd/r:8.1f} 1/s   "
              f"delta_T/R = {delta/r:.3f}")
        for n in (2, 4):
            w = omega(n, r)
            _, B = coeffs(n, t, False)
            g_qs, g_max = n - 1, r / delta
            s_solid = (n - 1) * rd / r
            s_qs = RHO_V / RHO_L * g_qs * rd / r
            s_max = RHO_V / RHO_L * g_max * rd / r
            recoil = 2 * RHO_V * rd ** 2 * g_max
            capill = (n - 1) * (n + 2) * SIGMA / r
            print(f"   n = {n}: omega_n = {w:9.3e} rad/s (period "
                  f"{2*np.pi/w*1e6:6.1f} us)  damping B/2 = {B/2:8.1f} 1/s")
            print(f"          Mullins-Sekerka rate, front slaved to the flux "
                  f"(solidification): {s_solid:8.1f} 1/s")
            print(f"          same for a vapour bubble, x rho_v/rho_l: "
                  f"{s_qs:6.2f} (quasi-static) .. {s_max:6.2f} 1/s (bound)  "
                  f"-> damping/growth >= {B/2/s_max:6.0f}")
            print(f"          recoil / capillary pressure perturbation "
                  f"<= {recoil/capill:.1e}")
        print()

    print("decay of a seed a(0) = a0 over the record (relative amplitude "
          "a/R against a0/R0):")
    for n in (2, 4):
        for bl in (False, True):
            t, a = integrate(n, bl)
            rel = a * R(0) / R(t)
            cycles = np.sum(np.diff(np.sign(a)) != 0) / 2.0
            row = "  ".join(
                f"{tt*1e3:.1f} ms: {envelope(t, rel, tt - 0.06e-3, tt):.3f}"
                for tt in (0.2e-3, 0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3))
            print(f"   n = {n}  {'with' if bl else 'no  '} BL correction: "
                  f"{cycles:5.1f} oscillations;  envelope of (a/R)/(a0/R0)  {row}")


if __name__ == "__main__":
    main()
