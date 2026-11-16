#!/usr/bin/env python3
"""Independent check of Section 4.6 (normalised donor value on a cubic cell).

Everything here is closed-form geometry on a unit cube (h = 1); nothing is read
from solver output.  It reproduces Table `tab:alpha_tilde` and tests the
sentences of Section 4.6 one by one.
"""
import itertools
import numpy as np

# ---------------------------------------------------------------- geometry ---
def _terms(nabs, d, power):
    """sum_v (-1)^|v| max(0, d - n.v)^power over the cube vertices v in {0,1}^k."""
    tot = 0.0
    k = len(nabs)
    for v in itertools.product((0, 1), repeat=k):
        arg = d - float(np.dot(nabs, v))
        if arg > 0:
            tot += (-1) ** sum(v) * arg ** power
    return tot


def vof_below(n, s):
    """Fraction of the unit cube centred at the origin with  n.x <= s.

    n is a unit vector; s the signed distance of the plane from the centre.
    Exact, from the standard cut-cube formula; degenerate (zero) components
    handled by dropping them (the cube is a product of intervals).
    """
    nabs = np.abs(np.asarray(n, float))
    keep = nabs > 1e-12
    a = nabs[keep]
    k = a.size
    d = s + a.sum() / 2.0            # distance measured from the lowest vertex
    if k == 0:
        return 1.0 if s >= 0 else 0.0
    if d <= 0:
        return 0.0
    if d >= a.sum():
        return 1.0
    from math import factorial
    return _terms(a, d, k) / (factorial(k) * np.prod(a))


def area_cross(n, s):
    """d(vof_below)/ds = cross-sectional area of the cube by the plane at s."""
    nabs = np.abs(np.asarray(n, float))
    keep = nabs > 1e-12
    a = nabs[keep]
    k = a.size
    d = s + a.sum() / 2.0
    if k == 0 or d <= 0 or d >= a.sum():
        return 0.0
    from math import factorial
    return _terms(a, d, k - 1) / (factorial(k - 1) * np.prod(a))


def alpha_of_zeta(n, zeta):
    """alpha(.) of Section 4.6: the volume fraction of the phase *behind* the
    front, for a cell whose centre sits at signed distance zeta from the plane,
    zeta > 0 meaning the centre is ahead of the front (less filled)."""
    return vof_below(n, -zeta)


def alpha_prime(n, zeta):
    return -area_cross(n, -zeta)     # d alpha / d zeta


# ---------------------------------------------------------- normalised donor --
def alpha_tilde_D(n, zeta, lam):
    aD = alpha_of_zeta(n, zeta)
    aU = alpha_of_zeta(n, zeta + lam)
    aA = alpha_of_zeta(n, zeta - lam)
    return (aD - aU) / (aA - aU), aU, aD, aA


def slope_at(n, lam, zeta=0.0, eps=1e-6):
    p, *_ = alpha_tilde_D(n, zeta + eps, lam)
    m, *_ = alpha_tilde_D(n, zeta - eps, lam)
    return (p - m) / (2 * eps)


def slope_formula(n, lam):
    """Closed form at zeta = 0:  (alpha'(0) - alpha'(lam)) / (1 - 2 alpha(lam))."""
    return (alpha_prime(n, 0.0) - alpha_prime(n, lam)) / (1 - 2 * alpha_of_zeta(n, lam))


ORIENT = {
    "axis":          np.array([1.0, 0.0, 0.0]),
    "face diagonal": np.array([1.0, 1.0, 0.0]) / np.sqrt(2),
    "body diagonal": np.array([1.0, 1.0, 1.0]) / np.sqrt(3),
}

print("=== Table tab:alpha_tilde, face along x, interface at the donor centre ===")
print(f"{'orientation':14s} {'|n|_1':>7s} {'W/h':>7s} {'lam/W':>7s} {'aU':>7s} {'aA':>7s} "
      f"{'slope(num)':>11s} {'slope(formula)':>15s} {'alpha_prime(0)':>15s} {'A_mid/h^2':>10s}")
for name, n in ORIENT.items():
    W = np.abs(n).sum()
    lam = abs(n[0])
    at, aU, aD, aA = alpha_tilde_D(n, 0.0, lam)
    print(f"{name:14s} {W:7.3f} {W:7.3f} {lam / W:7.3f} {aU:7.3f} {aA:7.3f} "
          f"{slope_at(n, lam):11.4f} {slope_formula(n, lam):15.4f} "
          f"{alpha_prime(n, 0.0):15.4f} {area_cross(n, 0.0):10.4f}")
print(f"exact values quoted in the caption: -1, -sqrt2 = {-np.sqrt(2):.4f}, "
      f"-15 sqrt3/23 = {-15 * np.sqrt(3) / 23:.4f}")
print(f"mid-plane cross-sections: 1, sqrt2 = {np.sqrt(2):.4f}, 3 sqrt3/4 = {3 * np.sqrt(3) / 4:.4f}")

print("\n=== denominator 1 - 2 alpha(lam) ('barely moves') ===")
for name, n in ORIENT.items():
    lam = abs(n[0])
    print(f"{name:14s} {1 - 2 * alpha_of_zeta(n, lam):.4f}")

print("\n=== profile shape in units of W: alpha at offsets -0.25W, -0.125W (text says 'sigmoidal') ===")
for name, n in ORIENT.items():
    W = np.abs(n).sum()
    print(f"{name:14s} alpha(0.25W) = {alpha_of_zeta(n, 0.25 * W):.4f}   "
          f"alpha(0.125W) = {alpha_of_zeta(n, 0.125 * W):.4f}   (linear ramp: 0.25, 0.375)")

print("\n=== counterfactual: a LINEAR ramp of the same width W, same sampling lam ===")
def lin_alpha(W, z):
    return float(np.clip(0.5 - z / W, 0.0, 1.0))
for name, n in ORIENT.items():
    W = np.abs(n).sum(); lam = abs(n[0]); eps = 1e-6
    def at(z):
        aD, aU, aA = lin_alpha(W, z), lin_alpha(W, z + lam), lin_alpha(W, z - lam)
        return (aD - aU) / (aA - aU)
    print(f"{name:14s} linear-ramp slope d(alpha_tilde)/d(zeta/h) = {(at(eps) - at(-eps)) / (2 * eps):8.4f}")

print("\n=== is the slope constant near the centre?  alpha_tilde_D(zeta) for |zeta| <= 0.3h ===")
for name, n in ORIENT.items():
    lam = abs(n[0])
    zs = np.linspace(-0.3, 0.3, 7)
    vals = [alpha_tilde_D(n, z, lam)[0] for z in zs]
    print(f"{name:14s} " + "  ".join(f"{v:.3f}" for v in vals))

# ------------------------------------------------ azimuthal scan, equatorial --
print("\n=== equatorial scan: per-face slope and its six-face aggregate vs azimuth ===")
phis = np.deg2rad(np.arange(0, 360, 2.5))
mean_slope, flux_slope, mean_at = [], [], []
for ph in phis:
    n = np.array([np.cos(ph), np.sin(ph), 0.0])
    faces = []
    for i in range(3):
        lam = abs(n[i])
        if lam < 1e-9:
            continue                       # theta_f = 90 deg, no normal flux, alpha_tilde undefined
        faces.append((lam, slope_at(n, lam), alpha_tilde_D(n, 0.15, lam)[0]))
    lam_arr = np.array([f[0] for f in faces]); sl = np.array([f[1] for f in faces]); at = np.array([f[2] for f in faces])
    mean_slope.append(sl.mean())
    flux_slope.append((lam_arr * sl).sum() / lam_arr.sum())
    mean_at.append(at.mean())

def fourier(phis, y, m):
    a = 2 * np.mean(y * np.cos(m * phis)); b = 2 * np.mean(y * np.sin(m * phis))
    amp = np.hypot(a, b) / np.mean(y); phase = np.rad2deg(np.arctan2(b, a) / m) % (360 / m)
    return amp, phase

for label, y in (("unweighted mean slope", np.array(mean_slope)),
                 ("|n_i|-weighted slope", np.array(flux_slope)),
                 ("mean alpha_tilde_D at zeta=0.15h", np.array(mean_at))):
    print(f"  {label:34s} mean = {np.mean(y):8.4f}   " +
          "   ".join(f"m={m}: {100 * fourier(phis, y, m)[0]:6.2f}% @ {fourier(phis, y, m)[1]:6.2f} deg"
                     for m in (2, 4, 6, 8)))
