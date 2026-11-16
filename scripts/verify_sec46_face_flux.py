#!/usr/bin/env python3
"""Does the orientation dependence of alpha_tilde_D produce an orientation-
dependent ERROR in the CICSAM face value, against the exact wetted fraction?

Planar interface with unit normal n, moving normal to itself.  Donor cell = unit
cube centred at the origin, face f at x = +1/2 (normal e_x), acceptor beyond it,
upwind behind it.  zeta = signed distance of the donor centre from the plane.
The exact face value is the fraction of the face square lying on the alpha side
of the plane.  CICSAM's face value is rebuilt from (alpha_U, alpha_D, alpha_A)
with the paper's Eq. (beta_f) and the limiter branches of App. A.1.
"""
import numpy as np
from verify_sec46_geometry import vof_below, alpha_of_zeta, ORIENT


def face_exact(n, zeta):
    """Wetted fraction of the +x face of the donor for the same phase convention
    as alpha_of_zeta: cell centre at distance zeta ahead of the front, so the
    plane is n.x = -zeta in the donor's frame; the face is x = 1/2."""
    s = -zeta - n[0] * 0.5                        # plane restricted to the face
    m = np.array([n[1], n[2]])
    if np.abs(m).sum() < 1e-12:
        return 1.0 if s >= 0 else 0.0
    return vof_below(m, s)


def cicsam_face(n, zeta, Co, gamma_override=None):
    lam = abs(n[0])
    # geometry: the front advances along +n; the upwind cell sits BEHIND the
    # donor (more alpha, closer to the front's origin), the acceptor AHEAD.
    # alpha_of_zeta decreases with zeta (zeta measured ahead of the front), so
    # the upwind centre is at zeta - lam and the acceptor at zeta + lam.
    aD = alpha_of_zeta(n, zeta)
    aU = alpha_of_zeta(n, zeta - lam)
    aA = alpha_of_zeta(n, zeta + lam)
    if abs(aA - aU) < 1e-12:
        return aD, np.nan, np.nan
    at = (aD - aU) / (aA - aU)
    if 0.0 <= at <= 1.0:
        hc = min(1.0, at / Co)
        uq = min((8 * Co * at + (1 - Co) * (6 * at + 3)) / 8.0, hc)
    else:
        hc = uq = at
    g = n[0] ** 2 if gamma_override is None else gamma_override
    af_t = g * hc + (1 - g) * uq
    af = aU + af_t * (aA - aU)
    beta = (af_t - at) / (1 - at) if at < 1 else 1.0
    return af, at, beta


print("=== face value: CICSAM (blend, Co) against the exact wetted fraction, +x face ===")
for Co in (0.2, 0.023):
    print(f"\n--- Co = {Co} ---")
    print(f"{'orientation':14s} {'zeta/h':>7s} {'a_tilde_D':>9s} {'beta_f':>7s} {'a_f CICSAM':>11s} {'a_f exact':>10s} {'error':>8s}")
    for name, n in ORIENT.items():
        for z in (-0.3, -0.15, 0.0, 0.15, 0.3):
            af, at, beta = cicsam_face(n, z, Co)
            ex = face_exact(n, z)
            print(f"{name:14s} {z:7.2f} {at:9.3f} {beta:7.3f} {af:11.4f} {ex:10.4f} {af - ex:8.4f}")

# ------------------------------------------------------------ azimuthal scan --
print("\n=== equatorial scan of the flux-weighted face-value error, all outflow faces ===")
phis = np.deg2rad(np.arange(0, 360, 2.5))


def fourier(phis, y, m):
    a = 2 * np.mean(y * np.cos(m * phis)); b = 2 * np.mean(y * np.sin(m * phis))
    return np.hypot(a, b), np.rad2deg(np.arctan2(b, a) / m) % (360 / m)


for Co in (0.2, 0.023):
    for gname, gov in (("blend gamma=n_i^2", None), ("gamma=1 (Hyper-C)", 1.0), ("gamma=0 (UQ)", 0.0)):
        for z0 in (-0.15, 0.0, 0.15):
            errs = []
            for ph in phis:
                n = np.array([np.cos(ph), np.sin(ph), 0.0])
                num = den = 0.0
                for i in range(3):
                    if abs(n[i]) < 1e-9:
                        continue
                    # rotate so that face i plays the role of +x, sign so that it is an outflow face
                    perm = [i, (i + 1) % 3, (i + 2) % 3]
                    nn = np.array([abs(n[perm[0]]), n[perm[1]], n[perm[2]]])
                    af, at, beta = cicsam_face(nn, z0, Co, gov)
                    ex = face_exact(nn, z0)
                    w = abs(n[i])                      # flux through face i ~ |n_i|
                    num += w * (af - ex); den += w
                errs.append(num / den)
            errs = np.array(errs)
            a4, p4 = fourier(phis, errs, 4); a2, _ = fourier(phis, errs, 2); a8, _ = fourier(phis, errs, 8)
            print(f"Co={Co:<6} {gname:20s} zeta={z0:+.2f}h  mean err={errs.mean():+.4f}  "
                  f"m4 amp={a4:.4f} @ {p4:5.1f} deg   m2={a2:.1e}  m8={a8:.4f}   "
                  f"err(0deg)={errs[0]:+.4f} err(45deg)={errs[18]:+.4f}")
