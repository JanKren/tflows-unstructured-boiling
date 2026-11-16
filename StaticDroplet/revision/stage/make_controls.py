#!/usr/bin/env python3
"""Controls of the long static-droplet runs for the CAF review (comment 2.9).

    python3 make_controls.py

Each control is the published one of ~/runs/Droplet3D/{hex64,poly64}/control
with only the listed keywords changed, so a diff against the original shows
exactly what a case varies.
"""
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = {"hex": "/home/user/runs/Droplet3D/hex64/control",
       "poly": "/home/user/runs/Droplet3D/poly64/control"}
L, N, SIG, SAFETY = 2.0, 64, 15.0e3, 0.16
H = L / N
T_END = 0.5                      # 3.1 viscous times rho_l R^2 / mu_l

#        name          mesh   cycles  rho(out, sphere)   mu(out, sphere)
CASES = [("m_hex_s12",  "hex",  12, (1.0, 1.0),    (1.0, 1.0)),
         ("m_poly_s4",  "poly",  4, (1.0, 1.0),    (1.0, 1.0)),
         ("m_hex_s4",   "hex",   4, (1.0, 1.0),    (1.0, 1.0)),
         ("m_poly_s12", "poly", 12, (1.0, 1.0),    (1.0, 1.0)),
         ("h_hex_s12",  "hex",  12, (1.0, 1.0e-3), (1.0, 1.0e-2)),
         ("h_poly_s4",  "poly",  4, (1.0, 1.0e-3), (1.0, 1.0e-2)),
         ("h_hex_s4",   "hex",   4, (1.0, 1.0e-3), (1.0, 1.0e-2))]


def sub(text, key, value):
    pat = re.compile(r"^(\s*%s\s+).*$" % re.escape(key), re.M)
    assert pat.search(text), key
    return pat.sub(lambda m: m.group(1) + value, text)


def main():
    os.makedirs(os.path.join(HERE, "controls"), exist_ok=True)
    for name, mesh, cyc, rho, mu in CASES:
        t = open(SRC[mesh]).read()
        dtmax = math.sqrt((rho[0] + rho[1]) * H**3 / (8.0 * SIG))
        dt = float("%.4e" % (SAFETY * dtmax))
        steps = int(math.ceil(T_END / dt / 100.0)) * 100
        t = sub(t, "NUMBER_OF_TIME_STEPS", str(steps))
        t = sub(t, "RESULTS_SAVE_INTERVAL", str(steps // 4))
        t = sub(t, "BACKUP_SAVE_INTERVAL", str(steps // 4))
        t = sub(t, "TIME_STEP", "%.4e" % dt)
        t = sub(t, "PHASE_DENSITIES", "%.4e  %.4e" % rho)
        t = sub(t, "PHASE_VISCOSITIES", "%.4e  %.4e" % mu)
        t = sub(t, "MAX_SMOOTHING_CYCLES_CURVATURE_VOF", str(cyc))
        t = re.sub(r"#   dt_max   = .*", "#   dt_max   = %.4e" % dtmax, t)
        t += ("\n#--- added for the long runs ---------------------------\n"
              " GeometricVOF            no\n"
              " WALL_TIME_MAX_HOURS     23.3\n"
              " SAVE_INITIAL_CONDITION  no\n")
        open(os.path.join(HERE, "controls", name + ".control"), "w").write(t)
        print("%-11s %-4s cycles %2d  rho %g/%g  mu %g/%g  dt %.4e "
              "(%.2f of %.4e)  steps %d  t_end %.4f s"
              % (name, mesh, cyc, rho[0], rho[1], mu[0], mu[1], dt,
                 dt / dtmax, dtmax, steps, steps * dt))


if __name__ == "__main__":
    main()
