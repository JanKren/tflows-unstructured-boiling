#!/usr/bin/env python3
"""Coarse preview, R/h = 6.4: the same four-cycle case on both mesh families,
matched and high density ratio, to 0.5 s.  The manuscript does not tabulate
this level (six cells per radius); the runs are here only because they finish
in minutes and show the time behaviour before the R/h = 12.8 jobs start.
"""
import math, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = {"hex": "/home/user/runs/Droplet3D/hex32_sm4/control",
       "poly": "/home/user/runs/Droplet3D/poly32/control"}
L, N, SIG, SAFETY, T_END = 2.0, 32, 15.0e3, 0.16, 0.5
H = L / N
CASES = [("c32m_hex_s4",  "hex",  (1.0, 1.0),    (1.0, 1.0)),
         ("c32m_poly_s4", "poly", (1.0, 1.0),    (1.0, 1.0)),
         ("c32h_hex_s4",  "hex",  (1.0, 1.0e-3), (1.0, 1.0e-2)),
         ("c32h_poly_s4", "poly", (1.0, 1.0e-3), (1.0, 1.0e-2))]
def sub(text, key, value):
    pat = re.compile(r"^(\s*%s\s+).*$" % re.escape(key), re.M)
    assert pat.search(text), key
    return pat.sub(lambda m: m.group(1) + value, text)
for name, mesh, rho, mu in CASES:
    t = open(SRC[mesh]).read()
    assert re.search(r"MAX_SMOOTHING_CYCLES_CURVATURE_VOF\s+4\b", t), name
    dtmax = math.sqrt((rho[0] + rho[1]) * H**3 / (8.0 * SIG))
    dt = float("%.4e" % (SAFETY * dtmax))
    steps = int(math.ceil(T_END / dt / 100.0)) * 100
    t = sub(t, "NUMBER_OF_TIME_STEPS", str(steps))
    t = sub(t, "RESULTS_SAVE_INTERVAL", str(steps))
    t = sub(t, "BACKUP_SAVE_INTERVAL", str(steps))
    t = sub(t, "TIME_STEP", "%.4e" % dt)
    t = sub(t, "PHASE_DENSITIES", "%.4e  %.4e" % rho)
    t = sub(t, "PHASE_VISCOSITIES", "%.4e  %.4e" % mu)
    t += ("\n#--- added for the long runs ---------------------------\n"
          " GeometricVOF            no\n WALL_TIME_MAX_HOURS     0.9\n"
          " SAVE_INITIAL_CONDITION  no\n")
    open(os.path.join(HERE, "controls", name + ".control"), "w").write(t)
    print("%-13s %-4s rho %g/%g dt %.4e (%.2f) steps %d" % (name, mesh, rho[0], rho[1], dt, dt/dtmax, steps))
