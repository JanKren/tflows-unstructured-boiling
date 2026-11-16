#!/usr/bin/env python3
"""Run the jobs of a file through run_stefan.py, keeping the number of this
directory's Process instances at or below a cap (default 4, the rank budget
of this work on the shared box).  A job is a line of run_stefan.py arguments."""
import os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 4

def mine():
    n = 0
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            if open(f"/proc/{p}/comm").read().strip() != "Process":
                continue
            if os.readlink(f"/proc/{p}/cwd").startswith(RUNS):
                n += 1
        except OSError:
            pass
    return n

jobs = [ln.split() for ln in open(sys.argv[1]) if ln.strip() and not ln.startswith("#")]
live = []
for j in jobs:
    while True:
        live = [p for p in live if p.poll() is None]
        if max(mine(), len(live)) < cap:
            break
        time.sleep(5)
    live.append(subprocess.Popen([sys.executable, os.path.join(HERE, "run_stefan.py")] + j, cwd=HERE))
    time.sleep(8)          # let Process appear before counting again
for p in live:
    p.wait()
print("dispatch done")
