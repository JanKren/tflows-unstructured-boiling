#!/usr/bin/env python3
"""Stage and run one Stefan case (one rank) in runs/<name>.

    run_stefan.py <name> <N> <dt> <nsteps> [KEY=VALUE ...]
                  [--mode M FLUX] [--binary PATH] [--template DIR] [--stage]

Build first: ../build_tools.sh, then ../build_case.sh . (see README.md).

N        cells along x (the paper's mesh is 100 x 5 x 5 on a 2 mm domain);
         0.5 mm, the initial interface, must fall on a cell face
KEY=VALUE  control keywords to replace (or append when absent); a value may
         hold several tokens, written with '+' for the blanks
--mode   writes caf_mode.ini (1 = uniform-temperature check, with the imposed
         interfacial mass flux in kg/m^2/s)
--cubic  cross-section 5 h x 5 h, so that the cells are cubes
--wallT  temperature of the heated wall (the paper's is 10.1)
--stage  only stage the directory, do not run
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def main():
    a = sys.argv[1:]
    if len(a) < 4:
        sys.exit(__doc__)
    name, n, dt, nsteps = a[0], int(a[1]), a[2], int(a[3])
    rest = a[4:]
    mode = None
    binary = os.path.join(HERE, "Process")
    tpl = HERE
    stage_only = False
    wall_t = None
    cubic = False
    over = {}
    i = 0
    while i < len(rest):
        if rest[i] == "--mode":
            mode = (rest[i + 1], rest[i + 2]); i += 3
        elif rest[i] == "--binary":
            binary = rest[i + 1]; i += 2
        elif rest[i] == "--template":
            tpl = rest[i + 1]; i += 2
        elif rest[i] == "--stage":
            stage_only = True; i += 1
        elif rest[i] == "--cubic":
            cubic = True; i += 1
        elif rest[i] == "--wallT":
            wall_t = rest[i + 1]; i += 2
        else:
            k, v = rest[i].split("=", 1)
            over[k] = v.replace("+", " "); i += 1

    d = os.path.join(HERE, "runs", name)
    os.makedirs(d, exist_ok=True)
    for f in os.listdir(d):
        if f.startswith("caf-") or f.endswith((".vtu", ".pvtu", ".backup")):
            os.remove(os.path.join(d, f))

    # mesh
    dom = open(os.path.join(tpl, "stefan.dom")).read()
    dom, k = re.subn(r"^\s*1\s+101\s+6\s+6\s*$", f"  1  {n + 1}    6    6",
                     dom, flags=re.M)
    assert k == 1, "block line not found in stefan.dom"
    half = 0.005
    if cubic:
        # cubic cells: the cross-section is 5 h x 5 h instead of 10 mm x 10 mm.
        # The problem is one-dimensional, so nothing physical changes; what
        # changes is the cell size V^(1/3) that the review branch's floor on
        # the interface distance (1 % of it) is measured against.
        half = 2.5 * 2.0e-3 / n
        dom = dom.replace("-0.005", f"{-half:.9e}").replace(
            " 0.005", f" {half:.9e}")
    open(os.path.join(d, "caf_area.txt"), "w").write(f"{(2 * half) ** 2:.12e}\n")
    open(os.path.join(d, "stefan.dom"), "w").write(dom)
    gen = os.path.join(ROOT, "Binaries", "Generate")
    with open(os.path.join(d, "generate.log"), "w") as lg:
        subprocess.run([gen], input="stefan\nskip\nskip\nskip\n\n", text=True,
                       cwd=d, stdout=lg, stderr=subprocess.STDOUT, check=True)

    # control
    over.setdefault("TIME_STEP", dt)
    over.setdefault("NUMBER_OF_TIME_STEPS", str(nsteps))
    over.setdefault("RESULTS_SAVE_INTERVAL", str(max(nsteps // 3, 1)))
    over.setdefault("BACKUP_SAVE_INTERVAL", "9999999")
    lines = open(os.path.join(tpl, "control")).read().split("\n")
    out, done = [], set()
    for ln in lines:
        s = ln.strip()
        key = s.split()[0] if s and not s.startswith("#") else None
        if key in over and key not in done:
            out.append(f"  {key:<38s} {over[key]}")
            done.add(key)
        else:
            out.append(ln)
    extra = [f"  {k:<38s} {v}" for k, v in over.items() if k not in done]
    if extra:
        # before the boundary conditions, which must stay last
        j = next(i for i, ln in enumerate(out)
                 if ln.strip().startswith("BOUNDARY_CONDITION"))
        out[j:j] = extra + [""]
    if wall_t is not None:
        # the temperature of the Dirichlet wall, last value of its VALUES line
        j = next(i for i, ln in enumerate(out) if ln.split()[:2] ==
                 ["BOUNDARY_CONDITION", "wall"])
        k = next(i for i in range(j, len(out)) if out[i].split()[:1] == ["VALUES"])
        tok = out[k].split()
        tok[-1] = wall_t
        out[k] = "    " + "   ".join(tok)
    open(os.path.join(d, "control"), "w").write("\n".join(out))

    shutil.copy(os.path.join(tpl, "box_parameters.ini"), d)
    shutil.copy(binary, os.path.join(d, "Process"))
    if mode:
        open(os.path.join(d, "caf_mode.ini"), "w").write(f"{mode[0]} {mode[1]}\n")
    elif os.path.exists(os.path.join(d, "caf_mode.ini")):
        os.remove(os.path.join(d, "caf_mode.ini"))
    if stage_only:
        return

    with open(os.path.join(d, "run.log"), "w") as lg:
        r = subprocess.run(["./Process"], cwd=d, stdout=lg,
                           stderr=subprocess.STDOUT)
    rows = sum(1 for _ in open(os.path.join(d, "caf-terms.dat"))) \
        if os.path.exists(os.path.join(d, "caf-terms.dat")) else 0
    print(f"{name}: exit {r.returncode}, {rows} rows of {nsteps}")


if __name__ == "__main__":
    main()
