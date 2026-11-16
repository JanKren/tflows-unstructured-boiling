#!/usr/bin/env python3
"""Source-level structural checks for the manuscript, for use without a compiler.

    python3 check_latex_structure.py

Written after an edit deleted a paragraph by searching forward to the next
blank line, which silently took the `\\end{table}` that sat between them. The
document still had balanced braces and even `$` counts; the first sign of
trouble was pdflatex reporting "Not in outer par mode" against a `\\begin{figure}`
thirty lines later. Environment balance is the check that catches that class of
damage, and it costs nothing to run before every commit.

Checks, in order of how badly they break a build:
  1. environment begin/end balance, per environment, per file
  2. brace balance and inline-math delimiter parity
  3. \\ref/\\eqref targets that have no \\label, and duplicated \\label names
  4. \\cite keys absent from references.bib
  5. \\includegraphics targets that do not exist on disk
  6. tabular/longtable rows whose cell count differs from the column spec

Exit status is 1 if anything in 1-5 fails, so it can gate a commit hook.
Row-count mismatches (6) are reported but do not set the status, because
\\multicolumn and p{} specs make a small false-positive rate unavoidable.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
FILES = ["Combined_revised.tex", "annular_flow_results_revised.tex"]
# environments worth tracking; \begin{document} and the like are not in these files
ENVS = ["table", "figure", "tabular", "longtable", "equation", "align",
        "subfigure", "center", "itemize", "enumerate", "abstract"]


def strip_comments(text):
    """Remove % comments without touching escaped \\% inside cells."""
    return re.sub(r"(?<!\\)%[^\n]*", "", text)


def main():
    fail = False
    sources = {}
    for name in FILES:
        path = os.path.join(ROOT, name)
        if not os.path.exists(path):
            print(f"MISSING SOURCE: {name}")
            return 1
        sources[name] = strip_comments(open(path, encoding="utf-8").read())

    print("1. environment balance")
    for name, s in sources.items():
        for env in ENVS:
            b = len(re.findall(r"\\begin\{" + env + r"\}", s))
            e = len(re.findall(r"\\end\{" + env + r"\}", s))
            if b != e:
                fail = True
                print(f"   UNBALANCED  {name}: {env} begin={b} end={e}")
    if not fail:
        print("   all balanced")

    print("2. braces and inline math")
    for name, s in sources.items():
        d = s.count("{") - s.count("}")
        # \$ is a literal dollar and must not count toward the parity
        n = len(re.findall(r"(?<!\\)\$", s))
        if d:
            fail = True
            print(f"   BRACE IMBALANCE  {name}: {d:+d}")
        if n % 2:
            fail = True
            print(f"   ODD $ COUNT      {name}: {n}")
    if not fail:
        print("   balanced")

    joined = "".join(sources.values())
    print("3. labels and references")
    labels = re.findall(r"\\label\{([^}]*)\}", joined)
    refs = set(re.findall(r"\\(?:ref|eqref)\{([^}]*)\}", joined))
    dangling = sorted(refs - set(labels))
    dup = sorted({l for l in labels if labels.count(l) > 1})
    if dangling:
        fail = True
        print(f"   DANGLING REFS: {dangling}")
    if dup:
        fail = True
        print(f"   DUPLICATE LABELS: {dup}")
    if not dangling and not dup:
        print(f"   {len(set(labels))} labels, {len(refs)} referenced, none dangling or duplicated")

    print("4. citations")
    bib = open(os.path.join(ROOT, "references.bib"), encoding="utf-8").read()
    keys = {k for grp in re.findall(r"\\cite\{([^}]*)\}", joined)
            for k in (x.strip() for x in grp.split(","))}
    missing = sorted(k for k in keys if "{" + k + "," not in bib)
    if missing:
        fail = True
        print(f"   NOT IN BIB: {missing}")
    else:
        print(f"   {len(keys)} keys, all present in references.bib")

    print("5. graphics")
    absent = []
    for g in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]*)\}", joined):
        if not any(os.path.exists(os.path.join(ROOT, g + ext))
                   for ext in ("", ".pdf", ".png", ".jpg")):
            absent.append(g)
    if absent:
        fail = True
        print(f"   MISSING FILES: {absent}")
    else:
        print("   every \\includegraphics target exists")

    print("6. table rows against column specs (advisory)")
    bad = 0
    for m in re.finditer(r"\\begin\{(tabular|longtable)\}\{((?:[^{}]|\{[^{}]*\})*)\}(.*?)"
                         r"\\end\{\1\}", joined, re.S):
        spec = re.sub(r"@\{[^}]*\}", "", m.group(2))
        ncol = len([c for c in spec if c in "lcrp"])
        for row in re.split(r"\\\\", m.group(3)):
            row = " ".join(row.split())
            if not row or "&" not in row or "multicolumn" in row:
                continue
            if row.count("&") + 1 != ncol:
                bad += 1
                print(f"   spec={ncol} got={row.count('&') + 1}: {row[:60]}")
    print(f"   {bad} row(s) differ from their column spec")

    print("\nFAIL" if fail else "\nOK")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
