#!/usr/bin/env python3
"""Report how far each float sits from the text that first cites it.

Run against the compiled PDF:

    python3 check_float_placement.py Combined_revised.pdf

A positive gap means the float lands that many pages *after* its first
mention, which is the failure this checks for -- the reader meets "Figure 7
shows..." several pages before Figure 7.

Appendix floats cited from the main text produce large gaps legitimately
(Table 13 is in an appendix and cited in Section 4), so those are reported
separately rather than flagged: the test is whether a float has drifted from
its own section, not whether it is far from every citation.

Captions render as "Figure N:" / "Table N:"; citations are the same token not
followed by a colon.  Both are read from `pdftotext -layout` output, so the
measurement is of the typeset result and not of the source order.
"""

import re
import subprocess
import sys

CAPTION = re.compile(r'\b(Figure|Table)s?\s+(\d+)\s*:')
CITATION = re.compile(r'\b(Figure|Table)s?\s+(\d+)(?!\s*:)')
CITATION_2ND = re.compile(r'\b(Figure|Table)s\s+\d+\s+and\s+(\d+)')
DRIFT = 2          # pages; at or above this a main-text float is flagged
APPENDIX_GAP = 10  # pages; above this the float is assumed to be back matter


def pages_of(pdf):
    out = subprocess.run(["pdftotext", "-layout", pdf, "-"],
                         capture_output=True, text=True, errors="ignore")
    if out.returncode:
        sys.exit(f"pdftotext failed on {pdf}")
    return out.stdout.split("\f")


def survey(pages):
    cap, ref = {}, {}
    for i, p in enumerate(pages, 1):
        for m in CAPTION.finditer(p):
            cap.setdefault((m.group(1), int(m.group(2))), i)
        for m in CITATION.finditer(p):
            ref.setdefault((m.group(1), int(m.group(2))), []).append(i)
        for m in CITATION_2ND.finditer(p):
            ref.setdefault((m.group(1), int(m.group(2))), []).append(i)
    return cap, ref


def main(pdf):
    cap, ref = survey(pages_of(pdf))
    drifted, backmatter = [], []
    print(f"{'float':>11} {'caption':>8} {'1st cite':>9} {'gap':>5}   cites")
    print("-" * 66)
    for key in sorted(cap, key=lambda k: (k[0], k[1])):
        c = cap[key]
        cites = sorted(set(ref.get(key, [])))
        first = cites[0] if cites else None
        gap = c - first if first else None
        note = ""
        if gap is None:
            note = "  <-- never cited"
        elif gap > APPENDIX_GAP:
            backmatter.append((key, gap)); note = "  (back matter)"
        elif abs(gap) >= DRIFT:
            drifted.append((key, gap)); note = f"  <-- {abs(gap)} pages"
        print(f"{key[0] + ' ' + str(key[1]):>11} {c:>8} {str(first):>9} "
              f"{str(gap):>5}   {cites[:6]}{note}")

    print(f"\n  drifted by {DRIFT}+ pages: {len(drifted)}")
    for k, g in drifted:
        print(f"     {k[0]} {k[1]:<3} {g:+d}")
    print(f"  back matter cited from the main text (expected): {len(backmatter)}")
    return 1 if drifted else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "Combined_revised.pdf"))
