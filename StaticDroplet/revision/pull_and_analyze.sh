#!/bin/bash
# Pull the static-droplet traces of the CAF review runs from Merlin7 (through
# the transfer node, small files only) and rebuild the table and the figure.
#
#   bash pull_and_analyze.sh            # production runs -> traces/<case>/
#   bash pull_and_analyze.sh pilots     # 400-step pilots -> traces_pilot/<case>/
#
# Safe to repeat while jobs run: the traces are appended every step.
set -euo pipefail
HERE=/home/user/runs/caf-review-2026-11/static-drop
REMOTE=user@service03.merlin7.psi.ch:/data/user/user/cases/caf_static_drop_2026-11
what=${1:-runs}
dest=$HERE/traces; [ "$what" = runs ] || dest=$HERE/traces_pilot
mkdir -p "$dest"
rsync -a --prune-empty-dirs -e "ssh -o BatchMode=yes -o ConnectTimeout=60" \
      --include='*/' --include='spurious.dat' --include='spurious_ext.dat' \
      --include='control' --include='run.info' --include='EXIT_STATUS' \
      --include='RANKS' --exclude='*' "$REMOTE/$what/" "$dest/" 2>&1 | grep -v "module" || true
ls "$dest"
if [ "$what" = runs ]; then
  cd /home/user/papers/UnstructuredBoiling
  python3 analyze_static_drop_long.py | tee "$HERE/analysis/table_$(date +%Y%m%d_%H%M).txt"
fi
