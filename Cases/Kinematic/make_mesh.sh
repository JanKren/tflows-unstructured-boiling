#!/bin/bash
# make_mesh.sh [NODES] [PARTS] : the 300 um cube of the Scriven case as a
# structured mesh in meshes/stNNN, divided into PARTS subdomains.
#
# NODES is the node count per edge of bubble.dom (default 75, which
# gives 74^3 = 405 224 cells: the "75^3" mesh of the paper); PARTS defaults to 8.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
N=${1:-75}; NP=${2:-8}
NAME=$(printf 'st%03d' "$N")
mkdir -p "$HERE/meshes/$NAME"; cd "$HERE/meshes/$NAME"
sed -E "s/^([[:space:]]*1[[:space:]]+)75([[:space:]]+)75([[:space:]]+)75[[:space:]]*$/\1$N\2$N\3$N/" \
    "$HERE/bubble.dom" > bubble.dom
grep -qE "^[[:space:]]*1[[:space:]]+$N[[:space:]]+$N[[:space:]]+$N" bubble.dom ||
  { echo "block line not found in bubble.dom"; exit 1; }
printf 'bubble\nskip\nskip\nskip\nskip\n' | "$ROOT/Binaries/Generate" > gen.log 2>&1
if [ "$NP" -gt 1 ]; then
  printf 'bubble\n%s\n' "$NP" | "$ROOT/Binaries/Divide" > div.log 2>&1
fi
find . -name "*.vtu" -type f -delete; find . -name "*.pvtu" -type f -delete
echo "meshes/$NAME: $(grep -i -m1 'number of cells' gen.log | tr -s ' ')"
