#!/bin/bash
# setup_expansion.sh NAME MESHDIR DT NSTEPS NSAVE [INNER] [INIT] [EXTRA]
#
# Creates the run folder NAME with links to the mesh in MESHDIR (made by
# make_mesh.sh) and its decomposition; the mesh is never copied.
#
#   INNER  0  liquid velocity continued through the interface (solenoidal)
#          1  no flux inside the exact sphere R_S(t)
#          2  flux scaled by the mean volume fraction of the two cells
#          3  flux kept only where that mean exceeds 1/2 (the velocity jump)
#   INIT   0  volume fraction of the production initialisation
#          1  sub-sampled hexahedra
#   EXTRA  one more line for the control file, for example "MAX_COURANT_VOF 1.0"
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
NAME=$1; MESH=$(cd "$2" && pwd); DT=$3; NSTEPS=$4; NSAVE=$5
INNER=${6:-0}; INIT=${7:-0}; EXTRA=${8:-}
D=$HERE/$NAME
mkdir -p "$D"; cd "$D"
ln -sf "$MESH/bubble.cfn" .; ln -sf "$MESH/bubble.dim" .
for s in "$MESH"/Sub/*; do
  [ -d "$s" ] || continue
  n=$(basename "$s"); mkdir -p "Sub/$n"
  ln -sf "$s/bubble.cfn" "Sub/$n/"; ln -sf "$s/bubble.dim" "Sub/$n/"
done
sed -e "s|@PROB@|bubble|" -e "s|@NSTEPS@|$NSTEPS|" -e "s|@NSAVE@|$NSAVE|" \
    -e "s|@DT@|$DT|" -e "s|@EXTRA@|$EXTRA|" "$HERE/control.expansion" > control
sed -e "s|@INNER@|$INNER|" -e "s|@INIT@|$INIT|" "$HERE/kinematic.in.expansion" > kinematic.in
echo "set up $D ($(ls Sub 2>/dev/null | wc -l) subdomains)"
