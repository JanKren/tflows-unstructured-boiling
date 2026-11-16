#!/bin/bash
# setup_expansion.sh NAME MESHDIR PROBLEM DT NSTEPS NSAVE [INNER] [INIT] [EXTRA]
# Creates run directory NAME with symlinks to an existing mesh and its
# decomposition (never copies the mesh).
set -e
ROOT=/home/user/runs/caf-review-2026-11/kinematic
NAME=$1; MESH=$2; PROB=$3; DT=$4; NSTEPS=$5; NSAVE=$6
INNER=${7:-0}; INIT=${8:-0}; EXTRA=${9:-}
D=$ROOT/$NAME
mkdir -p $D; cd $D
ln -sf $MESH/$PROB.cfn .; ln -sf $MESH/$PROB.dim .
for s in $MESH/Sub/*; do
  n=$(basename $s); mkdir -p Sub/$n
  ln -sf $s/$PROB.cfn Sub/$n/; ln -sf $s/$PROB.dim Sub/$n/
done
sed -e "s|@PROB@|$PROB|" -e "s|@NSTEPS@|$NSTEPS|" -e "s|@NSAVE@|$NSAVE|" \
    -e "s|@DT@|$DT|" \
    -e "s|@EXTRA@|$EXTRA|" $ROOT/control.expansion > control
sed -e "s|@INNER@|$INNER|" -e "s|@INIT@|$INIT|" $ROOT/kinematic.in.expansion > kinematic.in
echo "set up $D ($(ls Sub | wc -l) subdomains)"
