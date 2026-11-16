#!/bin/bash
# setup_adv.sh NAME MESHDIR PROB DT NSTEPS NSAVE COMAX "mode" "cx cy cz R" "p1 p2 p3 p4" INIT CYL
set -e
ROOT=/home/user/runs/caf-review-2026-11/kinematic/adv
NAME=$1; MESH=$2; PROB=$3; DT=$4; NSTEPS=$5; NSAVE=$6; COMAX=$7
D=$ROOT/$NAME; mkdir -p $D; cd $D
ln -sf $MESH/$PROB.cfn .; ln -sf $MESH/$PROB.dim .
if [ -d $MESH/Sub ]; then
  for s in $MESH/Sub/*; do n=$(basename $s); mkdir -p Sub/$n
    ln -sf $s/$PROB.cfn Sub/$n/; ln -sf $s/$PROB.dim Sub/$n/; done
fi
sed -e "s|@PROB@|$PROB|" -e "s|@NSTEPS@|$NSTEPS|" -e "s|@NSAVE@|$NSAVE|" \
    -e "s|@DT@|$DT|" -e "s|@COMAX@|$COMAX|" $ROOT/control.adv > control
cat > kinematic.in << KIN
$8                    mode: 1 expansion, 2 translation, 3 deformation
$9                    centre and radius
${10}                 mode 2: Ux Uy Uz t_reverse | mode 3: T
0                     inner (mode 1 only)
${11}                 init: 0 linear estimate, 1 sub-sampled hexahedra
${12}                 1: cylinder along z
KIN
