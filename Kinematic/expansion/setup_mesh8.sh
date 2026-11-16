#!/bin/bash
# setup_mesh8.sh NAME GLOBALDIR PROBLEM : an 8-way decomposition of an
# existing global mesh, in meshes/NAME (the global files are symlinked).
set -e
ROOT=/home/user/runs/caf-review-2026-11/kinematic
NAME=$1; SRC=$2; PROB=$3
mkdir -p $ROOT/meshes/$NAME; cd $ROOT/meshes/$NAME
ln -sf $SRC/$PROB.cfn .; ln -sf $SRC/$PROB.dim .
printf '%s\n8\n' "$PROB" | /home/user/src/T-Flows-trees/caf-kinematic/Binaries/Divide > div.log 2>&1
ls Sub | tr '\n' ' '; echo; tail -2 div.log
