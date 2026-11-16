#!/bin/bash
# Coarse Scriven (40^3, h = 7.5 um), paper solver, one control change per run.
#   ./setup_run.sh <name> <dt> <nsteps> <save_interval> [KEYWORD=value ...]
# Runs serially in the background; prints the PID.
set -e
TPL=/home/user/runs/caf-review-2026-11/stefan/scriven_tpl
BIN=/home/user/runs/caf-review-2026-11/stefan/Process_paper_scriven
GEN=/home/user/src/T-Flows-trees/caf-stefan-paper/Binaries/Generate
name=$1; dt=$2; n=$3; save=$4; shift 4
d=/home/user/runs/caf-review-2026-11/scriven-settings/$name
mkdir -p $d && cd $d
cp $TPL/bubble.dom $TPL/control $TPL/dT1.25_R50microns.txt $TPL/ellipsoid_parameters.ini .
sed -i 's/^  1  75   75   75/  1  41   41   41/' bubble.dom
grep -q "  1  41   41   41" bubble.dom
sed -i "s/^ TIME_STEP .*/ TIME_STEP               $dt/; s/^ NUMBER_OF_TIME_STEPS .*/ NUMBER_OF_TIME_STEPS   $n/; s/^ RESULTS_SAVE_INTERVAL .*/ RESULTS_SAVE_INTERVAL   $save/" control
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}
  if grep -q "^ *$k " control; then sed -i "s|^ *$k .*| $k    $v|" control; else echo " $k    $v" >> control; fi
done
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
printf "bubble\nskip\nskip\nskip\nskip\n" | $GEN > generate.log 2>&1
cp $BIN ./Process
ulimit -s unlimited
nohup ./Process > run.log 2>&1 &
echo "$name pid $!"
