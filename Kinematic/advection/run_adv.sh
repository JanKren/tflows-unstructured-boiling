#!/bin/bash
# run_adv.sh NAME NP [keep] : runs one advection case; volume output is
# removed afterwards unless "keep" is given
ROOT=/home/user/runs/caf-review-2026-11/kinematic/adv
NAME=$1; NP=${2:-1}; KEEP=${3:-no}
cd $ROOT/$NAME || exit 1
PROB=$(awk '/PROBLEM_NAME/{print $2}' control)
ulimit -s unlimited
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
rm -f bench-data.dat kin-flux.dat kin-final.dat kin-init.dat alert-dt-vof.dat
if [ "$NP" -gt 1 ]; then
  mpirun -np $NP $ROOT/../bin/Process_adv > run.log 2>&1
else
  $ROOT/../bin/Process_adv > run.log 2>&1
fi
echo "exit $?" > run.exit
grep "VOF_CORR" run.log | awk '{print $5, $10, $11}' | sort | uniq -c > vofcorr.summary
grep "VOF_CORR" run.log | awk '{split($10,a,"=");split($11,b,"=");lo+=a[2];hi+=b[2]} END{print "sum unbounded lo",lo,"hi",hi}' >> vofcorr.summary
find $ROOT/$NAME -maxdepth 1 -name "${PROB}-ts*.backup" -type f -delete
if [ "$KEEP" != "keep" ]; then
  find $ROOT/$NAME -name "${PROB}-ts*.vtu"  -type f -delete
  find $ROOT/$NAME -maxdepth 1 -name "${PROB}-ts*.pvtu" -type f -delete
fi
