#!/bin/bash
# run_case.sh NAME NP : runs one case and removes the volume output and the
# final backup of that case (the front files and the records are kept).
ROOT=/home/user/runs/caf-review-2026-11/kinematic
NAME=$1; NP=${2:-8}
cd $ROOT/$NAME || exit 1
PROB=$(awk '/PROBLEM_NAME/{print $2}' control)
ulimit -s unlimited
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
rm -f bench-data.dat kin-flux.dat kin-final.dat kin-init.dat
if [ "$NP" -gt 1 ]; then
  mpirun -np $NP $ROOT/bin/Process > run.log 2>&1
else
  $ROOT/bin/Process > run.log 2>&1
fi
echo "exit $?" > run.exit
# volume snapshots and backups of this case only
if [ ! -f keep_volume ]; then
  find $ROOT/$NAME -name "${PROB}-ts??????.vtu"    -type f -delete
  find $ROOT/$NAME -maxdepth 1 -name "${PROB}-ts??????.pvtu"   -type f -delete
fi
find $ROOT/$NAME -maxdepth 1 -name "${PROB}-ts??????.backup" -type f -delete
grep -c "VOF_CORR" run.log > vofcorr.count
grep "VOF_CORR" run.log | awk '{print $5, $10, $11}' | sort | uniq -c > vofcorr.summary
