#!/bin/bash
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
ulimit -s unlimited
P=/home/user/runs/caf-review-2026-11/static-drop/src-local/Binaries/Cpu/Process
cd /home/user/runs/caf-review-2026-11/static-drop/smoke
for m in 1 0; do ( cd exp$m && $P > solver.out 2> solver.err; echo "exp$m exit $?" ); done
echo EXP_DONE
