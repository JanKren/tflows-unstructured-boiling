#!/bin/bash
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
ulimit -s unlimited
P=/home/user/runs/caf-review-2026-11/static-drop/src-local/Binaries/Cpu/Process
cd /home/user/runs/caf-review-2026-11/static-drop/smoke
for m in 0 1 2; do ( cd ex$m && $P > solver.out 2> solver.err; echo "ex$m exit $?" ); done
echo EX_DONE
