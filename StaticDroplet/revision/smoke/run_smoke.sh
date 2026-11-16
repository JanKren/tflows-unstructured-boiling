#!/bin/bash
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
ulimit -s unlimited
P=/home/user/runs/caf-review-2026-11/static-drop/src-local/Binaries/Cpu/Process
cd /home/user/runs/caf-review-2026-11/static-drop/smoke
( cd ser_m_hex_s12 && $P > solver.out 2> solver.err; echo "ser_m exit $?" )
( cd ser_h_hex_s12 && $P > solver.out 2> solver.err; echo "ser_h exit $?" )
( cd par2_m_hex_s12 && mpirun -np 2 --bind-to none $P > solver.out 2> solver.err; echo "par2 exit $?" )
echo SMOKE_DONE
