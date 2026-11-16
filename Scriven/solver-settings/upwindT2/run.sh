#!/bin/bash
# The upwind-energy arm of the 40^3 Scriven settings check, with the keyword
# replaced in place (the arm "upwindT" kept minmod: its control has the keyword twice).
cd /home/user/runs/caf-review-2026-11/scriven-settings/upwindT2
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi/lib64:$LD_LIBRARY_PATH
printf "bubble\nskip\nskip\nskip\nskip\n" | /home/user/src/T-Flows-trees/caf-stefan-paper/Binaries/Generate > generate.log 2>&1
ulimit -s unlimited
./Process > run.log 2>&1
echo "exit $?"
wc -l bench-data.dat
