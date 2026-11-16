#!/bin/bash
# poly128 is declared finished at its current state, but RESULTS_SAVE_INTERVAL
# only fires at the end, so no field output exists to plot.  Load the newest
# complete backup and advance ONE step so the solver writes the .pvtu.
# One step at dt=1.27e-6 is physically inert; this is an extraction, not a run.
cd "$(dirname "$0")"
B=$(find . -maxdepth 1 -name "*.backup" -mmin +2 -printf "%T@ %f\n" | sort -n | tail -1 | cut -d" " -f2)
[ -z "$B" ] && { echo "no complete backup"; exit 1; }
cp control control.orig
sed -i "/^ *LOAD_BACKUP_NAME/d" control
sed -i "s/^ *NUMBER_OF_TIME_STEPS.*/ NUMBER_OF_TIME_STEPS  1/"  control
sed -i "s/^ *RESULTS_SAVE_INTERVAL.*/ RESULTS_SAVE_INTERVAL 1/" control
sed -i "/^ *BACKUP_SAVE_INTERVAL/a\\ LOAD_BACKUP_NAME      $B" control
echo "dumping field from $B"
sbatch --job-name=poly128_dump --partition=hourly --time=00:20:00 \
       --dependency=afterany:8075530 run.sh
