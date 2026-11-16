#!/bin/bash
# continue_case.sh <case> <partition> <hh:mm:ss> [wall_hours]
# Resume a run that stopped on its wall-time limit from its last backup.
# The traces are appended (the hook deletes them only at step 1).
set -euo pipefail
ROOT=/data/user/user/cases/caf_static_drop_2026-11
c=$1; part=$2; wall=$3; wh=${4:-}
d=$ROOT/runs/$c
cd "$d"
last=$(ls -1 *.backup | sort | tail -1)
[ -n "$last" ] || { echo "no backup in $d"; exit 2; }
n=$( (ls -d cont_* 2>/dev/null || true) | wc -l)
mkdir cont_$n && mv EXIT_STATUS run.info cont_$n/ 2>/dev/null || true
cp control cont_$n/control
grep -v "LOAD_BACKUP_NAME" cont_$n/control > control
echo " LOAD_BACKUP_NAME        $last" >> control
[ -z "$wh" ] || sed -i "s/^\( *WALL_TIME_MAX_HOURS\) .*/\1 $wh/" control
echo "resuming $c from $last ($(wc -l < spurious.dat) trace rows)"
cd $ROOT/logs
sbatch --job-name=sd_$c --partition=$part --time=$wall --ntasks=$(cat $d/RANKS) --mem=64G \
       --export=ALL,CASE_DIR=$d $ROOT/stage/run_case.sh
