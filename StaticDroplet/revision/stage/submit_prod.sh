#!/bin/bash
# submit_prod.sh <case> <hex64|poly64> <ranks> <hh:mm:ss>
# One daily job per case, sized to finish inside its wall; the solver stops by
# itself at WALL_TIME_MAX_HOURS (23.3 h) with a backup if it does not.
set -euo pipefail
ROOT=/data/user/user/cases/caf_static_drop_2026-11
S=$ROOT/stage
c=$1; mesh=$2; parts=$3; wall=$4
mkdir -p $ROOT/runs $ROOT/logs
bash $S/stage_case.sh $c $mesh $parts > /dev/null
cd $ROOT/logs
sbatch --job-name=sd_$c --partition=daily --time=$wall --ntasks=$parts --mem=64G \
       --export=ALL,CASE_DIR=$ROOT/runs/$c $S/run_case.sh
