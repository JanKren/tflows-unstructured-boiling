#!/bin/bash
# Coarse preview (R/h = 6.4, four cycles on both families), one hourly job each.
set -euo pipefail
ROOT=/data/user/user/cases/caf_static_drop_2026-11
S=$ROOT/stage
cd $ROOT/logs
while read c mesh; do
  bash $S/stage_case.sh $c $mesh 16 > /dev/null
  sbatch --job-name=sd_$c --partition=hourly --time=01:00:00 --ntasks=16 --mem=24G \
         --export=ALL,CASE_DIR=$ROOT/runs/$c $S/run_case.sh
done <<LIST
c32m_hex_s4 hex32
c32m_poly_s4 poly32
c32h_hex_s4 hex32
c32h_poly_s4 poly32
LIST
