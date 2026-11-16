#!/bin/bash
# Pilots on hourly: 400 steps each, to measure seconds per step and to compare
# the first rows of the matched hexahedral case with the published trace.
set -euo pipefail
ROOT=/data/user/user/cases/caf_static_drop_2026-11
S=$ROOT/stage
mkdir -p $ROOT/pilots $ROOT/logs
cd $ROOT/logs
while read c mesh parts; do
  bash $S/stage_case.sh $c $mesh $parts 400 > /dev/null
  sbatch --job-name=sdp_${c}_p${parts} --partition=hourly --time=00:20:00 \
         --ntasks=$parts --mem=48G \
         --export=ALL,CASE_DIR=$ROOT/pilots/${c}_p${parts} $S/run_case.sh
done <<LIST
m_hex_s12 hex64 32
m_hex_s4 hex64 32
m_poly_s4 poly64 32
h_hex_s4 hex64 32
h_hex_s12 hex64 32
h_poly_s4 poly64 32
h_poly_s4 poly64 64
LIST
