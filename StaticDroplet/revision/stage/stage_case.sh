#!/bin/bash
# stage_case.sh <case> <hex64|poly64> <parts> [steps]   -> runs/<case>[_pilot]
# Each case gets its own Sub (hard links of the mesh pieces only), so no case
# can overwrite another's fields.
set -euo pipefail
ROOT=/data/user/user/cases/caf_static_drop_2026-11
c=$1; mesh=$2; parts=$3; steps=${4:-}
d="$ROOT/runs/$c"; [ -z "$steps" ] || d="$ROOT/pilots/${c}_p${parts}"
[ ! -e "$d" ] || { echo "$d exists, refusing"; exit 2; }
mkdir -p "$d"
cp "$ROOT/stage/controls/$c.control" "$d/control"
cp "$ROOT/stage/ellipsoid_parameters.ini" "$ROOT/stage/sphere.stl" "$d/"
if [ -n "$steps" ]; then
  sed -i -e "s/^\( *NUMBER_OF_TIME_STEPS\) .*/\1 $steps/" \
         -e "s/^\( *RESULTS_SAVE_INTERVAL\) .*/\1 999999/" \
         -e "s/^\( *BACKUP_SAVE_INTERVAL\) .*/\1 999999/" "$d/control"
fi
cp -al "$ROOT/meshes/${mesh}_p${parts}/Sub" "$d/Sub"
cp "$ROOT/build/Binaries/Cpu/Process" "$d/Process"
echo "$parts" > "$d/RANKS"
grep -E "NUMBER_OF_TIME_STEPS|TIME_STEP|PHASE_|SMOOTHING_CYCLES_CURV|RESULTS_SAVE|BACKUP_SAVE|WALL_TIME" "$d/control"
echo "staged $d"
