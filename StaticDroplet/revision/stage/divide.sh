#!/bin/bash
# sbatch --export=ALL,MESH=<hex64|poly64>,NAME=<problem>,PARTS=<n> divide.sh
#SBATCH --partition=hourly
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --job-name=sd_divide
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
module purge
module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
ROOT=/data/user/user/cases/caf_static_drop_2026-11
D="$ROOT/meshes/${MESH}_p${PARTS}"
[ ! -e "$D" ] || { echo "$D exists, refusing"; exit 2; }
mkdir -p "$D"
cd "$D"
ln "$ROOT/stage/meshes/$MESH/$NAME.cfn" "$ROOT/stage/meshes/$MESH/$NAME.dim" .
cp "$ROOT/build/Binaries/Divide" ./Divide
printf '%s\n%s\n' "$NAME" "$PARTS" | ./Divide > divide.log 2>&1
test "$(ls Sub | wc -l)" -eq "$PARTS"
echo DIVIDE_DONE
