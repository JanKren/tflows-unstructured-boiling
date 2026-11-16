#!/bin/bash
# Build Process for one case: the solver with that case's User_Mod compiled in.
#
#   Cases/build_case.sh <case-dir> [output-binary]
#
# The binary is copied to <case-dir>/Process unless another path is given, so
# that a later build for another case cannot change it.
#
# T-Flows compiles the case's hooks (User_Mod/*.f90) into the solver.  Its
# makefile links them over the default hooks in Sources/Process/Cpu/User_Mod,
# and two things then go wrong when the case changes in the same tree:
# hooks of the previous case that the new one lacks stay linked, and a link
# keeps the date of its target, so make can take an old object for current.
# This script therefore restores every default hook from Cases/default_hooks,
# removes what is not a default hook, and touches the case's files and the
# modules that include hooks, before it builds.  Build serially: a parallel
# make races on the Fortran module files.
set -e
CASE=$(cd "${1:?case directory}" && pwd)
OUT=${2:-$CASE/Process}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/Sources/Process/Cpu
U=$SRC/User_Mod
DEFAULTS=$ROOT/Cases/default_hooks
[ -d "$CASE/User_Mod" ] || { echo "no User_Mod in $CASE"; exit 1; }

find "$U" -maxdepth 1 -type l -delete
for f in "$DEFAULTS"/*.f90; do cp -f "$f" "$U/$(basename "$f")"; done
for f in "$U"/*.f90; do
  [ -e "$DEFAULTS/$(basename "$f")" ] || rm -f "$f"
done
touch "$SRC/Vof_Mod.f90" "$SRC/User_Mod.f90" "$SRC/Process_Mod.f90"
touch "$CASE"/User_Mod/*.f90

LOG=$CASE/build_Process.log
( cd "$SRC" &&
  env -u FORTRAN -u PETSC_ARCH PETSC_DIR= \
    make FORTRAN=gnu MPI=yes DIR_CASE="$CASE" > "$LOG" 2>&1 ) ||
  { echo "build failed, see $LOG"; tail -20 "$LOG"; exit 1; }
cp "$ROOT/Binaries/Cpu/Process" "$OUT"
echo "ok: $OUT"
