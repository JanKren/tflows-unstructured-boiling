#!/bin/bash
# Build Process of one source tree against one case's User_Mod.
#
#   ./build_case.sh <tree> <case_dir> <out_binary> [marker string]
#
# <tree> is caf-stefan (T-Flows-jan-review) or caf-stefan-paper
# (extended-gradient-stencil) under ~/src/T-Flows-trees.  The ritual is the one
# of ~/runs/MassBalance/build_case.sh: clear the symlinks of the previous case,
# restore every pristine default hook, touch the modules that #include hooks,
# build serially, and copy the binary away so that a later build cannot change
# a running job.
set -e
TREE=${1:?tree}; CASE=$(readlink -f ${2:?case dir}); OUT=${3:?output binary}
MARK=${4:-caf-terms.dat}     # a string only this case's hooks put in the binary
HERE=$(dirname "$(readlink -f "$0")")
ROOT=/home/user/src/T-Flows-trees/$TREE
SRC=$ROOT/Sources/Process/Cpu
U=$SRC/User_Mod
PRISTINE=$HERE/.pristine_$TREE

[ -L "$U" ] && rm -f "$U"
mkdir -p "$U"
find "$U" -maxdepth 1 -type l -delete
for f in "$PRISTINE"/*.f90; do cp -f "$f" "$U/$(basename "$f")"; done
# files of a previous case that are not default hooks
for f in "$U"/*.f90; do
  [ -e "$PRISTINE/$(basename "$f")" ] || rm -f "$f"
done

touch "$SRC/Vof_Mod.f90" "$SRC/User_Mod.f90" "$SRC/Process_Mod.f90"
touch "$CASE"/User_Mod/*.f90

cd "$SRC"
env -u FORTRAN -u PETSC_ARCH PETSC_DIR= \
  make FORTRAN=gnu MPI=yes DIR_CASE="$CASE" > "$HERE/build_$TREE.log" 2>&1 \
  || { echo "BUILD FAILED, see $HERE/build_$TREE.log"; tail -30 "$HERE/build_$TREE.log"; exit 1; }
grep -iE "Error" "$HERE/build_$TREE.log" | head -5 || true

BIN=$ROOT/Binaries/Cpu/Process
[ -x "$BIN" ] || BIN=$ROOT/Binaries/Process
cp "$BIN" "$OUT"
strings "$OUT" | grep -q "$MARK" && echo "ok: $OUT carries $MARK" \
  || { echo "FAIL: $MARK absent from $OUT"; exit 1; }
