#!/bin/bash
# Build Process against one case's User_Mod, safely.
#
#   ./build_case.sh <case_dir>
#
# The manual sequence is easy to get wrong and did go wrong once here: hooks
# that are #included from Vof_Mod.f90 rather than User_Mod.f90 (notably
# End_Of_Compute_Vof.f90) are NOT tracked by the makefile's dependency list, so
# `make` leaves a stale Vof_Mod.o and silently carries the PREVIOUS case's hook
# into the new binary.  Stefan and Sucking were both run once with the static
# droplet's y-flattening hook compiled in before this was caught.
#
# The script therefore always: clears symlinks, restores every pristine default,
# touches Vof_Mod.f90 to force its rebuild, builds, and then VERIFIES the binary
# rather than trusting the make exit code.
set -e

CASE=${1:?usage: build_case.sh <case_dir>}
SRC=/home/user/src/T-Flows-trees/T-Flows-dfix/Sources/Process/Cpu
U=$SRC/User_Mod
PRISTINE=/home/user/runs/StaticDropletSym/.user_mod_pristine
BIN=/home/user/src/T-Flows-trees/T-Flows-dfix/Binaries/Cpu/Process

echo "== restoring pristine User_Mod =="
# The makefile may leave User_Mod as a SYMLINK to the case directory; a bare
# `find -type l -delete` then deletes the link itself and the tree loses the
# directory entirely.  Rebuild it as a real directory every time.
[ -L "$U" ] && rm -f "$U"
mkdir -p "$U"
find "$U" -maxdepth 1 -type l -delete
for f in "$PRISTINE"/*.f90; do cp -f "$f" "$U/$(basename "$f")"; done

echo "== forcing rebuild of the modules that #include user hooks =="
touch "$SRC/Vof_Mod.f90" "$SRC/User_Mod.f90"
touch "$CASE"/User_Mod/*.f90          # defeat the ln -s mtime trap

echo "== building against $CASE =="
cd "$SRC"
env -u FORTRAN -u PETSC_DIR -u PETSC_ARCH \
  make FORTRAN=gnu MPI=yes DIR_CASE="$CASE" 2>&1 | grep -iE "Error:|error #|Linking " | tail -4

echo "== verifying the binary =="
if strings "$BIN" | grep -q "VOF_FLATTEN"; then
  echo "  FAIL: binary still carries the static-droplet flattening hook"; exit 1
fi
strings "$BIN" | grep -q "mass-balance.dat" \
  && echo "  ok: mass-balance writer present" \
  || { echo "  FAIL: mass-balance writer absent"; exit 1; }

cp "$BIN" "$CASE/Process"
echo "  binary copied to $CASE/Process (private, so later builds cannot disturb a running job)"
