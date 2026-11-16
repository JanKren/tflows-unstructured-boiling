#!/bin/bash
# Build the pre-processors Generate, Divide and Convert (once per tree).
#
#   Cases/build_tools.sh
#
# The binaries land in Binaries/.  PETSc is not used; the variables that the
# makefiles would otherwise pick up from the environment are cleared.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
for p in Generate Divide Convert; do
  ( cd "$ROOT/Sources/$p" &&
    env -u FORTRAN -u PETSC_ARCH PETSC_DIR= make FORTRAN=gnu \
      > "$ROOT/Cases/build_$p.log" 2>&1 ) ||
    { echo "$p: build failed, see Cases/build_$p.log"; exit 1; }
  echo "$p: ok"
done
