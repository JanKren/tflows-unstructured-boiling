#!/bin/bash
#SBATCH --job-name=sd_build2
#SBATCH --partition=hourly
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=build-%j.out
#SBATCH --error=build-%j.err
# Builds T-Flows-jan-review (source.tar, one commit) with the static-droplet
# hooks, and Divide.  Serial make.  Refuses to overwrite a tree.
set -euo pipefail
module purge
module load gcc/14.2.0 openmpi/5.0.8
unset FORTRAN PETSC_DIR PETSC_ARCH
ulimit -s unlimited
ROOT=/data/user/user/cases/caf_static_drop_2026-11
cd "$ROOT"
( cd stage && sha256sum -c source.tar.sha256.remote )
[ ! -e build2 ] || { echo "build2 exists, refusing"; exit 2; }
mkdir build2
tar xf stage/source.tar -C build2 --strip-components=1
mkdir -p build2/Binaries/Cpu build2/case/User_Mod
cp stage/hooks/User_Mod/*.f90 build2/case/User_Mod/
( cd build2/Sources/Process/Cpu && make FORTRAN=gnu FCOMP=mpif90 MPI=yes OMP=no DIR_CASE="$ROOT/build2/case" ) > build2_process.log 2>&1
test -x build2/Binaries/Cpu/Process
( cd build2/Sources/Divide && make FORTRAN=gnu FCOMP=gfortran ) > build2_divide.log 2>&1
test -x build2/Binaries/Divide
sha256sum build2/Binaries/Cpu/Process build2/Binaries/Divide > BUILD2_SHA256
strings build2/Binaries/Cpu/Process | grep -c "EXACT_CURVATURE_TEST" > BUILD2_HOOK_PRESENT
echo BUILD_DONE
