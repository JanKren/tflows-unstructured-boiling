#!/bin/bash
#SBATCH --job-name=sd_build
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
[ ! -e build ] || { echo "build exists, refusing"; exit 2; }
mkdir build
tar xf stage/source.tar -C build --strip-components=1
mkdir -p build/Binaries/Cpu build/case/User_Mod
cp stage/hooks/User_Mod/*.f90 build/case/User_Mod/
( cd build/Sources/Process/Cpu && make FORTRAN=gnu FCOMP=mpif90 MPI=yes OMP=no DIR_CASE="$ROOT/build/case" ) > build_process.log 2>&1
test -x build/Binaries/Cpu/Process
( cd build/Sources/Divide && make FORTRAN=gnu FCOMP=gfortran ) > build_divide.log 2>&1
test -x build/Binaries/Divide
sha256sum build/Binaries/Cpu/Process build/Binaries/Divide > BUILD_SHA256
strings build/Binaries/Cpu/Process | grep -c "spurious_ext.dat" > BUILD_HOOK_PRESENT
echo BUILD_DONE
