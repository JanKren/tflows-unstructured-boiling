#!/bin/bash
# Exact-curvature test: sbatch --job-name=sd_ex_<base> --export=ALL,BASE=<m_hex_s4|...>,MESH=<hex64|poly64>,MODE=<1|2>,STEPS=<n> run_exact.sh
# Stages runs/ex<MODE>_<BASE> from the control of BASE with EXACT_CURVATURE_TEST
# added and the binary of build2, then runs it on 32 ranks.
#SBATCH --partition=hourly
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
module purge
module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GFORTRAN_UNBUFFERED_PRECONNECTED=y
ROOT=/data/user/user/cases/caf_static_drop_2026-11
d="$ROOT/runs/ex${MODE}_${BASE}"
[ ! -e "$d" ] || { echo "$d exists, refusing"; exit 2; }
mkdir -p "$d"
cp "$ROOT/stage/controls/$BASE.control" "$d/control"
cp "$ROOT/stage/ellipsoid_parameters.ini" "$d/"
sed -i -e "s/^\( *NUMBER_OF_TIME_STEPS\) .*/\1 $STEPS/" \
       -e "s/^\( *RESULTS_SAVE_INTERVAL\) .*/\1 $STEPS/" \
       -e "s/^\( *BACKUP_SAVE_INTERVAL\) .*/\1 999999/" "$d/control"
echo " EXACT_CURVATURE_TEST    $MODE" >> "$d/control"
cp -al "$ROOT/meshes/${MESH}_p32/Sub" "$d/Sub"
cp "$ROOT/build2/Binaries/Cpu/Process" "$d/Process"
echo 32 > "$d/RANKS"
cd "$d"
echo "start $(date -u +%FT%TZ) job $SLURM_JOB_ID process $(sha256sum Process | cut -c1-16)" > run.info
set +e
srun --ntasks=32 --cpu-bind=cores ./Process > solver.out 2> solver.err
code=$?
set -e
echo "$code" > EXIT_STATUS
echo "end $(date -u +%FT%TZ) exit $code" >> run.info
exit $code
