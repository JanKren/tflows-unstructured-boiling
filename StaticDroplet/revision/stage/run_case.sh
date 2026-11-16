#!/bin/bash
# sbatch --job-name=<case> --partition=<p> --time=<t> --ntasks=<n> --export=ALL,CASE_DIR=<dir> run_case.sh
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
module purge
module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 GFORTRAN_UNBUFFERED_PRECONNECTED=y
cd "$CASE_DIR"
RANKS=$(cat RANKS)
[ ! -e EXIT_STATUS ] || { echo "this case has run already: $CASE_DIR"; exit 2; }
echo "start $(date -u +%FT%TZ) job $SLURM_JOB_ID nodes $SLURM_JOB_NODELIST ranks $RANKS process $(sha256sum Process | cut -c1-16)" > run.info
set +e
srun --ntasks="$RANKS" --cpu-bind=cores ./Process > solver.out 2> solver.err
code=$?
set -e
echo "$code" > EXIT_STATUS
echo "end $(date -u +%FT%TZ) exit $code" >> run.info
exit $code
