#!/bin/bash
set -e
B=/data/scratch/shared/user/ScrivenPolyDfix
NEW=/data/user/user/T-Flows-extstencil/Binaries/Cpu/Process
declare -A WALL=( [75]="02:00:00" [100]="06:00:00" [125]="12:00:00" [150]="18:00:00" )
# keep the completed 20-rank poly-75 trace as a decomposition-sensitivity check
[ -s $B/run-75/bench-data.dat ] && cp -n $B/run-75/bench-data.dat $B/run-75/bench-data-20rank.dat
for L in 75 100 125 150; do
  d=$B/run-$L; cd $d
  cp -f $NEW ./Process
  rm -f bench-data.dat *.out *.pvtu *.backup
  find Sub -name "*-ts*.vtu" -delete 2>/dev/null || true
  [ -f control.orig ] && cp -f control.orig control
  sed -i -e "s/^\( *BACKUP_SAVE_INTERVAL *\).*/\1 100/" control
  cat > run.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dfix64-poly$L
#SBATCH --partition=daily --time=${WALL[$L]}
#SBATCH --nodes=1 --ntasks=64 --cpus-per-task=1 --mem=400G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
EOF
  printf "run-%-4s Sub=%s steps=%s backup=%s  " "$L" "$(ls Sub|wc -l)" \
    "$(grep -i NUMBER_OF_TIME_STEPS control|awk '{print $2}')" \
    "$(grep -i BACKUP_SAVE_INTERVAL control|awk '{print $2}')"
  sbatch run.sh
done
