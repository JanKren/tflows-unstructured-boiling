#!/bin/bash
set -e
B=/data/scratch/shared/user/ScrivenPolyDfix
NEW=/data/user/user/T-Flows-extstencil/Binaries/Cpu/Process
declare -A WALL=( [75]="02:00:00" [100]="04:00:00" [125]="06:00:00" [150]="08:00:00" )
for L in 75 100 125 150; do
  d=$B/run-$L; cd $d
  cp -f $NEW ./Process
  rm -f bench-data.dat *.out *.pvtu verify.sub
  find Sub -name "*-ts*.vtu" -delete 2>/dev/null || true
  rm -f *.backup
  # restore the full 600-step control if the verify run shortened it
  [ -f control.orig ] && cp -f control.orig control
  sed -i -e "s#^\( *#SBATCH --job-name=\).*#\1dfix-poly$L#" run.sh 2>/dev/null || true
  cat > run.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dfix-poly$L
#SBATCH --partition=daily --time=${WALL[$L]}
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=200G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
EOF
  printf "run-%-4s steps=%s  md5=%s  " "$L" "$(grep -i NUMBER_OF_TIME_STEPS control | awk '{print $2}')" "$(md5sum Process | cut -c1-8)"
  sbatch run.sh
done
