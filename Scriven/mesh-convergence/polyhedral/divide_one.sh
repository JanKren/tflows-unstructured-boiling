#!/bin/bash
#SBATCH --job-name=div64
#SBATCH --partition=hourly --time=01:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=240G
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
cd $RUNDIR
rm -rf Sub
printf "%s\n64\n" "$MESH" | /data/user/user/T-Flows-extstencil/Binaries/Divide
echo "EXIT=$?"
ls Sub | wc -l
