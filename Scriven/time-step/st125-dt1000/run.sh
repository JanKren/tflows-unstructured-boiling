#!/bin/bash
#SBATCH --job-name=st125-dt1000
#SBATCH --partition=daily --time=06:00:00
#SBATCH --nodes=1 --ntasks=80 --cpus-per-task=1 --mem=300G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
