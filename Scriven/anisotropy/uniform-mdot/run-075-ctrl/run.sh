#!/bin/bash
#SBATCH --job-name=u2-075-ctrl
#SBATCH --partition=daily --time=01:00:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=400G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
