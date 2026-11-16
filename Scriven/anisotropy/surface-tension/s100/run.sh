#!/bin/bash
#SBATCH --job-name=st-s100-ext
#SBATCH --partition=general --time=26:00:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=100G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
