#!/bin/bash
#SBATCH --job-name=st75-gamma-off
#SBATCH --partition=hourly --time=00:20:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=100G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
