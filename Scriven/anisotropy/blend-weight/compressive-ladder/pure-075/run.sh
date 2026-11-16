#!/bin/bash
#SBATCH --job-name=pure-075
#SBATCH --partition=daily --time=03:00:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=200G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
