#!/bin/bash
#SBATCH --job-name=d3d_poly64
#SBATCH --partition=daily --time=1-00:00:00
#SBATCH --nodes=1 --ntasks=32 --cpus-per-task=1 --mem=64G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/13.4.0 openmpi/5.0.8
ulimit -s unlimited
srun -n 32 ./Process
