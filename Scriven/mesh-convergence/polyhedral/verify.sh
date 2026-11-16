#!/bin/bash
set -e
B=/data/scratch/shared/user/ScrivenPolyDfix
NEW=/data/user/user/T-Flows-extstencil/Binaries/Cpu/Process
cd $B/run-75

# keep the decomposed mesh (Sub/*/*.cfn,*.dim); drop only previous OUTPUT
rm -f bench-data.dat *.pvtu *.out
find Sub -name "*-ts*.vtu" -delete 2>/dev/null || true

cp -f $NEW ./Process
[ -f control.orig ] || cp control control.orig
sed -e 's/^\( *NUMBER_OF_TIME_STEPS *\).*/\1 10/' \
    -e 's/^\( *RESULTS_SAVE_INTERVAL *\).*/\1 10/' control.orig > control

cat > verify.sub <<'EOF'
#!/bin/bash
#SBATCH --job-name=verify-poly75
#SBATCH --partition=hourly --time=00:20:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=200G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited
srun ./Process
EOF
sbatch verify.sub
