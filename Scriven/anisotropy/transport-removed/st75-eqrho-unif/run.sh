#!/bin/bash
#SBATCH --job-name=st75-eqrho-unif
#SBATCH --partition=hourly --time=01:00:00
#SBATCH --nodes=1 --ntasks=20 --cpus-per-task=1 --mem=100G
#SBATCH --hint=nomultithread
#SBATCH --output=%x-%j.out
module purge; module load gcc/14.2.0 openmpi/5.0.8
ulimit -s unlimited

TARGET=1000
MAX_ATTEMPTS=8
before=$(wc -l < bench-data.dat 2>/dev/null || echo 0)
attempt=$(cat .attempt 2>/dev/null || echo 0); attempt=$((attempt+1))
echo "$attempt" > .attempt
echo "=== attempt $attempt, starting from $before steps ==="

# --- restart wiring: point the control at the newest backup, if any ---------
bk=$(ls -1 *-ts*.backup 2>/dev/null | sed 's/.*-ts0*//;s/\.backup//' | sort -n | tail -1)
if [ -n "$bk" ] && [ "$bk" -gt 0 ] 2>/dev/null; then
  bkfile=$(ls -1 *-ts$(printf '%06d' "$bk").backup)
  grep -v '^[[:space:]]*LOAD_BACKUP_NAME' control > control.next
  printf '\n LOAD_BACKUP_NAME   %s\n' "$bkfile" >> control.next
  mv control.next control
  # bench-data appends, so drop any rows past the backup instant
  dt=$(awk '/^[[:space:]]*TIME_STEP[[:space:]]/{print $2}' control)
  if [ -f bench-data.dat ]; then
    cp -n bench-data.dat bench-data-orig.dat 2>/dev/null || true
    awk -v tc="$(awk -v s="$bk" -v d="$dt" 'BEGIN{printf "%.12e", s*d}')" \
        '$1 <= tc*1.000001' bench-data.dat > bench-data.trim && mv bench-data.trim bench-data.dat
  fi
  echo "restarting from $bkfile (step $bk)"
fi

srun ./Process

# --- chain: resubmit only while progress is being made ----------------------
after=$(wc -l < bench-data.dat 2>/dev/null || echo 0)
echo "=== reached $after of $TARGET steps ==="
if [ "$after" -ge "$TARGET" ]; then
  echo "TARGET REACHED"
elif [ "$attempt" -ge "$MAX_ATTEMPTS" ]; then
  echo "STOPPING: attempt cap reached at $after steps"
elif [ "$after" -le "$before" ]; then
  echo "STOPPING: no progress this attempt ($before -> $after)"
else
  echo "resubmitting"; sbatch run.sh
fi
