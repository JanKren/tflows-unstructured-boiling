#!/bin/bash
# run_case.sh NAME [RANKS] : runs one case with ./Process (built by
# ../build_case.sh .) and removes its volume output and backup; the front
# files and the records are kept.  RANKS must be the number of subdomains.
HERE=$(cd "$(dirname "$0")" && pwd)
NAME=$1; NP=${2:-8}
cd "$HERE/$NAME" || exit 1
ulimit -s unlimited
rm -f bench-data.dat kin-flux.dat kin-final.dat kin-init.dat
if [ "$NP" -gt 1 ]; then
  mpirun -np "$NP" "$HERE/Process" > run.log 2>&1
else
  "$HERE/Process" > run.log 2>&1
fi
echo "exit $?" > run.exit
if [ ! -f keep_volume ]; then
  find . -name "bubble-ts??????.vtu" -type f -delete
  find . -maxdepth 1 -name "bubble-ts??????.pvtu" -type f -delete
fi
find . -maxdepth 1 -name "bubble-ts??????.backup" -type f -delete
grep "VOF_CORR" run.log | awk '{print $5, $10, $11}' | sort | uniq -c > vofcorr.summary
cat run.exit
