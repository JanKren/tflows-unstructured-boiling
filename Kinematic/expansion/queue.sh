#!/bin/bash
# queue.sh NAME... : runs the listed cases one after another on 8 ranks
ROOT=/home/user/runs/caf-review-2026-11/kinematic
for n in "$@"; do
  echo "$(date +%H:%M:%S) start $n" >> $ROOT/queue.log
  $ROOT/run_case.sh $n 8
  echo "$(date +%H:%M:%S) done  $n  $(cat $ROOT/$n/run.exit)" >> $ROOT/queue.log
done
