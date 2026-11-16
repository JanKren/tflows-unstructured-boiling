#!/bin/bash
# queue3.sh WAITFOR : waits for case WAITFOR, then runs the cases listed in
# queue.list one after another (first line first; the list may be edited
# while the queue runs).  Lines:  exp NAME   |   adv NAME NP [keep]
# Stops when the list is empty.
ROOT=/home/user/runs/caf-review-2026-11/kinematic
W=$1
while [ ! -f $ROOT/$W/run.exit ]; do sleep 10; done
echo "$(date +%H:%M:%S) done  $W  $(cat $ROOT/$W/run.exit)" >> $ROOT/queue.log
while true; do
  line=$(head -1 $ROOT/queue.list)
  [ -z "$line" ] && break
  sed -i '1d' $ROOT/queue.list
  set -- $line
  echo "$(date +%H:%M:%S) start $line" >> $ROOT/queue.log
  if [ "$1" = exp ]; then
    $ROOT/run_case.sh $2 8
    echo "$(date +%H:%M:%S) done  $2  $(cat $ROOT/$2/run.exit)" >> $ROOT/queue.log
  else
    $ROOT/adv/run_adv.sh $2 $3 $4
    echo "$(date +%H:%M:%S) done  $2  $(cat $ROOT/adv/$2/run.exit)" >> $ROOT/queue.log
  fi
done
echo "$(date +%H:%M:%S) QUEUE FINISHED" >> $ROOT/queue.log
