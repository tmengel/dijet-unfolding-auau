#!/bin/bash
source setup_env.sh

for cent in 0 1 2 3; do
  RUN_UNCERTAINTIES=0 \
  RUN_SYSTEMATICS=1 \
  RUN_FLOW_SYS=1 \
    ./run_all_unfold_AA.sh 3 "$cent" \
    > "logs/rerun_systematics_cent${cent}.log" 2>&1
done
wait

for cent in 0 1 2 3; do
  ./run_everything_AA.sh 3 "$cent"
done