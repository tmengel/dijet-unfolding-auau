#!/bin/bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"
mkdir -p logs

echo "Running AA unfolding and systematics for cone=3, cent=0-3"
for cent in 0 1 2 3; do
  RUN_UNCERTAINTIES=1 \
  RUN_SYSTEMATICS=1 \
  RUN_FLOW_SYS=1 \
    ./run_all_unfold_AA.sh 3 "$cent" \
    > "logs/unfold_cent${cent}.log" 2>&1
done
wait

for cent in 0 1 2 3; do
  ./run_everything_AA.sh 3 "$cent" \
    > "logs/draw_cent${cent}.log" 2>&1
done

for cent in 0 1 2 3; do
  ./run_closure_tests_AA.sh 3 "$cent" > "logs/closure_AA_cent${cent}.log" 2>&1 &
done

./run_dphi_plots_AA.sh 3 > "logs/dphi_plots_AA.log" 2>&1

wait
