#!/bin/bash
set -euo pipefail

source "setup_env.sh"


pids=()
for cent in 0 1 2 3; do
    bash redo_sys.sh "$cent" \
    > "logs/new_unfold_cent${cent}.log" 2>&1 &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "$pid"
done
