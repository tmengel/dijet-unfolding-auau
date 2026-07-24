#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

cone_size="${1:-3}"
config="${2:-$script_dir/configs/binning_AA.config}"
plot_dir="$script_dir/dphi_plots"
log_dir="$script_dir/logs"

if [[ ! "$cone_size" =~ ^[0-9]+$ ]]; then
  echo "Cone size must be an integer: $cone_size" >&2
  exit 2
fi

if [[ ! -f "$config" ]]; then
  echo "Configuration file does not exist: $config" >&2
  exit 2
fi

mkdir -p "$plot_dir" "$log_dir"

# setup_env.sh defines ROOT and the paths consumed by read_binning.
source "$script_dir/setup_env.sh"
export AUAU_CONFIG="$config"

# Build the expensive simulation components only when their validated cache is
# absent.  Set FORCE_SIM_CACHE=1 after changing the simulation or reweights.
force_sim_cache=false
[[ "${FORCE_SIM_CACHE:-0}" == "1" ]] && force_sim_cache=true
sim_log="$log_dir/makeCOMBSimulationCache_AA.log"
echo "Preparing reusable simulation cache (log: $sim_log)"
root -l -b -q \
  "makeCOMBSimulationCache_AA.C+(${cone_size},-1,\"${config}\",${force_sim_cache})" \
  > "$sim_log" 2>&1
tail -n 4 "$sim_log"

for centrality in 0 1 2 3; do
  log_file="$log_dir/drawCOMBModulation_AA_v2_cent${centrality}.log"
  echo "Running centrality $centrality (log: $log_file)"
  root -l -b -q \
    "drawCOMBModulation_AA_v2.C+(${cone_size},${centrality},\"${config}\")" \
    > "$log_file" 2>&1
  tail -n 1 "$log_file"
  root -l -b -q \
    "drawCOMBModulationSimSimple_AA.C+(${cone_size},${centrality},\"${config}\")" \
    > "$log_dir/drawCOMBModulationSimSimple_AA_cent${centrality}.log" 2>&1
done

pdf_count="$(find "$plot_dir" -maxdepth 1 -type f -name '*.pdf' | wc -l)"

echo "Completed clean Au+Au dphi plotting run:"
echo "  PDF: $pdf_count"
