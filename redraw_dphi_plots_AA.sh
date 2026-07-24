#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

cone_size="${1:-3}"
config="${2:-binning_AA.config}"
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

mkdir -p "$log_dir"
source "$script_dir/setup_env.sh"

for centrality in 0 1 2 3; do
  radius="r$(printf '%02d' "$cone_size")"
  data_cache="$plot_dir/dphi_COMB_modulation_AA_cent_${centrality}_${radius}.root"
  sim_cache="$plot_dir/dphi_COMB_modulation_sim_combined_AA_cent_${centrality}_${radius}.root"
  for cache in "$data_cache" "$sim_cache"; do
    if [[ ! -f "$cache" ]]; then
      echo "Missing component cache: $cache" >&2
      echo "Run ./run_dphi_plots_AA.sh first." >&2
      exit 1
    fi
  done

  data_log="$log_dir/redrawCOMBModulationDataSimple_AA_cent${centrality}.log"
  echo "Redrawing data centrality $centrality (log: $data_log)"
  root -l -b -q \
    "drawCOMBModulationDataSimple_AA.C+(${cone_size},${centrality},\"${config}\")" \
    > "$data_log" 2>&1
  tail -n 1 "$data_log"

  sim_log="$log_dir/redrawCOMBModulationSimSimple_AA_cent${centrality}.log"
  echo "Redrawing simulation centrality $centrality (log: $sim_log)"
  root -l -b -q \
    "drawCOMBModulationSimSimple_AA.C+(${cone_size},${centrality},\"${config}\")" \
    > "$sim_log" 2>&1
  tail -n 1 "$sim_log"
done

echo "Redrew all data and combined-simulation COMB PDFs from cached components."
