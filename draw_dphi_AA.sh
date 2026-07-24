#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"
cone_size="3"
config="${AUAU_CONFIG}"

opts=(
  # data_raw
  # data_background
  # data_nominal
  # data_down
  # data_up
  # sim_raw
  # sim_background
  # sim_nominal
  # sim_down
  # sim_up
  data_pairs
  # data_eta
  data_fit
  data_background
  data_background_COMBDown
  data_background_COMBUp
  # data_nominal
  # data_down
  # data_up
  data_signal_region
  data_normalization_region
  # sim_pairs
  sim_eta
  sim_fit
  # sim_background
  # sim_background_COMBDown
  # sim_background_COMBUp
  #  "data_pairs", "data_eta", "data_fit", "data_background", "data_background_COMBDown", "data_background_COMBUp",
  #   "data_nominal", "data_down", "data_up",
  #   "data_signal_region", "data_normalization_region",
  #   "sim_pairs", "sim_eta", "sim_fit", "sim_background", "sim_background_COMBDown", "sim_background_COMBUp",
  #   "sim_nominal", "sim_down", "sim_up",
  #   "sim_signal_region", "sim_normalization_region",
  #   "sim_truth_pairs", "sim_truth_eta"
)
funcs=(
  data_flow_fit
  # data_flow_fit_nominal
  # data_flow_fit_COMBDown
  # data_flow_fit_COMBUp
  sim_flow_fit
  # sim_flow_fit_nominal
  # sim_flow_fit_COMBDown
  # sim_flow_fit_COMBUp
)
opts_str=$(IFS=,; echo "${opts[*]}")
funcs_str=$(IFS=,; echo "${funcs[*]}")



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


for centrality in 0 1 2 3; do

  radius="r$(printf '%02d' "$cone_size")"
  data_cache="$plot_dir/dphi_COMB_modulation_AA_cent_${centrality}_${radius}.root"
  sim_cache="$plot_dir/dphi_COMB_modulation_sim_combined_AA_cent_${centrality}_${radius}.root"

  for cache in "$data_cache" "$sim_cache"; do
    if [[ ! -f "$cache" ]]; then
      # echo "Missing component cache: $cache" >&2
      echo "Missing component cache: $cache"
      # echo "Run ./run_dphi_plots_AA.sh first." >&2
      exit 1
    fi
  done

  data_log="$log_dir/redrawCOMBModulationDataSim_AA_cent${centrality}.log"
  
  echo "Redrawing data centrality $centrality (log: $data_log)"
  root -l -b -q "drawCOMBDataSimSimple_AA.C+($cone_size,$centrality,\"$opts_str\",\"$funcs_str\", true,\"$config\")" 
  
done

echo "Redrew all data and combined-simulation COMB PDFs from cached components."
