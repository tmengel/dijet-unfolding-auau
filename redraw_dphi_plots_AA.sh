#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

cone_size="${1:-3}"
config="${2:-${AUAU_CONFIG}}"
plot_dir="$script_dir/dphi_plots"
log_dir="$DIJET_LOG_PATH"

export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_jet.root"
export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_jet.root"

export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)

# export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"

export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"

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
