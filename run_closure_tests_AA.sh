#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

cone_size="${1:-3}"
centrality="${2:-0}"
config="${3:-$DIJET_CONFIG_PATH/binning_AA.config}"
export AUAU_CONFIG="$config"

export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_jet.root"
export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_jet.root"
export AUAU_DATA_NAME="$(basename "$AUAU_DATA_FILE" .root)"
export AUAU_SIM_NAME="$(basename "$AUAU_SIM_FILE" .root)"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${cone_size}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${cone_size}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${cone_size}_jet30_${AUAU_SIM_NAME}.root"

if (( centrality < 0 || centrality > 3 )); then
    echo "Centrality index must be 0, 1, 2, or 3" >&2
    exit 1
fi

unset DIJET_MAX_MATCHED_EVENTS

echo "Building full-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx+(\"${config}\",2,10,${cone_size},${centrality},0)"

echo "Building reproducible half-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx+(\"${config}\",1,10,${cone_size},${centrality},0)"

echo "Drawing full and half closure tests for centrality ${centrality}"
root -l -b -q "drawClosureTests_AA.C(${cone_size},${centrality},1,\"${config}\")"
