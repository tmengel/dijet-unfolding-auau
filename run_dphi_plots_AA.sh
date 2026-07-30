#!/usr/bin/env bash
set -euo pipefail

source setup_env.sh

conesize="${1:-3}"


export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_jet.root"
export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_jet.root"

export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)

export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"


# Build the expensive simulation components only when their validated cache is
# absent.  Set FORCE_SIM_CACHE=1 after changing the simulation or reweights.
# force_sim_cache=true
# [[ "${FORCE_SIM_CACHE:-0}" == "1" ]] && force_sim_cache=true
# root -l -b -q  "makeCOMBSimulationCache_AA.C+(${conesize},-1,\"${AUAU_CONFIG}\",${force_sim_cache})" 

for centrality in 0 1 2 3; do
  root -l -b -q "drawCOMBModulation_AA_v2.C+(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
#   root -l -b -q "drawCOMBModulationSimSimple_AA.C+(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
#   root -l -b -q "drawCOMBModulationDataSimple_AA.C+(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
#   root -l -b -q "drawCOMBModulationSimSimple_AA.C+(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
#   root -l -b -q "drawCOMBModulation_AA.C+(${conesize}, ${centrality}, \"${AUAU_CONFIG}\")"
done


