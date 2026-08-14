#!/bin/bash
set -euo pipefail

source "setup_env.sh"

# conesize="${1:-3}"

# export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root"
# export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root"

# export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
# export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)

# export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"

# export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
# export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
# export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
# export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"

# dounfold=${DO_UNFOLD:-1}
# docahche=${DO_CACHE:-1}
# doplots=${DO_PLOTS:-1}

# pids=()
# for cent in 0 1; do
#   RUN_UNCERTAINTIES=1 \
#   RUN_SYSTEMATICS=1 \
#   RUN_FLOW_SYS=1 \
#     ./run_all_unfold_AA.sh 3 "$cent" \
#     > "logs/unfold_cent${cent}.log" 2>&1 &
#   pids+=("$!")
# done
pids=()
for cent in 0; do
  bash redo_sys.sh "$cent" \
    > "logs/unfold_cent${cent}.log" 2>&1 &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "$pid"
done

# force_sim_cache=false
# dphilogfile="logs/draw_dphi.log"
# [[ "${FORCE_SIM_CACHE:-0}" == "1" ]] && force_sim_cache=true
# root -l -b -q "makeCOMBSimulationCache_AA.C+(${conesize},-1,\"${AUAU_CONFIG}\",${force_sim_cache})" > "$dphilogfile" 2>&1
# wait

# pidsd=()
# for cent in 0; do
#   ./run_everything_AA.sh 3 "$cent" \
#     > "logs/draw_cent${cent}.log" 2>&1 &
#   pidsd+=("$!")
# done
# for pid in "${pidsd[@]}"; do
#   wait "$pid"
# done
# for cent in 0 1 2 3; do
  # root -l -q -b "drawFinalUnfold_AA.C(3,${cent}, \"${AUAU_CONFIG}\")"
  # root -l -q -b "drawFinalUnfold_AA_cent.C+(3,${cent}, \"${AUAU_CONFIG}\")"
  # root -l -q -b "drawFinalUnfold_AA_only.C+(3,${cent}, \"${AUAU_CONFIG}\", 2)"
# done 
