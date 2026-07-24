#!/usr/bin/env bash
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo"
source "$repo/setup_env.sh"

conesize=${1:-3}

TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_DATA_FILE

if [[ ! -f "$TNUPLE_DATA_FILE" ]]; then
    echo "Missing TNTUPLE data file: $TNUPLE_DATA_FILE"
    root -l -q -b "makeDataTreeAuAu.C(${conesize},1, \"${AUAU_DATA_FILE}\", \"${TNUPLE_DATA_FILE}\")" > "$DIJET_LOG_PATH/makeDataTreeAuAu_r0${conesize}.log" 2>&1
fi

DIJET_TNTUPLE_SIM_PATH="${DIJET_TNTUPLE_PATH}/${AUAU_SIM_NAME}"
mkdir -p "$DIJET_TNTUPLE_SIM_PATH"

TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_SIM_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_SIM_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_SIM_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET10
export TNUPLE_SIM_FILE_JET20
export TNUPLE_SIM_FILE_JET30

if [[ ! -f "$TNUPLE_SIM_FILE_JET10" || ! -f "$TNUPLE_SIM_FILE_JET20" || ! -f "$TNUPLE_SIM_FILE_JET30" ]]; then
    echo "Missing TNTUPLE simulation files: $TNUPLE_SIM_FILE_JET10, $TNUPLE_SIM_FILE_JET20, $TNUPLE_SIM_FILE_JET30"
    root -l -q -b "makeMatchedTreesInclusiveAuAu.C(${conesize}, \"${AUAU_SIM_FILE}\", \"${DIJET_TNTUPLE_SIM_PATH}\", \"${AUAU_SIM_NAME}\")" > "$DIJET_LOG_PATH/makeMatchedTreesInclusiveAuAu_r0${conesize}.log" 2>&1
fi

for cent in 0 1 2 3; do
    unfoldinghists_file="$DIJET_UNFOLDING_PATH/unfolding_hists/unfolding_hists_AA_cent_${cent}_r0${conesize}.root"
    if [[ -s "$unfoldinghists_file" && "${FORCE_UNFOLDING_HISTS:-0}" != "1" ]]; then
        echo "Using locally produced unfolding histograms $unfoldinghists_file"
    else
        root -l -q -b "makeUnfoldingHists.C(\"${AUAU_CONFIG}\", ${conesize}, ${cent})"
    fi
    background_file="$DIJET_UNFOLDING_PATH/unfolding_hists/background_hists_AA_cent_${cent}_r0${conesize}.root"
    if (( cent < 0 )); then
        echo "pp - skipping background"
    elif [[ -s "$background_file" && "${FORCE_BACKGROUND:-0}" != "1" ]]; then
        echo "Using locally produced background $background_file"
    else
        root -l -q -b "getBackground.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", \"${AUAU_DATA_FILE}\")"
    fi
done