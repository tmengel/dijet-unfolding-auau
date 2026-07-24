#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"


conesize=3
ROOTFILEBASE="/home/tmengel/PPG14/rootfiles/v001_20260720"
mkdir -p "$DIJET_UNFOLDING_PATH/rootfiles"

DATASET="run2auau_rho_jet"
SIMSET="hijing_rho_jet"

DATAFILE="${ROOTFILEBASE}/${DATASET}.root"
TNUPLE_DATAFILE="${DIJET_UNFOLDING_PATH}/rootfiles/TNTUPLE_DIJET_r0${conesize}_${DATASET}.root"
if [[ ! -f "$TNUPLE_DATAFILE" ]]; then
    echo "Missing TNTUPLE data file: $TNUPLE_DATAFILE"
    root -l -q -b "makeDataTreeAuAu.C(${conesize},1, \"${DATAFILE}\", \"${TNUPLE_DATAFILE}\")"
fi

SIMFILE="${ROOTFILEBASE}/${SIMSET}.root"
SIMFILE_TNTUPLE_DIR="${DIJET_UNFOLDING_PATH}/rootfiles/${SIMSET}"
mkdir -p "$SIMFILE_TNTUPLE_DIR"
SIMFILE_TNTUPLE_10="${SIMFILE_TNTUPLE_DIR}/TREE_MATCH_r0${conesize}_v15_10_new_ProdA_2024-00000030_sumeT.root"
SIMFILE_TNTUPLE_20="${SIMFILE_TNTUPLE_DIR}/TREE_MATCH_r0${conesize}_v15_20_new_ProdA_2024-00000030_sumeT.root"
SIMFILE_TNTUPLE_30="${SIMFILE_TNTUPLE_DIR}/TREE_MATCH_r0${conesize}_v15_30_new_ProdA_2024-00000030_sumeT.root"
if [[ ! -f "$SIMFILE_TNTUPLE_10" || ! -f "$SIMFILE_TNTUPLE_20" || ! -f "$SIMFILE_TNTUPLE_30" ]]; then
    echo "Missing TNTUPLE simulation files: $SIMFILE_TNTUPLE_10, $SIMFILE_TNTUPLE_20, $SIMFILE_TNTUPLE_30"
    root -l -q -b "makeMatchedTreesInclusiveAuAu.C(${conesize}, \"${SIMFILE}\", \"${SIMFILE_TNTUPLE_DIR}\")"
fi

for cent in 0 1 2 3; do
    unfoldinghists_file="$DIJET_UNFOLDING_PATH/unfolding_hists/unfolding_hists_AA_cent_${cent}_r0${conesize}.root"
    if [[ -s "$unfoldinghists_file" && "${FORCE_UNFOLDING_HISTS:-0}" != "1" ]]; then
        echo "Using locally produced unfolding histograms $unfoldinghists_file"
    else
        root -l -q -b "makeUnfoldingHists.C(\"binning_AA.config\", ${conesize}, ${cent})"
    fi
    background_file="$DIJET_UNFOLDING_PATH/unfolding_hists/background_hists_AA_cent_${cent}_r0${conesize}.root"
    if (( cent < 0 )); then
        echo "pp - skipping background"
    elif [[ -s "$background_file" && "${FORCE_BACKGROUND:-0}" != "1" ]]; then
        echo "Using locally produced background $background_file"
    else
        root -l -q -b "getBackground.C(${conesize}, ${cent}, \"binning_AA.config\", \"${AUAU_DATA_PATH}/run2auau_rho_jet.root\")"
    fi
done