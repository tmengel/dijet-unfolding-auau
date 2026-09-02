#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

conesize="${CONESIZE:-3}"

export THIS_FILE="/home/tmengel/PPG14/rootfiles/data_v004_20260821_calibrated_merged.root"
export DATA_NAME=$(basename "$THIS_FILE" .root)
export CONFIG="/home/tmengel/PPG14/version1/newdata/configs/binning_AA.config"
export TNTUPLEFILE="/home/tmengel/PPG14/version1/newdata/rootfiles/TNTUPLE_DIJET_r0${conesize}_${DATA_NAME}.root"

log() { echo "[cone ${conesize}  $(date +%H:%M:%S)] $*"; }
root -l -q -b "makeDataTreeAuAu.C(${conesize}, 1, \"${THIS_FILE}\", \"${TNTUPLEFILE}\", \"${CONFIG}\")"
# root -l -q -b "getProbability.C(${conesize}, 0, \"${CONFIG}\")"
for cent in 0 1 2 3; do
    log "Processing centrality bin ${cent}"
    root -l -q -b "makeUnfoldingHists.C(${conesize}, ${cent}, \"${CONFIG}\",\"${THIS_FILE}\" )"
    root -l -q -b "getBackground.C(${conesize}, ${cent}, \"${CONFIG}\", \"${THIS_FILE}\")"
done

log "Preprocessing complete for cone size ${conesize}"