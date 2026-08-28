#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

# Full qq / qg+gg flavor-composition cross-check for one centrality bin:
# builds both flavor-tagged responses (each its own primer1/primer2/nominal
# pass via run_flavor_sys_AA_exclusive.sh) and draws the unfolded x_J
# overlay against nominal (compareFlavorXj_AA.C). Not a systematic
# uncertainty -- see createResponse_exclusive_AA.cxx's FLAVOR handling and
# drawSys_AA.C, which deliberately does not include it.
#
# Same env vars run_full_cent_exclusive.sh and run_flavor_sys_AA_exclusive.sh
# use; set here too since compareFlavorXj_AA.C is called directly below, not
# just through the wrapper script.
conesize="${CONESIZE:-3}"
cent="${1:?Usage: $0 <centrality_bin 0-3>}"
niter="${NITER:-1}"

export AUAU_DATA_FILE="${AUAU_DATA_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
export AUAU_SIM_FILE="${AUAU_SIM_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root}"
export AUAU_DATA_NAME="$(basename "$AUAU_DATA_FILE" .root)"
export AUAU_SIM_NAME="$(basename "$AUAU_SIM_FILE" .root)"
export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"
export TNUPLE_DATA_FILE="${TNUPLE_DATA_FILE:-${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root}"

log() { echo "[flavor cent ${cent} $(date +%H:%M:%S)] $*"; }
log "=== run_flavor_workflow_AA.sh: conesize=${conesize}, cent=${cent} ==="

log "=== qq: primer1/primer2/nominal response + unfold ==="
bash run_flavor_sys_AA_exclusive.sh "${conesize}" "${cent}" qq

log "=== qg_gg: primer1/primer2/nominal response + unfold ==="
bash run_flavor_sys_AA_exclusive.sh "${conesize}" "${cent}" qg_gg

log "=== comparison plot: qq vs qg_gg vs nominal unfolded x_J ==="
root -l -q -b "compareFlavorXj_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", ${niter})"
root -l -q -b "drawPriorQA_xj_flavorCompare_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

log "=== done ==="
