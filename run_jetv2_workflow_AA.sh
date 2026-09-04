#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

# Full jet-v2 cross-check for one centrality bin: builds the v2-reweighted
# response/unfold via run_jetv2_sys_AA_exclusive.sh, then draws the unfolded
# x_J overlay against nominal (compareJetV2Xj_AA.C, whose `jetv2s` array
# lists exactly the sys_names built below -- add a call below and an entry
# in that array to compare another v2 value) and the applied per-pair weight
# (drawJetV2Weight_AA.C, reading h_jetv2_weight and, when present,
# h_jetv2_weight_vs_dpsi2_leg{1,2}/jetv2_amp straight out of the response
# file -- the weight-vs-dpsi2 profile against the analytic curve is the
# actual proof the v2 injection is doing what it claims, not just that it
# averages to 1). Mirrors run_flavor_workflow_AA.sh for the jet-v2
# "what if" cross-check.
# Not a systematic uncertainty -- see createResponse_exclusive_v2_AA.cxx's
# JETV2_SCALE handling and drawSys_AA.C, which deliberately does not include it.
#
# Same env vars run_full_cent_exclusive.sh and run_jetv2_sys_AA_exclusive.sh
# use; set here too since compareJetV2Xj_AA.C is called directly below, not
# just through the wrapper script.
conesize="${CONESIZE:-3}"
cent="${1:?Usage: $0 <centrality_bin 0-3> [v2]}"
jet_v2="${2:-0.03}"
niter="${NITER:-1}"

export AUAU_DATA_FILE="${AUAU_DATA_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
export AUAU_SIM_FILE="${AUAU_SIM_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root}"
export AUAU_DATA_NAME="$(basename "$AUAU_DATA_FILE" .root)"
export AUAU_SIM_NAME="$(basename "$AUAU_SIM_FILE" .root)"
export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"
export TNUPLE_DATA_FILE="${TNUPLE_DATA_FILE:-${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root}"

log() { echo "[jetv2 cent ${cent} $(date +%H:%M:%S)] $*"; }
log "=== run_jetv2_workflow_AA.sh: conesize=${conesize}, cent=${cent}, v2=${jet_v2} ==="

log "=== jet-v2: primer1/primer2/nominal response + unfold ==="
bash run_jetv2_sys_AA_exclusive.sh "${conesize}" "${cent}" "${jet_v2}"

# Matches run_jetv2_sys_AA_exclusive.sh's sys_name derivation so
# drawJetV2Weight_AA.C reads the response file this run just wrote.
pct="$(awk -v v="$jet_v2" 'BEGIN { printf "%d", (v * 100) + 0.5 }')"
sys_name="JETV2${pct}"

log "=== comparison plot: jet-v2 vs nominal unfolded x_J ==="
root -l -q -b "compareJetV2Xj_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", ${niter})"

log "=== weight QA plot: applied jet-v2 reweight distribution ==="
root -l -q -b "drawJetV2Weight_AA.C(${conesize}, ${cent}, \"${sys_name}\", \"${AUAU_CONFIG}\")"

log "=== done ==="
