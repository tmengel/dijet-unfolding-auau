#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

# Thin wrapper around run_all_sys_AA_exclusive_v2.sh for the jet-v2
# cross-check (see createResponse_exclusive_v2_AA.cxx's JETV2_SCALE block).
# HIJING+Pythia embeds the signal jets flat in azimuth relative to the event
# plane, so this asks "what if jets had a v2?" by reweighting each MC pair by
# 1 + 2*v2*cos(2*dpsi2) on BOTH truth legs. Because that shifts which part of
# the modulated UE each jet sits on, it moves the reco legs' background
# subtraction and so rebuilds the response -- it is not a flat scaling of the
# final spectrum, which is why the whole primer1/primer2/nominal chain is rerun
# here rather than just the last step.
#
# Like the flavor cross-check this is a "what if", NOT a measured up/down
# uncertainty, so it is deliberately not wired into drawSys_AA.C's total band.
# No new logic here: this just generates a config with the requested v2 and
# hands it to the normal exclusive-v2 chain.
#
# Usage:
#   run_jetv2_sys_AA_exclusive.sh conesize centrality [v2]
# v2 defaults to 0.03 (the 3% cross-check). sys_name is JETV2<percent>
# (JETV23 for 3%, JETV25 for 5%, ...), matching how the config's
# JETV2_SYSTEMATIC_NAME names the response/unfold/QA outputs, so several v2
# values can be run without colliding.
#
# Meant to be run standalone (like run_flavor_sys_AA_exclusive.sh, unlike
# run_all_sys_AA_exclusive_v2.sh which normally runs as a sub-step of
# run_full_cent_exclusive_v2.sh after that script has exported these), so it
# sets up the same env vars itself -- unfoldData_noempty_AA.cxx reads
# TNUPLE_DATA_FILE via an unguarded std::getenv() and throws std::logic_error
# if it is unset.
export AUAU_DATA_FILE="${AUAU_DATA_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
export AUAU_SIM_FILE="${AUAU_SIM_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root}"
export AUAU_DATA_NAME="$(basename "$AUAU_DATA_FILE" .root)"
export AUAU_SIM_NAME="$(basename "$AUAU_SIM_FILE" .root)"

conesize="${1:?Usage: $0 conesize centrality [v2]}"
cent="${2:?Usage: $0 conesize centrality [v2]}"
jet_v2="${3:-0.03}"

export TNUPLE_DATA_FILE="${TNUPLE_DATA_FILE:-${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root}"

# createResponse_exclusive_v2_AA.cxx rejects anything outside (0, 0.5): beyond
# that 1 + 2*v2*cos(2*dpsi2) goes negative. Fail here too so a typo does not
# cost a full response build first.
if ! awk -v v="$jet_v2" 'BEGIN { exit !(v > 0 && v < 0.5) }'; then
    echo "v2 must be a number in (0, 0.5), got: $jet_v2" >&2
    exit 1
fi

pct="$(awk -v v="$jet_v2" 'BEGIN { printf "%d", (v * 100) + 0.5 }')"
sys_name="JETV2${pct}"

# Not a preset per v2 value: generate a temp config off the committed
# binning_JETV2_AA.config template with JETV2_SCALE / JETV2_SYSTEMATIC_NAME
# substituted, the same way run_flavor_sys_AA_exclusive.sh handles MIX.
sysconfig="configs/.tmp_binning_${sys_name}_AA.config"
sed -e "s/^JETV2_SCALE:.*/JETV2_SCALE: ${jet_v2}/" \
    -e "s/^JETV2_SYSTEMATIC_NAME:.*/JETV2_SYSTEMATIC_NAME: ${sys_name}/" \
    configs/binning_JETV2_AA.config > "$sysconfig"

echo "Jet-v2 cross-check: v2 = ${jet_v2} -> sys_name = ${sys_name}, config = ${sysconfig}"

# Same exclusive input files as the nominal response -- the v2 is injected by
# reweighting, not by a different sample, so there is no JETV2 input directory.
export EXCLUSIVE_DIR="${EXCLUSIVE_DIR:-/home/tmengel/PPG14/rootfiles/dijet_match_08_31_2026/exclusive}"

bash run_all_sys_AA_exclusive_v2.sh "${conesize}" "${cent}" "${sysconfig}"

# QA plots for the v2-reweighted response/unfold -- same macros as every other
# systematic, pointed at this sys_name so they read
# response_matrix_..._JETV2<pct>.root instead of ..._nominal.root. No
# FULL_/HALF_ closure passes are run above, so drawResponse_AA.C skips those
# panels on its own.
export AUAU_CONFIG="${sysconfig}"
root -l -q -b "drawPrior_AA.C(${conesize}, ${cent}, \"${sys_name}\")"
root -l -q -b "drawResponse_AA.C(${conesize}, ${cent}, \"${sys_name}\")"
