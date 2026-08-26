#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

# Thin wrapper around run_all_sys_AA_exclusive.sh for the qq / qg+gg
# flavor-composition cross-check (see dijet_matching_flavor.C and
# createResponse_exclusive_AA.cxx's FLAVOR handling). Not a systematic
# uncertainty -- deliberately not wired into drawSys_AA.C's total-uncertainty
# band. No new logic here: this just points EXCLUSIVE_DIR at the
# flavor-tagged files and picks the matching binning_QQ_AA.config /
# binning_QGGG_AA.config.
#
# This is meant to be run standalone (unlike run_all_sys_AA_exclusive.sh,
# which normally runs as a sub-step of run_full_cent_exclusive.sh after that
# script has already exported these), so it sets up the same
# AUAU_DATA_FILE/AUAU_SIM_FILE/TNUPLE_DATA_FILE env vars itself --
# unfoldData_noempty_AA.cxx reads TNUPLE_DATA_FILE via an unguarded
# std::getenv() and crashes with std::logic_error if it's unset.
export AUAU_DATA_FILE="${AUAU_DATA_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
export AUAU_SIM_FILE="${AUAU_SIM_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root}"
export AUAU_DATA_NAME="$(basename "$AUAU_DATA_FILE" .root)"
export AUAU_SIM_NAME="$(basename "$AUAU_SIM_FILE" .root)"
conesize_for_tnuple="${1:-3}"
export TNUPLE_DATA_FILE="${TNUPLE_DATA_FILE:-${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize_for_tnuple}_${AUAU_DATA_NAME}.root}"

conesize="${1:?Usage: $0 conesize centrality qq|qg_gg}"
cent="${2:?Usage: $0 conesize centrality qq|qg_gg}"
flavor="${3:?Usage: $0 conesize centrality qq|qg_gg}"

case "$flavor" in
    qq)    sysconfig="configs/binning_QQ_AA.config"; sys_name="QQ" ;;
    qg_gg) sysconfig="configs/binning_QGGG_AA.config"; sys_name="QGGG" ;;
    *)     echo "Usage: $0 conesize centrality qq|qg_gg" >&2; exit 1 ;;
esac

export EXCLUSIVE_DIR="${EXCLUSIVE_DIR:-/home/tmengel/PPG14/rootfiles/out/flavor_v2}"

bash run_all_sys_AA_exclusive_v2.sh "${conesize}" "${cent}" "${sysconfig}"

# QA plots for the flavor-tagged response/unfold -- same macros used for
# every other systematic, pointed at this flavor's sys_name so they read
# response_matrix_..._QQ.root / ..._QGGG.root instead of ..._nominal.root.
# No FULL_/HALF_ closure passes exist for a flavor (only primer1/primer2/
# nominal are run above), so drawResponse_AA.C skips those panels on its own.
export AUAU_CONFIG="${sysconfig}"
root -l -q -b "drawPrior_AA.C(${conesize}, ${cent}, \"${sys_name}\")"
root -l -q -b "drawResponse_AA.C(${conesize}, ${cent}, \"${sys_name}\")"
