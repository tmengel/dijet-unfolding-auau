#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

OPTION_A_SETTINGS="${OPTION_A_SETTINGS:-false}"
export OPTION_A_SETTINGS
echo "OPTION_A_SETTINGS=${OPTION_A_SETTINGS}"

conesize=$1
cent=$2
sysconfig=${3:-}
if [[ -z "$cent" || -z "$conesize" || -z "$sysconfig" ]]; then
    echo "Usage: $0 conesize centrality config"
    exit 1
fi

system="AA_cent_${cent}"
case "$(basename "$sysconfig")" in
    binning_COMBDown_AA.config) sys_name="COMBDown"; response_invariant=1 ;;
    binning_COMBUp_AA.config)   sys_name="COMBUp";   response_invariant=1 ;;
    binning_Inclusive_AA.config) sys_name="INCLUSIVE"; response_invariant=1 ;;
    *) sys_name=""; response_invariant=0 ;;
esac


root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 1)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"

root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${sysconfig}\")"

# root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 2 )"

# root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 2)"

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 0 )"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
