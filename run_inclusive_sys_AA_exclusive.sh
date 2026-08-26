#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

OPTION_A_SETTINGS="${OPTION_A_SETTINGS:-false}"
export OPTION_A_SETTINGS
echo "OPTION_A_SETTINGS=${OPTION_A_SETTINGS}"

# createResponse_exclusive_AA_inclusive.cxx variant of run_inclusive_sys_AA.sh:
# same three calls (primer1, primer2, nominal), pointed at the per-leg
# "leading paired with every other truth jet" TTree written by
# dijet_matching_inclusive.C instead of the AllPairs (C(n,2)) tn_match-style
# TNtuple written by makeMatchedTreesAllPairsTaggedAuAu.C that
# createResponse_noempty_AA_inclusive.cxx reads. Standalone, like the
# original -- not wired into run_full_cent_exclusive.sh, same as
# run_inclusive_sys_AA.sh isn't wired into run_full_cent.sh.
inclusive_dir="${INCLUSIVE_DIR:-/home/tmengel/PPG14/rootfiles/out/inclusive}"

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


root -l -q -b "createResponse_exclusive_AA_inclusive.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 1, \"${inclusive_dir}\")"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"
root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${sysconfig}\")"

root -l -q -b "createResponse_exclusive_AA_inclusive.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 2, \"${inclusive_dir}\")"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 2)"

root -l -q -b "createResponse_exclusive_AA_inclusive.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 0, \"${inclusive_dir}\")"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
