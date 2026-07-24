#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"


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

if (( response_invariant )); then
    # These variations alter the measured signal definition, not detector
    # smearing.  Reuse the current-run nominal primer response and event
    # weights, then rebuild the final response with the varied unfolded prior.
    cp "response_matrices/response_matrix_${system}_r0${conesize}_PRIMER1_nominal.root" \
       "response_matrices/response_matrix_${system}_r0${conesize}_PRIMER1_${sys_name}.root"
    root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"
    for directory in vertex centrality sumeT; do
        cp "${directory}/${directory}_reweight_${system}_r0${conesize}_nominal.root" \
           "${directory}/${directory}_reweight_${system}_r0${conesize}_${sys_name}.root"
    done
    root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent})"
    root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
    exit 0
fi

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 1)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"

root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${sysconfig}\")"

# root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 2)"

# root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 2)"

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent})"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
