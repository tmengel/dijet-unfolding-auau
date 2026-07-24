#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"


conesize=$1
cent=$2
if [[ -z "$cent" || -z "$conesize" ]]; then
    echo "Usage: $0 conesize centrality"
    exit 1
fi

run_uncertainties=${RUN_UNCERTAINTIES:-1}
run_systematics=${RUN_SYSTEMATICS:-1}
run_flow_sys=${RUN_FLOW_SYS:-1}

if [[ -n "${DIJET_PROOF_EVENTS:-}" ]]; then
    export DIJET_MAX_DATA_EVENTS="${DIJET_MAX_DATA_EVENTS:-$DIJET_PROOF_EVENTS}"
    export DIJET_MAX_MATCHED_EVENTS="${DIJET_MAX_MATCHED_EVENTS:-100000}"
    run_uncertainties=${RUN_UNCERTAINTIES:-0}
    run_systematics=${RUN_SYSTEMATICS:-0}
    echo "Proof/sample mode: DIJET_MAX_DATA_EVENTS=${DIJET_MAX_DATA_EVENTS}, DIJET_MAX_MATCHED_EVENTS=${DIJET_MAX_MATCHED_EVENTS}"
fi

probability_file="$DIJET_UNFOLDING_PATH/unfolding_hists/probability_hists_AA_r0${conesize}.root"
if (( cent >= 0 )) && [[ ! -f "$probability_file" ]]; then
    echo "Missing required probability-correction file: $probability_file" >&2
    echo "Provide this analysis prerequisite before running the unfolding workflow." >&2
    exit 1
fi

#root -l -q -b "unfoldDataUncertainties_noempty.cxx(\"binning.config\"${conesize})"
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

root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent}, 1)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent}, 1)"

root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"binning_AA.config\")"

root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent}, 2)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent}, 2)"

root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent})"

root -l -q -b "validateReweighting_AA.C(${conesize}, ${cent}, \"binning_AA.config\")"

root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent})"

ntuple_file="$DIJET_UNFOLDING_PATH/rootfiles/unfolded_ntuple_AA_cent_${cent}_r0${conesize}.root"
# root -l -q -b "makeDataTreeAuAu.C(${conesize}, ${cent})"

if (( run_uncertainties )); then
    root -l -q -b "unfoldDataUncertainties_noempty_AA.cxx(10, ${conesize}, ${cent})"
else
    echo "Skipping unfolding uncertainties (RUN_UNCERTAINTIES=0)"
fi

if (( run_systematics )); then
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_negJES_AA.config
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_posJES_AA.config
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_negJER_AA.config
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_posJER_AA.config

    if (( run_flow_sys )); then
        root -l -b -q "makeFlowModulationPreload_AA.C(${conesize}, ${cent})"
        bash run_all_sys_AA.sh ${conesize} ${cent} binning_COMBDown_AA.config
        bash run_all_sys_AA.sh ${conesize} ${cent} binning_COMBUp_AA.config
    else
        echo "Skipping flow-modulation systematic variation (RUN_FLOW_SYS=0)"
    fi
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_prior_AA.config
    bash run_all_sys_AA.sh ${conesize} ${cent} binning_Inclusive_AA.config
else
    echo "Skipping systematic variations (RUN_SYSTEMATICS=0)"
fi
