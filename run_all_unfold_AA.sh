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
run_jes_sys=${RUN_JES_SYS:-1}
run_jer_sys=${RUN_JER_SYS:-1}
run_prior_sys=${RUN_PRIOR_SYS:-1}
run_inclusive_sys=${RUN_INCLUSIVE_SYS:-1}

probability_file="$DIJET_UNFOLDING_PATH/unfolding_hists/probability_hists_AA_r0${conesize}.root"
if (( cent >= 0 )) && [[ ! -f "$probability_file" ]]; then
    echo "Missing required probability-correction file: $probability_file" >&2
    echo "Provide this analysis prerequisite before running the unfolding workflow." >&2
    exit 1
fi

export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_jet.root"
export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_jet.root"

export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)

export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"

export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"

force_reprocess_tntuple=${FORCE_REPROCESS_TNTUPLE:-0}
echo "Running unfolding workflow for conesize=$conesize, centrality=$cent"
if [[ ! -f "$TNUPLE_DATA_FILE" || "$force_reprocess_tntuple" == "1" ]]; then
    echo "Missing TNTUPLE data file: $TNUPLE_DATA_FILE"
    root -l -q -b "makeDataTreeAuAu.C(${conesize}, 1, \"${AUAU_DATA_FILE}\", \"${TNUPLE_DATA_FILE}\")"
fi
if [[ ! -f "$TNUPLE_SIM_FILE_JET10" || ! -f "$TNUPLE_SIM_FILE_JET20" || ! -f "$TNUPLE_SIM_FILE_JET30" || "$force_reprocess_tntuple" == "1" ]]; then
    echo "Missing TNTUPLE simulation files: $TNUPLE_SIM_FILE_JET10, $TNUPLE_SIM_FILE_JET20, $TNUPLE_SIM_FILE_JET30"
    root -l -q -b "makeMatchedTreesInclusiveAuAu.C(${conesize}, \"${AUAU_SIM_FILE}\", \"${DIJET_TNTUPLE_PATH}\", \"${AUAU_CONFIG}\", \"${AUAU_SIM_NAME}\")" 
fi

unfoldinghists_file="$DIJET_UNFOLDING_PATH/unfolding_hists/unfolding_hists_AA_cent_${cent}_r0${conesize}.root"
if [[ -s "$unfoldinghists_file" && "${FORCE_UNFOLDING_HISTS:-0}" != "1" ]]; then
    echo "Using locally produced unfolding histograms $unfoldinghists_file"
else
    echo "Missing unfolding histograms file: $unfoldinghists_file"
    root -l -q -b "makeUnfoldingHists.C(\"${AUAU_CONFIG}\", ${conesize}, ${cent})"
fi

background_file="$DIJET_UNFOLDING_PATH/unfolding_hists/background_hists_AA_cent_${cent}_r0${conesize}.root"
if (( cent < 0 )); then
    echo "pp - skipping background"
elif [[ -s "$background_file" && "${FORCE_BACKGROUND:-0}" != "1" ]]; then
    echo "Using locally produced background $background_file"
else
    echo "This will be ignored"
    # root -l -q -b "getBackground.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", \"${AUAU_DATA_FILE}\")"
fi

root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 1)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 1)"

root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 2)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 2)"

root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent})"

root -l -q -b "validateReweighting_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent})"

if (( run_uncertainties )); then
    root -l -q -b "unfoldDataUncertainties_noempty_AA.cxx(10, ${conesize}, ${cent})"
else
    echo "Skipping unfolding uncertainties (RUN_UNCERTAINTIES=0)"
fi

#function for running systematic variations, we will call run_all_sys_AA.sh with the appropriate config file
# invariant_sys(){
#     local sysconfig=$1
#     local sys_name=$2
#     cp "${DIJET_UNFOLDING_PATH}/response_matrices/response_matrix_AA_cent_${cent}_r0${conesize}_PRIMER1_nominal.root" \
#        "${DIJET_UNFOLDING_PATH}/response_matrices/response_matrix_AA_cent_${cent}_r0${conesize}_PRIMER1_${sys_name}.root"
#     root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"
#     for directory in vertex centrality sumeT; do
#         cp "${DIJET_UNFOLDING_PATH}/${directory}/${directory}_reweight_AA_cent_${cent}_r0${conesize}_nominal.root" \
#            "${DIJET_UNFOLDING_PATH}/${directory}/${directory}_reweight_AA_cent_${cent}_r0${conesize}_${sys_name}.root"
#     done
#     root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent})"
#     root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
# }
# regular_sys(){
#     local sysconfig=$1
#     local sys_name=$2
#     echo "Running systematic variation: $sys_name"
#     root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent},1)"
#     root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent},1)"
#     root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${sysconfig}\")"
#     root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent})"
#     root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"
# }

if (( run_systematics )); then

    if (( run_jes_sys )); then
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_negJES_AA.config
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_posJES_AA.config
        # regular_sys "${DIJET_CONFIG_PATH}/binning_negJES_AA.config" "negJES"
        # regular_sys "${DIJET_CONFIG_PATH}/binning_posJES_AA.config" "posJES"
    else
        echo "Skipping JES systematic variations (RUN_JES_SYS=0)"
    fi
    if (( run_jer_sys )); then
        # regular_sys "${DIJET_CONFIG_PATH}/binning_negJER_AA.config" "negJER"
        # regular_sys "${DIJET_CONFIG_PATH}/binning_posJER_AA.config" "posJER"
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_negJER_AA.config
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_posJER_AA.config
    else
        echo "Skipping JER systematic variations (RUN_JER_SYS=0)"
    fi

    if (( run_flow_sys )); then

        root -l -b -q "makeFlowModulationPreload_AA.C(${conesize}, ${cent} \
            , \"${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config\" \
            , \"${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config\")"
        # invariant_sys "${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config" "COMBDown"
        # invariant_sys "${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config" "COMBUp"
        # invariant_sys "${DIJET_CONFIG_PATH}/binning_Inclusive_AA.config" "INCLUSIVE"
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config
    else
        echo "Skipping flow-modulation systematic variation (RUN_FLOW_SYS=0)"
    fi

    if (( run_prior_sys )); then
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_prior_AA.config
        # regular_sys "${DIJET_CONFIG_PATH}/binning_prior_AA.config" "PRIOR"
    else
        echo "Skipping prior systematic variation (RUN_PRIOR_SYS=0)"
    fi

    if (( run_inclusive_sys )); then
        # invariant_sys "${DIJET_CONFIG_PATH}/binning_Inclusive_AA.config" "INCLUSIVE"
        bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_Inclusive_AA.config
    else
        echo "Skipping inclusive systematic variation (RUN_INCLUSIVE_SYS=0)"
    fi
else
    echo "Skipping systematic variations (RUN_SYSTEMATICS=0)"
fi



# do_closure_tests=${DO_CLOSURE_TESTS:-1}
# if (( do_closure_tests )); then
#     bash run_closure_tests_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_AA.config
# else
#     echo "Skipping closure tests (DO_CLOSURE_TESTS=0)"
# fi

