#!/bin/bash
set -euo pipefail

source "setup_env.sh"

conesize=3
cent=0


export AUAU_DATA_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root"
export AUAU_SIM_FILE="/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root"

export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)

export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"

export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"

# root -l -q -b "makeDataTreeAuAu.C(${conesize}, 1, \"${AUAU_DATA_FILE}\", \"${TNUPLE_DATA_FILE}\")"
# root -l -q -b "makeMatchedTreesInclusiveAuAu.C(${conesize}, \"${AUAU_SIM_FILE}\", \"${DIJET_TNTUPLE_PATH}\", \"${AUAU_CONFIG}\", \"${AUAU_SIM_NAME}\")" 

# root -l -q -b "makeUnfoldingHists.C(\"${AUAU_CONFIG}\", ${conesize}, ${cent})"
# root -l -q -b "getBackground.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", \"${AUAU_DATA_FILE}\")"

OPTION_A_SETTINGS=false
export OPTION_A_SETTINGS

# root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 1, ${OPTION_A_SETTINGS})"

# root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 1)"

# root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

# root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 2, ${OPTION_A_SETTINGS})"

# root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 2)"

# root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 0, ${OPTION_A_SETTINGS})"

# root -l -q -b "validateReweighting_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

# root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent})"

# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_negJES_AA.config
# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_posJES_AA.config
# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_negJER_AA.config
# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_posJER_AA.config

# root -l -b -q "makeFlowModulationPreload_AA.C(${conesize}, ${cent} \
#     , \"${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config\" \
#     , \"${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config\")"
# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config
# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config

bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_prior_AA.config

# bash run_all_sys_AA.sh ${conesize} ${cent} ${DIJET_CONFIG_PATH}/binning_Inclusive_AA.config

root -l -b -q "drawSys_AA.C(${conesize}, ${cent} )"

# root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",2,10,${conesize},${cent},0 , ${OPTION_A_SETTINGS})"
# root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",1,10,${conesize},${cent},0 , ${OPTION_A_SETTINGS})"

# root -l -b -q "drawClosureTests_AA.C(${conesize},${cent},1,\"${AUAU_CONFIG}\")"

# root -l -q -b "makeIterationPlot_AA.C(3,${cent},0)" 

root -l -q -b "drawFinalUnfold_AA.C(3,${cent}, \"${AUAU_CONFIG}\")"
root -l -b -q "drawFinalUnfold_AA_only.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"
root -l -q -b "drawFinalUnfold_AA_cent.C(3,${cent}, \"${AUAU_CONFIG}\")"

# root -l -q -b "drawPrior_AA.C(${conesize}, ${cent} )"