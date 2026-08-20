#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

conesize="${CONESIZE:-3}"
cent="${1:?Usage: $0 <centrality_bin 0-3>}"

# Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2, matching the
# prior_iteration constant in createResponse_noempty_AA.cxx. The iteration scan
# gives the same optimum for every centrality bin.
nominal_iter=1

export AUAU_DATA_FILE="${AUAU_DATA_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
export AUAU_SIM_FILE="${AUAU_SIM_FILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_alljet_rho_jet.root}"
export AUAU_DATA_NAME=$(basename "$AUAU_DATA_FILE" .root)
export AUAU_SIM_NAME=$(basename "$AUAU_SIM_FILE" .root)
export AUAU_CONFIG="${DIJET_CONFIG_PATH}/binning_AA.config"
export TNUPLE_DATA_FILE="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_r0${conesize}_${AUAU_DATA_NAME}.root"
export TNUPLE_SIM_FILE_JET10="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet10_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET20="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet20_${AUAU_SIM_NAME}.root"
export TNUPLE_SIM_FILE_JET30="${DIJET_TNTUPLE_PATH}/TNTUPLE_DIJET_SIM_r0${conesize}_jet30_${AUAU_SIM_NAME}.root"

# Reference (previous analysis version) used by compareSystematics_AA.C.
V0_REF="${V0_REF:-/home/tmengel/PPG14/version0/v001_20260715}"
compare_range=1

probability_file="${DIJET_UNFOLDING_PATH}/unfolding_hists/probability_hists_AA_r0${conesize}.root"
if [[ -f "$probability_file" ]]; then
    echo "[cent ${cent}] Missing shared probability file: $probability_file" >&2
    echo "[cent ${cent}] Build it once, serially, before launching concurrent jobs:" >&2
    echo "  root -l -q -b \"getProbability.C(${conesize}, 0, \\\"\${AUAU_CONFIG}\\\")\"" >&2
    exit 1
fi

log() { echo "[cent ${cent} $(date +%H:%M:%S)] $*"; }

unfolding_preload_file="${DIJET_UNFOLDING_PATH}/unfolding_hists/unfolding_preload_AA_r0${conesize}.root"
if [[ ! -s "$unfolding_preload_file" ]]; then
    log "=== make unfolding preload ==="
    root -l -q -b "makeUnfoldingHists.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\",\"${AUAU_DATA_FILE}\" )"
fi
background_file="${DIJET_UNFOLDING_PATH}/unfolding_hists/background_hists_AA_cent${cent}_r0${conesize}.root"
if [[ ! -s "$background_file" ]]; then
    log "=== get background ==="
    root -l -b -q "getBackground.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", \"${AUAU_DATA_FILE}\")"
fi


log "=== primer1: build response, unfold, centrality reweight ==="
root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 1)"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 1)"
root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

log "=== primer2: build response, unfold, centrality reweight ==="
root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 2)"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent}, 2)"

log "=== nominal: build response, validate reweighting, unfold, stat uncertainties ==="
root -l -q -b "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\", 0, 10, ${conesize}, ${cent}, 0)"
root -l -q -b "validateReweighting_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"
root -l -q -b "unfoldData_noempty_AA.cxx(\"${AUAU_CONFIG}\", 10, ${conesize}, ${cent})"
root -l -q -b "unfoldDataUncertainties_noempty_AA.cxx(10, ${conesize}, ${cent})"

log "=== JES / JER systematics ==="
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_negJES_AA.config"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_posJES_AA.config"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_negJER_AA.config"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_posJER_AA.config"

log "=== combinatoric (COMB) flow-modulation systematic ==="
root -l -b -q "makeFlowModulationPreload_AA.C(${conesize}, ${cent} \
    , \"${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config\" \
    , \"${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config\")"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_COMBDown_AA.config"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_COMBUp_AA.config"

log "=== prior / inclusive systematics ==="
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_prior_AA.config"
bash run_all_sys_AA.sh "${conesize}" "${cent}" "${DIJET_CONFIG_PATH}/binning_Inclusive_AA.config"

log "=== full + half closure responses and tests ==="
root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",2,10,${conesize},${cent},0)"
root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",1,10,${conesize},${cent},0)"
root -l -b -q "drawClosureTests_AA.C(${conesize},${cent},${nominal_iter},\"${AUAU_CONFIG}\")"

log "=== plots: systematics, iteration scan, final unfold, prior, dphi/COMB modulation ==="
root -l -b -q "drawSys_AA.C(${conesize}, ${cent})"
root -l -q -b "makeIterationPlot_AA.C(${conesize},${cent},0)"
root -l -q -b "compareSystematics_AA.C(${compare_range}, ${nominal_iter} \
    , \"${DIJET_UNFOLDING_PATH}/uncertainties/systematics_AA_cent_${cent}_r0${conesize}.root\" \
    , \"${V0_REF}/uncertainties/systematics_AA_cent_${cent}_r0${conesize}.root\" \
    , \"${V0_REF}/uncertainties/uncertainties_AA_cent_${cent}_r0${conesize}_nominal.root\" \
    , \"${V0_REF}/unfolding_hists\" \
    , \"${DIJET_UNFOLDING_PATH}/comparison_plots\" \
    , ${cent})"
root -l -b -q "drawFinalUnfold_AA.C(${conesize},${cent}, \"${AUAU_CONFIG}\")"
root -l -b -q "drawFinalUnfold_AA_only.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\", ${nominal_iter})"
root -l -b -q "drawFinalUnfold_AA_cent.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"
root -l -q -b "drawPrior_AA.C(${conesize}, ${cent})"

root -l -b -q "drawResponse_AA.C(${conesize}, ${cent})"
root -l -b -q "drawFullClosure_AA.C(${conesize}, ${cent})"
root -l -b -q "drawHalfClosure_AA.C(${conesize}, ${cent})"

log "=== COMB modulation plots ==="

force_sim_cache=false
[[ "${FORCE_SIM_CACHE:-0}" == "1" ]] && force_sim_cache=true
root -l -b -q "makeCOMBSimulationCache_AA.C(${conesize},-1,\"${AUAU_CONFIG}\",${force_sim_cache})" 

root -l -b -q "drawCOMBModulation_AA_v2.C(${conesize},${cent},\"${AUAU_CONFIG}\")"
root -l -b -q "drawCOMBModulationSimSimple_AA.C(${conesize},${cent},\"${AUAU_CONFIG}\")"
root -l -b -q "drawCOMBModulationDataSimple_AA.C(${conesize},${cent},\"${AUAU_CONFIG}\")"
root -l -b -q "drawCOMBModulation_AA.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"
root -l -b -q "drawCOMBDataSimSimple_AA_v2.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"
root -l -b -q "drawCOMBModulation_AA_v2.C(${conesize}, ${cent}, \"${AUAU_CONFIG}\")"

log "=== done ==="
