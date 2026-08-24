#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

source "$script_dir/setup_env.sh"


# prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/final_hists_NocentCut_jes4_jerv1_rho_comb_norm0prior0_sumeTZ"
# prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/preliminary_results"
# prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/final_hists_centCut_jes4_jerv1_rho_comb_norm05prior1_sumeTZ_tbin0"
# prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/final_hists_Nocentcut_jes4_jerv1_rho_comb_norm1prior05_sumeTZ_tbin0"
# prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/final_hists_centCut_jes4_jerv1_rho_comb_norm1prior1_sumeTZ_tbin0"
prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/final_hists_NocentCut_jes2"
output_dir="${DIJET_UNFOLDING_PATH}/compare_prelim"
mkdir -p "$output_dir"
conesize=3
for cent in 0 1 2 3 ; do
    prelim_file="${prelim_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    updated_file="${DIJET_UNFOLDING_PATH}/final_hists/final_hists_AA_cent_${cent}_r0${conesize}.root"
    # Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2, matching the
    # prior_iteration constant in createResponse_noempty_AA.cxx.
    iter=1
    for range in 1; do
        output_file="${output_dir}/compare_final_xj_AA_cent_${cent}_r0${conesize}_iter${iter}_range${range}.pdf"
        root -l -q -b \
            "compare_final_xj.C(\"${prelim_file}\",\"${updated_file}\",${range},${iter},${cent},\"Preliminary\",\"Updated\",\"${output_file}\",false)"
    done
done

