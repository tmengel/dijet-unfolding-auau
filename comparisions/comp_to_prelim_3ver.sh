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
v1_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/preliminary_results"
v1_label="Preliminary"

v2_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/results/no_cent_cut_2026-08-27_v01/final_hists"
v2_label="No Cent Cut"

v3_result_dir="${DIJET_UNFOLDING_PATH}/final_hists"
v3_label="Cent Cut"

output_dir="${DIJET_UNFOLDING_PATH}/compare_prelim_3ver"
mkdir -p "$output_dir"
conesize=3
for cent in 0 1 2 3 ; do
    v1_file="${v1_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    v2_file="${v2_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    v3_file="${v3_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    # Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2, matching the
    # prior_iteration constant in createResponse_noempty_AA.cxx.
    iter=1
    for range in 1; do
        output_file="${output_dir}/compare_final_xj_3ver_AA_cent_${cent}_r0${conesize}_iter${iter}_range${range}.pdf"
        root -l -q -b \
            "compare_final_xj_3ver.C(\"${v1_file}\",\"${v2_file}\",\"${v3_file}\",${range},${iter},${cent},\"${v1_label}\",\"${v2_label}\",\"${v3_label}\",\"${output_file}\",false)"
    done
done
