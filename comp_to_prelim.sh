#!/usr/bin/env bash
set -euo pipefail

# script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# cd "$script_dir"

# source "$script_dir/setup_env.sh"


prelim_result_dir="/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/preliminary_results"
output_dir="/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/compare_prelim"
mkdir -p "$output_dir"
conesize=3
for cent in 0 1 2 3 ; do
    prelim_file="${prelim_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    updated_file="final_hists/final_hists_AA_cent_${cent}_r0${conesize}.root"
    # Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2, matching the
    # prior_iteration constant in createResponse_noempty_AA.cxx.
    iter=1
    for range in 1; do
        output_file="${output_dir}/compare_final_xj_AA_cent_${cent}_r0${conesize}_iter${iter}_range${range}.pdf"
        root -l -q -b \
            "compare_final_xj.C(\"${prelim_file}\",\"${updated_file}\",${range},${iter},${cent},\"Preliminary\",\"Updated\",\"${output_file}\",false)"
    done
done

