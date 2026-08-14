#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

source "$script_dir/setup_env.sh"


prelim_result_dir="/home/tmengel/PPG14/version1/dijet-unfolding-final/preliminary_results"
output_dir="${DIJET_UNFOLDING_PATH}/compare_prelim"
mkdir -p "$output_dir"
conesize=3
for cent in 0 ; do
    prelim_file="${prelim_result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
    updated_file="${DIJET_UNFOLDING_PATH}/final_hists/final_hists_AA_cent_${cent}_r0${conesize}.root"
    for iter in 1; do
        for range in 1; do
            output_file="${output_dir}/compare_final_xj_AA_cent_${cent}_r0${conesize}_iter${iter}_range${range}.pdf"
            root -l -q -b \
                "compare_final_xj.C(\"${prelim_file}\",\"${updated_file}\",${range},${iter},${cent},\"Preliminary\",\"Updated\",\"${output_file}\",false)"
        done
    done
done

