#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

source "$script_dir/setup_env.sh"

# Different result directories were saved in different formats -- some only
# kept the raw "final_hists_*.root" (every Bayesian iteration), others only
# kept the "final_plots/.../final_xj_*_AAonly.root" (already-chosen final
# iteration). Set each version's *_source to "hists" or "plots" to match
# what that directory actually has; *_iter is only used when source=hists.
v2_result_dir="${script_dir}/results/final_2026-08-18_v01"
v2_label="IAN"
v2_source="plots"   # "hists" or "plots"
v2_iter=1            # only used when v1_source=hists

v1_result_dir="${script_dir}/preliminary_results"
v1_label="Preliminary"
v1_source="hists"
v1_iter=1

v3_result_dir="${DIJET_UNFOLDING_PATH}"
v3_label="Current"
v3_source="plots"
v3_iter=1

output_dir="${DIJET_UNFOLDING_PATH}/compare_prelim_3ver_finalplots"
mkdir -p "$output_dir"
conesize=3

# The final_plots/ and final_hists/ trees are inconsistently nested across
# runs -- e.g. final_plots/final_plots/final_xj_*.root vs a flat
# final_plots/final_xj_*.root, and similarly for final_hists/ -- so try both
# a nested and a flat layout, plus the file living directly in result_dir.
resolve_version_file() {
    local result_dir="$1"
    local source="$2"
    local cent="$3"
    local candidates=()
    case "$source" in
        plots)
            candidates=(
                "${result_dir}/final_plots/final_plots/final_xj_AA_cent_${cent}_r0${conesize}_AAonly.root"
                "${result_dir}/final_plots/final_xj_AA_cent_${cent}_r0${conesize}_AAonly.root"
            )
            ;;
        hists)
            candidates=(
                "${result_dir}/final_hists/final_hists/final_hists_AA_cent_${cent}_r0${conesize}.root"
                "${result_dir}/final_hists/final_hists_AA_cent_${cent}_r0${conesize}.root"
                "${result_dir}/final_hists_AA_cent_${cent}_r0${conesize}.root"
            )
            ;;
        *)
            echo "ERROR: unknown source type '${source}' (expected 'hists' or 'plots')" >&2
            return 1
            ;;
    esac
    for candidate in "${candidates[@]}"; do
        if [[ -f "$candidate" ]]; then
            echo "$candidate"
            return 0
        fi
    done
    echo "ERROR: could not find a '${source}' x_J file for cent ${cent} under ${result_dir}" >&2
    return 1
}

for cent in 0 1 2 3; do
    v1_file="$(resolve_version_file "$v1_result_dir" "$v1_source" "$cent")"
    v2_file="$(resolve_version_file "$v2_result_dir" "$v2_source" "$cent")"
    v3_file="$(resolve_version_file "$v3_result_dir" "$v3_source" "$cent")"

    # compare_final_xj_3ver_finalplots.C reads a version at a fixed iteration
    # (iter_vN >= 0) only when that version's source is "hists"; pass -1 for
    # "plots" sources, whose files already store just the final iteration.
    [[ "$v1_source" == "hists" ]] && iv1="$v1_iter" || iv1=-1
    [[ "$v2_source" == "hists" ]] && iv2="$v2_iter" || iv2=-1
    [[ "$v3_source" == "hists" ]] && iv3="$v3_iter" || iv3=-1

    for range in 1; do
        output_file="${output_dir}/compare_final_xj_3ver_finalplots_AA_cent_${cent}_r0${conesize}_range${range}.pdf"
        root -l -q -b \
            "compare_final_xj_3ver_finalplots.C(\"${v1_file}\",\"${v2_file}\",\"${v3_file}\",${range},${cent},\"${v1_label}\",\"${v2_label}\",\"${v3_label}\",\"${output_file}\",${iv1},${iv2},${iv3})"
    done
done
