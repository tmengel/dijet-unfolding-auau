#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

conesize="${1:-3}"
cent="${2:-0}"
config="${3:-binning_AA.config}"

mkdir -p logs final_plots systematic_plots uncertainties

system="AA_cent_${cent}"
required=(
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_nominal.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_posJES.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_negJES.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_posJER.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_negJER.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_COMBDown.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_COMBUp.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_PRIOR.root"
    "unfolding_hists/unfolding_hists_${system}_r0${conesize}_INCLUSIVE.root"
    "response_matrices/response_matrix_${system}_r0${conesize}_nominal.root"
    "uncertainties/uncertainties_${system}_r0${conesize}_nominal.root"
)

missing=0
for path in "${required[@]}"; do
    if [[ ! -f "$path" ]]; then
        echo "Missing required draw input: $path" >&2
        missing=1
    fi
done

if (( missing )); then
    echo "Run the AA unfolding/systematics first, or set RUN_AA_RECOMPUTE=1 to recompute from this script." >&2
    if [[ "${RUN_AA_RECOMPUTE:-0}" != "1" ]]; then
        exit 1
    fi
fi

if [[ "${RUN_AA_RECOMPUTE:-0}" == "1" ]]; then
    bash run_all_unfold_AA.sh "$conesize" "$cent"
fi

if [[ "${RUN_AA_LEGACY_DRAWSYS:-1}" == "1" ]]; then
    echo "Drawing AA systematics with coherent COMBDown/COMBUp flow variations"
    root -l -b -q "drawSys_AA.C(${conesize}, ${cent})" > "logs/drawSys_AA_cent${cent}.log" 2>&1
fi

echo "Drawing aggregate delta-phi COMB modulation diagnostic"
root -l -b -q "drawCOMBModulation_AA.C(${conesize}, ${cent})" > "logs/drawCOMBModulation_AA_cent${cent}.log" 2>&1

echo "Drawing AA-only final plots for cone=${conesize}, cent=${cent}"
root -l -b -q "drawFinalUnfold_AA_only.C(${conesize}, ${cent}, \"${config}\")" > "logs/drawFinalUnfold_AA_cent${cent}.log" 2>&1

echo "Done. Key outputs:"
find final_plots systematic_plots uncertainties -maxdepth 1 -type f \
    \( -name "*AA_cent_${cent}_r0${conesize}*" -o -name "systematics_AA_cent_${cent}_r0${conesize}.root" \) \
    -printf "%p\n" | sort | tail -40
