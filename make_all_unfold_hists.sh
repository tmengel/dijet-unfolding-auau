#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

conesize="${CONESIZE:-3}"
config="${CONFIG:-binning_AA.config}"
infile="${INFILE:-/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root}"
jobs="${JOBS:-2}"
force="${FORCE:-0}"

if [[ $# -gt 0 ]]; then
    conesize="$1"
fi

mkdir -p logs unfolding_hists

if [[ "${CENTS:-}" != "" ]]; then
    read -r -a cents <<< "$CENTS"
else
    cents=(0 1 2 3)
fi

run_cent() {
    local cent="$1"
    local output="unfolding_hists/unfolding_hists_preload_AA_cent_${cent}_r0${conesize}_nominal.root"
    local log="logs/makeUnfoldingHists_AA_cent${cent}_r0${conesize}.log"

    if [[ "$force" != "1" && -s "$output" ]]; then
        echo "cent ${cent}: exists, skipping $output"
        return 0
    fi

    echo "cent ${cent}: starting -> $log"
    root -l -b -q "makeUnfoldingHists.C(\"${config}\", ${conesize}, ${cent}, 0, \"${infile}\")" > "$log" 2>&1
    echo "cent ${cent}: done -> $output"
}

active_jobs=0
for cent in "${cents[@]}"; do
    run_cent "$cent" &
    ((++active_jobs))
    if (( active_jobs >= jobs )); then
        wait -n
        ((--active_jobs))
    fi
done

wait

echo "Done. Outputs:"
for cent in "${cents[@]}"; do
    output="unfolding_hists/unfolding_hists_preload_AA_cent_${cent}_r0${conesize}_nominal.root"
    if [[ -f "$output" ]]; then
        ls -lh "$output"
    else
        echo "missing: $output"
    fi
done
