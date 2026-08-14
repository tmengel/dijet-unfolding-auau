#!/bin/bash

today=$(date +%Y-%m-%d)
name="${1:-output_${today}}"

storage_dir="$(pwd)/results"
input_dir="$(pwd)"

this_save="${storage_dir}/${name}_${today}"
version_tag=1
if [[ -d "$this_save" ]]; then
  #iterative version until we find a free one
  while [[ -d "${this_save}_v$(printf '%02d' "$version_tag")" ]]; do
    version_tag=$((version_tag + 1))
  done
fi
this_save="${this_save}_v$(printf '%02d' "$version_tag")"

mkdir -p "$this_save"
output_dirs=(
  final_plots
  dphi_plots
  logs
  closure_results
  closure_plots
  configs
  response_matrices
  systematic_plots
  truth_hists
  uncertainties
  unfolding_hists
  unfolding_plots
  vertex
  centrality
  sumeT
)

for dir in "${output_dirs[@]}"; do
  mkdir -p "$this_save/$dir"
  cp -r "$input_dir/$dir" "$this_save/$dir"
done
