#/bin/bash


today=$(date +%Y-%m-%d)
name="${1:-${today}}"
storage_dir="$(pwd)/tmp"
input_dir="$(pwd)"

# storage_dir="/home/tmengel/PPG14/dijet-unfolding-final/results/temp"
mkdir -p "$storage_dir"
# input_dir="/home/tmengel/PPG14/dijet-unfolding-final/"

this_save="${storage_dir}/${name}_${today}"
version_tag=1
this_save="${this_save}_v$(printf '%02d' "$version_tag")"


output_dirs=(
  final_plots
  dphi_plots
  logs
  closure_results
  closure_plots
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

mkdir -p "$this_save"
for dir in "${output_dirs[@]}"; do
  mkdir -p "$this_save/$dir"
  cp -r "$input_dir/$dir" "$this_save/$dir"
done

for dir in "${output_dirs[@]}"; do
  rm -rf "$input_dir/$dir"
done