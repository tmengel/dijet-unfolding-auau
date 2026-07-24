#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

cone_size="${1:-3}"
centrality="${2:-0}"
config="${3:-binning_AA.config}"

if (( centrality < 0 || centrality > 3 )); then
    echo "Centrality index must be 0, 1, 2, or 3" >&2
    exit 1
fi

unset DIJET_MAX_MATCHED_EVENTS

echo "Building full-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx(\"${config}\",2,10,${cone_size},${centrality},0)"

echo "Building reproducible half-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx(\"${config}\",1,10,${cone_size},${centrality},0)"

echo "Drawing full and half closure tests for centrality ${centrality}"
root -l -b -q "drawClosureTests_AA.C(${cone_size},${centrality},1,\"${config}\")"
