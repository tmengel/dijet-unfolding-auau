#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"
source "$script_dir/setup_env.sh"

conesize=$1
centrality=$2

root -l -b -q "drawSys_AA.C(${conesize}, ${centrality} )"

root -l -b -q "drawFinalUnfold_AA_only.C(${conesize}, ${centrality}, \"${AUAU_CONFIG}\")"

root -l -b -q "drawCOMBModulation_AA_v2.C(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
root -l -b -q "drawCOMBModulationSimSimple_AA.C(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
root -l -b -q "drawCOMBModulationDataSimple_AA.C(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
root -l -b -q "drawCOMBModulationSimSimple_AA.C(${conesize},${centrality},\"${AUAU_CONFIG}\")" 
root -l -b -q "drawCOMBModulation_AA.C(${conesize}, ${centrality}, \"${AUAU_CONFIG}\")"

echo "Building full-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",2,10,${conesize},${centrality},0)"

echo "Building reproducible half-sample closure response for centrality ${centrality}"
root -l -b -q "createResponse_noempty_AA.cxx(\"${AUAU_CONFIG}\",1,10,${conesize},${centrality},0)"

echo "Drawing full and half closure tests for centrality ${centrality}"
root -l -b -q "drawClosureTests_AA.C(${conesize},${centrality},1,\"${AUAU_CONFIG}\")"

root -l -q -b "makeIterationPlot_AA.C(3,${centrality},0)" 

root -l -q -b "drawFinalUnfold_AA.C(3,${centrality}, \"${AUAU_CONFIG}\")"

# root -l -q -b "drawFinalUnfold_AA_cent.C(3,${centrality}, \"${AUAU_CONFIG}\")"

# root -l -q -b "drawFinalUnfold_AA_only.C(3,${centrality}, \"${AUAU_CONFIG}\", 2)"


