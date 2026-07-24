#!/bin/bash
source setup_env.sh
for centrality in 0 1 2 3; do
    root -l -q -b "makeIterationPlot_AA.C+(3,${centrality},0)" 
    root -l -q -b "drawFinalUnfold_AA.C+(3,${centrality}, \"${AUAU_CONFIG}\")"
    root -l -q -b "drawFinalUnfold_AA_cent.C+(3,${centrality}, \"${AUAU_CONFIG}\")"
    root -l -q -b "drawFinalUnfold_AA_only.C+(3,${centrality}, \"${AUAU_CONFIG}\", 2)"
    # root -l -q -b "drawPrior_AA.C+(3,${centrality})"

#   root -l -q -b "drawResponse_AA.C+(3,${centrality})"

done