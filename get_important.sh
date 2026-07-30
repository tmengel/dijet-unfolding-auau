#!/bin/bash


indir=/home/tmengel/PPG14/dijet-unfolding-final
subdirs=(
    closure_plots
    final_plots
    systematic_plots
)
outdir=new_plots
mkdir -p $outdir
for subdir in "${subdirs[@]}"; do
    echo "Copying important files from $subdir"
    mkdir -p "$outdir/$subdir"
    cp "$indir/$subdir"/*range_1*.pdf "$outdir/$subdir"
done

