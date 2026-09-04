#!/bin/bash
config="/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config"
for c in 0 1 2 3; do root -l -b -q "drawPriorQA_xj_flavorCompare_AA.C(3, ${c}, \"${config}\")" ; done