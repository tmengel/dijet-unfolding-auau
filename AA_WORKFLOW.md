# Au+Au unfolding workflow

This directory is configured for the 2026-07-20 R=0.3 inputs:

- data: `/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root`
- simulation: `/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root`
- pair probability: `/home/tmengel/PPG14/rootfiles/probs.root`
- JES/JER response calibration: `jer/jer_smear_functions.root`

The simulation adapter treats every event as minimum-bias, uses the R=0.3
truth collection for both the truth dijet and `maxpttruth`, and splits the
single input tree by `jetid` 10, 20, and 30.  It deliberately preserves the
legacy one-to-one, first-match jet-matching behavior.

## Production order

```bash
source setup_env.sh
make clean && make -j"$(nproc)"

root -l -b -q 'prepareProbabilityCorrections_AA.C+(3,"/home/tmengel/PPG14/rootfiles/probs.root","binning_AA.config")'
root -l -b -q 'makeDataTreeAuAu.C+(3,1,"/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root","rootfiles/TNTUPLE_DIJET_AA_r03.root")'
root -l -b -q 'makeMatchedTreesInclusiveAuAu.C+(3,"/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root","rootfiles")'

FORCE=1 JOBS=4 ./make_all_unfold_hists.sh 3
for cent in 0 1 2 3; do
  root -l -b -q "getBackground.C+(3,${cent},\"binning_AA.config\",\"/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root\")"
done

for cent in 0 1 2 3; do
  ./run_all_unfold_AA.sh 3 "$cent" > "logs/full_unfold_AA_cent${cent}.log" 2>&1 &
done
wait

for cent in 0 1 2 3; do
  ./run_everything_AA.sh 3 "$cent"
done
```

`makeUnfoldingHists.C` is the current equivalent of the historical
`makeAuAuHistos` stage.  `validateReweighting_AA.C` produces data versus
simulation comparisons before and after the event weights for vertex,
centrality, and sum-ET.

## Systematics

- JES: `binning_negJES_AA.config`, `binning_posJES_AA.config`
- JER: `binning_negJER_AA.config`, `binning_posJER_AA.config`
- Prior: `binning_prior_AA.config`
- Inclusive: `binning_Inclusive_AA.config`
- COMB: coherent 0.7 and 1.3 scaling of both fitted `v22` and `v33`, configured
  by `binning_COMBDown_AA.config` and `binning_COMBUp_AA.config`

`makeFlowModulationPreload_AA.C` writes both COMB signal inputs and the fit,
coefficient, and signal diagnostics. `drawSys_AA.C` compares COMBDown and
COMBUp after projection and normalization, then selects the larger absolute
shift in each final xJ bin for combination with the other sources.
Bayesian iterations are retained as unfolding outputs, but variation of the
iteration count is not included as a systematic uncertainty.
The former zero-yield-at-minimum systematic is not used.

The preload fit explicitly releases both harmonic parameters after a sparse
bin.  This is a minimal correctness fix to the legacy reused-`TF1` state: the
old loop released `v22` but could leave `v33` fixed to zero for every later
bin after the first sparse fit.

## Production outputs

- `final_plots/`: twelve Au+Au-only final spectra (three leading-jet ranges
  for each of four centralities), as PDF, PNG, and per-centrality ROOT files.
- `systematic_plots/`: individual systematic comparisons and the multipage
  flow-fit/COMB diagnostics.
- `uncertainties/`: statistical-toy and combined/component systematic ROOT
  files.
- `unfolding_plots/`: response diagnostics and vertex, centrality, and sum-ET
  comparisons before and after reweighting, plus one aggregate COMB-modulated
  delta-phi plot per centrality from `drawCOMBModulation_AA.C`.
- `logs/full_unfold_AA_cent*.log` and `logs/full_draw_AA_cent*.log`: complete
  production records.
- `obsolete/`: retired files retained for provenance; these are not used by
  any active entry point.

## Proof mode

For a fast plumbing test, set limits before running individual stages:

```bash
export DIJET_MAX_DATA_EVENTS=1000000
export DIJET_MAX_SIM_EVENTS=100000
export DIJET_MAX_MATCHED_EVENTS=100000
```

Unset all three variables for production.

## Closure tests

The closure driver uses the current matched Au+Au simulation files and nominal
binning. Full closure trains and tests on the complete simulation sample. Half
closure uses a fixed-seed 50/50 split, with disjoint training and pseudo-data
events. The plotted result uses two Bayesian iterations, matching the final
Au+Au result; additional iteration curves are diagnostic only and are not a
systematic uncertainty.

```bash
for cent in 0 1 2 3; do
  ./run_closure_tests_AA.sh 3 "$cent" > "logs/closure_AA_cent${cent}.log" 2>&1 &
done
wait
```

Outputs are written to `closure_plots/` and `closure_results/`.
