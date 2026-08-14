#include "histo_opps.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Usage (after setup_env.sh):
//   root -l -b -q 'priorClosureTest_AA.C("response_matrices/response_matrix_AA_cent_0_r03_nominal.root",
//                                        "unfolding_hists/unfolding_hists_AA_cent_0_r03_nominal.root", 1)'
// niter_index = Bayes iterations - 1 (1 -> the analysis default of 2 iterations).
//
// Data-driven reweighted-closure estimate of the prior/unfolding systematic:
// reweight the MC truth fully (fraction 1.0) to the unfolded data shape, fold it
// through the nominal response, unfold with the nominal (unreweighted-prior)
// response at the analysis iteration count, and take the non-closure vs the
// reweighted truth. No arbitrary amplitude enters: the target is the data itself.
void priorClosureTest_AA(const char* respfile, const char* unffile, int niter_index = 1)
{
  const int nbins = 15;
  float ipt[16] = {7, 8.39655, 10.0717, 12.0811, 14.4914, 17.3825, 20.8505, 25.0103, 30, 35.9852, 43.1645, 51.7762, 62.1059, 74.4965, 89.3591, 107.187};
  float ixj[16] = {0.0653065, 0.0783356, 0.0939641, 0.112711, 0.135197, 0.16217, 0.194524, 0.233333, 0.279885, 0.335724, 0.402703, 0.483046, 0.579417, 0.695015, 0.833676, 1};
  TFile *fr = TFile::Open(respfile, "READ");
  TFile *fu = TFile::Open(unffile, "READ");
  auto *resp = (RooUnfoldResponse*) fr->Get("response_noempty");
  auto *h_truth_skim = (TH1D*) fr->Get("h_flat_truth_skim");
  auto *h_tmap = (TH1D*) fr->Get("h_flat_truth_mapping");
  auto *h_truth_full = (TH1D*) fr->Get("h_truth_flat_pt1pt2");     // full truth (test), complete after fix
  auto *h_unf = (TH1D*) fu->Get(Form("h_flat_unfold_pt1pt2_%d", niter_index)); // data unfold, full space
  if (!resp||!h_truth_skim||!h_tmap||!h_truth_full||!h_unf) { printf("missing input\n"); return; }

  // full-space weights w(k) = normalized unfold / normalized truth
  double su = h_unf->Integral(), st = h_truth_full->Integral();
  // reweighted truth in skim space
  TH1D *h_tprime_skim = (TH1D*) h_truth_skim->Clone("h_tprime_skim");
  for (int k = 0; k < nbins*nbins; k++)
  {
    int skim = (int) h_tmap->GetBinContent(k+1);
    if (skim <= 0) continue;
    double t = h_truth_full->GetBinContent(k+1)/st;
    double u = h_unf->GetBinContent(k+1)/su;
    double w = (t > 0) ? u/t : 1.0;
    h_tprime_skim->SetBinContent(skim, h_truth_skim->GetBinContent(skim)*w);
  }
  // fold and unfold
  TH1D *h_dprime = (TH1D*) resp->ApplyToTruth(h_tprime_skim, "h_dprime");
  RooUnfoldBayes unfold(resp, h_dprime, niter_index + 1);
  TH1D *h_uprime_skim = (TH1D*) unfold.Hunfold();
  // back to full space
  TH1D *h_tprime_full = (TH1D*) h_truth_full->Clone("h_tprime_full"); h_tprime_full->Reset();
  TH1D *h_uprime_full = (TH1D*) h_truth_full->Clone("h_uprime_full"); h_uprime_full->Reset();
  for (int k = 0; k < nbins*nbins; k++)
  {
    int skim = (int) h_tmap->GetBinContent(k+1);
    if (skim <= 0) continue;
    h_tprime_full->SetBinContent(k+1, h_tprime_skim->GetBinContent(skim));
    h_uprime_full->SetBinContent(k+1, h_uprime_skim->GetBinContent(skim));
  }
  // project both to xJ, range 1 (30-43.2), and range edges for context
  for (int irange = 0; irange < 3; irange++)
  {
    int lo[3] = {6, 8, 10}, hi[3] = {8, 10, 12};
    TH2D *h2t = new TH2D(Form("h2t%d",irange), "", nbins, ipt, nbins, ipt);
    TH2D *h2u = new TH2D(Form("h2u%d",irange), "", nbins, ipt, nbins, ipt);
    histo_opps::make_sym_pt1pt2(h_tprime_full, h2t, nbins);
    histo_opps::make_sym_pt1pt2(h_uprime_full, h2u, nbins);
    TH1D *hxt = new TH1D(Form("hxt%d",irange), "", nbins, ixj);
    TH1D *hxu = new TH1D(Form("hxu%d",irange), "", nbins, ixj);
    histo_opps::project_xj(h2t, hxt, nbins, lo[irange], hi[irange], 2, nbins-2);
    histo_opps::project_xj(h2u, hxu, nbins, lo[irange], hi[irange], 2, nbins-2);
    histo_opps::normalize_histo(hxt, nbins);
    histo_opps::normalize_histo(hxu, nbins);
    printf("range %d non-closure (U'-T')/T':", irange);
    for (int b = 8; b <= nbins; b++)
    {
      double t = hxt->GetBinContent(b);
      printf("  %.3f:%+6.3f", hxt->GetBinCenter(b), t>0 ? (hxu->GetBinContent(b)-t)/t : 0);
    }
    printf("\n");
  }
}
