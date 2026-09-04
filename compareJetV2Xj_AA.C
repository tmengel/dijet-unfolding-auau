#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

#include "dlUtility.h"
#include "histo_opps.h"
#include "PlotUtils.h"
#include "read_binning.h"

// Standalone cross-check plot: overlays the unfolded x_J spectrum for the
// jet-v2 reweighted response(s) (see createResponse_exclusive_v2_AA.cxx's
// JETV2_SCALE block and run_jetv2_sys_AA_exclusive.sh) against the nominal
// (v2 = 0, flat-in-azimuth) one, at a single Bayesian iteration. Mirrors
// compareFlavorXj_AA.C for the jet-v2 "what if" cross-check.
//
// Reads unfolding_hists_<system>_r0<cone>_<sys>.root for sys in {nominal,
// JETV2, ...} (written by unfoldData_noempty_AA.cxx) and re-projects the
// saved flat (pT1,pT2) unfolded histogram to x_J with the same histo_opps
// helpers unfoldData_noempty_AA.cxx and drawSys_AA.C already use, since the
// x_J projection itself isn't saved to file.
//
// The default `jetv2s` entry below, sys = "JETV2", is the single cross-check
// point already produced from configs/binning_JETV2_AA.config's default
// JETV2_SCALE: 0.03 (run directly, not through the wrapper). To compare
// additional v2 values, run
//   run_jetv2_sys_AA_exclusive.sh <conesize> <cent> <v2>
// which names its output sys_name JETV2<pct> (e.g. JETV25 for v2 = 0.05),
// then add a matching entry to the `jetv2s` array below.
//
// Like the flavor cross-check this is a "what if", NOT a measured up/down
// uncertainty, so it is deliberately not wired into drawSys_AA.C's total
// band.
void compareJetV2Xj_AA(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string configfile = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config",
  const int niter = 1
)
{
  PlotUtils::set_sphenix_style();
  gStyle->SetOptStat(0);

  read_binning rb(configfile.c_str());
  const std::string system_string = rb.get_system_string(centrality_bin);
  const int cone = cone_size;

  const int nbins = rb.get_nbins();
  float ipt_bins[nbins + 1];
  float ixj_bins[nbins + 1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  const int measure_leading_bin = rb.get_measure_leading_bin();
  const int measure_subleading_bin = rb.get_measure_subleading_bin();

  auto load_xj = [&](const std::string &sys_name, const int color, const int marker) -> TH1D *
  {
    const TString path = Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_%s.root",
                               rb.get_code_location().c_str(), system_string.c_str(), cone, sys_name.c_str());
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie())
    {
      std::cerr << "Missing " << path << std::endl;
      return nullptr;
    }

    TH1D *hflat = (TH1D *)f->Get(Form("h_flat_unfold_pt1pt2_%d", niter));
    if (!hflat)
    {
      std::cerr << "Missing h_flat_unfold_pt1pt2_" << niter << " in " << path << std::endl;
      f->Close();
      return nullptr;
    }
    hflat = (TH1D *)hflat->Clone(Form("h_flat_unfold_pt1pt2_%s", sys_name.c_str()));
    hflat->SetDirectory(nullptr);
    f->Close();

    TH2D *h2 = new TH2D(Form("h_pt1pt2_unfold_%s", sys_name.c_str()), "", nbins, ipt_bins, nbins, ipt_bins);
    histo_opps::make_sym_pt1pt2(hflat, h2, nbins);

    TH1D *hxj = new TH1D(Form("h_xj_unfold_%s", sys_name.c_str()), ";x_{J};#frac{1}{N}#frac{dN}{dx_{J}}", nbins, ixj_bins);
    histo_opps::project_xj(h2, hxj, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    histo_opps::normalize_histo(hxj, nbins);

    dlutility::SetLineAtt(hxj, color, 2, 1);
    dlutility::SetMarkerAtt(hxj, color, 1.1, marker);
    dlutility::SetFont(hxj, 42, 0.05, 0.045, 0.045, 0.045);

    return hxj;
  };

  TH1D *h_nom = load_xj("nominal", kBlack, 20);

  struct JetV2 { std::string sys; std::string label; int color; int marker; };
  const std::array<JetV2, 1> jetv2s = {
    JetV2{"JETV2", "v_{2} = 0.03 reweight", kMagenta + 1, 21}
  };

  if (!h_nom)
  {
    std::cerr << "compareJetV2Xj_AA: missing required unfolded x_J input nominal, aborting." << std::endl;
    return;
  }
  std::vector<TH1D *> h_jetv2s;
  for (const auto &j : jetv2s)
  {
    TH1D *h = load_xj(j.sys, j.color, j.marker);
    if (!h)
    {
      std::cerr << "compareJetV2Xj_AA: missing required unfolded x_J input " << j.sys
                << ", aborting." << std::endl;
      return;
    }
    h_jetv2s.push_back(h);
  }

  TCanvas *c = new TCanvas("cJetV2Xj", "cJetV2Xj", 650, 750);
  TPad *p1 = new TPad("p1", "p1", 0, 0.35, 1, 1);
  p1->SetBottomMargin(0.02);
  p1->Draw();
  TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.35);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.32);
  p2->Draw();

  p1->cd();
  h_nom->SetTitle(Form(";;#frac{1}{N}#frac{dN}{dx_{J}}   (%s, R=0.%d, iter %d)", system_string.c_str(), cone, niter));
  h_nom->GetXaxis()->SetRangeUser(ixj_bins[0], ixj_bins[nbins]);
  h_nom->SetMinimum(0.0);
  h_nom->GetXaxis()->SetLabelSize(0);
  h_nom->Draw("p");
  for (TH1D *h : h_jetv2s) { h->Draw("p same"); }

  TLegend *leg = new TLegend(0.62, 0.68, 0.85, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_nom, "Nominal (v_{2} = 0)", "p");
  for (std::size_t i = 0; i < jetv2s.size(); i++) { leg->AddEntry(h_jetv2s[i], jetv2s[i].label.c_str(), "p"); }
  leg->Draw();

  p2->cd();
  std::vector<TH1D *> r_jetv2s;
  for (std::size_t i = 0; i < jetv2s.size(); i++)
  {
    TH1D *r = (TH1D *)h_jetv2s[i]->Clone(Form("r_%s", jetv2s[i].sys.c_str()));
    r->Divide(h_nom);
    r_jetv2s.push_back(r);
  }

  r_jetv2s[0]->SetTitle(";x_{J};Jet-v2 / Nominal");
  r_jetv2s[0]->SetMinimum(0.0);
  r_jetv2s[0]->SetMaximum(2.0);
  r_jetv2s[0]->GetXaxis()->SetRangeUser(ixj_bins[0], ixj_bins[nbins]);
  dlutility::SetFont(r_jetv2s[0], 42, 0.09, 0.08, 0.08, 0.08);
  r_jetv2s[0]->Draw("p");
  for (std::size_t i = 1; i < r_jetv2s.size(); i++) { r_jetv2s[i]->Draw("p same"); }

  TLine *line = new TLine(ixj_bins[0], 1.0, ixj_bins[nbins], 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->Draw("same");

  const TString outpath = Form("%s/unfolding_plots/jetv2_compare_xj_%s_r%02d_iter%d.pdf",
                                rb.get_code_location().c_str(), system_string.c_str(), cone, niter);
  c->Print(outpath);
  std::cout << "Wrote " << outpath << std::endl;
}
