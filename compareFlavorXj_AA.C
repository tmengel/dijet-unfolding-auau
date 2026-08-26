#include <iostream>
#include <string>

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
// qq and qg/gg flavor-tagged responses (see createResponse_exclusive_AA.cxx's
// FLAVOR handling, dijet_matching_flavor.C) against the nominal
// (flavor-inclusive) one, at a single Bayesian iteration. Deliberately not
// part of drawSys_AA.C's systematic-uncertainty machinery -- this is a
// physics robustness check, not an up/down uncertainty.
//
// Reads unfolding_hists_<system>_r0<cone>_<sys>.root for sys in
// {nominal, QQ, QGGG} (written by unfoldData_noempty_AA.cxx) and re-projects
// the saved flat (pT1,pT2) unfolded histogram to x_J with the same
// histo_opps helpers unfoldData_noempty_AA.cxx and drawSys_AA.C already use,
// since the x_J projection itself isn't saved to file.
void compareFlavorXj_AA(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string configfile = "configs/binning_AA.config",
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
  TH1D *h_qq = load_xj("QQ", kRed + 1, 21);
  TH1D *h_qggg = load_xj("QGGG", kAzure - 6, 22);

  if (!h_nom || !h_qq || !h_qggg)
  {
    std::cerr << "compareFlavorXj_AA: missing at least one required unfolded x_J input, aborting." << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("cFlavorXj", "cFlavorXj", 650, 750);
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
  h_qq->Draw("p same");
  h_qggg->Draw("p same");

  TLegend *leg = new TLegend(0.66, 0.72, 0.85, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_nom, "Nominal (inclusive)", "p");
  leg->AddEntry(h_qq, "qq dijets", "p");
  leg->AddEntry(h_qggg, "qg /gg dijets", "p");
  leg->Draw();

  p2->cd();
  TH1D *r_qq = (TH1D *)h_qq->Clone("r_qq");
  r_qq->Divide(h_nom);
  TH1D *r_qggg = (TH1D *)h_qggg->Clone("r_qggg");
  r_qggg->Divide(h_nom);

  r_qq->SetTitle(";x_{J};Flavor / Nominal");
  r_qq->SetMinimum(0.0);
  r_qq->SetMaximum(2.0);
  r_qq->GetXaxis()->SetRangeUser(ixj_bins[0], ixj_bins[nbins]);
  dlutility::SetFont(r_qq, 42, 0.09, 0.08, 0.08, 0.08);
  r_qq->Draw("p");
  r_qggg->Draw("p same");

  TLine *line = new TLine(ixj_bins[0], 1.0, ixj_bins[nbins], 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->Draw("same");

  const TString outpath = Form("%s/unfolding_plots/flavor_compare_xj_%s_r%02d_iter%d.pdf",
                                rb.get_code_location().c_str(), system_string.c_str(), cone, niter);
  c->Print(outpath);
  std::cout << "Wrote " << outpath << std::endl;
}
