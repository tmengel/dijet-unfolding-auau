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
// qq, qg/gg, and qq/qg+gg-mix flavor-tagged responses (see
// createResponse_exclusive_v2_AA.cxx's FLAVOR handling, dijet_matching_flavor.C)
// against the nominal (flavor-inclusive) one, at a single Bayesian iteration.
// Deliberately not part of drawSys_AA.C's systematic-uncertainty machinery --
// this is a physics robustness check, not an up/down uncertainty.
//
// Reads unfolding_hists_<system>_r0<cone>_<sys>.root for sys in
// {nominal, QQ, QGGG, MIX66, MIX80} (written by unfoldData_noempty_AA.cxx)
// and re-projects the saved flat (pT1,pT2) unfolded histogram to x_J with
// the same histo_opps helpers unfoldData_noempty_AA.cxx and drawSys_AA.C
// already use, since the x_J projection itself isn't saved to file. The
// MIX66/MIX80 points are the qq/qg+gg mix cross-check
// (run_flavor_sys_AA_exclusive.sh's "mix" mode) at 66% and 80% QQ share --
// add/remove entries in the `flavors` array below to change which mix
// points (or flavor variants) this overlay shows.
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

  struct Flavor { std::string sys; std::string label; int color; int marker; };
  const std::array<Flavor, 4> flavors = {
    Flavor{"QQ",    "qq dijets",     kRed + 1,    21},
    Flavor{"QGGG",  "qg /gg dijets", kAzure - 6,  22},
    Flavor{"MIX66", "66% qq mix",    kGreen + 2,  23},
    Flavor{"MIX80", "80% qq mix",    kOrange + 7, 33}
  };

  if (!h_nom)
  {
    std::cerr << "compareFlavorXj_AA: missing required unfolded x_J input nominal, aborting." << std::endl;
    return;
  }
  std::vector<TH1D *> h_flavors;
  for (const auto &f : flavors)
  {
    TH1D *h = load_xj(f.sys, f.color, f.marker);
    if (!h)
    {
      std::cerr << "compareFlavorXj_AA: missing required unfolded x_J input " << f.sys
                << ", aborting." << std::endl;
      return;
    }
    h_flavors.push_back(h);
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
  for (TH1D *h : h_flavors) { h->Draw("p same"); }

  TLegend *leg = new TLegend(0.62, 0.60, 0.85, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_nom, "Nominal (inclusive)", "p");
  for (std::size_t i = 0; i < flavors.size(); i++) { leg->AddEntry(h_flavors[i], flavors[i].label.c_str(), "p"); }
  leg->Draw();

  p2->cd();
  std::vector<TH1D *> r_flavors;
  for (std::size_t i = 0; i < flavors.size(); i++)
  {
    TH1D *r = (TH1D *)h_flavors[i]->Clone(Form("r_%s", flavors[i].sys.c_str()));
    r->Divide(h_nom);
    r_flavors.push_back(r);
  }

  r_flavors[0]->SetTitle(";x_{J};Flavor / Nominal");
  r_flavors[0]->SetMinimum(0.0);
  r_flavors[0]->SetMaximum(2.0);
  r_flavors[0]->GetXaxis()->SetRangeUser(ixj_bins[0], ixj_bins[nbins]);
  dlutility::SetFont(r_flavors[0], 42, 0.09, 0.08, 0.08, 0.08);
  r_flavors[0]->Draw("p");
  for (std::size_t i = 1; i < r_flavors.size(); i++) { r_flavors[i]->Draw("p same"); }

  TLine *line = new TLine(ixj_bins[0], 1.0, ixj_bins[nbins], 1.0);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->Draw("same");

  const TString outpath = Form("%s/unfolding_plots/flavor_compare_xj_%s_r%02d_iter%d.pdf",
                                rb.get_code_location().c_str(), system_string.c_str(), cone, niter);
  c->Print(outpath);
  std::cout << "Wrote " << outpath << std::endl;
}
