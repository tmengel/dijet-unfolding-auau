#include <array>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"

#include "dlUtility.h"
#include "histo_opps.h"
#include "PlotUtils.h"
#include "priorReweightQA.h"
#include "read_binning.h"

// Standalone comparison macro: reproduces compareFlavorXj_AA.C's inclusive
// x_J overlay of the nominal/qq/qg-gg unfolded data (top pad, data unfolded
// each with its own flavor-tagged response; bottom pad, flavor/nominal
// ratio) and additionally overlays, for qq and qg/gg, the prior_qa MC truth
// (unreweighted) and prior-reweighted truth built to match the nominal
// unfolded data -- the same reweighting target used by
// createResponse_exclusive_v2_AA.cxx (prior_data_sys_name), since the data
// itself carries no flavor tag. Reads the same cached
// response_matrices/unfolding_hists ROOT files that
// createResponse_exclusive_v2_AA.cxx / unfoldData_noempty_AA.cxx /
// compareFlavorXj_AA.C already use (primer-less, i.e. the actual production
// unfold, not the primer1/primer2 QA passes) and reruns only the (cheap)
// prior_qa blend/apply/project steps -- it does not touch or rerun any
// existing macro. Read-only relative to everything already on disk.
namespace
{
// Mirrors the "project" lambda inside prior_qa::draw_xj / compareFlavorXj_AA.C's
// load_xj so the shapes here match what those plots already show.
TH1D *projectToXj(const TH1D *flat, const int nbins, const float *ipt_bins,
                  const float *ixj_bins, const int lead_lo, const int lead_hi,
                  const int sub_lo, const int sub_hi, const char *name)
{
  TH2D h2(Form("h2_%s", name), "", nbins, ipt_bins, nbins, ipt_bins);
  h2.SetDirectory(nullptr);
  histo_opps::make_sym_pt1pt2(const_cast<TH1D *>(flat), &h2, nbins);

  TH1D *hxj = new TH1D(name, "", nbins, ixj_bins);
  hxj->SetDirectory(nullptr);
  histo_opps::project_xj(&h2, hxj, nbins, lead_lo, lead_hi, sub_lo, sub_hi);
  histo_opps::normalize_histo(hxj, nbins);
  return hxj;
}

// A small constant horizontal shift (as a fraction of each bin's width) so
// that data series which land on nearly the same y at a given x_J bin --
// e.g. nominal and a flavor's own data -- don't fully occlude one another.
TGraphErrors *offsetGraph(const TH1D *h, const double xOffsetFrac, const char *name)
{
  const int n = h->GetNbinsX();
  TGraphErrors *g = new TGraphErrors(n);
  g->SetName(name);
  for (int b = 1; b <= n; b++)
    {
      g->SetPoint(b - 1, h->GetBinCenter(b) + xOffsetFrac*h->GetBinWidth(b), h->GetBinContent(b));
      g->SetPointError(b - 1, 0.0, h->GetBinError(b));
    }
  return g;
}

TH1D *loadUnfoldedFlat(const TString &location, const std::string &system,
                       const int cone_size, const std::string &sys_name,
                       const int niter)
{
  const TString path = Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_%s.root",
                            location.Data(), system.c_str(), cone_size, sys_name.c_str());
  std::unique_ptr<TFile> f(TFile::Open(path, "READ"));
  if (!f || f->IsZombie())
    {
      std::cerr << "Cannot open " << path << std::endl;
      return nullptr;
    }
  TH1D *h = dynamic_cast<TH1D *>(f->Get(Form("h_flat_unfold_pt1pt2_%d", niter)));
  if (!h)
    {
      std::cerr << "Missing h_flat_unfold_pt1pt2_" << niter << " in " << path << std::endl;
      return nullptr;
    }
  h = (TH1D *) h->Clone(("h_flat_unfold_pt1pt2_" + sys_name).c_str());
  h->SetDirectory(nullptr);
  return h;
}
}

void drawPriorQA_xj_flavorCompare_AA(
  const int cone_size = 3, const int centrality_bin = 0,
  const std::string config = "binning_AA.config", const int niter = 1)
{
  PlotUtils::set_sphenix_style();
  gStyle->SetOptStat(0);

  read_binning rb(config);
  float centbins[4];
  rb.get_centrality_bins(centbins);
  const float centLo = centbins[centrality_bin];
  const float centHi = centbins[centrality_bin + 1];
  std::string cent_str = Form("%.0f-%.0f%%", centLo, centHi);
  
  const std::string system = rb.get_system_string(centrality_bin);
  const int nbins = rb.get_nbins();
  std::unique_ptr<float[]> ptBins(new float[nbins + 1]);
  rb.get_pt_bins(ptBins.get());
  std::unique_ptr<float[]> xjBins(new float[nbins + 1]);
  rb.get_xj_bins(xjBins.get());
  const int leadLo = rb.get_measure_leading_bin();
  const int leadHi = nbins - 2;
  const int subLo = rb.get_measure_subleading_bin();
  const int subHi = nbins - 2;

  const int measureBins = rb.get_measure_bins();
  std::vector<prior_qa::lead_group> leadGroups;
  for (int i = 0; i < measureBins; i++)
    leadGroups.push_back({rb.get_measure_region(i), rb.get_measure_region(i + 1)});

  const TString location = rb.get_code_location();

  // fraction = 1.0 (prior_norm in createResponse_exclusive_v2_AA.cxx) applies
  // to both QQ and QGGG: PRIOR is its own systematic (prior_sys true,
  // fraction 0.5) and neither QQ nor QGGG is it.
  const double fraction = 1.0;

  // The actual measurement, unfolded once with each flavor's own response --
  // exactly what compareFlavorXj_AA.C plots. The nominal curve additionally
  // serves as the prior_qa reweighting target below, since the data itself
  // carries no flavor tag.
  std::unique_ptr<TH1D> unfoldNominal(loadUnfoldedFlat(location, system, cone_size, "nominal", niter));
  std::unique_ptr<TH1D> unfoldQQ(loadUnfoldedFlat(location, system, cone_size, "QQ", niter));
  std::unique_ptr<TH1D> unfoldQGGG(loadUnfoldedFlat(location, system, cone_size, "QGGG", niter));
  if (!unfoldNominal || !unfoldQQ || !unfoldQGGG)
    {
      std::cerr << "drawPriorQA_xj_flavorCompare_AA: missing at least one required "
                << "unfolded x_J input, aborting." << std::endl;
      return;
    }

  std::unique_ptr<TH1D> hNom(projectToXj(unfoldNominal.get(), nbins, ptBins.get(), xjBins.get(),
                                         leadLo, leadHi, subLo, subHi, "h_xj_unfold_nominal"));
  std::unique_ptr<TH1D> hQQ(projectToXj(unfoldQQ.get(), nbins, ptBins.get(), xjBins.get(),
                                        leadLo, leadHi, subLo, subHi, "h_xj_unfold_QQ"));
  std::unique_ptr<TH1D> hQGGG(projectToXj(unfoldQGGG.get(), nbins, ptBins.get(), xjBins.get(),
                                          leadLo, leadHi, subLo, subHi, "h_xj_unfold_QGGG"));

  dlutility::SetMarkerAtt(hNom.get(), kBlack, 1.1, 20);
  dlutility::SetMarkerAtt(hQQ.get(), kRed + 1, 1.1, 21);
  dlutility::SetMarkerAtt(hQGGG.get(), kAzure - 6, 1.1, 22);
  dlutility::SetFont(hNom.get(), 42, 0.05, 0.045, 0.045, 0.045);

  // dataOffset shifts each flavor's data series by a fraction of the bin
  // width, away from nominal, so the two points that can otherwise coincide
  // at a given x_J bin fan out into a readable pair instead of stacking on
  // top of each other.
  struct Selection
  {
    std::string sys; std::string label; TH1D *data; int color;
    int dataMarker; double dataOffset;
  };
  const std::array<Selection, 2> selections = {
    Selection{"QQ",   "qq dijets",    hQQ.get(),   kRed + 1,   21, -0.15},
    Selection{"QGGG", "qg/gg dijets", hQGGG.get(), kAzure - 6, 22, +0.15}
  };

  std::vector<std::unique_ptr<TH1D>> truthHistograms;
  std::vector<std::unique_ptr<TH1D>> reweightedHistograms;
  std::vector<std::unique_ptr<TGraphErrors>> dataGraphs;

  for (const Selection &selection : selections)
    {
      // PRIMER2, not the primer-less production pass used for the unfolded
      // data curves above: the primer-less h_truth_flat_pt1pt2 for QQ/QGGG
      // turns out to be the same shape as nominal's, just rescaled (no real
      // flavor selection baked into that truth-level histogram) -- only the
      // PRIMER1/PRIMER2 QA passes carry a genuinely flavor-selected truth,
      // which is what the "MC unreweighted" curve is supposed to show.
      const TString truthPath = Form(
        "%s/response_matrices/response_matrix_%s_r%02d_PRIMER2_%s.root",
        location.Data(), system.c_str(), cone_size, selection.sys.c_str());
      std::unique_ptr<TFile> truthFile(TFile::Open(truthPath, "READ"));
      if (!truthFile || truthFile->IsZombie())
        {
          std::cerr << "Cannot open " << truthPath << std::endl;
          continue;
        }
      TH1D *truthFlat = dynamic_cast<TH1D*>(truthFile->Get("h_truth_flat_pt1pt2"));
      if (!truthFlat)
        {
          std::cerr << "Missing h_truth_flat_pt1pt2 in " << truthPath << std::endl;
          continue;
        }

      const prior_qa::weights weights = prior_qa::build_xj(
        truthFlat, unfoldNominal.get(), fraction, nbins, ptBins.get(), xjBins.get(),
        leadGroups, subLo, subHi, "_" + selection.sys);
      if (!weights.truth || !weights.weight)
        {
          std::cerr << "Prior reweighting failed for " << selection.sys << std::endl;
          continue;
        }
      std::unique_ptr<TH1D> reweightedFlat(prior_qa::apply(
        weights.truth, weights.weight, nbins,
        "h_prior_reweighted_truth_" + selection.sys));

      // Same lead/sub range as hNom/hQQ/hQGGG above so the overlay is a
      // like-for-like comparison to the plotted data.
      std::unique_ptr<TH1D> truthXj(projectToXj(
        weights.truth, nbins, ptBins.get(), xjBins.get(), leadLo, leadHi,
        subLo, subHi, Form("h_xj_truth_%s", selection.sys.c_str())));
      std::unique_ptr<TH1D> reweightedXj(projectToXj(
        reweightedFlat.get(), nbins, ptBins.get(), xjBins.get(), leadLo, leadHi,
        subLo, subHi, Form("h_xj_reweighted_%s", selection.sys.c_str())));

      // Truth dashed, reweighted solid -- same color, same line style family
      // as compareFlavorXj_AA.C's dashed step curves, distinguished from each
      // other by dash pattern alone since they're drawn as points.
      dlutility::SetLineAtt(truthXj.get(), selection.color, 2, 2);
      dlutility::SetLineAtt(reweightedXj.get(), selection.color, 2, 1);

      std::unique_ptr<TGraphErrors> dataGraph(offsetGraph(
        selection.data, selection.dataOffset, Form("g_data_%s", selection.sys.c_str())));
      dlutility::SetMarkerAtt(dataGraph.get(), selection.color, 1.1, selection.dataMarker);
      dlutility::SetLineAtt(dataGraph.get(), selection.color, 2, 1);

      truthHistograms.push_back(std::move(truthXj));
      reweightedHistograms.push_back(std::move(reweightedXj));
      dataGraphs.push_back(std::move(dataGraph));
    }

  TCanvas *c = new TCanvas("cFlavorXjPriorQA", "cFlavorXjPriorQA", 650, 750);
  TPad *p1 = new TPad("p1", "p1", 0, 0.35, 1, 1);
  p1->SetBottomMargin(0.02);
  p1->Draw();
  TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.35);
  p2->SetTopMargin(0.03);
  p2->SetBottomMargin(0.32);
  p2->Draw();

  p1->cd();
  double maximum = std::max({hNom->GetMaximum(), hQQ->GetMaximum(), hQGGG->GetMaximum()});
  for (const auto &histogram : truthHistograms)
    maximum = std::max(maximum, histogram->GetMaximum());
  for (const auto &histogram : reweightedHistograms)
    maximum = std::max(maximum, histogram->GetMaximum());

  // The lowest x_J bins carry essentially no signal (see the earlier full-range
  // plot) and just waste frame space, so start the visible window where the
  // distribution actually turns on.
  const double xjPlotMin = 0.2;

  hNom->SetTitle(Form(";;#frac{1}{N}#frac{dN}{dx_{J}}"));
  hNom->GetXaxis()->SetRangeUser(xjPlotMin, xjBins[nbins]);
  // Zoomed in: just enough headroom above the tallest curve for the legend.
  hNom->GetYaxis()->SetRangeUser(0.0, 1.2*maximum);
  hNom->GetXaxis()->SetLabelSize(0);
  hNom->Draw("axis");
  for (const auto &histogram : truthHistograms)
    histogram->Draw("hist same");
  for (const auto &histogram : reweightedHistograms)
    histogram->Draw("hist same");

  std::unique_ptr<TGraphErrors> nomGraph(offsetGraph(hNom.get(), 0.0, "g_data_nominal"));
  dlutility::SetMarkerAtt(nomGraph.get(), kBlack, 0.8, 20);
  nomGraph->Draw("P SAME");
  for (const auto &graph : dataGraphs)
    graph->Draw("P SAME");

  TLegend legNom(0.18, 0.79, 0.4, 0.89);
  legNom.SetBorderSize(0);
  legNom.SetFillStyle(0);
  legNom.SetTextSize(0.032);
  legNom.AddEntry(nomGraph.get(), "Nominal", "pe");
  legNom.AddEntry(dataGraphs[0].get(), "qq response", "pe");
  legNom.AddEntry(dataGraphs[1].get(), "qg & gg response", "pe");
  legNom.Draw();

  TLegend legQQ(0.45, 0.79, 0.6, 0.89);
  legQQ.SetBorderSize(0);
  legQQ.SetFillStyle(0);
  legQQ.SetTextSize(0.032);
  legQQ.SetHeader("qq response");

  TLegend legQGGG(0.7, 0.79, 0.85, 0.89);
  legQGGG.SetBorderSize(0);
  legQGGG.SetFillStyle(0);
  legQGGG.SetTextSize(0.032);
  legQGGG.SetHeader("qg/gg response");



  std::array<TLegend *, 2> flavorLegends = {&legQQ, &legQGGG};
  for (std::size_t i = 0; i < truthHistograms.size(); i++)
    {
      // flavorLegends[i]->AddEntry(dataGraphs[i].get(), "Unfolded w/ response", "pe");
      flavorLegends[i]->AddEntry(truthHistograms[i].get(), "Unweighted Truth", "l");
      flavorLegends[i]->AddEntry(reweightedHistograms[i].get(), "Reweighted Truth", "l");
    }
  legQQ.Draw();
  legQGGG.Draw();

  // std::string cent_str = "";
  // centrality_bin
  // if ( centrality 
  PlotUtils::myText(0.55, 0.1,kBlack, Form("%s, R=0.%d, iter %d", cent_str.c_str(), cone_size, niter), 0.04);

  p2->cd();
  std::unique_ptr<TH1D> rQQ(dynamic_cast<TH1D*>(hQQ->Clone("r_qq")));
  rQQ->Divide(hNom.get());
  std::unique_ptr<TH1D> rQGGG(dynamic_cast<TH1D*>(hQGGG->Clone("r_qggg")));
  rQGGG->Divide(hNom.get());

  rQQ->SetTitle(";x_{J};Unfolded w/ Response / Nominal");
  rQQ->GetYaxis()->SetTitleOffset(1.2);
  rQQ->SetMinimum(0.5);
  rQQ->SetMaximum(1.5);
  rQQ->GetXaxis()->SetRangeUser(xjPlotMin, xjBins[nbins]);
  dlutility::SetFont(rQQ.get(), 42, 0.09, 0.08, 0.08, 0.08);
  rQQ->Draw("axis");

  // TH1's "P" draws a horizontal tick spanning the full bin width regardless
  // of the "X0" option, so use graphs (explicit zero x error) instead, same
  // as the top pad.
  std::unique_ptr<TGraphErrors> ratioQQ(offsetGraph(rQQ.get(),  -0.15, "g_ratio_qq"));
  dlutility::SetMarkerAtt(ratioQQ.get(), kRed, 1.1, 21);
  dlutility::SetLineAtt(ratioQQ.get(), kRed, 2, 1);
  std::unique_ptr<TGraphErrors> ratioQGGG(offsetGraph(rQGGG.get(),  0.15, "g_ratio_qggg"));
  dlutility::SetMarkerAtt(ratioQGGG.get(), kAzure - 6, 1.1, 22);
  dlutility::SetLineAtt(ratioQGGG.get(), kAzure - 6, 2, 1);
  ratioQQ->Draw("P SAME");
  ratioQGGG->Draw("P SAME");

  TLine line(xjPlotMin, 1.0, xjBins[nbins], 1.0);
  line.SetLineStyle(2);
  line.SetLineColor(kGray + 2);
  line.Draw("same");

  const TString outputStem = Form(
    "%s/unfolding_plots/priorQA_xj_flavorCompare_%s_r%02d_iter%d",
    location.Data(), system.c_str(), cone_size, niter);
  c->Print(outputStem + ".pdf");
  std::cout << "Wrote " << outputStem << ".pdf" << std::endl;
}
