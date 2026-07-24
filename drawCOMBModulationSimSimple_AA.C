#include <algorithm>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"

#include "dlUtility.h"
#include "read_binning.h"

// Redraw the combined-simulation COMB plot from cached histograms.  This macro
// does no tree reading, event reweighting, fitting, or background calculation.
void drawCOMBModulationSimSimple_AA(
  const int cone_size = 3, const int centrality_bin = 1,
  const std::string config = "binning_AA.config")
{
  read_binning rb(config);
  std::unique_ptr<float[]> ptBins(new float[rb.get_nbins() + 1]);
  rb.get_pt_bins(ptBins.get());
  const TString stem = Form(
    "%s/dphi_plots/dphi_COMB_modulation_sim_combined_AA_cent_%d_r%02d",
    rb.get_code_location().c_str(), centrality_bin, cone_size);
  std::unique_ptr<TFile> input(TFile::Open(stem + ".root", "READ"));
  if (!input || input->IsZombie())
    {
      std::cerr << "Cannot open component cache " << stem << ".root" << std::endl;
      return;
    }

  auto get = [&](const char *name)
    {
      TH1D *histogram = dynamic_cast<TH1D*>(input->Get(name));
      if (!histogram) std::cerr << "Missing " << name << " in cache" << std::endl;
      return histogram;
    };
  TH1D *raw = get("h_dphi_pairs");
  TH1D *background = get("h_flow_background_nominal");
  TH1D *nominal = get("h_dphi_nominal_subtracted_sim");
  TH1D *down = get("h_dphi_down_subtracted_sim");
  TH1D *up = get("h_dphi_up_subtracted_sim");
  TH1D *signalRegion = get("h_signal_region");
  TH1D *normalizationRegion = get("h_flow_normalization_region_sim");
  if (!raw || !background || !nominal || !down || !up || !signalRegion ||
      !normalizationRegion) return;

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  raw->SetLineColor(kBlack);
  raw->SetLineWidth(2);
  background->SetLineColor(kGray + 2);
  background->SetLineWidth(2);
  background->SetLineStyle(7);
  nominal->SetLineColor(kRed + 1);
  nominal->SetLineWidth(2);
  down->SetLineColor(kBlue + 1);
  down->SetLineWidth(2);
  down->SetLineStyle(2);
  up->SetLineColor(kMagenta + 2);
  up->SetLineWidth(2);
  up->SetLineStyle(3);
  normalizationRegion->SetFillColorAlpha(kYellow, 0.28);
  normalizationRegion->SetLineColor(kOrange + 7);
  signalRegion->SetFillColorAlpha(kGreen + 1, 0.22);
  signalRegion->SetLineColor(kGreen + 3);
  double minimum = 0;
  for (const TH1D *histogram : {nominal, down, up})
    minimum = std::min(minimum, histogram->GetMinimum());
  raw->SetMinimum(std::min(-0.08*raw->GetMaximum(), 1.15*minimum));
  raw->SetMaximum(1.55*raw->GetMaximum());

  TCanvas canvas("c_dphi_sim_comb_simple", "c_dphi_sim_comb_simple", 650, 560);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.04);
  canvas.SetBottomMargin(0.12);
  raw->Draw("hist");
  normalizationRegion->Draw("hist same");
  signalRegion->Draw("hist same");
  raw->Draw("hist same");
  background->Draw("hist same");
  nominal->Draw("hist same");
  down->Draw("hist same");
  up->Draw("hist same");

  float centralityBins[5] = {0};
  rb.get_centrality_bins(centralityBins);
  dlutility::DrawSPHENIX(0.18, 0.91, 0.040);
  dlutility::drawText("Combined reweighted simulation", 0.18, 0.82,
                      0, kBlack, 0.034);
  dlutility::drawText(Form("p_{T,1} > %.1f GeV", rb.get_reco_leading_cut()),
                      0.18, 0.76, 0, kBlack, 0.034);
  dlutility::drawText(Form("p_{T,2} > %.1f GeV", rb.get_reco_subleading_cut()),
                      0.18, 0.71, 0, kBlack, 0.034);
  dlutility::drawText(Form("%.0f - %.0f %%", centralityBins[centrality_bin],
                           centralityBins[centrality_bin + 1]),
                      0.18, 0.66, 0, kBlack, 0.034);
  TLegend legend(0.49, 0.52, 0.94, 0.93);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.029);
  legend.AddEntry(raw, "Raw reweighted simulation", "l");
  legend.AddEntry(background, "Nominal fitted modulation", "l");
  legend.AddEntry(nominal, "Nominal modulation sub.", "l");
  legend.AddEntry(down, "COMBDown: 0.7(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(up, "COMBUp: 1.3(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(signalRegion, "Signal region", "f");
  legend.AddEntry(normalizationRegion, "Flow normalization region", "f");
  legend.Draw();
  canvas.SaveAs(stem + ".pdf");
  std::cout << "Redrew " << stem << ".pdf from cached components" << std::endl;
}
