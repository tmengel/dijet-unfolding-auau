#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"

#include "dlUtility.h"
#include "read_binning.h"

namespace
{
constexpr double kFlowFitLow = 0.0;
constexpr double kFlowFitHigh = 2.5;
constexpr double kNormalizationLow = 0.8;
constexpr double kNormalizationHigh = 2.5;

void addNormalizedBackground(TH1D *sum, TH1D *pairDphi,
                             const TF1 &fit, const double v22Scale,
                             const double v33Scale)
{
  TF1 normalized(fit);
  normalized.SetParameter(1, v22Scale*fit.GetParameter(1));
  normalized.SetParameter(2, v33Scale*fit.GetParameter(2));

  const int firstBin = pairDphi->FindBin(kNormalizationLow);
  const int lastBin = pairDphi->FindBin(kNormalizationHigh);
  const int binCount = std::max(1, lastBin - firstBin);
  const double dataLevel = pairDphi->Integral(firstBin, lastBin)/binCount;
  const double fitLevel = normalized.Integral(kNormalizationLow,
                                               kNormalizationHigh)
    /(kNormalizationHigh - kNormalizationLow);
  normalized.SetParameter(0, normalized.GetParameter(0) + dataLevel - fitLevel);
  normalized.SetRange(0, TMath::Pi());

  for (int bin = 1; bin <= sum->GetNbinsX(); ++bin)
    sum->AddBinContent(bin, normalized.Eval(sum->GetBinCenter(bin)));
}

std::unique_ptr<TH1D> subtractBackground(const TH1D *raw,
                                         const TH1D *background,
                                         const char *name)
{
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(raw->Clone(name)));
  result->SetDirectory(nullptr);
  result->Add(background, -1.0);
  return result;
}

std::unique_ptr<TH1D> makeRegion(const TH1D *raw, const char *name,
                                 const double low, const double high)
{
  auto region = std::unique_ptr<TH1D>(static_cast<TH1D*>(raw->Clone(name)));
  region->SetDirectory(nullptr);
  for (int bin = 0; bin <= region->GetNbinsX() + 1; ++bin)
    {
      const double center = region->GetBinCenter(bin);
      if (center < low || center >= high)
        {
          region->SetBinContent(bin, 0);
          region->SetBinError(bin, 0);
        }
    }
  return region;
}
}

void drawCOMBModulation_AA(const int cone_size = 3,
                           const int centrality_bin = 0,
                           const std::string config = "binning_AA.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(config);
  const int nbins = rb.get_nbins();
  std::unique_ptr<float[]> ptBins(new float[nbins + 1]);
  rb.get_pt_bins(ptBins.get());
  const int leadingBin = rb.get_reco_leading_bin();
  const int subleadingBin = rb.get_reco_subleading_bin();
  const double leadingCut = rb.get_reco_leading_cut();
  const double subleadingCut = rb.get_reco_subleading_cut();
  const double signalCut = rb.get_dphicut();

  read_binning downBinning("binning_COMBDown_AA.config");
  read_binning upBinning("binning_COMBUp_AA.config");
  const std::array<double, 3> v22Scales = {
    1.0, downBinning.get_flow_sys(), upBinning.get_flow_sys()};
  const std::array<double, 3> v33Scales = {
    1.0, downBinning.get_flow_v33_sys(), upBinning.get_flow_v33_sys()};

  const std::string system = "AA_cent_" + std::to_string(centrality_bin);
  const TString inputPath = Form(
    "%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.c_str(), cone_size);
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  if (!input || input->IsZombie())
    {
      std::cerr << "Cannot open " << inputPath << std::endl;
      return;
    }

  auto raw = std::make_unique<TH1D>("h_dphi_pairs", ";#Delta#phi;Counts",
                                    32, 0, TMath::Pi());
  auto etaSeparated = std::make_unique<TH1D>(
    "h_dphi_eta_separated", ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  std::array<std::unique_ptr<TH1D>, 3> backgrounds;
  for (int variation = 0; variation < 3; ++variation)
    backgrounds[variation] = std::make_unique<TH1D>(
      Form("h_flow_background_%d", variation), ";#Delta#phi;Counts",
      32, 0, TMath::Pi());
  raw->Sumw2();
  etaSeparated->Sumw2();

  int fittedBins = 0;
  for (int i = leadingBin; i < nbins; ++i)
    for (int j = subleadingBin; j <= i; ++j)
      {
        TH1D *pairDphi = static_cast<TH1D*>(input->Get(Form(
          "h_dphi_exclusive_%d_%d", i, j)));
        TH1D *fitSource = static_cast<TH1D*>(input->Get(Form(
          "h_dphi_eta_inclusive_%d_%d", i, j)));
        if (!pairDphi || !fitSource) continue;
        raw->Add(pairDphi);
        etaSeparated->Add(fitSource);

        std::unique_ptr<TH1D> fitHistogram(static_cast<TH1D*>(fitSource->Clone(
          Form("h_fit_work_%d_%d", i, j))));
        fitHistogram->SetDirectory(nullptr);
        TF1 fit(Form("flow_fit_%d_%d", i, j),
                "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))",
                kFlowFitLow, kFlowFitHigh);
        fit.SetParLimits(0, 0, 10000);
        fit.SetParLimits(1, 0, 0.5);
        fit.SetParLimits(2, 0, 0.5);
        const int firstFitBin = fitHistogram->FindBin(kFlowFitLow);
        const int lastFitBin = fitHistogram->FindBin(kFlowFitHigh);
        const int fitBinCount = std::max(1, lastFitBin - firstFitBin);
        const double fitCounts = fitHistogram->Integral(firstFitBin, lastFitBin);
        if (fitCounts < fitBinCount)
          {
            fit.SetParameter(0, fitCounts/fitBinCount);
            fit.FixParameter(1, 0);
            fit.FixParameter(2, 0);
          }
        else
          fit.SetParameters(fitCounts/fitBinCount, 0.01, 0.01);
        if (fitCounts > 0) fitHistogram->Fit(&fit, "0RlQ");
        if (fitCounts > 0 || pairDphi->Integral() > 0) ++fittedBins;
        for (int variation = 0; variation < 3; ++variation)
          addNormalizedBackground(backgrounds[variation].get(), pairDphi, fit,
                                  v22Scales[variation], v33Scales[variation]);
      }

  if (raw->Integral() <= 0)
    {
      std::cerr << "No selected #Delta#phi pairs for " << system << std::endl;
      return;
    }

  auto nominalSubtracted = subtractBackground(
    raw.get(), backgrounds[0].get(), "h_dphi_nominal_subtracted");
  auto downSubtracted = subtractBackground(
    raw.get(), backgrounds[1].get(), "h_dphi_COMBDown_subtracted");
  auto upSubtracted = subtractBackground(
    raw.get(), backgrounds[2].get(), "h_dphi_COMBUp_subtracted");
  auto normalizationRegion = makeRegion(
    raw.get(), "h_flow_normalization_region", kNormalizationLow,
    kNormalizationHigh);
  auto signalRegion = makeRegion(raw.get(), "h_signal_region", signalCut,
                                 TMath::Pi() + 1e-6);

  raw->SetLineColor(kBlack);
  raw->SetLineWidth(2);
  nominalSubtracted->SetLineColor(kRed + 1);
  nominalSubtracted->SetLineWidth(2);
  downSubtracted->SetLineColor(kBlue + 1);
  downSubtracted->SetLineWidth(2);
  downSubtracted->SetLineStyle(2);
  upSubtracted->SetLineColor(kMagenta + 2);
  upSubtracted->SetLineWidth(2);
  upSubtracted->SetLineStyle(3);
  backgrounds[0]->SetLineColor(kGray + 2);
  backgrounds[0]->SetLineWidth(2);
  backgrounds[0]->SetLineStyle(7);
  normalizationRegion->SetFillColorAlpha(kYellow, 0.28);
  normalizationRegion->SetLineColor(kOrange + 7);
  normalizationRegion->SetLineStyle(3);
  signalRegion->SetFillColorAlpha(kGreen + 1, 0.22);
  signalRegion->SetLineColor(kGreen + 3);
  signalRegion->SetLineStyle(3);

  double minimum = 0;
  for (const TH1D *histogram : {nominalSubtracted.get(), downSubtracted.get(),
                                upSubtracted.get()})
    minimum = std::min(minimum, histogram->GetMinimum());
  raw->SetMinimum(std::min(-0.08*raw->GetMaximum(), 1.15*minimum));
  raw->SetMaximum(1.55*raw->GetMaximum());
  raw->GetXaxis()->SetTitleOffset(1.05);
  raw->GetYaxis()->SetTitleOffset(1.25);

  gSystem->mkdir(Form("%s/unfolding_plots", rb.get_code_location().c_str()), true);
  TCanvas canvas("c_dphi_comb", "c_dphi_comb", 650, 560);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.04);
  canvas.SetBottomMargin(0.12);
  raw->Draw("hist");
  normalizationRegion->Draw("hist same");
  signalRegion->Draw("hist same");
  raw->Draw("hist same");
  backgrounds[0]->Draw("hist same");
  nominalSubtracted->Draw("hist same");
  downSubtracted->Draw("hist same");
  upSubtracted->Draw("hist same");

  float centralityBins[5] = {0};
  rb.get_centrality_bins(centralityBins);
  dlutility::DrawSPHENIX(0.18, 0.91, 0.040);
  dlutility::drawText(Form("#it{p}_{T,1} #geq %.1f GeV", leadingCut),
                      0.18, 0.80, 0, kBlack, 0.038);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %.1f GeV", subleadingCut),
                      0.18, 0.75, 0, kBlack, 0.038);
  dlutility::drawText(Form("%.0f - %.0f %%", centralityBins[centrality_bin],
                           centralityBins[centrality_bin + 1]),
                      0.18, 0.70, 0, kBlack, 0.038);

  TLegend legend(0.50, 0.55, 0.94, 0.93);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(raw.get(), "Raw #Delta#phi pairs", "l");
  legend.AddEntry(backgrounds[0].get(), "Nominal flow modulation", "l");
  legend.AddEntry(nominalSubtracted.get(), "Nominal modulation sub.", "l");
  legend.AddEntry(downSubtracted.get(), "COMBDown: 0.7(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(upSubtracted.get(), "COMBUp: 1.3(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(signalRegion.get(), "Signal region", "f");
  legend.AddEntry(normalizationRegion.get(), "Flow normalization region", "f");
  legend.Draw();

  const TString outputStem = Form(
    "%s/unfolding_plots/dphi_COMB_modulation_%s_r%02d",
    rb.get_code_location().c_str(), system.c_str(), cone_size);
  canvas.SaveAs(outputStem + ".pdf");
  canvas.SaveAs(outputStem + ".png");

  std::unique_ptr<TFile> output(TFile::Open(outputStem + ".root", "RECREATE"));
  if (output && !output->IsZombie())
    {
      raw->Write();
      etaSeparated->Write();
      backgrounds[0]->Write("h_flow_background_nominal");
      backgrounds[1]->Write("h_flow_background_COMBDown");
      backgrounds[2]->Write("h_flow_background_COMBUp");
      nominalSubtracted->Write();
      downSubtracted->Write();
      upSubtracted->Write();
      signalRegion->Write();
      normalizationRegion->Write();
      output->Write();
    }
  std::cout << "Wrote aggregate COMB modulation plot for " << system
            << " using " << fittedBins << " populated pT-bin fits to "
            << outputStem << ".{pdf,png,root}" << std::endl;
}
