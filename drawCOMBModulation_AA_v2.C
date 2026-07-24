#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TSystem.h"
#include "TTree.h"

#include "dlUtility.h"
#include "read_binning.h"

namespace
{
constexpr double kFlowFitLow = 0.0;
constexpr double kFlowFitHigh = 2.5;
constexpr double kNormalizationLow = 0.8;
constexpr double kNormalizationHigh = 2.5;
constexpr bool kSavePng = false;
constexpr bool kSaveRoot = false;

std::unique_ptr<TF1> addNormalizedBackground(TH1D *sum, TH1D *pairDphi,
                             const TF1 &fit, const double v22Scale,
                             const double v33Scale)
{
  // TF1 normalized(fit);
  std::unique_ptr<TF1> normalized(static_cast<TF1*>(fit.Clone("normalized")));
  normalized->SetParameter(1, v22Scale*fit.GetParameter(1));
  normalized->SetParameter(2, v33Scale*fit.GetParameter(2));

  const int firstBin = pairDphi->FindBin(kNormalizationLow);
  const int lastBin = pairDphi->FindBin(
    std::nextafter(kNormalizationHigh, kNormalizationLow));
  const double dataCounts = pairDphi->Integral(firstBin, lastBin);
  normalized->SetParameter(0, 1.0);
  double unitShapeCounts = 0;
  for (int bin = firstBin; bin <= lastBin; ++bin)
    unitShapeCounts += normalized->Eval(pairDphi->GetBinCenter(bin));
  normalized->SetParameter(0, unitShapeCounts > 0
    ? dataCounts/unitShapeCounts : 0.0);
  normalized->SetRange(0, TMath::Pi());

  for (int bin = 1; bin <= sum->GetNbinsX(); ++bin)
    sum->AddBinContent(bin, normalized->Eval(sum->GetBinCenter(bin)));

  return normalized;
}

void addFitShape(TH1D *sum, const TF1 &fit)
{
  for (int bin = 1; bin <= sum->GetNbinsX(); ++bin)
    sum->AddBinContent(bin, fit.Eval(sum->GetBinCenter(bin)));
}

std::unique_ptr<TH1D> getSimulationDphi(
  const TString &inputPath, const TString &name, read_binning &rb,
  int coneSize, int centralityBin, int sampleIndex,
  double centralityLow, double centralityHigh,
  double leadingCut, double subleadingCut, bool useTruth = false
);
std::unique_ptr<TH1D> getSimulationDphi_EtaSeperated(
  const TString &inputPath, const TString &name, read_binning &rb,
  int coneSize, int centralityBin, int sampleIndex,
  double centralityLow, double centralityHigh,
  double leadingCut, double subleadingCut, bool useTruth = false,
  bool requireEtaGap = true, bool inclusivePairs = true
);
std::unique_ptr<TH1D> getCombinedSimulationDphi(
  const TString &name, read_binning &rb, int coneSize, int centralityBin,
  double centralityLow, double centralityHigh,
  double leadingCut, double subleadingCut, bool useTruth = false
);
std::unique_ptr<TH1D> getCombinedSimulationDphi_EtaSeperated(
  const TString &name, read_binning &rb, const int coneSize,
  const int centralityBin, const double centralityLow,
  const double centralityHigh, const double leadingCut,
  const double subleadingCut, const bool useTruth = false,
  const bool requireEtaGap = true, const bool inclusivePairs = true
);
std::unique_ptr<TH1D> subtractBackground(  const TH1D *raw, const TH1D *background, const char *name);
std::unique_ptr<TH1D> makeRegion(  const TH1D *raw, const char *name, double low, double high);

void drawPtBinFit(const TH1D *pairs, const TH1D *fitInput, const TF1 &sourceFit,
                  const int pt1Bin, const int pt2Bin, const float *ptBins,
                  const TString &system, const int coneSize,
                  const TString &outputDirectory)
{
  if (pairs->Integral() <= 0 && fitInput->Integral() <= 0) return;
  auto raw = std::unique_ptr<TH1D>(static_cast<TH1D*>(pairs->Clone(
    Form("h_dphi_pairs_pt1_%d_pt2_%d", pt1Bin, pt2Bin))));
  auto fitted = std::unique_ptr<TH1D>(static_cast<TH1D*>(fitInput->Clone(
    Form("h_dphi_fit_input_pt1_%d_pt2_%d", pt1Bin, pt2Bin))));
  raw->SetDirectory(nullptr);
  fitted->SetDirectory(nullptr);
  raw->SetLineColor(kBlack);
  raw->SetLineWidth(2);
  fitted->SetMarkerColor(kBlue + 1);
  fitted->SetLineColor(kBlue + 1);
  fitted->SetMarkerStyle(20);
  fitted->SetMarkerSize(0.75);
  TF1 fit(sourceFit);
  fit.SetName(Form("f_flow_fit_pt1_%d_pt2_%d", pt1Bin, pt2Bin));
  fit.SetLineColor(kRed + 1);
  fit.SetLineWidth(2);
  raw->SetMaximum(1.35*std::max(raw->GetMaximum(), fitted->GetMaximum()));

  TCanvas canvas(Form("c_dphi_pt1_%d_pt2_%d", pt1Bin, pt2Bin),
                 "c_dphi_pt_bins", 650, 560);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.04);
  canvas.SetBottomMargin(0.12);
  raw->Draw("hist");
  fitted->Draw("E1 same");
  fit.Draw("same");
  dlutility::DrawSPHENIX(0.18, 0.91, 0.040);
  dlutility::drawText(Form("%.1f < p_{T,1} < %.1f GeV", ptBins[pt1Bin],
                           ptBins[pt1Bin + 1]), 0.18, 0.79, 0, kBlack, 0.034);
  dlutility::drawText(Form("%.1f < p_{T,2} < %.1f GeV", ptBins[pt2Bin],
                           ptBins[pt2Bin + 1]), 0.18, 0.74, 0, kBlack, 0.034);
  TLegend legend(0.48, 0.65, 0.94, 0.92);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(raw.get(), "Exclusive pairs", "l");
  legend.AddEntry(fitted.get(), "#Delta#eta-separated fit input", "ep");
  legend.AddEntry(&fit, Form("fit: v_{2,2}=%.4g, v_{3,3}=%.4g",
                             fit.GetParameter(1), fit.GetParameter(2)), "l");
  legend.Draw();
  const TString stem = Form("%s/dphi_fit_%s_r%02d_pt1_%02d_pt2_%02d",
                            outputDirectory.Data(), system.Data(), coneSize,
                            pt1Bin, pt2Bin);
  canvas.SaveAs(stem + ".pdf");
  // canvas.SaveAs(stem + ".png");
  std::unique_ptr<TFile> output(kSaveRoot
    ? TFile::Open(stem + ".root", "RECREATE") : nullptr);
  if (output && !output->IsZombie())
    {
      raw->Write();
      fitted->Write();
      fit.Write();
      output->Write();
    }
}
struct Jet
{
  int id = -1;
  float pt = 0.0F;
  float eta = 0.0F;
  float phi = 0.0F;
};
float deltaPhi(const float first, const float second)
{
  float difference = first - second;
  while (difference > TMath::Pi()) difference -= 2.0F*TMath::Pi();
  while (difference < -TMath::Pi()) difference += 2.0F*TMath::Pi();
  return std::fabs(difference);
}

void drawSimulationDphi(const TString &inputPath, const TString &sample,
                        const int coneSize, const int centralityBin,
                        const double centralityLow, const double centralityHigh,
                        const double leadingCut, const double subleadingCut,
                        const TString &outputDirectory, read_binning &rb,
                        const int sampleIndex)
{
  const int nbins = rb.get_nbins();
  std::unique_ptr<float[]> ptBins(new float[nbins + 1]);
  rb.get_pt_bins(ptBins.get());
  const double signalCut = rb.get_dphicut();

  read_binning downBinning("binning_COMBDown_AA.config");
  read_binning upBinning("binning_COMBUp_AA.config");
  const std::array<double, 3> v22Scales = { 1.0, downBinning.get_flow_sys(), upBinning.get_flow_sys()};
  const std::array<double, 3> v33Scales = { 1.0, downBinning.get_flow_v33_sys(), upBinning.get_flow_v33_sys()};

  const std::string system = "AA_cent_" + std::to_string(centralityBin);
  const TString plotDirectory = Form("%s/dphi_plots", rb.get_code_location().c_str());
  gSystem->mkdir(plotDirectory, true);
  std::cout << "Drawing COMB modulation for " << system << std::endl;

  auto raw = getCombinedSimulationDphi_EtaSeperated("h_dphi_pairs_sim", 
    rb, coneSize,
    centralityBin, 
    centralityLow, centralityHigh, 
    leadingCut, subleadingCut, false, false, false
  );
  if (!raw)
  {
    std::cerr << "Failed to load raw histogram for " << system << std::endl;
    return;
  }
  auto etaSeparated = getCombinedSimulationDphi_EtaSeperated("h_dphi_eta_separated_sim", 
    rb, coneSize,
    centralityBin, 
    centralityLow, centralityHigh, leadingCut, subleadingCut,
    false, true, true
  );
  if (!etaSeparated)
  {
    std::cerr << "Failed to load  eta histogram for " << system << std::endl;
    return;
  }
  auto fittedEtaSeparated = std::make_unique<TH1D>( "h_dphi_eta_separated_fit_sim", ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  std::array<std::unique_ptr<TH1D>, 3> backgrounds;
  
  auto truth_raw = getCombinedSimulationDphi_EtaSeperated("h_dphi_pairs_truth", rb, coneSize,
    centralityBin, 
    centralityLow, centralityHigh, leadingCut, subleadingCut,
    true, false, false
  );
  if (!truth_raw)
  {
    std::cerr << "Failed to load histogram for " << system << std::endl;
    return;
  }
  auto truth_etaSeparated = getCombinedSimulationDphi_EtaSeperated("h_dphi_eta_separated_truth", rb, coneSize,
    centralityBin, 
    centralityLow, centralityHigh, leadingCut, subleadingCut,
    true, true, true
  );
  if (!truth_etaSeparated)
  {
    std::cerr << "Failed to load histogram for " << system << std::endl;
    return;
  }

  if (!raw || !etaSeparated || !truth_raw || !truth_etaSeparated)
  {
    
    std::cerr << "Failed to load simulation histograms for " << system
              << std::endl;
    return;
  }
  if (raw->Integral() <= 0 || etaSeparated->Integral() <= 0)
  {
    std::cerr << "No selected #Delta#phi pairs for " << system << std::endl;
    // return;
  }
  
  for (int variation = 0; variation < 3; ++variation)
  {
    backgrounds[variation] = std::make_unique<TH1D>( Form("h_flow_background_%d_sim", variation), ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  }
  TF1 fit("flow_fit_sim", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))", kFlowFitLow, kFlowFitHigh );
  fit.SetParLimits(0, 0, 10000);
  fit.SetParLimits(1, 0, 0.5);
  fit.SetParLimits(2, 0, 0.5);
  const int firstFitBin = etaSeparated->FindBin(kFlowFitLow);
  const int lastFitBin = etaSeparated->FindBin(kFlowFitHigh);
  const int fitBinCount = std::max(1, lastFitBin - firstFitBin + 1);
  const double fitCounts = etaSeparated->Integral(firstFitBin, lastFitBin);
  if (fitCounts < fitBinCount)
  {
    fit.SetParameter(0, fitCounts/fitBinCount);
    fit.FixParameter(1, 0);
    fit.FixParameter(2, 0);
  }
  else fit.SetParameters(fitCounts/fitBinCount, 0.01, 0.01);
  
  if (fitCounts > 0) etaSeparated->Fit(&fit, "0RlQ");
  addFitShape(fittedEtaSeparated.get(), fit);
  // TF1 * fitNorms[3] = {nullptr, nullptr, nullptr};
  // std::vector<std::unique_ptr<TF1>> fitNorms(3);
  std::array<std::unique_ptr<TF1>, 3> fitNorms;
  for (int variation = 0; variation < 3; ++variation)
  {
    // addNormalizedBackground(
    fitNorms[variation] = addNormalizedBackground(
      backgrounds[variation].get(), 
      raw.get() , fit, v22Scales[variation], v33Scales[variation]
    );
    fitNorms[variation]->SetName(Form("f_flow_fit_%d_sim", variation));
  }

  auto nominalSubtracted = subtractBackground( raw.get(), backgrounds[0].get(), "h_dphi_nominal_subtracted_sim");
  auto downSubtracted = subtractBackground( raw.get(), backgrounds[1].get(), "h_dphi_down_subtracted_sim");
  auto upSubtracted = subtractBackground( raw.get(), backgrounds[2].get(), "h_dphi_up_subtracted_sim");
  auto normalizationRegion = makeRegion(
    raw.get(), "h_flow_normalization_region_sim", kNormalizationLow,
    kNormalizationHigh
  );
  auto signalRegion = makeRegion(raw.get(), "h_signal_region", signalCut, TMath::Pi() + 1e-6);

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
  etaSeparated->SetLineColor(kCyan + 2);
  etaSeparated->SetLineWidth(2);
  etaSeparated->SetLineStyle(9);
  fittedEtaSeparated->SetLineColor(kOrange + 7);
  fittedEtaSeparated->SetLineWidth(3);
  normalizationRegion->SetFillColorAlpha(kYellow, 0.28);
  normalizationRegion->SetLineColor(kOrange + 7);
  normalizationRegion->SetLineStyle(3);
  signalRegion->SetFillColorAlpha(kGreen + 1, 0.22);
  signalRegion->SetLineColor(kGreen + 3);
  signalRegion->SetLineStyle(3);

  double minimum = 0;
  for (
    const TH1D *histogram : {
      nominalSubtracted.get(), downSubtracted.get(), upSubtracted.get()
    }
  )
  {
    minimum = std::min(minimum, histogram->GetMinimum());
  }
  raw->SetMinimum(std::min(-0.08*raw->GetMaximum(), 1.15*minimum));
  raw->SetMaximum(1.55*raw->GetMaximum());
  raw->GetXaxis()->SetTitleOffset(1.05);
  raw->GetYaxis()->SetTitleOffset(1.25);

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
  etaSeparated->Draw("hist same");
  fittedEtaSeparated->Draw("hist same");
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
  dlutility::drawText(Form("%.0f - %.0f %%", centralityBins[centralityBin],
                           centralityBins[centralityBin + 1]),
                      0.18, 0.70, 0, kBlack, 0.038);

  TLegend legend(0.47, 0.49, 0.94, 0.93);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(raw.get(), "Raw #Delta#phi pairs", "l");
  legend.AddEntry(backgrounds[0].get(), "Nominal flow modulation", "l");
  legend.AddEntry(etaSeparated.get(), "#Delta#eta-separated fit input", "l");
  legend.AddEntry(fittedEtaSeparated.get(), "Fitted #Delta#eta histogram", "l");
  legend.AddEntry(nominalSubtracted.get(), "Nominal modulation sub.", "l");
  legend.AddEntry(downSubtracted.get(), "COMBDown: 0.7(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(upSubtracted.get(), "COMBUp: 1.3(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(signalRegion.get(), "Signal region", "f");
  legend.AddEntry(normalizationRegion.get(), "Flow normalization region", "f");
  legend.Draw();

  const TString outputStem = Form(
    "%s/dphi_COMB_modulation_sim_combined_AA_cent_%d_r%02d",
    plotDirectory.Data(), centralityBin, coneSize);
  canvas.SaveAs(outputStem + ".pdf");
  if (kSavePng) canvas.SaveAs(outputStem + ".png");

  // Cache the aggregate data components for fast redraws, matching the
  // combined-simulation COMB cache.
  std::unique_ptr<TFile> output(TFile::Open(outputStem + ".root", "RECREATE"));
  if (output && !output->IsZombie())
    {
      raw->Write("h_dphi_pairs");
      etaSeparated->Write("h_dphi_eta_separated");
      fittedEtaSeparated->Write("h_dphi_eta_separated_fit");
      backgrounds[0]->Write("h_flow_background_nominal");
      backgrounds[1]->Write("h_flow_background_COMBDown");
      backgrounds[2]->Write("h_flow_background_COMBUp");
      nominalSubtracted->Write();
      downSubtracted->Write();
      upSubtracted->Write();
      signalRegion->Write();
      normalizationRegion->Write();
      truth_raw->Write("h_dphi_pairs_truth");
      truth_etaSeparated->Write("h_dphi_eta_separated_truth");
      fit.Write("f_flow_fit");
      for (int variation = 0; variation < 3; ++variation)
        fitNorms[variation]->Write();
      output->Write();
    }

  std::cout << "Finished drawing COMB modulation for " << system
            << std::endl;
}

std::unique_ptr<TH1D> getSimulationDphi(
  const TString &inputPath, const TString &name, read_binning &rb,
  const int coneSize, const int centralityBin, const int sampleIndex,
  const double centralityLow, const double centralityHigh,
  const double leadingCut, const double subleadingCut, const bool useTruth)
{
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  TTree *tree = input && !input->IsZombie()
    ? dynamic_cast<TTree*>(input->Get("tn_match")) : nullptr;
  if (!tree) return nullptr;
  auto histogram = std::make_unique<TH1D>(name, ";#Delta#phi;;#frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}",
                                          32, 0, TMath::Pi());
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();

  auto loadWeight = [&](const TString &path, const char *histogramName)
  {
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    TH1D *source = file && !file->IsZombie()
      ? dynamic_cast<TH1D*>(file->Get(histogramName)) : nullptr;
    if (!source) return std::unique_ptr<TH1D>();
    auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(source->Clone()));
    result->SetDirectory(nullptr);
    return result;
  };
  const TString system = Form("AA_cent_%d", centralityBin);
  auto centralityWeight = loadWeight(Form(
    "%s/centrality/centrality_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_centrality_reweight");
  auto vertexWeight = loadWeight(Form(
    "%s/vertex/vertex_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_mbd_reweight");
  auto sumETWeight = loadWeight(Form(
    "%s/sumeT/sumeT_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_sumeT_reweight");
  if (!centralityWeight || !vertexWeight || !sumETWeight)
  {
    std::cerr << "Missing nominal simulation reweighting histogram for "
              << system << std::endl;
    return nullptr;
  }
  auto weightAt = [](const TH1D *weights, const double value)
  {
    const int bin = weights->FindFixBin(value);
    return bin >= 1 && bin <= weights->GetNbinsX()
      ? weights->GetBinContent(bin) : 1.0;
  };

  float maxPtTruth = 0;
  float pt1Reco = 0;
  float pt2Reco = 0;
  float dphiReco = 0;
  float centrality = 0;
  float vertex = 0;
  float sumET = 0;
  tree->SetBranchAddress("maxpttruth", &maxPtTruth);
  tree->SetBranchAddress(useTruth ? "pt1_truth" : "pt1_reco", &pt1Reco);
  tree->SetBranchAddress(useTruth ? "pt2_truth" : "pt2_reco", &pt2Reco);
  tree->SetBranchAddress(useTruth ? "dphi_truth" : "dphi_reco", &dphiReco);
  tree->SetBranchAddress("centrality", &centrality);
  tree->SetBranchAddress("mbd_vertex", &vertex);
  tree->SetBranchAddress("sumeT", &sumET);
  const double sampleLow = rb.get_sample_boundary(sampleIndex);
  const double sampleHigh = rb.get_sample_boundary(sampleIndex + 1);
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree->GetEntry(entry);
      if (maxPtTruth < sampleLow || maxPtTruth >= sampleHigh) continue;
      if (centrality < centralityLow || centrality >= centralityHigh) continue;
      if (pt1Reco < leadingCut || pt2Reco < subleadingCut) continue;
      const double eventWeight = weightAt(centralityWeight.get(), centrality)
        * weightAt(vertexWeight.get(), vertex) * weightAt(sumETWeight.get(), sumET);
      histogram->Fill(dphiReco, eventWeight);
    }
  return histogram;
}
std::unique_ptr<TH1D> getSimulationDphi_EtaSeperated(
  const TString &inputPath, const TString &name, read_binning &rb,
  const int coneSize, 
  const int centralityBin, 
  const int sampleIndex,
  const double centralityLow, 
  const double centralityHigh,
  const double leadingCut, 
  const double subleadingCut,
  const bool useTruth,
  const bool requireEtaGap,
  const bool inclusivePairs
)
{
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  TTree *tree = input && !input->IsZombie() ? dynamic_cast<TTree*>(input->Get("T")) : nullptr;
  
  if (!tree) return nullptr;

  auto histogram = std::make_unique<TH1D>(
    name, ";#Delta#phi;;#frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}",
    32, 0, TMath::Pi()
  );
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();

  auto loadWeight = [&](const TString &path, const char *histogramName)
  {
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    TH1D *source = file && !file->IsZombie()
      ? dynamic_cast<TH1D*>(file->Get(histogramName)) : nullptr;
    if (!source) return std::unique_ptr<TH1D>();
    auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(source->Clone()));
    result->SetDirectory(nullptr);
    return result;
  };
  const TString system = Form("AA_cent_%d", centralityBin);
  auto centralityWeight = loadWeight(Form(
    "%s/centrality/centrality_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_centrality_reweight"
  );
  auto vertexWeight = loadWeight(Form(
    "%s/vertex/vertex_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_mbd_reweight"
  );
  auto sumETWeight = loadWeight(Form(
    "%s/sumeT/sumeT_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_sumeT_reweight"
  );
  if (!centralityWeight || !vertexWeight || !sumETWeight)
  {
    std::cerr << "Missing nominal simulation reweighting histogram for "
              << system << std::endl;
    return nullptr;
  }
  auto weightAt = [](const TH1D *weights, const double value)
  {
    const int bin = weights->FindFixBin(value);
    return bin >= 1 && bin <= weights->GetNbinsX()
      ? weights->GetBinContent(bin) : 1.0;
  };

  tree->SetBranchStatus("*", 0);
  std::vector<std::string> required = {
    "jetid", "zvrtx", "cent", "sumeT",
    "truth_jet_pT", "truth_jet_eta", "truth_jet_phi",
    "jet_pT", "jet_E", "jet_unsub_E", "jet_eta", "jet_phi"
  };
  for (const auto& branch : required) tree->SetBranchStatus(branch.c_str(), 1);
  tree->SetCacheSize(256*1024*1024);
  tree->AddBranchToCache("*", true);

  int jetId = 0;
  int centrality = 0;
  float vertex = 0.0F;
  float sumET = 0.0F;
  std::vector<float>* truthPt = nullptr;
  std::vector<float>* truthEta = nullptr;
  std::vector<float>* truthPhi = nullptr;
  std::vector<float>* recoPt = nullptr;
  std::vector<float>* recoE = nullptr;
  std::vector<float>* recoUnsubE = nullptr;
  std::vector<float>* recoEta = nullptr;
  std::vector<float>* recoPhi = nullptr;
  
  tree->SetBranchAddress("jetid", &jetId);
  tree->SetBranchAddress("zvrtx", &vertex);
  tree->SetBranchAddress("cent", &centrality);
  tree->SetBranchAddress("sumeT", &sumET);

  tree->SetBranchAddress("truth_jet_pT", &truthPt);
  tree->SetBranchAddress("truth_jet_eta", &truthEta);
  tree->SetBranchAddress("truth_jet_phi", &truthPhi);
  
  tree->SetBranchAddress("jet_pT", &recoPt);
  tree->SetBranchAddress("jet_E", &recoE);
  tree->SetBranchAddress("jet_unsub_E", &recoUnsubE);
  tree->SetBranchAddress("jet_eta", &recoEta);
  tree->SetBranchAddress("jet_phi", &recoPhi);

  const double sampleLow = rb.get_sample_boundary(sampleIndex);
  const double sampleHigh = rb.get_sample_boundary(sampleIndex + 1);

  std::vector<Jet> truthJets;
  std::vector<Jet> recoJets;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
      tree->GetEntry(entry);
      if (jetId != 10*(sampleIndex + 1)) continue;
      if (std::fabs(vertex) > 60.0F) continue;
      if (centrality < centralityLow || centrality >= centralityHigh) continue;

      if (!truthPt || truthPt->empty() || !truthEta || !truthPhi ||
          !recoPt || !recoE || !recoUnsubE || !recoEta || !recoPhi) continue;
      float maxPtTruth = *std::max_element(truthPt->begin(), truthPt->end());
      if (maxPtTruth < sampleLow || maxPtTruth >= sampleHigh) continue;
      
      truthJets.clear();
      recoJets.clear();
      for (std::size_t jet = 0; jet < truthPt->size(); ++jet)
      {
        if (truthPt->at(jet) >= 3.0F && std::fabs(truthEta->at(jet)) <= 0.8F )
        {
          truthJets.push_back({static_cast<int>(jet), truthPt->at(jet), truthEta->at(jet), truthPhi->at(jet)});
        }
      }
      for (std::size_t jet = 0; jet < recoPt->size(); ++jet)
      {
        if( recoPt->at(jet) >= 3.0F && recoE->at(jet) >= 0.0F &&
            recoUnsubE->at(jet) >= 0.0F && std::fabs(recoEta->at(jet)) <= 0.8F
        )
        {
          recoJets.push_back({static_cast<int>(jet), recoPt->at(jet), recoEta->at(jet), recoPhi->at(jet)});
        }
      }

      const auto descendingPt = [](const Jet& first, const Jet& second){ return first.pt > second.pt; };
      std::sort(truthJets.begin(), truthJets.end(), descendingPt);
      std::sort(recoJets.begin(), recoJets.end(), descendingPt);
      if (truthJets.size() < 2 || recoJets.size() < 2) continue;

      float pt1Reco = recoJets[0].pt;
      float pt2Reco = recoJets[1].pt;
      if (pt1Reco < leadingCut || pt2Reco < subleadingCut) continue;
    
      const double eventWeight =  (
        weightAt(centralityWeight.get(), centrality)
        * weightAt(vertexWeight.get(), vertex) 
        * weightAt(sumETWeight.get(), sumET)
      );
      if ( useTruth )
      {
        const std::size_t pairCount = inclusivePairs ? truthJets.size() : 2;
        for (std::size_t i = 1 ; i < pairCount; ++i)
        {
          if (truthJets[i].pt < subleadingCut) continue;
          float dphiTruth = deltaPhi(truthJets[0].phi, truthJets[i].phi);
          
          float dEtaTruth = truthJets[0].eta - truthJets[i].eta;
          if (requireEtaGap && std::fabs(dEtaTruth) < 0.8F) continue;

          histogram->Fill(dphiTruth, eventWeight);
        }
      }
      else
      {
        const std::size_t pairCount = inclusivePairs ? recoJets.size() : 2;
        for (std::size_t i = 1 ; i < pairCount; ++i)
        {
          if (recoJets[i].pt < subleadingCut) continue;
          float dphiReco = deltaPhi(recoJets[0].phi, recoJets[i].phi);
          float dEtaReco = recoJets[0].eta - recoJets[i].eta;
          if (requireEtaGap && std::fabs(dEtaReco) < 0.8F) continue;

          histogram->Fill(dphiReco, eventWeight);
        }
      }


  }
  input->Close();
  return histogram;
}
std::unique_ptr<TH1D> getSimulationDphi_EtaSeperated_legacy(
  const TString &inputPath, const TString &name, read_binning &rb,
  const int coneSize, 
  const int centralityBin, 
  const int sampleIndex,
  const double centralityLow, 
  const double centralityHigh,
  const double leadingCut, 
  const double subleadingCut,
  const bool useTruth
)
{
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  TTree *tree = input && !input->IsZombie() ? dynamic_cast<TTree*>(input->Get("T")) : nullptr;
  
  if (!tree) return nullptr;

  auto histogram = std::make_unique<TH1D>(
    name, ";#Delta#phi;;#frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}",
    32, 0, TMath::Pi()
  );
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();

  auto loadWeight = [&](const TString &path, const char *histogramName)
  {
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    TH1D *source = file && !file->IsZombie()
      ? dynamic_cast<TH1D*>(file->Get(histogramName)) : nullptr;
    if (!source) return std::unique_ptr<TH1D>();
    auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(source->Clone()));
    result->SetDirectory(nullptr);
    return result;
  };
  const TString system = Form("AA_cent_%d", centralityBin);
  auto centralityWeight = loadWeight(Form(
    "%s/centrality/centrality_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_centrality_reweight"
  );
  auto vertexWeight = loadWeight(Form(
    "%s/vertex/vertex_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_mbd_reweight"
  );
  auto sumETWeight = loadWeight(Form(
    "%s/sumeT/sumeT_reweight_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.Data(), coneSize),
    "h_sumeT_reweight"
  );
  if (!centralityWeight || !vertexWeight || !sumETWeight)
  {
    std::cerr << "Missing nominal simulation reweighting histogram for "
              << system << std::endl;
    return nullptr;
  }
  auto weightAt = [](const TH1D *weights, const double value)
  {
    const int bin = weights->FindFixBin(value);
    return bin >= 1 && bin <= weights->GetNbinsX()
      ? weights->GetBinContent(bin) : 1.0;
  };

  tree->SetBranchStatus("*", 0);
  std::vector<std::string> required = {
    "jetid", "zvrtx", "cent", "sumeT",
    "truth_jet_pT", "truth_jet_eta", "truth_jet_phi",
    "jet_pT", "jet_E", "jet_unsub_E", "jet_eta", "jet_phi"
  };
  for (const auto& branch : required) tree->SetBranchStatus(branch.c_str(), 1);
  tree->SetCacheSize(256*1024*1024);
  tree->AddBranchToCache("*", true);

  int jetId = 0;
  int centrality = 0;
  float vertex = 0.0F;
  float sumET = 0.0F;
  std::vector<float>* truthPt = nullptr;
  std::vector<float>* truthEta = nullptr;
  std::vector<float>* truthPhi = nullptr;
  std::vector<float>* recoPt = nullptr;
  std::vector<float>* recoE = nullptr;
  std::vector<float>* recoUnsubE = nullptr;
  std::vector<float>* recoEta = nullptr;
  std::vector<float>* recoPhi = nullptr;
  
  tree->SetBranchAddress("jetid", &jetId);
  tree->SetBranchAddress("zvrtx", &vertex);
  tree->SetBranchAddress("cent", &centrality);
  tree->SetBranchAddress("sumeT", &sumET);

  tree->SetBranchAddress("truth_jet_pT", &truthPt);
  tree->SetBranchAddress("truth_jet_eta", &truthEta);
  tree->SetBranchAddress("truth_jet_phi", &truthPhi);
  
  tree->SetBranchAddress("jet_pT", &recoPt);
  tree->SetBranchAddress("jet_E", &recoE);
  tree->SetBranchAddress("jet_unsub_E", &recoUnsubE);
  tree->SetBranchAddress("jet_eta", &recoEta);
  tree->SetBranchAddress("jet_phi", &recoPhi);

  const double sampleLow = rb.get_sample_boundary(sampleIndex);
  const double sampleHigh = rb.get_sample_boundary(sampleIndex + 1);

  std::vector<Jet> truthJets;
  std::vector<Jet> recoJets;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
      tree->GetEntry(entry);
      if ( jetId != sampleIndex ) continue;
      if (std::fabs(vertex) > 60.0F) continue;
      if (centrality < centralityLow || centrality >= centralityHigh) continue;

      if ( !truthPt || !truthEta || !truthPhi || !recoPt || !recoE  || !recoUnsubE || !recoEta || !recoPhi ) continue;
      float maxPtTruth = *std::max_element(truthPt->begin(), truthPt->end());
      if (maxPtTruth < sampleLow || maxPtTruth >= sampleHigh) continue;
      
      truthJets.clear();
      recoJets.clear();
      for (std::size_t jet = 0; jet < truthPt->size(); ++jet)
      {
        if (truthPt->at(jet) >= 3.0F && std::fabs(truthEta->at(jet)) <= 0.8F )
        {
          truthJets.push_back({static_cast<int>(jet), truthPt->at(jet), truthEta->at(jet), truthPhi->at(jet)});
        }
      }
      for (std::size_t jet = 0; jet < recoPt->size(); ++jet)
      {
        if( recoPt->at(jet) >= 3.0F && recoE->at(jet) >= 0.0F &&
            recoUnsubE->at(jet) >= 0.0F && std::fabs(recoEta->at(jet)) <= 0.8F
        )
        {
          recoJets.push_back({static_cast<int>(jet), recoPt->at(jet), recoEta->at(jet), recoPhi->at(jet)});
        }
      }

      const auto descendingPt = [](const Jet& first, const Jet& second){ return first.pt > second.pt; };
      std::sort(truthJets.begin(), truthJets.end(), descendingPt);
      std::sort(recoJets.begin(), recoJets.end(), descendingPt);
      if (truthJets.size() < 2 || recoJets.size() < 2) continue;

      float pt1Reco = recoJets[0].pt;
      float pt2Reco = recoJets[1].pt;
      if (pt1Reco < leadingCut || pt2Reco < subleadingCut) continue;
    
      const double eventWeight =  (
        weightAt(centralityWeight.get(), centrality)
        * weightAt(vertexWeight.get(), vertex) 
        * weightAt(sumETWeight.get(), sumET)
      );
      if ( useTruth )
      {
        for (std::size_t i = 1 ; i < truthJets.size(); ++i)
        {
          float dphiTruth = deltaPhi(truthJets[0].phi, truthJets[i].phi);
          
          float dEtaTruth = truthJets[0].eta - truthJets[i].eta;
          if ( std::fabs(dEtaTruth) < 0.8F ) continue; // skip jets that are too close in eta

          histogram->Fill(dphiTruth, eventWeight);
        }
      }
      else
      {
        for (std::size_t i = 1 ; i < recoJets.size(); ++i)
        {
          float dphiReco = deltaPhi(recoJets[0].phi, recoJets[i].phi);
          float dEtaReco = recoJets[0].eta - recoJets[i].eta;
          if ( std::fabs(dEtaReco) < 0.8F ) continue;

          histogram->Fill(dphiReco, eventWeight);
        }
      }


  }
  input->Close();
  return histogram;
}
std::unique_ptr<TH1D> getCombinedSimulationDphi(
  const TString &name, read_binning &rb, const int coneSize,
  const int centralityBin, const double centralityLow,
  const double centralityHigh, const double leadingCut,
  const double subleadingCut, const bool useTruth)
{
  const std::array<int, 3> samples = {10, 20, 30};
  const std::array<double, 3> crossSections = {
    2.889e-6, 5.4067742e-8, 2.505e-9};
  std::array<double, 3> eventCounts = {0, 0, 0};
  std::array<TString, 3> paths;
  for (std::size_t index = 0; index < samples.size(); ++index)
    {
      paths[index] = Form(
        "%s/rootfiles/TREE_MATCH_r%02d_v15_%d_new_ProdA_2024-00000030_sumeT.root",
        rb.get_code_location().c_str(), coneSize, samples[index]);
      std::unique_ptr<TFile> input(TFile::Open(paths[index], "READ"));
      TTree *stats = input && !input->IsZombie()
        ? dynamic_cast<TTree*>(input->Get("tn_stats")) : nullptr;
      float events = 0;
      if (!stats)
        {
          std::cerr << "Missing tn_stats in " << paths[index] << std::endl;
          return nullptr;
        }
      stats->SetBranchAddress("nevents", &events);
      stats->GetEntry(0);
      eventCounts[index] = events;
    }
  auto combined = std::make_unique<TH1D>(name, ";#Delta#phi;Weighted dijets",
                                         32, 0, TMath::Pi());
  combined->SetDirectory(nullptr);
  combined->Sumw2();
  for (std::size_t index = 0; index < samples.size(); ++index)
    {
      auto component = getSimulationDphi(
        paths[index], Form("%s_component_%d", name.Data(), samples[index]),
        rb, coneSize, centralityBin, index, centralityLow, centralityHigh,
        leadingCut, subleadingCut, useTruth);
      if (!component) return nullptr;
      const double scale = index == 2 ? 1.0
        : (eventCounts[2]/eventCounts[index])
          *(crossSections[index]/crossSections[2]);
        if (index ==0) continue; // skip the first sample for now, as it is not used in the final combination
      combined->Add(component.get(), scale);
    }
  return combined;
}

std::unique_ptr<TH1D> getCombinedSimulationDphi_EtaSeperated(
  const TString &name, read_binning &rb, const int coneSize,
  const int centralityBin, const double centralityLow,
  const double centralityHigh, const double leadingCut,
  const double subleadingCut, const bool useTruth,
  const bool requireEtaGap, const bool inclusivePairs)
{
  const std::array<int, 3> samples = {10, 20, 30};
  	// 2.502e+03 
    // 0.000003997
    // 6.217999999e-8
    //2.502E-9
  const std::array<double, 3> crossSections = { 0.000003997, 5.4067742e-8, 2.505e-9};
  std::array<double, 3> eventCounts = {0, 0, 0};
  std::array<TString, 3> paths;
  std::cout << "Loading simulation histograms for cone size " << coneSize
            << ", centrality bin " << centralityBin
            << ", leading cut " << leadingCut
            << ", subleading cut " << subleadingCut
            << ", and truth mode " << (useTruth ? "enabled" : "disabled")
  << std::endl;

  for (std::size_t index = 0; index < samples.size(); ++index)
  {
    paths[index] = Form("%s/rootfiles/TREE_MATCH_r%02d_v15_%d_new_ProdA_2024-00000030_sumeT.root", rb.get_code_location().c_str(), coneSize, samples[index]);
    std::unique_ptr<TFile> input(TFile::Open(paths[index], "READ"));
    TTree *stats = input && !input->IsZombie() ? dynamic_cast<TTree*>(input->Get("tn_stats")) : nullptr;
    float events = 0;
    if (!stats)
    {
      std::cerr << "Missing tn_stats in " << paths[index] << std::endl;
      return nullptr;
    }
    stats->SetBranchAddress("nevents", &events);
    stats->GetEntry(0);
    eventCounts[index] = events;
    input->Close();
    
  }
  std::cout << "Event counts: ";
  for (const auto& count : eventCounts) std::cout << count << " ";
  std::cout << std::endl;
  
  auto combined = std::make_unique<TH1D>(name, ";#Delta#phi;Weighted dijets", 32, 0, TMath::Pi());
  combined->SetDirectory(nullptr);
  combined->Sumw2();
  for (std::size_t index = 1; index < samples.size(); ++index)
  {

    auto component = getSimulationDphi_EtaSeperated(
      "/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root",
      Form("%s_component_%d", name.Data(), samples[index]),
      rb, coneSize, centralityBin, index, centralityLow, centralityHigh,
      leadingCut, subleadingCut, useTruth, requireEtaGap, inclusivePairs
    );
    if (!component) return nullptr;
    const double scale = index == 2 ? 1.0
      : (eventCounts[2]/eventCounts[index])*(crossSections[index]/crossSections[2]);
    if (index ==0) continue; // skip the first sample for now, as it is not used in the final combination

    combined->Add(component.get(), scale);
  }
  return combined;
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

} // namespace

void drawCOMBModulation_AA_v2(const int cone_size = 3,
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
  std::cout << "Using leading cut " << leadingCut
            << ", subleading cut " << subleadingCut
            << ", and signal cut " << signalCut << std::endl;

  read_binning downBinning("binning_COMBDown_AA.config");
  read_binning upBinning("binning_COMBUp_AA.config");
  const std::array<double, 3> v22Scales = { 1.0, downBinning.get_flow_sys(), upBinning.get_flow_sys()};
  const std::array<double, 3> v33Scales = { 1.0, downBinning.get_flow_v33_sys(), upBinning.get_flow_v33_sys()};
  std::cout << "Using v2 scales: " << v22Scales[0] << ", " << v22Scales[1] << ", " << v22Scales[2] << std::endl;
  std::cout << "Using v3 scales: " << v33Scales[0] << ", " << v33Scales[1] << ", " << v33Scales[2] << std::endl;

  const std::string system = "AA_cent_" + std::to_string(centrality_bin);
  const TString plotDirectory = Form("%s/dphi_plots", rb.get_code_location().c_str());
  gSystem->mkdir(plotDirectory, true);
  const TString inputPath = Form(
    "%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root", rb.get_code_location().c_str(), system.c_str(), cone_size
  );
  
  std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << "Cannot open " << inputPath << std::endl;
    return;
  }

  auto raw = std::make_unique<TH1D>("h_dphi_pairs", ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  auto etaSeparated = std::make_unique<TH1D>("h_dphi_eta_separated", ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  auto fittedEtaSeparated = std::make_unique<TH1D>( "h_dphi_eta_separated_fit", ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  std::array<std::unique_ptr<TH1D>, 3> backgrounds;
  for (int variation = 0; variation < 3; ++variation)
  {
    backgrounds[variation] = std::make_unique<TH1D>( Form("h_flow_background_%d", variation), ";#Delta#phi;Counts", 32, 0, TMath::Pi());
  }
  
  raw->Sumw2();
  etaSeparated->Sumw2();

  int fittedBins = 0;
  for (int i = leadingBin; i < nbins; ++i)
    for (int j = subleadingBin; j <= i; ++j)
      {
        TH1D *pairDphi = static_cast<TH1D*>(input->Get(Form( "h_dphi_exclusive_%d_%d", i, j)));
        TH1D *fitSource = static_cast<TH1D*>(input->Get(Form( "h_dphi_eta_inclusive_%d_%d", i, j)));
       
        if (!pairDphi || !fitSource) continue;
        raw->Add(pairDphi);
        etaSeparated->Add(fitSource);

        // std::unique_ptr<TH1D> fitHistogram(static_cast<TH1D*>(fitSource->Clone( Form("h_fit_work_%d_%d", i, j))));
        
        // fitHistogram->SetDirectory(nullptr);
        
        
        // TF1 fit(Form("flow_fit_%d_%d", i, j),
        //         "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))",
        //         kFlowFitLow, kFlowFitHigh
        // );
        
        // fit.SetParLimits(0, 0, 10000);
        // fit.SetParLimits(1, 0, 0.5);
        // fit.SetParLimits(2, 0, 0.5);
        // const int firstFitBin = fitHistogram->FindBin(kFlowFitLow);
        // const int lastFitBin = fitHistogram->FindBin(kFlowFitHigh);
        // const int fitBinCount = std::max(1, lastFitBin - firstFitBin);
        // const double fitCounts = fitHistogram->Integral(firstFitBin, lastFitBin);
        // if (fitCounts < fitBinCount)
        //   {
        //     fit.SetParameter(0, fitCounts/fitBinCount);
        //     fit.FixParameter(1, 0);
        //     fit.FixParameter(2, 0);
        //   }
        // else
        //   fit.SetParameters(fitCounts/fitBinCount, 0.01, 0.01);
        // if (fitCounts > 0) fitHistogram->Fit(&fit, "0RlQ");
        // if (fitCounts > 0 || pairDphi->Integral() > 0) ++fittedBins;
        // addFitShape(fittedEtaSeparated.get(), fit);
        // for (int variation = 0; variation < 3; ++variation)
        //   addNormalizedBackground(backgrounds[variation].get(), pairDphi, fit,
        //                           v22Scales[variation], v33Scales[variation]);
      }

      // add global fit
      TF1 globalFit("global_flow_fit", "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))", kFlowFitLow, kFlowFitHigh);
      globalFit.SetParLimits(0, 0, 10000);
      globalFit.SetParLimits(1, 0, 0.5);
      globalFit.SetParLimits(2, 0, 0.5);

      const int firstFitBin = etaSeparated->FindBin(kFlowFitLow);
      const int lastFitBin = etaSeparated->FindBin(kFlowFitHigh);
      const int fitBinCount = std::max(1, lastFitBin - firstFitBin);
      const double fitCounts = etaSeparated->Integral(firstFitBin, lastFitBin);

      if (fitCounts < fitBinCount)
      {
        globalFit.SetParameter(0, fitCounts/fitBinCount);
        globalFit.FixParameter(1, 0);
        globalFit.FixParameter(2, 0);
      }
      else
        globalFit.SetParameters(fitCounts/fitBinCount, 0.01, 0.01);
      if (fitCounts > 0) etaSeparated->Fit(&globalFit, "0RlQ");
      // if (fitCounts > 0 || raw->Integral() > 0) ++fittedBins;
      addFitShape(fittedEtaSeparated.get(), globalFit);
      // TF1 * fitNorms[3];
      std::array<std::unique_ptr<TF1>, 3> fitNorms;
      for (int variation = 0; variation < 3; ++variation)
      {
        fitNorms[variation] = addNormalizedBackground(backgrounds[variation].get(), raw.get(), globalFit,
                                v22Scales[variation], v33Scales[variation]);  
        fitNorms[variation]->SetName(Form("f_flow_fit_%d", variation));
      }

  if (raw->Integral() <= 0)
  {
    std::cerr << "No selected #Delta#phi pairs for " << system << std::endl;
    return;
  }

  auto nominalSubtracted = subtractBackground( raw.get(), backgrounds[0].get(), "h_dphi_nominal_subtracted");
  auto downSubtracted = subtractBackground( raw.get(), backgrounds[1].get(), "h_dphi_COMBDown_subtracted");
  auto upSubtracted = subtractBackground( raw.get(), backgrounds[2].get(), "h_dphi_COMBUp_subtracted");
  auto normalizationRegion = makeRegion( raw.get(), "h_flow_normalization_region", kNormalizationLow, kNormalizationHigh);
  auto signalRegion = makeRegion(raw.get(), "h_signal_region", signalCut, TMath::Pi() + 1e-6);

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
  etaSeparated->SetLineColor(kCyan + 2);
  etaSeparated->SetLineWidth(2);
  etaSeparated->SetLineStyle(9);
  fittedEtaSeparated->SetLineColor(kOrange + 7);
  fittedEtaSeparated->SetLineWidth(3);
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
  etaSeparated->Draw("hist same");
  fittedEtaSeparated->Draw("hist same");
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

  TLegend legend(0.47, 0.49, 0.94, 0.93);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(raw.get(), "Raw #Delta#phi pairs", "l");
  legend.AddEntry(backgrounds[0].get(), "Nominal flow modulation", "l");
  legend.AddEntry(etaSeparated.get(), "#Delta#eta-separated fit input", "l");
  legend.AddEntry(fittedEtaSeparated.get(), "Fitted #Delta#eta histogram", "l");
  legend.AddEntry(nominalSubtracted.get(), "Nominal modulation sub.", "l");
  legend.AddEntry(downSubtracted.get(), "COMBDown: 0.7(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(upSubtracted.get(), "COMBUp: 1.3(v_{2,2},v_{3,3})", "l");
  legend.AddEntry(signalRegion.get(), "Signal region", "f");
  legend.AddEntry(normalizationRegion.get(), "Flow normalization region", "f");
  legend.Draw();

  const TString outputStem = Form(
    "%s/dphi_COMB_modulation_%s_r%02d",
    plotDirectory.Data(), system.c_str(), cone_size);
  canvas.SaveAs(outputStem + ".pdf");
  if (kSavePng) canvas.SaveAs(outputStem + ".png");

  // Cache the aggregate data components for fast redraws, matching the
  // combined-simulation COMB cache.
  std::unique_ptr<TFile> output(TFile::Open(outputStem + ".root", "RECREATE"));
  if (output && !output->IsZombie())
  {
    raw->Write();
    etaSeparated->Write();
    fittedEtaSeparated->Write();
    backgrounds[0]->Write("h_flow_background_nominal");
    backgrounds[1]->Write("h_flow_background_COMBDown");
    backgrounds[2]->Write("h_flow_background_COMBUp");
    nominalSubtracted->Write();
    downSubtracted->Write();
    upSubtracted->Write();
    signalRegion->Write();
    normalizationRegion->Write();
    output->Write();
    globalFit.Write("global_flow_fit");
    for (int variation = 0; variation < 3; ++variation)
      fitNorms[variation]->Write();
  }

  // drawSimulationDphi("", "combined", 
  //   cone_size, centrality_bin,
  //   centralityBins[centrality_bin],
  //   centralityBins[centrality_bin + 1],
  //   leadingCut, subleadingCut, plotDirectory, rb, 0
  // );
  std::cout << "Wrote aggregate COMB modulation plot for " << system
            << " using " << fittedBins << " populated pT-bin fits to "
            << outputStem << ".pdf" << std::endl;
}
