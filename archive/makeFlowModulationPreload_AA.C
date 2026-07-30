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
#include "TH2D.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"

#include "read_binning.h"

namespace
{
struct FlowVariation
{
  std::string name;
  double v22Scale = 1.0;
  double v33Scale = 1.0;
  std::unique_ptr<TH2D> exclusive;
  std::unique_ptr<TH2D> inclusive;
  std::unique_ptr<TH1D> flatExclusive;
  std::unique_ptr<TH1D> flatInclusive;
};

void detach(TH1 *histogram)
{
  if (histogram) histogram->SetDirectory(nullptr);
}

void setSymmetric(TH2D *histogram, const int i, const int j,
                  const double value, const double error)
{
  for (const int bin : {histogram->GetBin(i + 1, j + 1),
                        histogram->GetBin(j + 1, i + 1)})
    {
      histogram->SetBinContent(bin, value);
      histogram->SetBinError(bin, error);
    }
}

void setSymmetric(TH1D *histogram, const int nbins, const int i, const int j,
                  const double value, const double error)
{
  for (const int bin : {1 + i + nbins*j, 1 + j + nbins*i})
    {
      histogram->SetBinContent(bin, value);
      histogram->SetBinError(bin, error);
    }
}

std::pair<double, double> subtractFlow(
  TH1D *dphi, const TF1 &fittedShape, const double v22Scale,
  const double v33Scale, const double dphiCut, const double probability)
{
  constexpr double normalizationLow = 0.8;
  constexpr double normalizationHigh = 2.5;
  TF1 normalized(fittedShape);
  normalized.SetParameter(1, v22Scale*fittedShape.GetParameter(1));
  normalized.SetParameter(2, v33Scale*fittedShape.GetParameter(2));

  const int firstNormBin = dphi->FindBin(normalizationLow);
  const int lastNormBin = dphi->FindBin(normalizationHigh);
  const int normBinCount = std::max(1, lastNormBin - firstNormBin);
  const double dataLevel = dphi->Integral(firstNormBin, lastNormBin)/normBinCount;
  const double fitLevel = normalized.Integral(normalizationLow, normalizationHigh)
    /(normalizationHigh - normalizationLow);
  normalized.SetParameter(0, normalized.GetParameter(0) + dataLevel - fitLevel);
  normalized.SetRange(0, TMath::Pi());

  const double signalBinCount = (TMath::Pi() - dphiCut)/(TMath::Pi()/32.0);
  const double background = signalBinCount*normalized.Integral(dphiCut, TMath::Pi())
    /(TMath::Pi() - dphiCut);
  const double signalRegion = dphi->Integral(dphi->FindBin(dphiCut), -1);
  double signal = std::max(0.0, signalRegion - background);
  double error = std::sqrt(std::max(0.0, signalRegion));
  if (probability > 0)
    {
      signal /= probability;
      error /= probability;
    }
  else
    signal = error = 0;
  return {signal, error};
}

void copyInputAndWriteVariation(TFile *input, const TString &outputPath,
                                FlowVariation &variation)
{
  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie())
    {
      std::cerr << "Cannot create " << outputPath << std::endl;
      return;
    }
  TIter next(input->GetListOfKeys());
  while (TKey *key = static_cast<TKey*>(next()))
    {
      std::unique_ptr<TObject> object(key->ReadObj());
      output->cd();
      object->Write(key->GetName());
    }
  output->cd();
  variation.exclusive->Write("h_pt1_pt2_signal_exclusive", TObject::kOverwrite);
  variation.inclusive->Write("h_pt1_pt2_signal_inclusive", TObject::kOverwrite);
  variation.flatExclusive->Write("h_data_flat_pt1pt2_signal_exclusive", TObject::kOverwrite);
  variation.flatInclusive->Write("h_data_flat_pt1pt2_signal_inclusive", TObject::kOverwrite);
  output->Write();
  std::cout << "Created " << outputPath << " with coherent (v22,v33) scales ("
            << variation.v22Scale << "," << variation.v33Scale << ")" << std::endl;
}
}

void makeFlowModulationPreload_AA(
  const int cone_size = 3, const int centrality_bin = 0,
  const std::string down_config = "binning_COMBDown_AA.config",
  const std::string up_config = "binning_COMBUp_AA.config")
{
  read_binning nominalBinning(std::getenv("AUAU_CONFIG"));
  read_binning downBinning(down_config);
  read_binning upBinning(up_config);
  const int nbins = nominalBinning.get_nbins();
  const double dphiCut = nominalBinning.get_dphicut();
  const std::string system = "AA_cent_" + std::to_string(centrality_bin);
  const TString nominalPath = Form(
    "%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root",
    nominalBinning.get_code_location().c_str(), system.c_str(), cone_size);
  std::unique_ptr<TFile> input(TFile::Open(nominalPath, "READ"));
  if (!input || input->IsZombie())
    {
      std::cerr << "Cannot open nominal preload " << nominalPath << std::endl;
      return;
    }

  std::unique_ptr<TFile> probabilityFile(TFile::Open(Form(
    "%s/unfolding_hists/probability_hists_AA_r%02d.root",
    nominalBinning.get_code_location().c_str(), cone_size), "READ"));
  TH1D *pt2Probability = probabilityFile
    ? static_cast<TH1D*>(probabilityFile->Get(Form(
        "h_pt2_bin_log_correction_%d", centrality_bin))) : nullptr;
  TH2D *template2D = static_cast<TH2D*>(input->Get("h_pt1_pt2_signal_exclusive"));
  TH1D *templateFlat = static_cast<TH1D*>(input->Get(
    "h_data_flat_pt1pt2_signal_exclusive"));
  if (!pt2Probability || !template2D || !templateFlat)
    {
      std::cerr << "Missing probability correction or signal templates for "
                << system << std::endl;
      return;
    }

  std::array<FlowVariation, 2> variations = {{
    {downBinning.get_flow_systematic_name(), downBinning.get_flow_sys(),
     downBinning.get_flow_v33_sys()},
    {upBinning.get_flow_systematic_name(), upBinning.get_flow_sys(),
     upBinning.get_flow_v33_sys()}
  }};
  for (auto &variation : variations)
    {
      variation.exclusive.reset(static_cast<TH2D*>(template2D->Clone()));
      variation.inclusive.reset(static_cast<TH2D*>(template2D->Clone()));
      variation.flatExclusive.reset(static_cast<TH1D*>(templateFlat->Clone()));
      variation.flatInclusive.reset(static_cast<TH1D*>(templateFlat->Clone()));
      detach(variation.exclusive.get());
      detach(variation.inclusive.get());
      detach(variation.flatExclusive.get());
      detach(variation.flatInclusive.get());
      variation.exclusive->Reset("ICES");
      variation.inclusive->Reset("ICES");
      variation.flatExclusive->Reset("ICES");
      variation.flatInclusive->Reset("ICES");
    }

  auto makeMap = [&](const char *name, const char *title)
    {
      auto result = std::unique_ptr<TH2D>(static_cast<TH2D*>(template2D->Clone(name)));
      detach(result.get());
      result->Reset("ICES");
      result->SetTitle(title);
      return result;
    };
  auto v22Nominal = makeMap("h_v22_nominal", ";p_{T,1} bin;p_{T,2} bin;v_{2,2}");
  auto v22Down = makeMap("h_v22_COMBDown", ";p_{T,1} bin;p_{T,2} bin;0.7v_{2,2}");
  auto v22Up = makeMap("h_v22_COMBUp", ";p_{T,1} bin;p_{T,2} bin;1.3v_{2,2}");
  auto v33Nominal = makeMap("h_v33_nominal", ";p_{T,1} bin;p_{T,2} bin;v_{3,3}");
  auto v33Down = makeMap("h_v33_COMBDown", ";p_{T,1} bin;p_{T,2} bin;0.7v_{3,3}");
  auto v33Up = makeMap("h_v33_COMBUp", ";p_{T,1} bin;p_{T,2} bin;1.3v_{3,3}");
  auto signalNominal = makeMap("h_signal_nominal", ";p_{T,1};p_{T,2};signal");

  gSystem->mkdir("systematic_plots", true);
  const TString diagnosticPdf = Form(
    "%s/systematic_plots/flow_modulation_diagnostics_%s_r%02d.pdf",
    nominalBinning.get_code_location().c_str(), system.c_str(), cone_size);
  const TString diagnosticRoot = Form(
    "%s/systematic_plots/flow_modulation_diagnostics_%s_r%02d.root",
    nominalBinning.get_code_location().c_str(), system.c_str(), cone_size);
  TCanvas canvas("flow_diagnostic", "flow_diagnostic", 900, 850);
  canvas.Print(diagnosticPdf + "[");
  int diagnosticPages = 0;

  for (int i = 0; i < nbins; ++i)
    for (int j = 0; j <= i; ++j)
      {
        TH1D *fitHistogram = static_cast<TH1D*>(input->Get(Form(
          "h_dphi_eta_inclusive_%d_%d", i, j)));
        TH1D *exclusive = static_cast<TH1D*>(input->Get(Form(
          "h_dphi_exclusive_%d_%d", i, j)));
        TH1D *inclusive = static_cast<TH1D*>(input->Get(Form(
          "h_dphi_inclusive_%d_%d", i, j)));
        if (!fitHistogram || !exclusive || !inclusive) continue;

        TF1 fit(Form("flow_fit_nominal_%d_%d", i, j),
                "[0]*(1+2*[1]*cos(2*x)+2*[2]*cos(3*x))", 0.0, 2.5);
        fit.SetParLimits(0, 0, 10000);
        fit.SetParLimits(1, 0, 0.5);
        fit.SetParLimits(2, 0, 0.5);
        const int firstFitBin = fitHistogram->FindBin(0.0);
        const int lastFitBin = fitHistogram->FindBin(2.5);
        const int fitBins = std::max(1, lastFitBin - firstFitBin);
        const double fitCounts = fitHistogram->Integral(firstFitBin, lastFitBin);
        if (fitCounts < fitBins)
          {
            fit.SetParameter(0, fitCounts/fitBins);
            fit.FixParameter(1, 0);
            fit.FixParameter(2, 0);
          }
        else
          {
            fit.SetParameters(fitCounts/fitBins, 0.01, 0.01);
          }
        if (fitCounts > 0) fitHistogram->Fit(&fit, "0RlQ");
        const double nominalV22 = fit.GetParameter(1);
        const double nominalV33 = fit.GetParameter(2);
        setSymmetric(v22Nominal.get(), i, j, nominalV22, fit.GetParError(1));
        setSymmetric(v22Down.get(), i, j, variations[0].v22Scale*nominalV22, 0);
        setSymmetric(v22Up.get(), i, j, variations[1].v22Scale*nominalV22, 0);
        setSymmetric(v33Nominal.get(), i, j, nominalV33, fit.GetParError(2));
        setSymmetric(v33Down.get(), i, j, variations[0].v33Scale*nominalV33, 0);
        setSymmetric(v33Up.get(), i, j, variations[1].v33Scale*nominalV33, 0);

        const double probability = pt2Probability->GetBinContent(std::min(i, j) + 1);
        const auto nominalExclusive = subtractFlow(
          exclusive, fit, 1.0, 1.0, dphiCut, probability);
        setSymmetric(signalNominal.get(), i, j,
                     nominalExclusive.first, nominalExclusive.second);
        std::array<double, 2> variedSignals = {0, 0};
        for (std::size_t index = 0; index < variations.size(); ++index)
          {
            auto &variation = variations[index];
            const auto variedExclusive = subtractFlow(
              exclusive, fit, variation.v22Scale, variation.v33Scale,
              dphiCut, probability);
            const auto variedInclusive = subtractFlow(
              inclusive, fit, variation.v22Scale, variation.v33Scale,
              dphiCut, 1.0);
            variedSignals[index] = variedExclusive.first;
            setSymmetric(variation.exclusive.get(), i, j,
                         variedExclusive.first, variedExclusive.second);
            setSymmetric(variation.inclusive.get(), i, j,
                         variedInclusive.first, variedInclusive.second);
            setSymmetric(variation.flatExclusive.get(), nbins, i, j,
                         variedExclusive.first, variedExclusive.second);
            setSymmetric(variation.flatInclusive.get(), nbins, i, j,
                         variedInclusive.first, variedInclusive.second);
          }

        if (fitCounts > 0 || exclusive->Integral() > 0)
          {
            ++diagnosticPages;
            canvas.Clear();
            canvas.Divide(1, 2);
            canvas.cd(1);
            fitHistogram->SetTitle(Form(
              "%s, p_{T} bins (%d,%d);#Delta#phi;counts", system.c_str(), i, j));
            fitHistogram->SetMarkerStyle(20);
            fitHistogram->Draw("E");
            TF1 nominalFit(fit);
            TF1 downFit(fit);
            TF1 upFit(fit);
            nominalFit.SetLineColor(kBlack);
            downFit.SetParameter(1, variations[0].v22Scale*nominalV22);
            downFit.SetParameter(2, variations[0].v33Scale*nominalV33);
            downFit.SetLineColor(kBlue + 1);
            downFit.SetLineStyle(2);
            upFit.SetParameter(1, variations[1].v22Scale*nominalV22);
            upFit.SetParameter(2, variations[1].v33Scale*nominalV33);
            upFit.SetLineColor(kRed + 1);
            upFit.SetLineStyle(2);
            nominalFit.Draw("same");
            downFit.Draw("same");
            upFit.Draw("same");
            TLegend legend(0.42, 0.61, 0.89, 0.89);
            legend.SetBorderSize(0);
            legend.AddEntry(&nominalFit, Form("nominal: v_{2,2}=%.4g, v_{3,3}=%.4g",
                                             nominalV22, nominalV33), "l");
            legend.AddEntry(&downFit, Form("COMBDown: %.4g, %.4g",
                                          variations[0].v22Scale*nominalV22,
                                          variations[0].v33Scale*nominalV33), "l");
            legend.AddEntry(&upFit, Form("COMBUp: %.4g, %.4g",
                                        variations[1].v22Scale*nominalV22,
                                        variations[1].v33Scale*nominalV33), "l");
            legend.Draw();
            canvas.cd(2);
            TH1D signalComparison("signal_comparison", ";variation;extracted exclusive signal", 3, 0, 3);
            signalComparison.SetDirectory(nullptr);
            signalComparison.GetXaxis()->SetBinLabel(1, "nominal");
            signalComparison.GetXaxis()->SetBinLabel(2, "COMBDown");
            signalComparison.GetXaxis()->SetBinLabel(3, "COMBUp");
            signalComparison.SetBinContent(1, nominalExclusive.first);
            signalComparison.SetBinError(1, nominalExclusive.second);
            signalComparison.SetBinContent(2, variedSignals[0]);
            signalComparison.SetBinContent(3, variedSignals[1]);
            signalComparison.SetMarkerStyle(20);
            signalComparison.Draw("E1");
            canvas.Print(diagnosticPdf);
          }
      }
  canvas.Print(diagnosticPdf + "]");

  std::unique_ptr<TFile> diagnostics(TFile::Open(diagnosticRoot, "RECREATE"));
  if (diagnostics && !diagnostics->IsZombie())
    {
      v22Nominal->Write();
      v22Down->Write();
      v22Up->Write();
      v33Nominal->Write();
      v33Down->Write();
      v33Up->Write();
      signalNominal->Write();
      variations[0].exclusive->Write("h_signal_COMBDown");
      variations[1].exclusive->Write("h_signal_COMBUp");
      diagnostics->Write();
    }

  for (auto &variation : variations)
    copyInputAndWriteVariation(input.get(), Form(
      "%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_%s.root",
      nominalBinning.get_code_location().c_str(), system.c_str(), cone_size,
      variation.name.c_str()), variation);
  std::cout << "Wrote " << diagnosticPages << " fit/signal diagnostic pages to "
            << diagnosticPdf << std::endl;
}
