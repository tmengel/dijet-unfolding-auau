#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TSystem.h"

#include "dlUtility.h"
#include "histo_opps.h"
#include "read_binning.h"

namespace
{
TH1D *requiredHistogram(TFile *file, const TString &name)
{
  TH1D *histogram = file ? dynamic_cast<TH1D*>(file->Get(name)) : nullptr;
  if (!histogram)
    std::cerr << "Missing " << name << " in "
              << (file ? file->GetName() : "null file") << std::endl;
  return histogram;
}

std::unique_ptr<TH1D> makeFinalXj(
  TH1D *flat, const TString &name, const int nbins, const float *ptBins,
  const float *xjBins, const double *doubleXjBins, const int firstBin,
  const double firstXj, const int leadingLow, const int leadingHigh,
  const int subleadingLow, const int subleadingHigh)
{
  auto pt1pt2 = std::make_unique<TH2D>(name + "_pt1pt2", "", nbins, ptBins,
                                       nbins, ptBins);
  auto xj = std::make_unique<TH1D>(name + "_xj", "", nbins, xjBins);
  auto finalXj = std::make_unique<TH1D>(name + "_final", "", nbins, xjBins);
  pt1pt2->SetDirectory(nullptr);
  xj->SetDirectory(nullptr);
  finalXj->SetDirectory(nullptr);
  histo_opps::make_sym_pt1pt2(flat, pt1pt2.get(), nbins);
  histo_opps::project_xj(pt1pt2.get(), xj.get(), nbins, leadingLow,
                         leadingHigh, subleadingLow, subleadingHigh);
  histo_opps::normalize_histo(xj.get(), nbins);
  histo_opps::finalize_xj(xj.get(), finalXj.get(), nbins, firstXj);
  TH1D *rebinned = static_cast<TH1D*>(finalXj->Rebin(
    nbins - firstBin, name, &doubleXjBins[firstBin]));
  rebinned->SetDirectory(nullptr);
  return std::unique_ptr<TH1D>(rebinned);
}

std::unique_ptr<TH1D> makeResidual(const TH1D *unfolded,
                                   const TH1D *truth,
                                   const TString &name)
{
  auto residual = std::unique_ptr<TH1D>(
    static_cast<TH1D*>(truth->Clone(name)));
  residual->SetDirectory(nullptr);
  residual->Reset("ICES");
  residual->SetTitle(";x_{J};(Unfolded - Truth) / Truth");
  for (int bin = 1; bin <= truth->GetNbinsX(); ++bin)
    {
      const double t = truth->GetBinContent(bin);
      const double u = unfolded->GetBinContent(bin);
      if (t == 0) continue;
      const double et = truth->GetBinError(bin);
      const double eu = unfolded->GetBinError(bin);
      residual->SetBinContent(bin, u/t - 1.0);
      residual->SetBinError(bin, std::sqrt(
        std::pow(eu/t, 2) + std::pow(u*et/(t*t), 2)));
    }
  return residual;
}

std::pair<double, double> closureMetrics(const TH1D *residual,
                                         const TH1D *truth)
{
  double maximum = 0;
  double sum = 0;
  int populated = 0;
  for (int bin = 1; bin <= residual->GetNbinsX(); ++bin)
    {
      if (truth->GetBinContent(bin) == 0) continue;
      const double value = std::abs(residual->GetBinContent(bin));
      if (!std::isfinite(value)) continue;
      maximum = std::max(maximum, value);
      sum += value;
      ++populated;
    }
  return {maximum, populated ? sum/populated : 0.0};
}
}

void drawClosureTests_AA(const int cone_size = 3,
                         const int centrality_bin = 0,
                         const int selected_iteration = 1,
                         const std::string config = "binning_AA.config")
{
  constexpr int niterations = 10;
  if (selected_iteration < 0 || selected_iteration >= niterations)
    {
      std::cerr << "selected_iteration must be in [0,9]" << std::endl;
      return;
    }
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(config);
  const int nbins = rb.get_nbins();
  const int ranges = rb.get_measure_bins();
  std::vector<float> ptBins(nbins + 1);
  std::vector<float> xjBins(nbins + 1);
  std::vector<double> doubleXjBins(nbins + 1);
  rb.get_pt_bins(ptBins.data());
  rb.get_xj_bins(xjBins.data());
  const double firstXj = rb.get_first_xj();
  int firstBin = 0;
  for (int bin = 0; bin <= nbins; ++bin)
    {
      doubleXjBins[bin] = xjBins[bin];
      if (xjBins[bin] >= firstXj && firstBin == 0) firstBin = bin;
    }
  std::vector<int> measureBins(ranges + 1);
  for (int range = 0; range <= ranges; ++range)
    measureBins[range] = rb.get_measure_region(range);
  const int subleadingBin = rb.get_measure_subleading_bin();
  float centralityBins[5] = {0};
  rb.get_centrality_bins(centralityBins);
  const std::string system = "AA_cent_" + std::to_string(centrality_bin);
  const std::array<std::string, 2> modes = {"FULL", "HALF"};

  gSystem->mkdir(Form("%s/closure_plots", rb.get_code_location().c_str()), true);
  gSystem->mkdir(Form("%s/closure_results", rb.get_code_location().c_str()), true);

  for (const std::string &mode : modes)
    {
      const TString inputPath = Form(
        "%s/response_matrices/response_matrix_%s_r%02d_%s_nominal.root",
        rb.get_code_location().c_str(), system.c_str(), cone_size, mode.c_str());
      std::unique_ptr<TFile> input(TFile::Open(inputPath, "READ"));
      if (!input || input->IsZombie())
        {
          std::cerr << "Cannot open closure response " << inputPath << std::endl;
          continue;
        }
      TH1D *flatTruth = requiredHistogram(input.get(), "h_truth_flat_pt1pt2");
      TH1D *flatReco = requiredHistogram(input.get(), "h_reco_flat_pt1pt2");
      std::array<TH1D*, niterations> flatUnfold{};
      for (int iteration = 0; iteration < niterations; ++iteration)
        flatUnfold[iteration] = requiredHistogram(
          input.get(), Form("h_flat_unfold_pt1pt2_%d", iteration));
      if (!flatTruth || !flatReco || !flatUnfold[selected_iteration] ||
          flatTruth->Integral() == 0 || flatReco->Integral() == 0)
        {
          std::cerr << "Empty or incomplete " << mode << " closure input for "
                    << system << std::endl;
          continue;
        }

      const std::string modeLower = mode == "FULL" ? "full" : "half";
      std::unique_ptr<TFile> output(TFile::Open(Form(
        "%s/closure_results/%s_closure_%s_r%02d.root",
        rb.get_code_location().c_str(), modeLower.c_str(), system.c_str(),
        cone_size), "RECREATE"));

      for (int range = 0; range < ranges; ++range)
        {
          auto truth = makeFinalXj(
            flatTruth, Form("h_%s_truth_range_%d", modeLower.c_str(), range),
            nbins, ptBins.data(), xjBins.data(), doubleXjBins.data(), firstBin,
            firstXj, measureBins[range], measureBins[range + 1],
            subleadingBin, nbins - 2);
          auto reco = makeFinalXj(
            flatReco, Form("h_%s_reco_range_%d", modeLower.c_str(), range),
            nbins, ptBins.data(), xjBins.data(), doubleXjBins.data(), firstBin,
            firstXj, measureBins[range], measureBins[range + 1],
            subleadingBin, nbins - 2);
          std::array<std::unique_ptr<TH1D>, niterations> unfolded;
          std::array<std::unique_ptr<TH1D>, niterations> residual;
          for (int iteration = 0; iteration < niterations; ++iteration)
            {
              unfolded[iteration] = makeFinalXj(
                flatUnfold[iteration], Form("h_%s_unfold_range_%d_iter_%d",
                modeLower.c_str(), range, iteration), nbins, ptBins.data(),
                xjBins.data(), doubleXjBins.data(), firstBin, firstXj,
                measureBins[range], measureBins[range + 1], subleadingBin,
                nbins - 2);
              residual[iteration] = makeResidual(
                unfolded[iteration].get(), truth.get(), Form(
                "h_%s_residual_range_%d_iter_%d", modeLower.c_str(), range,
                iteration));
            }

          const auto metrics = closureMetrics(
            residual[selected_iteration].get(), truth.get());
          std::cout << "CLOSURE_METRIC mode=" << modeLower
                    << " cent=" << centrality_bin << " range=" << range
                    << " iteration=" << selected_iteration + 1
                    << " max_abs_fractional=" << metrics.first
                    << " mean_abs_fractional=" << metrics.second << std::endl;

          TCanvas *comparison = new TCanvas(
            Form("c_%s_closure_%d", modeLower.c_str(), range), "", 600, 720);
          dlutility::ratioPanelCanvas(comparison);
          comparison->cd(1);
          truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
          truth->SetLineColor(kRed + 1);
          truth->SetMarkerColor(kRed + 1);
          truth->SetMarkerStyle(20);
          reco->SetLineColor(kGray + 2);
          reco->SetMarkerColor(kGray + 2);
          reco->SetMarkerStyle(24);
          unfolded[selected_iteration]->SetLineColor(kAzure - 6);
          unfolded[selected_iteration]->SetMarkerColor(kAzure - 6);
          unfolded[selected_iteration]->SetMarkerStyle(21);
          const double maximum = std::max({truth->GetMaximum(), reco->GetMaximum(),
            unfolded[selected_iteration]->GetMaximum()});
          truth->SetMinimum(0);
          truth->SetMaximum(1.45*maximum);
          dlutility::SetFont(truth.get(), 42, 0.05);
          truth->Draw("E1 p");
          reco->Draw("E1 p same");
          unfolded[selected_iteration]->Draw("E1 p same");
          dlutility::DrawSPHENIX(0.20, 0.86);
          dlutility::drawText(Form("%s closure", mode.c_str()), 0.20, 0.76);
          dlutility::drawText(Form("%d Bayesian iterations",
            selected_iteration + 1), 0.20, 0.71);
          dlutility::drawText(Form("%.1f #leq #it{p}_{T,1} < %.1f GeV",
            ptBins[measureBins[range]], ptBins[measureBins[range + 1]]),
            0.20, 0.66);
          dlutility::drawText(Form("#it{p}_{T,2} #geq %.1f GeV",
            ptBins[subleadingBin]), 0.20, 0.61);
          dlutility::drawText(Form("%.0f - %.0f %%", centralityBins[centrality_bin],
            centralityBins[centrality_bin + 1]), 0.20, 0.56);
          TLegend legend(0.57, 0.69, 0.88, 0.88);
          legend.SetBorderSize(0);
          legend.AddEntry(truth.get(), "Pseudo-data truth", "p");
          legend.AddEntry(reco.get(), "Pseudo-data reco", "p");
          legend.AddEntry(unfolded[selected_iteration].get(), "Unfolded", "p");
          legend.Draw();

          comparison->cd(2);
          double residualMaximum = 0.25;
          for (int bin = 1; bin <= residual[selected_iteration]->GetNbinsX(); ++bin)
            residualMaximum = std::max(residualMaximum,
              std::abs(residual[selected_iteration]->GetBinContent(bin)) +
              residual[selected_iteration]->GetBinError(bin));
          residual[selected_iteration]->SetMinimum(-1.15*residualMaximum);
          residual[selected_iteration]->SetMaximum(1.15*residualMaximum);
          dlutility::SetFont(residual[selected_iteration].get(), 42, 0.10,
                             0.07, 0.07, 0.07);
          residual[selected_iteration]->SetLineColor(kBlack);
          residual[selected_iteration]->SetMarkerColor(kBlack);
          residual[selected_iteration]->SetMarkerStyle(20);
          residual[selected_iteration]->Draw("E1 p");
          TLine zero(residual[selected_iteration]->GetXaxis()->GetXmin(), 0,
                     residual[selected_iteration]->GetXaxis()->GetXmax(), 0);
          zero.SetLineColor(kRed + 2);
          zero.SetLineStyle(2);
          zero.Draw();
          const TString comparisonStem = Form(
            "%s/closure_plots/%s_closure_%s_r%02d_range_%d",
            rb.get_code_location().c_str(), modeLower.c_str(), system.c_str(),
            cone_size, range);
          comparison->SaveAs(comparisonStem + ".pdf");
          comparison->SaveAs(comparisonStem + ".png");
          delete comparison;

          TCanvas iterations(Form("c_%s_iterations_%d", modeLower.c_str(), range),
                             "", 650, 550);
          const std::array<int, 5> shown = {0, 1, 3, 5, 7};
          const std::array<int, 5> colors = {kBlue + 1, kCyan + 2, kGreen + 2,
                                             kViolet + 1, kRed + 1};
          double iterationMaximum = 0.25;
          for (const int iteration : shown)
            for (int bin = 1; bin <= residual[iteration]->GetNbinsX(); ++bin)
              iterationMaximum = std::max(iterationMaximum,
                std::abs(residual[iteration]->GetBinContent(bin)) +
                residual[iteration]->GetBinError(bin));
          TLegend iterationLegend(0.66, 0.66, 0.88, 0.88);
          iterationLegend.SetBorderSize(0);
          for (std::size_t index = 0; index < shown.size(); ++index)
            {
              TH1D *histogram = residual[shown[index]].get();
              histogram->SetLineColor(colors[index]);
              histogram->SetMarkerColor(colors[index]);
              histogram->SetMarkerStyle(20 + index);
              histogram->SetMinimum(-1.15*iterationMaximum);
              histogram->SetMaximum(1.15*iterationMaximum);
              dlutility::SetFont(histogram, 42, 0.045);
              histogram->Draw(index == 0 ? "E1 p" : "E1 p same");
              iterationLegend.AddEntry(histogram, Form("%d iterations",
                shown[index] + 1), "p");
            }
          TLine iterationZero(residual[0]->GetXaxis()->GetXmin(), 0,
                              residual[0]->GetXaxis()->GetXmax(), 0);
          iterationZero.SetLineStyle(2);
          iterationZero.Draw();
          dlutility::DrawSPHENIX(0.16, 0.90);
          dlutility::drawText(Form("%s closure, %.0f - %.0f %%", mode.c_str(),
            centralityBins[centrality_bin], centralityBins[centrality_bin + 1]),
            0.16, 0.79);
          dlutility::drawText(Form("%.1f #leq #it{p}_{T,1} < %.1f GeV",
            ptBins[measureBins[range]], ptBins[measureBins[range + 1]]),
            0.16, 0.74);
          iterationLegend.Draw();
          const TString iterationStem = Form(
            "%s/closure_plots/%s_closure_iterations_%s_r%02d_range_%d",
            rb.get_code_location().c_str(), modeLower.c_str(), system.c_str(),
            cone_size, range);
          iterations.SaveAs(iterationStem + ".pdf");
          iterations.SaveAs(iterationStem + ".png");

          if (output && !output->IsZombie())
            {
              output->cd();
              truth->Write();
              reco->Write();
              for (int iteration = 0; iteration < niterations; ++iteration)
                {
                  unfolded[iteration]->Write();
                  residual[iteration]->Write();
                }
              TParameter<double>(Form("max_abs_fractional_range_%d_iter_%d",
                range, selected_iteration), metrics.first).Write();
              TParameter<double>(Form("mean_abs_fractional_range_%d_iter_%d",
                range, selected_iteration), metrics.second).Write();
            }
        }
      if (output && !output->IsZombie())
        {
          output->Close();
        }
    }
}
