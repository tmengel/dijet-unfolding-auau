#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"

#include "dlUtility.h"
#include "read_binning.h"

namespace
{
  // Fill one probability histogram from raw jet counts.  This is equivalent to
  // the calculation in getProbability_v0.C, but keeps the Poisson mean in
  // count/event form so no intermediate unfolding histogram is needed.
  void fillProbabilityHistogram(TH1D* probability, const TH1D* jetCounts,
                                const double nEvents, const double area1,
                                const double area2)
  {
    if (!probability || !jetCounts || nEvents <= 0.0)
    {
      return;
    }

    const int nBins = jetCounts->GetNbinsX();
    for (int bin = 1; bin <= probability->GetNbinsX(); ++bin)
    {
      const double pt = probability->GetBinCenter(bin);
      const int countBin = jetCounts->GetXaxis()->FindBin(pt);
      const double nCurrent = jetCounts->GetBinContent(countBin);
      // Include ROOT's overflow bin: it represents jets above the last pT bin.
      const double nAbove = jetCounts->Integral(countBin + 1, nBins + 1);

      const double lambda = (area1 * nCurrent + area2 * nAbove) / nEvents;
      const double probabilityValue = std::exp(-lambda);
      const double lambdaVariance =
        (area1 * area1 * nCurrent + area2 * area2 * nAbove) /
        (nEvents * nEvents);

      probability->SetBinContent(bin, probabilityValue);
      probability->SetBinError(bin,
                               probabilityValue * std::sqrt(lambdaVariance));
    }
  }
}

void getProbability(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string& configfile = "binning_AA.config",
  const std::string& infile = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root")
{
  (void) centrality_bin;  // All configured centrality bins are produced.

  read_binning rb(configfile.c_str());
  const std::string unfoldingHistsPath = rb.get_unfolding_hists_path();
  const int nCentralityBins = rb.get_number_centrality_bins();
  const int nPtBins = rb.get_nbins();

  if (nCentralityBins <= 0 || nPtBins <= 0)
  {
    std::cerr << "Invalid centrality or pT binning in " << configfile << std::endl;
    return;
  }

  std::vector<float> centralityEdges(nCentralityBins + 1);
  std::vector<float> ptEdges(nPtBins + 1);
  rb.get_centrality_bins(centralityEdges.data());
  rb.get_pt_bins(ptEdges.data());

  const double jetRadius = cone_size * 0.1;
  const double jetArea = TMath::Pi() * jetRadius * jetRadius;
  const double totalArea = 2.0 * TMath::Pi() * 1.6;
  const double awaySideArea = 2.0 * (TMath::Pi() - rb.get_dphicut()) * 1.6;
  const double area1 = (totalArea - awaySideArea - jetArea) / totalArea;
  const double area2 = (totalArea - jetArea) / totalArea;

  if (area1 < 0.0 || area1 > 1.0 || area2 < 0.0 || area2 > 1.0)
  {
    std::cerr << "Invalid geometric acceptance factors: A1=" << area1
              << ", A2=" << area2 << std::endl;
    return;
  }

  const float vertexCut = static_cast<float>(rb.get_vtx_cut());
  const float etaCut = rb.get_abs_eta_acceptance(jetRadius);
  const float recoPtCut = static_cast<float>(rb.get_reco_pt_min_cut());

  auto* input = TFile::Open(infile.c_str(), "READ");
  if (!input || input->IsZombie())
  {
    std::cerr << "Could not open " << infile << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(input->Get("T"));
  if (!tree)
  {
    std::cerr << "Could not find tree T in " << infile << std::endl;
    input->Close();
    return;
  }

  tree->SetBranchStatus("*", 0);
  for (const char* branch : {"cent", "zvrtx", "jet_pT", "jet_unsub_pT",
                             "jet_E", "jet_unsub_E", "jet_eta"})
  {
    tree->SetBranchStatus(branch, 1);
  }

  int centrality = -1;
  float vertexZ = 0.0;
  std::vector<float>* jetPt = nullptr;
  std::vector<float>* jetUnsubPt = nullptr;
  std::vector<float>* jetEnergy = nullptr;
  std::vector<float>* jetUnsubEnergy = nullptr;
  std::vector<float>* jetEta = nullptr;
  tree->SetBranchAddress("cent", &centrality);
  tree->SetBranchAddress("zvrtx", &vertexZ);
  tree->SetBranchAddress("jet_pT", &jetPt);
  tree->SetBranchAddress("jet_unsub_pT", &jetUnsubPt);
  tree->SetBranchAddress("jet_E", &jetEnergy);
  tree->SetBranchAddress("jet_unsub_E", &jetUnsubEnergy);
  tree->SetBranchAddress("jet_eta", &jetEta);

  constexpr int kFineCentralityBins = 100;
  std::vector<double> eventsByCentrality(kFineCentralityBins, 0.0);
  std::vector<TH1D*> jetCountsByCentrality(kFineCentralityBins, nullptr);
  for (int cent = 0; cent < kFineCentralityBins; ++cent)
  {
    jetCountsByCentrality[cent] = new TH1D(
      Form("h_jet_spectra_etacut_%d", cent),
      ";#it{p}_{T} [GeV];Jets", 45, 5.0, 50.0);
    jetCountsByCentrality[cent]->Sumw2();
  }

  const Long64_t entries = tree->GetEntries();
  const Long64_t progressStep = std::max<Long64_t>(1, entries / 10);
  TF1 backgroundCut("fcut", "[0]+[1]*TMath::Exp(-[2]*x)", 0.0, 100.0);
  backgroundCut.SetParameters(2.5, 36.2, 0.035);

  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    tree->GetEntry(entry);
    if (entry > 0 && entry % progressStep == 0)
    {
      std::cout << "Event " << entry << " / " << entries << std::endl;
    }

    if (std::abs(vertexZ) > vertexCut || centrality < 0 ||
        centrality >= kFineCentralityBins)
    {
      continue;
    }

    // The denominator is the accepted minimum-bias event count, as in the
    // v0 reweighting procedure and the direct example in genProbs.C.
    eventsByCentrality[centrality] += 1.0;

    if (!jetPt || !jetUnsubPt || !jetEnergy || !jetUnsubEnergy || !jetEta ||
        jetPt->size() != jetUnsubPt->size() ||
        jetPt->size() != jetEnergy->size() ||
        jetPt->size() != jetUnsubEnergy->size() ||
        jetPt->size() != jetEta->size())
    {
      continue;
    }

    const float cutValue = backgroundCut.Eval(centrality);
    for (std::size_t jet = 0; jet < jetPt->size(); ++jet)
    {
      const float pt = jetPt->at(jet);
      if (!std::isfinite(pt) || !std::isfinite(jetEta->at(jet)) ||
          pt < recoPtCut || jetEnergy->at(jet) < 0.0 ||
          jetUnsubEnergy->at(jet) < 0.0 || std::abs(jetEta->at(jet)) > etaCut ||
          jetUnsubPt->at(jet) - pt > cutValue)
      {
        continue;
      }
      // Do not upper-cut pT: jets over 50 GeV are retained in the overflow bin.
      jetCountsByCentrality[centrality]->Fill(pt);
    }
  }

  std::vector<TH1D*> probabilityFine(kFineCentralityBins, nullptr);
  for (int cent = 0; cent < kFineCentralityBins; ++cent)
  {
    probabilityFine[cent] = new TH1D(
      Form("h_pt2_correction_%d", cent),
      ";#it{p}_{T} [GeV];Subleading Efficiency", 45, 5.0, 50.0);
    fillProbabilityHistogram(probabilityFine[cent], jetCountsByCentrality[cent],
                             eventsByCentrality[cent], area1, area2);
  }

  std::vector<TH1D*> jetCounts(nCentralityBins, nullptr);
  std::vector<TH1D*> jetSpectra(nCentralityBins, nullptr);
  std::vector<TH1D*> probability(nCentralityBins, nullptr);
  std::vector<TH1D*> probabilityLog(nCentralityBins, nullptr);
  std::vector<double> events(nCentralityBins, 0.0);

  for (int cent = 0; cent < nCentralityBins; ++cent)
  {
    jetCounts[cent] = new TH1D(Form("h_jet_counts_%d", cent),
      ";#it{p}_{T} [GeV];Jets", 45, 5.0, 50.0);
    jetCounts[cent]->Sumw2();
    probability[cent] = new TH1D(Form("h_pt2_bin_correction_%d", cent),
      ";#it{p}_{T} [GeV];Subleading Efficiency", 45, 5.0, 50.0);
    probabilityLog[cent] = new TH1D(Form("h_pt2_bin_log_correction_%d", cent),
      ";#it{p}_{T} [GeV];Subleading Efficiency", nPtBins, ptEdges.data());

    for (int fineCent = 0; fineCent < kFineCentralityBins; ++fineCent)
    {
      if (fineCent >= centralityEdges[cent] && fineCent < centralityEdges[cent + 1])
      {
        events[cent] += eventsByCentrality[fineCent];
        jetCounts[cent]->Add(jetCountsByCentrality[fineCent]);
      }
    }

    fillProbabilityHistogram(probability[cent], jetCounts[cent], events[cent], area1, area2);
    fillProbabilityHistogram(probabilityLog[cent], jetCounts[cent], events[cent], area1, area2);
    jetSpectra[cent] = dynamic_cast<TH1D*>(jetCounts[cent]->Clone(
      Form("h_jet_spectra_meas_%d", cent)));
    if (events[cent] > 0.0)
    {
      jetSpectra[cent]->Scale(1.0 / events[cent], "width");
    }
  }

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  // Show both probability binnings.  The third panel is intentionally drawn
  // from h_pt2_bin_log_correction rather than a rebinned display histogram,
  // so it is a direct check of the correction consumed by the unfolding.
  auto* canvas = new TCanvas("c", "c", 1500, 450);
  canvas->Divide(3, 1);
  const int colors[] = {kRed, kBlue, kGreen + 2, kOrange + 7};
  const int nDraw = std::min(4, nCentralityBins);
  auto* legend = new TLegend(0.62, 0.62, 0.88, 0.88);
  legend->SetLineWidth(0);
  legend->SetTextSize(0.04);

  canvas->cd(1);
  gPad->SetLogy();
  for (int cent = 0; cent < nDraw; ++cent)
  {
    dlutility::SetLineAtt(jetSpectra[cent], colors[cent], 2, 1);
    jetSpectra[cent]->SetTitle(";#it{p}_{T} [GeV];#frac{1}{N_{evt}}#frac{dN_{jet}}{d#it{p}_{T}}");
    jetSpectra[cent]->GetXaxis()->SetRangeUser(5, 30);
    jetSpectra[cent]->Draw(cent == 0 ? "hist" : "hist same");
    legend->AddEntry(jetSpectra[cent], Form("%.0f-%.0f%%", centralityEdges[cent], centralityEdges[cent + 1]), "l");
  }
  legend->Draw();

  canvas->cd(2);
  for (int cent = 0; cent < nDraw; ++cent)
  {
    dlutility::SetLineAtt(probability[cent], colors[cent], 2, 1);
    probability[cent]->SetMinimum(0.0);
    probability[cent]->SetMaximum(1.05);
    probability[cent]->GetXaxis()->SetRangeUser(5, 30);
    probability[cent]->Draw(cent == 0 ? "hist" : "hist same");
  }
  auto* unity = new TLine(5, 1, 30, 1);
  unity->SetLineStyle(4);
  unity->Draw();
  legend->Draw();

  canvas->cd(3);
  auto* logLegend = new TLegend(0.62, 0.62, 0.88, 0.88);
  logLegend->SetLineWidth(0);
  logLegend->SetTextSize(0.04);
  for (int cent = 0; cent < nDraw; ++cent)
  {
    dlutility::SetLineAtt(probabilityLog[cent], colors[cent], 2, 1);
    probabilityLog[cent]->SetMinimum(0.0);
    probabilityLog[cent]->SetMaximum(1.05);
    probabilityLog[cent]->SetTitle(";#it{p}_{T} [GeV];Subleading Efficiency");
    probabilityLog[cent]->GetXaxis()->SetRangeUser(5, 30);
    probabilityLog[cent]->Draw(cent == 0 ? "hist" : "hist same");
    logLegend->AddEntry(probabilityLog[cent],
      Form("%.0f-%.0f%%", centralityEdges[cent], centralityEdges[cent + 1]), "l");
  }
  auto* logUnity = new TLine(5, 1, 30, 1);
  logUnity->SetLineStyle(4);
  logUnity->Draw();
  logLegend->Draw();
  canvas->Print(Form("probabilities_AA_r%02d.pdf", cone_size));

  const std::string output = unfoldingHistsPath +
    "/probability_hists_AA_r0" + std::to_string(cone_size) + ".root";
  auto* outputFile = TFile::Open(output.c_str(), "RECREATE");
  if (!outputFile || outputFile->IsZombie())
  {
    std::cerr << "Could not create " << output << std::endl;
    input->Close();
    return;
  }

  for (int cent = 0; cent < kFineCentralityBins; ++cent)
  {
    probabilityFine[cent]->Write();
  }
  for (int cent = 0; cent < nCentralityBins; ++cent)
  {
    probability[cent]->Write();
    probabilityLog[cent]->Write();
    jetSpectra[cent]->Write();
  }
  outputFile->Close();
  input->Close();
}
