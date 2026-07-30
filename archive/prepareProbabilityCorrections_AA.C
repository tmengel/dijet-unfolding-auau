#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"

#include "read_binning.h"

namespace
{
double overlapIntegral(const TH1D& source, const double low, const double high)
{
  double sum = 0.0;
  for (int bin = 1; bin <= source.GetNbinsX(); ++bin)
    {
      const double sourceLow = source.GetXaxis()->GetBinLowEdge(bin);
      const double sourceHigh = source.GetXaxis()->GetBinUpEdge(bin);
      const double overlap = std::max(
        0.0, std::min(high, sourceHigh) - std::max(low, sourceLow));
      const double width = sourceHigh - sourceLow;
      if (width > 0.0) sum += source.GetBinContent(bin)*overlap/width;
    }
  return sum;
}

double overlapWeightedProbability(const TH1D& probability,
                                  const TH1D& counts,
                                  const double low,
                                  const double high)
{
  double sum = 0.0;
  for (int bin = 1; bin <= counts.GetNbinsX(); ++bin)
    {
      const double sourceLow = counts.GetXaxis()->GetBinLowEdge(bin);
      const double sourceHigh = counts.GetXaxis()->GetBinUpEdge(bin);
      const double overlap = std::max(
        0.0, std::min(high, sourceHigh) - std::max(low, sourceLow));
      const double width = sourceHigh - sourceLow;
      if (width > 0.0)
        sum += counts.GetBinContent(bin)*probability.GetBinContent(bin)
          *overlap/width;
    }
  return sum;
}
}

void prepareProbabilityCorrections_AA(
  const int cone_size = 3,
  const std::string input_path = "/home/tmengel/PPG14/rootfiles/probs.root",
  const std::string configfile = "binning_AA.config")
{
  read_binning rb(configfile);
  const int nbins = rb.get_nbins();
  std::vector<float> float_edges(static_cast<std::size_t>(nbins + 1));
  rb.get_pt_bins(float_edges.data());
  std::vector<double> edges(float_edges.begin(), float_edges.end());

  std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
  if (!input || input->IsZombie())
    {
      std::cerr << "Cannot open probability input " << input_path << std::endl;
      return;
    }

  const std::string output_path = Form(
    "%s/unfolding_hists/probability_hists_AA_r0%d.root",
    rb.get_code_location().c_str(), cone_size);
  std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "RECREATE"));
  if (!output || output->IsZombie())
    {
      std::cerr << "Cannot create " << output_path << std::endl;
      return;
    }

  for (int centrality = 0; centrality < 4; ++centrality)
    {
      auto* counts = dynamic_cast<TH1D*>(
        input->Get(Form("h1_pt_counts_cent%d", centrality)));
      auto* reconstructed = dynamic_cast<TH1D*>(
        input->Get(Form("h1_pt_reco_cent%d", centrality)));
      auto* supplied_probability = dynamic_cast<TH1D*>(
        input->Get(Form("h1_prob_cent%d", centrality)));
      if (!counts || !reconstructed || !supplied_probability)
        {
          std::cerr << "Missing probability inputs for centrality "
                    << centrality << std::endl;
          return;
        }

      TH1D correction(
        Form("h_pt2_bin_log_correction_%d", centrality),
        ";p_{T,2} [GeV];reconstruction probability",
        nbins, edges.data());
      correction.SetDirectory(nullptr);
      for (int bin = 1; bin <= nbins; ++bin)
        {
          const double denominator = overlapIntegral(
            *counts, edges[static_cast<std::size_t>(bin - 1)],
            edges[static_cast<std::size_t>(bin)]);
          // h1_pt_reco is a normalized spectrum, not reconstructed raw
          // counts. Preserve the supplied probability definition by taking
          // its count-weighted average in each analysis pT interval.
          const double numerator = overlapWeightedProbability(
            *supplied_probability, *counts,
            edges[static_cast<std::size_t>(bin - 1)],
            edges[static_cast<std::size_t>(bin)]);
          const double probability = denominator > 0.0
            ? numerator/denominator
            : supplied_probability->GetBinContent(
                supplied_probability->GetNbinsX());
          correction.SetBinContent(bin, probability);
          correction.SetBinError(bin, denominator > 0.0
            ? std::sqrt(std::max(0.0, probability*(1.0 - probability)/denominator))
            : 0.0);
          std::cout << "cent " << centrality << " bin " << bin << " ["
                    << edges[static_cast<std::size_t>(bin - 1)] << ", "
                    << edges[static_cast<std::size_t>(bin)] << "): "
                    << numerator << "/" << denominator << " = "
                    << probability << std::endl;
        }

      output->cd();
      correction.Write();
      auto* counts_copy = static_cast<TH1D*>(counts->Clone(
        Form("h_probability_source_counts_cent%d", centrality)));
      auto* reco_copy = static_cast<TH1D*>(reconstructed->Clone(
        Form("h_probability_source_reco_cent%d", centrality)));
      auto* probability_copy = static_cast<TH1D*>(supplied_probability->Clone(
        Form("h_probability_source_ratio_cent%d", centrality)));
      counts_copy->Write();
      reco_copy->Write();
      probability_copy->Write();

      TCanvas canvas(Form("c_probability_cent%d", centrality), "", 650, 550);
      supplied_probability->SetLineColor(kGray + 2);
      supplied_probability->SetMarkerColor(kGray + 2);
      supplied_probability->SetMarkerStyle(24);
      supplied_probability->SetMinimum(0.0);
      supplied_probability->SetMaximum(1.1);
      supplied_probability->Draw("E1");
      correction.SetLineColor(kRed + 1);
      correction.SetMarkerColor(kRed + 1);
      correction.SetMarkerStyle(20);
      correction.Draw("E1 SAME");
      TLegend legend(0.52, 0.20, 0.88, 0.36);
      legend.SetBorderSize(0);
      legend.AddEntry(supplied_probability, "Supplied 1 GeV bins", "lep");
      legend.AddEntry(&correction, "Analysis p_{T} bins", "lep");
      legend.Draw();
      canvas.Print(Form("%s/unfolding_plots/probability_AA_cent_%d_r0%d.pdf",
                        rb.get_code_location().c_str(), centrality, cone_size));
      canvas.Print(Form("%s/unfolding_plots/probability_AA_cent_%d_r0%d.png",
                        rb.get_code_location().c_str(), centrality, cone_size));
      delete counts_copy;
      delete reco_copy;
      delete probability_copy;
    }
  output->Close();
  std::cout << "Wrote " << output_path << std::endl;
}
