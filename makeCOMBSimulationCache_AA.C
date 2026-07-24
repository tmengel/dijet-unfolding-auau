#include <array>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"
#include "TString.h"
#include "TSystem.h"

// Reuse the validated simulation selection, event weights, fitting, and cache
// writer from the full diagnostic macro.  This file only provides a small,
// simulation-only entry point.
#include "drawCOMBModulation_AA_v2.C"

namespace
{
bool validCOMBSimulationCache(const TString &path)
{
  std::unique_ptr<TFile> input(TFile::Open(path, "READ"));
  if (!input || input->IsZombie()) return false;

  const std::array<const char *, 17> required = {
    "h_dphi_pairs",
    "h_dphi_eta_separated",
    "h_dphi_eta_separated_fit",
    "h_flow_background_nominal",
    "h_flow_background_COMBDown",
    "h_flow_background_COMBUp",
    "h_dphi_nominal_subtracted_sim",
    "h_dphi_down_subtracted_sim",
    "h_dphi_up_subtracted_sim",
    "h_signal_region",
    "h_flow_normalization_region_sim",
    "h_dphi_pairs_truth",
    "h_dphi_eta_separated_truth",
    "f_flow_fit",
    "f_flow_fit_0_sim",
    "f_flow_fit_1_sim",
    "f_flow_fit_2_sim"
  };
  for (const char *name : required)
    if (!input->Get(name)) return false;
  return true;
}
} // namespace

// centrality_bin = -1 creates/reuses caches for all four centralities.
// Set overwrite = true only when the simulation input, selection, or event
// reweighting files have changed.  Normal plot redraws should leave it false.
void makeCOMBSimulationCache_AA(
  const int cone_size = 3, const int centrality_bin = -1,
  const std::string config = "binning_AA.config",
  const bool overwrite = false)
{
  if (centrality_bin < -1 || centrality_bin > 3)
    {
      std::cerr << "centrality_bin must be -1 or 0 through 3" << std::endl;
      return;
    }

  read_binning rb(config);
  std::unique_ptr<float[]> ptBins(new float[rb.get_nbins() + 1]);
  rb.get_pt_bins(ptBins.get());
  const double leadingCut = rb.get_reco_leading_cut();
  const double subleadingCut = rb.get_reco_subleading_cut();
  float centralityBins[5] = {0};
  rb.get_centrality_bins(centralityBins);

  const TString plotDirectory = Form(
    "%s/dphi_plots", rb.get_code_location().c_str());
  gSystem->mkdir(plotDirectory, true);

  const int firstCentrality = centrality_bin < 0 ? 0 : centrality_bin;
  const int lastCentrality = centrality_bin < 0 ? 3 : centrality_bin;
  for (int centrality = firstCentrality;
       centrality <= lastCentrality; ++centrality)
    {
      const TString cache = Form(
        "%s/dphi_COMB_modulation_sim_combined_AA_cent_%d_r%02d.root",
        plotDirectory.Data(), centrality, cone_size);
      if (!overwrite && validCOMBSimulationCache(cache))
        {
          std::cout << "Reusing simulation cache " << cache << std::endl;
          continue;
        }

      std::cout << "Producing simulation cache for centrality "
                << centrality << std::endl;
      drawSimulationDphi(
        "", "combined", cone_size, centrality,
        centralityBins[centrality], centralityBins[centrality + 1],
        leadingCut, subleadingCut, plotDirectory, rb, 0);

      if (!validCOMBSimulationCache(cache))
        {
          std::cerr << "Simulation cache validation failed: " << cache
                    << std::endl;
          return;
        }
      std::cout << "Saved and validated " << cache << std::endl;
    }
}
