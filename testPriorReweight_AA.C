// Standalone verification of the prior reweighting used by
// createResponse_noempty_AA.cxx.
//
//   root -l -b -q 'testPriorReweight_AA.C(3, 0)'
//   root -l -b -q 'testPriorReweight_AA.C(3, 0, "nominal", 1, "0,0.25,0.5,1.0")'
//
// It reads the PRIMER1 unfolded data and the PRIMER1 MC truth, builds the prior
// weights with prior_qa::build()/build_xj() and applies them with
// prior_qa::apply() -- the same code the event loop uses -- for a list of prior
// fractions and for both blend spaces (per (pt1,pt2) bin, and in xJ with the
// ratio mapped back onto the pt1pt2 bins).  Applying the
// weights bin-by-bin to the flat truth is exactly what the event loop does
// (every event in truth bin k is scaled by w(k)), so no ntuple pass is needed.
//
// Expected result when the implementation is correct:
//   kBlendPerBin : reweighted/target = 1 in every FLAT bin, for every fraction
//   kBlendXJ     : reweighted/target = 1 in every xJ bin; the flat-space ratio
//                  is deliberately not 1, since the xJ blend does not chase
//                  per-bin structure
// Any other structure in those ratios is a bug in the reweighting, not physics.
//
// The macro also re-runs the whole thing with bin_offset = -1, reproducing the
// historical GetBinContent(k) lookup, so the size of that off-by-one is visible
// side by side with the fixed lookup.

#include <string>
#include <sstream>
#include <vector>

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
#include "priorReweightQA.h"
#include "PlotUtils.h"

namespace
{
  std::vector<double> parse_fractions(const std::string &csv)
  {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ','))
    {
      if (item.empty()) continue;
      out.push_back(std::atof(item.c_str()));
    }
    if (out.empty()) out.push_back(0.5);
    return out;
  }
}

void testPriorReweight_AA(const int cone_size = 3,
                          const int centrality_bin = 0,
                          const std::string sys_name = "nominal",
                          const int prior_iteration = 1,
                          const std::string fraction_list = "0.0,0.25,0.5,1.0",
                          // -1 = both blend spaces, 0 = per (pt1,pt2) bin only,
                          // 1 = blend in xJ and map the ratio back onto pt1pt2.
                          const int blend_space = -1,
                          const bool also_broken_lookup = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  PlotUtils::set_sphenix_style();

  read_binning rb("/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config");

  const int nbins = rb.get_nbins();
  float ipt_bins[nbins + 1];
  float ixj_bins[nbins + 1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  const int mbins = rb.get_measure_bins();
  int measure_bins[10] = {0};
  for (int ir = 0; ir < mbins + 1; ir++) { measure_bins[ir] = rb.get_measure_region(ir); }
  const int measure_subleading_bin = rb.get_measure_subleading_bin();

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins + 1];
  rb.get_centrality_bins(icentrality_bins);

  const std::string system_string = rb.get_system_string(centrality_bin);
  const std::string cent_label = Form("%d - %d %%", (int) icentrality_bins[centrality_bin],
                                                    (int) icentrality_bins[centrality_bin + 1]);

  // ---- inputs: exactly the two histograms createResponse reads -------------
  const std::string unfold_path = Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER2_%s.root",
                                       system_string.c_str(), cone_size, sys_name.c_str());
  const std::string resp_path = Form("response_matrices/response_matrix_%s_r%02d_PRIMER2_%s.root",
                                     system_string.c_str(), cone_size, sys_name.c_str());

  TFile *fun = TFile::Open(unfold_path.c_str(), "READ");
  if (!fun || fun->IsZombie()) { std::cerr << "Missing " << unfold_path << std::endl; return; }
  TH1D *h_unfold_flat = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
  if (!h_unfold_flat) { std::cerr << "Missing h_flat_unfold_pt1pt2_" << prior_iteration << std::endl; return; }
  h_unfold_flat = (TH1D*) h_unfold_flat->Clone("h_unfold_flat_in");
  h_unfold_flat->SetDirectory(nullptr);
  fun->Close();

  TFile *ftr = TFile::Open(resp_path.c_str(), "READ");
  if (!ftr || ftr->IsZombie()) { std::cerr << "Missing " << resp_path << std::endl; return; }
  TH1D *h_truth_flat = (TH1D*) ftr->Get("h_truth_flat_pt1pt2");
  if (!h_truth_flat) { std::cerr << "Missing h_truth_flat_pt1pt2" << std::endl; return; }
  h_truth_flat = (TH1D*) h_truth_flat->Clone("h_truth_flat_in");
  h_truth_flat->SetDirectory(nullptr);
  ftr->Close();

  std::cout << "=== prior reweighting test: " << system_string << " r" << cone_size
            << " " << sys_name << ", iteration index " << prior_iteration << " ===" << std::endl;

  const std::vector<double> fractions = parse_fractions(fraction_list);

  std::vector<int> spaces;
  if (blend_space < 0) { spaces = {prior_qa::kBlendPerBin, prior_qa::kBlendXJ}; }
  else                 { spaces = {blend_space}; }

  for (const int space : spaces)
  for (const double fraction : fractions)
  {
    const std::string tag = Form("_%s_f%03d", prior_qa::blend_name(space),
                                 (int) std::lround(100*fraction));
    std::vector<prior_qa::lead_group> lead_groups;
    for (int irange = 0; irange < mbins; irange++)
    {
      lead_groups.push_back({measure_bins[irange], measure_bins[irange + 1]});
    }

    prior_qa::weights w = (space == prior_qa::kBlendXJ)
      ? prior_qa::build_xj(h_truth_flat, h_unfold_flat, fraction, nbins,
                           ipt_bins, ixj_bins, lead_groups,
                           measure_subleading_bin, nbins - 2, tag)
      : prior_qa::build(h_truth_flat, h_unfold_flat, fraction, nbins, tag);
    if (!w.weight) { std::cerr << "weight build failed for fraction " << fraction << std::endl; continue; }

    const int nmodes = also_broken_lookup ? 2 : 1;
    for (int mode = 0; mode < nmodes; mode++)
    {
      const int bin_offset = (mode == 0) ? 0 : -1;
      const std::string mode_name = (mode == 0) ? "fixed" : "offbyone";
      const std::string mode_text = (mode == 0) ? "lookup: GetBinContent(k+1)"
                                                : "lookup: GetBinContent(k) (off-by-one)";

      TH1D *reweighted = prior_qa::apply(w.truth, w.weight, nbins,
                                         "h_prior_reweighted" + tag + "_" + mode_name, bin_offset);

      const std::string label = Form("blend=%s f=%.2f %s", prior_qa::blend_name(space),
                                     fraction, mode_name.c_str());
      prior_qa::print_summary(w, reweighted, nbins, label);
      for (int irange = 0; irange < mbins; irange++)
      {
        prior_qa::print_xj_summary(w, reweighted, nbins, ipt_bins, ixj_bins,
                                   measure_bins[irange], measure_bins[irange + 1],
                                   measure_subleading_bin, nbins - 2,
                                   Form("%s range%d", label.c_str(), irange));
      }

      std::vector<std::string> captions = {
        Form("Prior fraction %.2f", fraction),
        Form("blended in %s", space == prior_qa::kBlendXJ
                              ? "#it{x}_{J}" : "(#it{p}_{T,1}, #it{p}_{T,2})"),
        cent_label
      };
      if (also_broken_lookup) { captions.push_back(mode_text); }

      prior_qa::draw_flat(w, reweighted, captions,
        Form("prior_qa_plots/priorQA_flat_%s_r%02d_%s_blend%s_f%03d_%s.pdf",
             system_string.c_str(), cone_size, sys_name.c_str(),
             prior_qa::blend_name(space),
             (int) std::lround(100*fraction), mode_name.c_str()));

      for (int irange = 1; irange < 2; irange++)
      {
        std::vector<std::string> xj_captions = captions;
        xj_captions.push_back(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV",
                                   ipt_bins[measure_bins[irange]],
                                   ipt_bins[measure_bins[irange + 1]]));
        prior_qa::draw_xj(w, reweighted, nbins, ipt_bins, ixj_bins,
                          measure_bins[irange], measure_bins[irange + 1],
                          measure_subleading_bin, nbins - 2, xj_captions,
          Form("prior_qa_plots/priorQA_xj_%s_r%02d_%s_blend%s_f%03d_%s_range%d.pdf",
               system_string.c_str(), cone_size, sys_name.c_str(),
               prior_qa::blend_name(space),
               (int) std::lround(100*fraction), mode_name.c_str(), irange));
      }

      // Mean xJ before/after, the number the systematic is ultimately built on.
      for (int irange = 1; irange < 2; irange++)
      {
        auto project = [&](TH1D *flat, const char *tag2)
        {
          TH2D *h2 = new TH2D(Form("h2_mean_%s_%s%s_%d", tag2, mode_name.c_str(),
                                   tag.c_str(), irange),
                              "", nbins, ipt_bins, nbins, ipt_bins);
          h2->SetDirectory(nullptr);
          histo_opps::make_sym_pt1pt2(flat, h2, nbins);
          TH1D *hxj = new TH1D(Form("hxj_mean_%s_%s%s_%d", tag2, mode_name.c_str(),
                                    tag.c_str(), irange),
                               "", nbins, ixj_bins);
          hxj->SetDirectory(nullptr);
          histo_opps::project_xj(h2, hxj, nbins, measure_bins[irange], measure_bins[irange + 1],
                                 measure_subleading_bin, nbins - 2);
          histo_opps::normalize_histo(hxj, nbins);
          delete h2;
          return hxj;
        };
        TH1D *xj_target = project(w.target, "target");
        TH1D *xj_truth  = project(w.truth,  "truth");
        TH1D *xj_rw     = project(reweighted, "rw");
        printf("   range %d (%4.1f - %4.1f GeV)  <xJ>: truth %.4f  target %.4f  reweighted %.4f\n",
               irange, ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange + 1]],
               histo_opps::get_average_xj(xj_truth),
               histo_opps::get_average_xj(xj_target),
               histo_opps::get_average_xj(xj_rw));
        delete xj_target; delete xj_truth; delete xj_rw;
      }
    }
  }

  std::cout << "=== plots written to prior_qa_plots/priorQA_* ===" << std::endl;
}
