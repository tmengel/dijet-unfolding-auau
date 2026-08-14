#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TString.h"

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "histo_opps.h"
#include "read_binning.h"


int createResponse(
  const std::string & configfile = "binning.config",
  const int full_or_half = 0,
  const int niterations = 10,
  const int cone_size = 4,
  const int centrality_bin = 0,
  const int primer = 0
)
{
  if (full_or_half < 0 || full_or_half > 2)
  {
    std::cerr << "Closure mode must be 0 (analysis), 1 (half closure), or 2 (full closure)" << std::endl;
    return 1;
  }

  const bool half_closure = (full_or_half == 1);
  const bool full_closure = (full_or_half == 2);
  const bool full_or_half_closure = (full_or_half != 0);

  const bool ispp = (centrality_bin < 0);
  read_binning rb(configfile.c_str());

  const std::string system_string = rb.get_system_string(centrality_bin);

  const std::string j10_file = std::getenv("TNUPLE_SIM_FILE_JET10");
  const std::string j20_file = std::getenv("TNUPLE_SIM_FILE_JET20");
  const std::string j30_file = std::getenv("TNUPLE_SIM_FILE_JET30");
  std::cout << "Using matched simulation files:" << std::endl;
  std::cout << "  Jet10: " << j10_file << std::endl;
  std::cout << "  Jet20: " << j20_file << std::endl;
  std::cout << "  Jet30: " << j30_file << std::endl;

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float dphi_reco[3];
  float match[3];
  float mbd_vertex[3];
  float centrality[3];
  float sumeT[3] = {0, 0, 0};

  float n_events[3];
  float b_n_events = 0;
  float event_weight[3] = {1, 1, 1};
  bool has_event_weight[3] = {false, false, false};
  bool has_sumeT[3] = {false, false, false};

  TFile *fin[3];
  fin[0] = new TFile(j10_file.c_str(), "r");
  fin[1] = new TFile(j20_file.c_str(), "r");
  fin[2] = new TFile(j30_file.c_str(), "r");
  TNtuple *tn[3];
  for (int i = 0; i < 3; i++)
  {
    if (!fin[i] || fin[i]->IsZombie())
    {
      std::cerr << "Cannot open matched simulation file " << fin[i]->GetName() << std::endl;
      return 1;
    }

    tn[i] = (TNtuple*) fin[i]->Get("tn_match");
    if (!tn[i])
    {
      std::cerr << "Missing tn_match in " << fin[i]->GetName() << std::endl;
      return 1;
    }

    tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
    tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
    tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
    tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
    tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
    tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
    tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
    tn[i]->SetBranchAddress("matched", &match[i]);
    tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
    tn[i]->SetBranchAddress("centrality", &centrality[i]);
    has_sumeT[i] = tn[i]->GetBranch("sumeT") != nullptr;
    if (has_sumeT[i]) { tn[i]->SetBranchAddress("sumeT", &sumeT[i]); }
    else if (!ispp) { std::cerr << "Warning: missing sumeT in " << fin[i]->GetName() << "; sumeT reweighting is disabled for this sample." << std::endl; }

    has_event_weight[i] = tn[i]->GetBranch("weight") != nullptr;
    if (has_event_weight[i]) { tn[i]->SetBranchAddress("weight", &event_weight[i]); }

    TNtuple *tn_stats = (TNtuple*) fin[i]->Get("tn_stats");
    if (!tn_stats)
    {
      std::cerr << "Missing tn_stats in " << fin[i]->GetName() << std::endl;
      return 1;
    }
    tn_stats->SetBranchAddress("nevents", &b_n_events);
    tn_stats->GetEntry(0);
    n_events[i] = b_n_events;
  }

  // pt-hat cross sections used to combine the jet10/jet20/jet30 samples into
  // one luminosity-normalized spectrum
  const float cs_10 = 0.000003997;
  const float cs_20 = 6.218e-8;
  const float cs_30 = 2.505e-9;

  float scale_factor[3];
  scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
  scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30;
  scale_factor[2] = 1;

  const int minentries = rb.get_minentries();
  const int nbins = rb.get_nbins();

  const double dphicut = rb.get_dphicut();
  const double dphicuttruth = dphicut;

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);

  const int zyam_sys = rb.get_zyam_sys();
  const double flow_v22_scale = rb.get_flow_sys();
  const double flow_v33_scale = rb.get_flow_v33_sys();
  const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
  const int inclusive_sys = rb.get_inclusive_sys();
  const double JES_sys = rb.get_jes_sys();
  const double JER_sys = rb.get_jer_sys();
  const int prior_sys = rb.get_prior_sys();

  TF1 *f_smear = nullptr;
  if (!ispp)
  {
    f_smear = (TF1*) rb.get_smear_function(centrality_bin);
    if (!f_smear)
    {
      std::cout << " NO SMEAR " << std::endl;
      exit(-1);
    }
  }

  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  const float width = 0.8 + JER_sys;
  fgaus->SetParameters(1, 0, width);

  std::string sys_name = "nominal";
  if (prior_sys)      { sys_name = "PRIOR"; }
  if (zyam_sys)        { sys_name = "ZYAM"; }
  if (flow_sys)         { sys_name = rb.get_flow_systematic_name(); }
  if (inclusive_sys)    { sys_name = "INCLUSIVE"; }
  if (JER_sys != 0)
  {
    sys_name = (JER_sys < 0) ? "negJER" : "posJER";
    std::cout << "Calculating JER extra = " << JER_sys << std::endl;
  }
  if (JES_sys != 0)
  {
    sys_name = (JES_sys < 0) ? "negJES" : "posJES";
    std::cout << "Calculating JES extra = " << JES_sys << std::endl;
  }

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0; i < nbins+1; i++)
  {
    std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
  }

  const int max_reco_bin = rb.get_maximum_reco_bin();

  const int centrality_bins = rb.get_number_centrality_bins();
  float icentrality_bins[centrality_bins+1];
  rb.get_centrality_bins(icentrality_bins);

  const int prior_iteration = 1;

  // Vertex / sumeT / centrality reweighting: each factor file holds a 1D
  // histogram whose bin content is the reweight factor for events below that
  // bin's upper edge; load once into (upper_edge, factor) pairs.
  auto load_reweight = [](const std::string & path, const char * histname, std::vector<std::pair<float, float>> & out) -> bool
  {
    TFile f(path.c_str(), "READ");
    TH1D * h = (!f.IsZombie()) ? (TH1D*) f.Get(histname) : nullptr;
    if (!h) { return false; }
    for (int ib = 0; ib < h->GetNbinsX(); ib++)
    {
      out.push_back(std::make_pair(h->GetBinLowEdge(ib+1) + h->GetBinWidth(ib+1), h->GetBinContent(ib+1)));
    }
    return true;
  };

  std::vector<std::pair<float, float>> centrality_scales;
  std::vector<std::pair<float, float>> vertex_scales;
  std::vector<std::pair<float, float>> sumeT_scales;

  if (primer != 1 && !full_or_half)
  {
    if (!load_reweight(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "h_mbd_reweight", vertex_scales))
    {
      std::cerr << "Missing required vertex reweight for " << system_string << " " << sys_name << std::endl;
      return 1;
    }

    if (!ispp)
    {
      if (!load_reweight(Form("centrality/centrality_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "h_centrality_reweight", centrality_scales))
      {
        std::cerr << "Missing required centrality reweight for " << system_string << " " << sys_name << std::endl;
        return 1;
      }
      if (!load_reweight(Form("sumeT/sumeT_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "h_sumeT_reweight", sumeT_scales))
      {
        std::cerr << "Missing required sumeT reweight file/hist for " << system_string << " r0" << cone_size << " " << sys_name << std::endl;
        return 1;
      }
    }
  }

  const float truth_subleading_cut = rb.get_truth_subleading_cut();
  const float reco_leading_cut = rb.get_reco_leading_cut();
  const float reco_subleading_cut = rb.get_reco_subleading_cut();
  const float measure_leading_cut = rb.get_measure_leading_cut();
  const float measure_subleading_cut = rb.get_measure_subleading_cut();
  const int measure_leading_bin = rb.get_measure_leading_bin();
  const int measure_subleading_bin = rb.get_measure_subleading_bin();

  float sample_boundary[4] = {0};
  for (int ib = 0; ib < 4; ib++)
  {
    sample_boundary[ib] = rb.get_sample_boundary(ib);
    std::cout << sample_boundary[ib] << std::endl;
  }
  std::cout << "Max reco bin: " << max_reco_bin << std::endl;
  std::cout << "Reco 1: " << reco_leading_cut << std::endl;
  std::cout << "Meas 1: " << measure_leading_cut << std::endl;
  std::cout << "Reco 2: " << reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " << measure_subleading_cut << std::endl;

  TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);

  TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2", ";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

  // Prior-reweight bootstrap: a primer==1 pass produces an unfolded truth
  // spectrum with a flat prior; this (primer==0) pass reads that back and
  // reweights truth-level fills toward it, so the response's prior is closer
  // to data. PRIOR systematic = blend only half the events toward it instead
  // of all of them (prior_fraction), same "scale away from nominal" idiom
  // preprocess.C uses for its systematics.
  TH1D *h_flatreweight_pt1pt2 = new TH1D("h_unfold_flat_pt1pt2", ";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  const double prior_fraction = prior_sys ? 0.5 : 1.0;

  if (!full_or_half && !primer)
  {
    std::cout << "Building prior reweight matrix with fraction " << prior_fraction << std::endl;

    TFile *fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
    TFile *ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");

    TH1D *h_unfold_flat = fun ? (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration)) : nullptr;
    TH1D *h_truth_flat = ftr ? (TH1D*) ftr->Get("h_truth_flat_to_response_pt1pt2") : nullptr;

    if (!h_unfold_flat || !h_truth_flat)
    {
      std::cerr << "Missing PRIMER1 inputs for prior reweighting" << std::endl;
      return 1;
    }

    if (h_unfold_flat->Integral() != 0) { h_unfold_flat->Scale(1./h_unfold_flat->Integral()); }
    if (h_truth_flat->Integral() != 0)  { h_truth_flat->Scale(1./h_truth_flat->Integral()); }

    for (int ibin = 0; ibin < nbins*nbins; ibin++)
    {
      const float v = h_unfold_flat->GetBinContent(ibin+1);
      const float b = h_truth_flat->GetBinContent(ibin+1);
      h_flatreweight_pt1pt2->SetBinContent(ibin+1, b > 0 ? v/b : 1);
    }
  }

  const int nbin_response = nbins*nbins;
  RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

  // A fixed seed makes the half-closure train/test assignment reproducible.
  TRandom rng(314159);

  auto get_pt_bin = [&ipt_bins, nbins](float pt) -> int
  {
    for (int ib = 0; ib < nbins; ib++)
    {
      if (pt >= ipt_bins[ib] && pt < ipt_bins[ib+1]) { return ib; }
    }
    return nbins;
  };

  for (int isample = 0; isample < 3; isample++)
  {
    std::cout << "Sample " << isample << std::endl;
    const int entries = tn[isample]->GetEntries();

    for (int i = 0; i < entries; i++)
    {
      tn[isample]->GetEntry(i);
      double event_scale = has_event_weight[isample] ? event_weight[isample] : scale_factor[isample];

      if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) { continue; }
      if (!ispp && (centrality[isample] < icentrality_bins[centrality_bin] || centrality[isample] >= icentrality_bins[centrality_bin+1])) { continue; }

      bool fill_response = false;
      if (full_or_half_closure && rng.Uniform() >= 0.5) { fill_response = true; }

      // fill_response selects which half of a closure sample trains the
      // response; the other half is reserved to test it
      const bool use_for_response = !half_closure || fill_response;
      const bool use_for_test = full_closure || (half_closure && !fill_response);

      if (primer != 1 && !full_or_half)
      {
        for (size_t ib = 0; ib < vertex_scales.size(); ib++)
        {
          if (mbd_vertex[isample] < vertex_scales[ib].first)
          {
            event_scale *= vertex_scales[ib].second;
            break;
          }
        }

        if (!ispp && has_sumeT[isample])
        {
          for (size_t ib = 0; ib < sumeT_scales.size(); ib++)
          {
            if (sumeT[isample] < sumeT_scales[ib].first)
            {
              event_scale *= sumeT_scales[ib].second;
              break;
            }
          }
        }

        if (!ispp)
        {
          for (size_t ib = 0; ib < centrality_scales.size(); ib++)
          {
            if (centrality[isample] < centrality_scales[ib].first)
            {
              event_scale *= centrality_scales[ib].second;
              break;
            }
          }
        }
      }

      float max_truth, min_truth, max_reco, min_reco;
      if (pt1_truth[isample] >= pt2_truth[isample])
      {
        max_truth = pt1_truth[isample]; max_reco = pt1_reco[isample];
        min_truth = pt2_truth[isample]; min_reco = pt2_reco[isample];
      }
      else
      {
        max_truth = pt2_truth[isample]; max_reco = pt2_reco[isample];
        min_truth = pt1_truth[isample]; min_reco = pt1_reco[isample];
      }

      const float e1 = max_truth;
      const float e2 = min_truth;
      float es1 = max_reco;
      float es2 = min_reco;

      double smear1 = 0;
      double smear2 = 0;
      if (ispp && !full_or_half)
      {
        smear1 = fgaus->GetRandom();
        smear2 = fgaus->GetRandom();
      }
      else if (!full_or_half)
      {
        fgaus->SetParameter(2, f_smear->Eval(e1));
        smear1 = fgaus->GetRandom();
        fgaus->SetParameter(2, f_smear->Eval(e2));
        smear2 = fgaus->GetRandom();
      }

      if (JES_sys != 0)
      {
        es1 += (JES_sys + smear1)*e1;
        es2 += (JES_sys + smear2)*e2;
      }
      else
      {
        es1 += smear1*e1;
        es2 += smear2*e2;
      }

      const float maxi = std::max(es1, es2);
      const float mini = std::min(es1, es2);
      const float maxit = std::max(e1, e2);
      const float minit = std::min(e1, e2);

      if (maxi >= ipt_bins[max_reco_bin]) { continue; }

      const int pt1_truth_bin = get_pt_bin(e1);
      const int pt2_truth_bin = get_pt_bin(e2);
      const int pt1_reco_bin = get_pt_bin(es1);
      const int pt2_reco_bin = get_pt_bin(es2);

      const bool truth_good = (maxit >= sample_boundary[1] && minit >= truth_subleading_cut && dphi_truth[isample] >= dphicuttruth);
      const bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] >= dphicut);

      if (!truth_good && !reco_good) { continue; }

      if (!primer && !full_or_half)
      {
        const int recorrectbin = pt1_truth_bin*nbins + pt2_truth_bin;
        if (rng.Uniform() <= prior_fraction)
        {
          event_scale *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
        }
      }

      if (truth_good && !reco_good)
      {
        if (use_for_response)
        {
          h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
          h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
          h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);

          rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
        }

        if (use_for_test)
        {
          h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
        }

        continue;
      }

      if (match[isample] && reco_good && truth_good)
      {
        if (use_for_response)
        {
          h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
          h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
          h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
          h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
          rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
          h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
          h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
          h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
          h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
        }
        else if (use_for_test)
        {
          h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
          h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
          h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
          h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
        }

        h_mbd_vertex->Fill(mbd_vertex[isample], event_scale);
        if (!ispp && has_sumeT[isample]) { h_sumeT->Fill(sumeT[isample], event_scale); }
        if (!ispp) { h_centrality->Fill(centrality[isample], event_scale); }
      }
    } // end of entry loop
  } // end of sample loop

  // fills happen in symmetric (a,b)+(b,a) pairs throughout, so every flat/2D
  // histogram above is exactly double-counted
  h_flat_reco_pt1pt2->Scale(0.5);
  h_flat_truth_pt1pt2->Scale(0.5);
  h_flat_reco_to_response_pt1pt2->Scale(0.5);
  h_flat_truth_to_response_pt1pt2->Scale(0.5);
  h_flat_response_pt1pt2->Scale(0.5);

  // Skim down to bins with at least `minentries` counts on both truth and
  // reco sides, and remap the response/truth/reco histograms onto the
  // resulting smaller (non-empty) binning before unfolding.
  std::cout << "Minimum entries == " << minentries << std::endl;

  TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone("h_flat_truth_mapping");
  TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone("h_flat_reco_mapping");
  h_flat_truth_mapping->Reset();
  h_flat_reco_mapping->Reset();

  int nempty_reco = 0;
  int nempty_truth = 0;
  std::vector<int> binnumbers_reco;
  std::vector<int> binnumbers_truth;
  int nrecobins = 0;
  int ntruthbins = 0;
  for (int ib = 0; ib < nbins*nbins; ib++)
  {
    if (h_count_flat_reco_pt1pt2->GetBinContent(ib+1) < minentries)
    {
      nempty_reco++;
      binnumbers_reco.push_back(0);
      h_flat_reco_mapping->SetBinContent(ib+1, 0);
    }
    else
    {
      nrecobins++;
      binnumbers_reco.push_back(nrecobins);
      h_flat_reco_mapping->SetBinContent(ib+1, nrecobins);
    }
    if (h_count_flat_truth_pt1pt2->GetBinContent(ib+1) < minentries)
    {
      nempty_truth++;
      binnumbers_truth.push_back(0);
      h_flat_truth_mapping->SetBinContent(ib+1, 0);
    }
    else
    {
      ntruthbins++;
      binnumbers_truth.push_back(ntruthbins);
      h_flat_truth_mapping->SetBinContent(ib+1, ntruthbins);
    }
  }

  TH1D *h_flat_truth_skim = new TH1D("h_flat_truth_skim", "", ntruthbins, 0, ntruthbins);
  TH1D *h_flat_reco_skim = new TH1D("h_flat_reco_skim", "", nrecobins, 0, nrecobins);
  TH1D *h_flat_truth_to_unfold_skim = new TH1D("h_flat_truth_to_unfold_skim", "", ntruthbins, 0, ntruthbins);
  TH1D *h_flat_reco_to_unfold_skim = new TH1D("h_flat_reco_to_unfold_skim", "", nrecobins, 0, nrecobins);
  TH2D *h_flat_response_skim = new TH2D("h_flat_response_skim", "", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

  for (int ib = 0; ib < nbins*nbins; ib++)
  {
    const int bintruth = binnumbers_truth[ib];
    const int binreco = binnumbers_reco[ib];

    if (binreco)
    {
      h_flat_reco_skim->SetBinContent(binreco, h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1));
      h_flat_reco_skim->SetBinError(binreco, h_flat_reco_to_response_pt1pt2->GetBinError(ib+1));
      h_flat_reco_to_unfold_skim->SetBinContent(binreco, h_flat_reco_pt1pt2->GetBinContent(ib+1));
      h_flat_reco_to_unfold_skim->SetBinError(binreco, h_flat_reco_pt1pt2->GetBinError(ib+1));
    }

    if (bintruth)
    {
      h_flat_truth_skim->SetBinContent(bintruth, h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1));
      h_flat_truth_skim->SetBinError(bintruth, h_flat_truth_to_response_pt1pt2->GetBinError(ib+1));
      h_flat_truth_to_unfold_skim->SetBinContent(bintruth, h_flat_truth_pt1pt2->GetBinContent(ib+1));
      h_flat_truth_to_unfold_skim->SetBinError(bintruth, h_flat_truth_pt1pt2->GetBinError(ib+1));

      for (int ibr = 0; ibr < nbins*nbins; ibr++)
      {
        const int binreco2 = binnumbers_reco[ibr];
        if (binreco2)
        {
          const int rbin = h_flat_response_skim->GetBin(binreco2, bintruth);
          const int tbin = h_flat_response_pt1pt2->GetBin(ibr+1, ib+1);
          h_flat_response_skim->SetBinContent(rbin, h_flat_response_pt1pt2->GetBinContent(tbin));
          h_flat_response_skim->SetBinError(rbin, h_flat_response_pt1pt2->GetBinError(tbin));
        }
      }
    }
  }

  std::cout << "Nbins skim reco = " << h_flat_reco_skim->GetNbinsX() << std::endl;
  std::cout << "Nbins skim truth = " << h_flat_truth_skim->GetNbinsX() << std::endl;
  std::cout << "Nbins skim response = " << h_flat_response_skim->GetXaxis()->GetNbins()
             << " -- " << h_flat_response_skim->GetYaxis()->GetNbins() << std::endl;
  std::cout << "Reco empty: " << nempty_reco << std::endl;
  std::cout << "Truth empty: " << nempty_truth << std::endl;

  RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);
  rooResponsehist.SetName("response_noempty");

  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
  {
    RooUnfoldBayes unfold(&rooResponsehist, h_flat_reco_to_unfold_skim, iter + 1);
    TH1D *h_flat_unfold_skim = (TH1D*) unfold.Hunfold();

    h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone(Form("h_flat_unfold_pt1pt2_%d", iter));
    h_flat_unfold_pt1pt2[iter]->Reset();
    for (int ib = 0; ib < nbins*nbins; ib++)
    {
      const int bin = binnumbers_truth[ib];
      if (bin)
      {
        h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim->GetBinContent(bin));
        h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim->GetBinError(bin));
      }
    }
  }

  // x_J projections of the truth/reco/unfolded pt1-pt2 matrices, reusing the
  // same histo_opps helpers the draw macros expect these histograms to have
  // been built with
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  for (int iter = 0; iter < niterations; iter++)
  {
    h_pt1pt2_unfold[iter] = new TH2D(Form("h_pt1pt2_unfold_iter%d", iter), ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
    histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
  }

  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];

  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
  {
    h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};", nbins, ixj_bins);
    histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
  }

  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
  {
    histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
  }

  std::cout << "writing" << std::endl;
  TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
  if (primer > 0)
  {
    responsepath += "_PRIMER" + std::to_string(primer);
  }
  if (half_closure)
  {
    responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_HALF";
  }
  else if (full_closure)
  {
    responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_FULL";
  }
  responsepath += "_" + sys_name;
  responsepath += ".root";

  TFile *fr = new TFile(responsepath.Data(), "recreate");
  rooResponsehist.Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_flat_truth_to_response_pt1pt2->Write();
  h_flat_response_pt1pt2->Write();
  h_flat_truth_mapping->Write();
  h_flat_reco_mapping->Write();
  h_flat_truth_skim->Write();
  h_flat_reco_skim->Write();
  h_flat_response_skim->Write();
  h_mbd_vertex->Write();
  h_centrality->Write();
  h_sumeT->Write();
  h_xj_truth->Write();
  h_xj_reco->Write();
  for (int iter = 0; iter < niterations; iter++)
  {
    h_flat_unfold_pt1pt2[iter]->Write();
    h_xj_unfold[iter]->Write();
  }
  fr->Close();

  return 0;
}
