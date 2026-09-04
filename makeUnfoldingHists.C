#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

#include "read_binning.h"


struct jet
{
  int id;
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
  int ptbin = -1;

  void print()
  {
    std::cout << "jet " << id << " : " << pt << " ( " << ptbin << " ) " << std::endl;
  }
};


void makeUnfoldingHists(
  const int cone_size = 3, 
  const int centrality_bin = 0,
  const std::string & configfile = "binning_AA.config",  
  const std::string & infile = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root"
)
{
  const bool ispp = ( centrality_bin < 0 );


  read_binning rb(configfile.c_str());

  std::string system_string = rb.get_system_string(centrality_bin);
  std::cout << "System string: " << system_string << std::endl;

  const double dphicut = rb.get_dphicut();
  const int nbins = rb.get_nbins();

  const float vertex_cut = static_cast<float>( rb.get_vtx_cut() );
  // Per-jet fiducial |eta| acceptance now comes from the jet_accept_eta flag
  // in the input tree rather than a hard cut computed from the cone size.
  const float etacut_bkg = rb.get_eta_cut_bkg();
  const float reco_cut = static_cast<float>( rb.get_reco_pt_min_cut() );

  int zyam_sys = rb.get_zyam_sys();

  double flow_v22_scale = rb.get_flow_sys();
  double flow_v33_scale = rb.get_flow_v33_sys();
  const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
  
  int prior_sys = rb.get_prior_sys();

  double JES_sys = rb.get_jes_sys();
  double JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  std::string sys_name = "nominal";
  
  if (prior_sys) sys_name = "PRIOR";
  if (zyam_sys) sys_name = "ZYAM";
  if (flow_sys) sys_name = rb.get_flow_systematic_name();
  if (JER_sys < 0) sys_name = "negJER";
  if (JER_sys > 0) sys_name = "posJER";
  if (JES_sys < 0) sys_name = "negJES";
  if (JES_sys > 0) sys_name = "posJES";
    
  float low_fit = rb.get_flow_fit_low();
  float high_fit = rb.get_flow_fit_high();
  float zyam_integral_low = rb.get_zyam_region_low();
  float zyam_integral_high = rb.get_zyam_region_high();
  
  float dphicut_z_low = rb.get_zyam_low();
  float dphicut_z_high = rb.get_zyam_high();

  float nbins_diff = (TMath::Pi() - dphicut)/(TMath::Pi()/32.);
  float nbins_zyam = (zyam_integral_high - zyam_integral_low)/(TMath::Pi()/32.);
 
  float ZYAM_scale = (TMath::Pi() - dphicut)/(dphicut_z_high - dphicut_z_low);

  const int n_centrality_bins = rb.get_number_centrality_bins();  
  float icentrality_bins[n_centrality_bins+1];
  rb.get_centrality_bins(icentrality_bins);
  for (int i = 0; i < n_centrality_bins + 1; i++)
  {
    std::cout << "Centrality bin " << i << ": " << icentrality_bins[i] << std::endl;
  }
  
  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
  {
    std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
  }

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();

  float sample_boundary[4] = {0};
  float low_trigger[3] = {0};
  for (int ib = 0; ib < 4; ib++)
  {
    sample_boundary[ib] = rb.get_sample_boundary(ib);
    std::cout <<  sample_boundary[ib] << std::endl;
  }

  int max_reco_bin = rb.get_maximum_reco_bin();
  std::cout << "Max reco: " << max_reco_bin << " - ( " << ipt_bins[max_reco_bin] << ")" << std::endl;
  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  double scale_to_signal[nbins][nbins];
  for (int i = 0; i < nbins; i++)
  {
    for (int j = 0; j < nbins; j++)
    {
      scale_to_signal[i][j] = 1;
    }
  }

  TH1D *h_pt2_prob = nullptr;
  
  if (!ispp)
  {

    TString probpath = rb.get_code_location() + "/unfolding_hists/probability_hists_AA_r0" + std::to_string(cone_size) + ".root";

    TFile *fprob = new TFile(probpath.Data(),"read");
    if (!fprob || fprob->IsZombie())
    {
      std::cerr << "Missing probability-correction file: " << probpath << std::endl;
      return;
    }
    h_pt2_prob = (TH1D*) fprob->Get(Form("h_pt2_bin_log_correction_%d", centrality_bin));
    if (!h_pt2_prob)
    {
      std::cerr << "Missing correction histogram for centrality bin " << centrality_bin
                << " in " << probpath << std::endl;
      return;
    }

  }

  TFile *f = new TFile(infile.c_str(), "r");
  TTree *t = (TTree*) f->Get("T");
  if (!t)
  {
    std::cerr << "Could not find tree T in " << infile << std::endl;
    return;
  }
  if (cone_size != 3)
  {
    std::cerr << "This input file contains R=0.3 jets only; use cone_size = 3." << std::endl;
    return;
  }
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("cent", 1);
  t->SetBranchStatus("zvrtx", 1);
  t->SetBranchStatus("sumeT", 1);
  t->SetBranchStatus("jet_pT", 1);
  bool has_jet_comp_pT = (t->GetBranch("jet_comp_pT") != nullptr);
  if ( has_jet_comp_pT ) t->SetBranchStatus("jet_comp_pT", 1);
  else t->SetBranchStatus("jet_unsub_pT", 1);

  t->SetBranchStatus("jet_E", 1);
  t->SetBranchStatus("jet_unsub_E", 1);
  t->SetBranchStatus("jet_eta", 1);
  t->SetBranchStatus("jet_phi", 1);
  const bool has_jet_accept_eta = (t->GetBranch("jet_accept_eta") != nullptr);
  if ( has_jet_accept_eta ) t->SetBranchStatus("jet_accept_eta", 1);
  std::cout << "Has jet_accept_eta branch: " << has_jet_accept_eta << std::endl;
  t->SetCacheSize(256 * 1024 * 1024);
  t->AddBranchToCache("*", kTRUE);
  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_unsub = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;
  std::vector<int>   *reco_jet_accept_eta = 0;
  int centrality;
  float mbd_vertex_z;
  float sumeT = 0;
  t->SetBranchAddress("cent", &centrality);
  t->SetBranchAddress("jet_pT", &reco_jet_pt);
  // t->SetBranchAddress("jet_comp_pT", &reco_jet_pt_unsub); // renamed: v004 background-quality discriminant is jet_comp_pT, not jet_unsub_pT
  if ( has_jet_comp_pT ) t->SetBranchAddress("jet_comp_pT", &reco_jet_pt_unsub);
  else t->SetBranchAddress("jet_unsub_pT", &reco_jet_pt_unsub);

  t->SetBranchAddress("jet_E", &reco_jet_e);
  t->SetBranchAddress("jet_unsub_E", &reco_jet_e_unsub);
  t->SetBranchAddress("jet_eta", &reco_jet_eta);
  t->SetBranchAddress("jet_phi", &reco_jet_phi);
  // t->SetBranchAddress("jet_accept_eta", &reco_jet_accept_eta);
  if ( has_jet_accept_eta ) t->SetBranchAddress("jet_accept_eta", &reco_jet_accept_eta);
  t->SetBranchAddress("zvrtx", &mbd_vertex_z);
  t->SetBranchAddress("sumeT", &sumeT);

  // This reduced tree has no minbias, gl1_scaled, or mbd_time_zero branches.
  // mbd_time_zero was unused.  The min-bias and GL1-trigger cuts below are
  // therefore treated as already applied when this input format is used.
  std::cout << "Input tree has no minbias or GL1 trigger branch; assuming the input is preselected." << std::endl;

  std::vector<struct jet> myrecojets_inclusive;
  std::vector<struct jet> myrecojets;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  
  int entries = t->GetEntries();
  

  TH1D *h_dphi_exclusive_all = new TH1D("h_dphi_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_exclusive_all = new TH1D("h_dphi_eta_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_inclusive_all = new TH1D("h_dphi_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_inclusive_all = new TH1D("h_dphi_eta_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_exclusive[nbins][nbins];
  TH1D *h_dphi_eta_exclusive[nbins][nbins];

  TH1D *h_dphi_inclusive[nbins][nbins];
  TH1D *h_dphi_eta_inclusive[nbins][nbins];
  
  for (int i = 0; i < nbins; i++)
  {
    for (int j = 0; j < nbins; j++)
    {
      h_dphi_exclusive[i][j] = new TH1D(Form("h_dphi_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      h_dphi_eta_exclusive[i][j] = new TH1D(Form("h_dphi_eta_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      
      h_dphi_inclusive[i][j] = new TH1D(Form("h_dphi_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      h_dphi_eta_inclusive[i][j] = new TH1D(Form("h_dphi_eta_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
    }
  }

  TH2D *h_pt1_pt2_signal_zyam_exclusive = new TH2D("h_pt1_pt2_signal_zyam_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1_pt2_signal_exclusive = new TH2D("h_pt1_pt2_signal_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1_pt2_zyam_exclusive = new TH2D("h_pt1_pt2_zyam_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);

  TH2D *h_pt1_pt2_signal_zyam_inclusive = new TH2D("h_pt1_pt2_signal_zyam_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1_pt2_signal_inclusive = new TH2D("h_pt1_pt2_signal_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1_pt2_zyam_inclusive = new TH2D("h_pt1_pt2_zyam_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
  TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);

  TH1D *h_flat_data_pt1pt2_signal_zyam_exclusive = new TH1D("h_data_flat_pt1pt2_signal_zyam_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_data_pt1pt2_signal_exclusive = new TH1D("h_data_flat_pt1pt2_signal_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_data_pt1pt2_zyam_exclusive = new TH1D("h_data_flat_pt1pt2_zyam_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  TH1D *h_flat_data_pt1pt2_signal_zyam_inclusive = new TH1D("h_data_flat_pt1pt2_signal_zyam_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_data_pt1pt2_signal_inclusive = new TH1D("h_data_flat_pt1pt2_signal_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_data_pt1pt2_zyam_inclusive = new TH1D("h_data_flat_pt1pt2_zyam_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  
  std::pair<int, float> id_leaders[2];

  TF1 *fcut = new TF1("fcut","[0]+[1]*TMath::Exp(-[2]*x)",0.0,100.0);
  if ( has_jet_accept_eta ) 
  {
    // new verison of file
    fcut -> SetParameters( -2.29, 41.7,0.029 );
    // fcut -> SetParameters( 7.95, 34.3, 0.047 );
  }
  else 
  {
    fcut -> SetParameters( 0, 40.0, 0.035 );
  }
  std::cout << "fcut parameters: " << fcut->GetParameter(0) << ", " << fcut->GetParameter(1) << ", " << fcut->GetParameter(2) << std::endl;

  for (int i = 0; i < entries; i++)
  {
    t->GetEntry( i);

    if (i % int(entries/10) == 0) std::cout << "Event: " << i << " / " << entries << "\r" << std::flush;

    if (fabs(mbd_vertex_z) > vertex_cut) continue;
    if (centrality < icentrality_bins[centrality_bin] || centrality >= icentrality_bins[centrality_bin+1]) continue;
    
    int nrecojets = reco_jet_pt->size();
    
    myrecojets.clear();
    myrecojets_inclusive.clear();
    
    bool found_reco_dijet = false;
    
    id_leaders[0] = std::make_pair(0, 0);
    id_leaders[1] = std::make_pair(0, 0);

    int njet_good = 0;
    int bad_event = 0;
    float cut_value = fcut->Eval(centrality);

    for (int j = 0; j < nrecojets;j++)
	  {
      if (reco_jet_pt->at(j) < reco_cut) continue;
      if (reco_jet_e->at(j) < 0) continue;
      if (reco_jet_e_unsub->at(j) < 0) continue;
      // if (!reco_jet_accept_eta->at(j)) continue;
      bool eta_good = has_jet_accept_eta ? (reco_jet_accept_eta->at(j) == 1)  : (fabs(reco_jet_eta->at(j)) < 0.8f);
      if (!eta_good) continue;

	    float pt_unsub = reco_jet_pt_unsub->at(j) - reco_jet_pt->at(j);

      if (pt_unsub > cut_value)
      {
        bad_event = 1;
        continue;
      }


	    njet_good++;

      if (reco_jet_pt->at(j) > id_leaders[0].second)
      {
        id_leaders[1] = id_leaders[0];
        id_leaders[0] = std::make_pair(j, reco_jet_pt->at(j));
      }
      else if (reco_jet_pt->at(j) > id_leaders[1].second)
      {
        id_leaders[1] = std::make_pair(j, reco_jet_pt->at(j));
      } 	  

      struct jet tempjet;
      tempjet.pt = reco_jet_pt->at(j);
      tempjet.eta = reco_jet_eta->at(j);
      tempjet.phi = reco_jet_phi->at(j);
      tempjet.id = j;
      int pt_reco_bin = nbins;

	    float es = tempjet.pt;

      if (es >= ipt_bins[max_reco_bin] )
      {
        bad_event = 1;
        continue;
      }
	    
      for (int ib = 0; ib < nbins; ib++)
	    {

	      if ( es < ipt_bins[ib+1] && es >= ipt_bins[ib])
        {
          pt_reco_bin = ib;
        }
	    }
      tempjet.ptbin = pt_reco_bin;
      myrecojets_inclusive.push_back(tempjet);
	  
	  }

    if (id_leaders[0].second < reco_leading_cut)
    {
      continue;
    }
    if (bad_event) continue;
    if (njet_good < 2) continue;

    h_centrality->Fill(centrality);
    h_mbd_vertex->Fill(mbd_vertex_z);
    h_sumeT->Fill(sumeT);
    struct jet tempjet1;
    tempjet1.pt = reco_jet_pt->at(id_leaders[0].first);
    tempjet1.eta = reco_jet_eta->at(id_leaders[0].first);
    tempjet1.phi = reco_jet_phi->at(id_leaders[0].first);
    tempjet1.id = id_leaders[0].first;

    int pt_reco_bin1 = nbins;

    float es1 = tempjet1.pt;

    for (int ib = 0; ib < nbins; ib++)
	  {

      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
      {
        pt_reco_bin1 = ib;
      }
    }
    tempjet1.ptbin = pt_reco_bin1;
    
    myrecojets.push_back(tempjet1);
    struct jet tempjet2;
    tempjet2.pt = reco_jet_pt->at(id_leaders[1].first);
    tempjet2.eta = reco_jet_eta->at(id_leaders[1].first);
    tempjet2.phi = reco_jet_phi->at(id_leaders[1].first);
    tempjet2.id = id_leaders[1].first;;

    int pt_reco_bin2 = nbins;

    float es2 = tempjet2.pt;

    for (int ib = 0; ib < nbins; ib++)
    {

      if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
      {
        pt_reco_bin2 = ib;
      }
    }

    tempjet2.ptbin = pt_reco_bin2;
    
    myrecojets.push_back(tempjet2);

    double dphir = 0;
    double detar = 0;
    auto leading_iter = myrecojets.begin();


    auto subleading_iter = myrecojets.begin() + 1;
    if (subleading_iter != myrecojets.end())
	  {
      detar = fabs(leading_iter->eta - subleading_iter->eta);
      dphir = fabs(leading_iter->phi - subleading_iter->phi);
      if (dphir > TMath::Pi())
      {
        dphir = 2*TMath::Pi() - dphir;
      }

      float pt1_reco_bin = leading_iter->ptbin;;
      float pt2_reco_bin = subleading_iter->ptbin;

      bool fill_good = (leading_iter->pt >= reco_leading_cut && subleading_iter->pt >= reco_subleading_cut);
      bool signal_good = (leading_iter->pt >= reco_leading_cut && subleading_iter->pt >= reco_subleading_cut && dphir >= dphicut);
      bool ZYAM_good = (leading_iter->pt >= reco_leading_cut && subleading_iter->pt >= reco_subleading_cut && dphir >= dphicut_z_low && dphir < dphicut_z_high);

      if (fill_good)
      {
        h_dphi_exclusive_all->Fill(dphir);
        h_dphi_exclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);
        if (detar > etacut_bkg)
        {
          h_dphi_eta_exclusive_all->Fill(dphir);
          h_dphi_eta_exclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);
        }
      }
      if (signal_good)
      {
        h_flat_data_pt1pt2_signal_zyam_exclusive->Fill(leading_iter->ptbin + nbins*subleading_iter->ptbin);
        h_flat_data_pt1pt2_signal_zyam_exclusive->Fill(subleading_iter->ptbin + nbins*leading_iter->ptbin);
        h_pt1_pt2_signal_zyam_exclusive->Fill(leading_iter->pt, subleading_iter->pt);
        h_pt1_pt2_signal_zyam_exclusive->Fill(subleading_iter->pt, leading_iter->pt);
      }

	    if (ZYAM_good)
	    {
	      h_flat_data_pt1pt2_zyam_exclusive->Fill(leading_iter->ptbin + nbins*subleading_iter->ptbin);
	      h_flat_data_pt1pt2_zyam_exclusive->Fill(subleading_iter->ptbin + nbins*leading_iter->ptbin);
	      h_pt1_pt2_zyam_exclusive->Fill(leading_iter->pt, subleading_iter->pt);
	      h_pt1_pt2_zyam_exclusive->Fill(subleading_iter->pt, leading_iter->pt);
	    }
	  }
      
    for (auto sub_iter :  myrecojets_inclusive)
	  {
      if (sub_iter.id == leading_iter->id) continue;

      detar = fabs(leading_iter->eta - sub_iter.eta);
      dphir = fabs(leading_iter->phi - sub_iter.phi);
      if (dphir > TMath::Pi())
      {
        dphir = 2*TMath::Pi() - dphir;
      }
      
      bool fill_good = (leading_iter->pt >= reco_leading_cut && sub_iter.pt >= reco_subleading_cut);
      bool signal_good = (leading_iter->pt >= reco_leading_cut && sub_iter.pt >= reco_subleading_cut && dphir >= dphicut);
      bool ZYAM_good = (leading_iter->pt >= reco_leading_cut && sub_iter.pt >= reco_subleading_cut && dphir >= dphicut_z_low && dphir < dphicut_z_high);

      if (fill_good)
      {
        h_dphi_inclusive_all->Fill(dphir);
        h_dphi_inclusive[leading_iter->ptbin][sub_iter.ptbin]->Fill(dphir);
        if (detar > etacut_bkg)
        {
          h_dphi_eta_inclusive_all->Fill(dphir);
          h_dphi_eta_inclusive[leading_iter->ptbin][sub_iter.ptbin]->Fill(dphir);
        }
      }
      if (signal_good)
      {
        h_flat_data_pt1pt2_signal_zyam_inclusive->Fill(leading_iter->ptbin + nbins*sub_iter.ptbin);
        h_flat_data_pt1pt2_signal_zyam_inclusive->Fill(sub_iter.ptbin + nbins*leading_iter->ptbin);
        h_pt1_pt2_signal_zyam_inclusive->Fill(leading_iter->pt, sub_iter.pt);
        h_pt1_pt2_signal_zyam_inclusive->Fill(sub_iter.pt, leading_iter->pt);
      }

      if (ZYAM_good)
      {
        h_flat_data_pt1pt2_zyam_inclusive->Fill(leading_iter->ptbin + nbins*sub_iter.ptbin);
        h_flat_data_pt1pt2_zyam_inclusive->Fill(sub_iter.ptbin + nbins*leading_iter->ptbin);
        h_pt1_pt2_zyam_inclusive->Fill(leading_iter->pt, sub_iter.pt);
        h_pt1_pt2_zyam_inclusive->Fill(sub_iter.pt, leading_iter->pt);
      }
	  }
	  
  } // end of event loop


  TF1 *ffit = new TF1("ffit", "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", low_fit, high_fit);

  TF1 *fits_ex[nbins][nbins];
  TF1 *fits_in[nbins][nbins];

  ffit->SetParLimits(0, 0, 10000);
  ffit->SetParLimits(1, 0, 0.5);
  ffit->SetParLimits(2, 0, 0.5);
  // iterate through leading bins
  for (int i = 0; i < nbins; i++)
  {
    // iterate through subleading bins
    for (int j = 0; j <= i; j++)
	  {
      std::cout << " ---------- " << i << " - " << j << " -----------" << std::endl;
      
      int nbins_over = h_dphi_eta_inclusive[i][j]->FindBin(high_fit) - h_dphi_eta_inclusive[i][j]->FindBin(low_fit);
      int ncounts = h_dphi_eta_inclusive[i][j]->Integral(h_dphi_eta_inclusive[i][j]->FindBin(low_fit), h_dphi_eta_inclusive[i][j]->FindBin(high_fit));

      // if below this count - use linear fit of fit region
      std::cout << ncounts << std::endl;	  
      if (ncounts < nbins_over)
      {
        ffit->SetParameters(0, 0, 0);
        ffit->FixParameter(1, 0);
        ffit->FixParameter(2, 0);
      }
      else
      {
        ffit->ReleaseParameter(1);
        ffit->ReleaseParameter(2);
        ffit->ReleaseParameter(0);
        ffit->SetParLimits(0, 0, 10000);
        ffit->SetParLimits(1, 0, 0.5);
        ffit->SetParLimits(2, 0, 0.5);
        ffit->SetParameters(0, 0, 0);
      }


	    h_dphi_eta_inclusive[i][j]->Fit("ffit", "0RlQ");
      if (flow_sys)
      {
        std::cout << sys_name << " flow systematic: scaling fitted (v22,v33) from ("
                  << ffit->GetParameter(1) << "," << ffit->GetParameter(2) << ") to ("
                  << flow_v22_scale * ffit->GetParameter(1) << ","
                  << flow_v33_scale * ffit->GetParameter(2) << ")"
                  << " for bin " << i << "," << j << std::endl;
        ffit->SetParameter(1, flow_v22_scale * ffit->GetParameter(1));
        ffit->SetParameter(2, flow_v33_scale * ffit->GetParameter(2));
      }

      fits_ex[i][j] = (TF1*) ffit->Clone();
      fits_ex[i][j]->SetName(Form("ffit_ex_%d_%d", i, j));

      fits_ex[i][j]->SetParameters(ffit->GetParameter(0),ffit->GetParameter(1),ffit->GetParameter(2));

      // get first and second bin
      
      int first_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_low);
      int second_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_high);
      
      float bin_diff = second_bin - first_bin;
      
      //calculate the difference and scale up the histogram
      float dphi_zyam_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin)/bin_diff;
      float dphi_flow_integral = ffit->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high -  zyam_integral_low);
      
      float par1 = ffit->GetParameter(0);
      
      float integral_diff = dphi_zyam_integral - dphi_flow_integral;

      par1 += integral_diff;

      fits_ex[i][j]->SetParameter(0,par1);
      fits_ex[i][j]->SetRange(0, TMath::Pi());

      // now get the signal
          
      first_bin = h_dphi_exclusive[i][j]->FindBin(dphicut);
      second_bin = -1;
      
      float flow_signal = nbins_diff*fits_ex[i][j]->Integral(dphicut, TMath::Pi())/(TMath::Pi() - dphicut);
      float dphi_signal_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin);
      
      float signal_from_fit = dphi_signal_integral - flow_signal;
      if (signal_from_fit < 0) signal_from_fit = 0;

      int sbin = std::min(i, j);
      signal_from_fit /= h_pt2_prob->GetBinContent(sbin + 1);
      float nerr = sqrt(dphi_signal_integral)/h_pt2_prob->GetBinContent(sbin + 1);

      int gbin1 = h_pt1_pt2_signal_exclusive->GetBin(i+1, j+1);
      int gbin2 = h_pt1_pt2_signal_exclusive->GetBin(j+1, i+1);
      h_pt1_pt2_signal_exclusive->SetBinContent(gbin1, signal_from_fit);
      h_pt1_pt2_signal_exclusive->SetBinError(gbin1, nerr);
      h_pt1_pt2_signal_exclusive->SetBinContent(gbin2, signal_from_fit);
      h_pt1_pt2_signal_exclusive->SetBinError(gbin2, nerr);
	  
      int ggbin1 = 1 + i + nbins*j;
      int ggbin2 = 1 + j + nbins*i;
      h_flat_data_pt1pt2_signal_exclusive->SetBinContent(ggbin1, signal_from_fit);
      h_flat_data_pt1pt2_signal_exclusive->SetBinContent(ggbin2, signal_from_fit);
      h_flat_data_pt1pt2_signal_exclusive->SetBinError(ggbin1, nerr);
      h_flat_data_pt1pt2_signal_exclusive->SetBinError(ggbin2, nerr);
      
      // now inclusive
      fits_in[i][j] = (TF1*) ffit->Clone();
      fits_in[i][j]->SetName(Form("ffit_ex_%d_%d", i, j));

      fits_in[i][j]->SetParameters(ffit->GetParameter(0), ffit->GetParameter(1), ffit->GetParameter(2));

      // get first and second bin
      
      first_bin = h_dphi_inclusive[i][j]->FindBin(zyam_integral_low);
      second_bin = h_dphi_inclusive[i][j]->FindBin(zyam_integral_high);
      
      bin_diff = second_bin - first_bin;
      
      //calculate the difference and scale up the histogram
      dphi_zyam_integral = h_dphi_inclusive[i][j]->Integral(first_bin, second_bin)/bin_diff;
      dphi_flow_integral = ffit->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high -  zyam_integral_low);
      
      par1 = ffit->GetParameter(0);
      
      integral_diff = dphi_zyam_integral - dphi_flow_integral;
      
      par1 += integral_diff;
      
      fits_in[i][j]->SetParameter(0,par1);
      fits_in[i][j]->SetRange(0, TMath::Pi());
      
      // now get the signal
      
      first_bin = h_dphi_inclusive[i][j]->FindBin(dphicut);
      second_bin = -1;
      
      flow_signal = nbins_diff*fits_in[i][j]->Integral(dphicut, TMath::Pi())/(TMath::Pi() - dphicut);
      dphi_signal_integral = h_dphi_inclusive[i][j]->Integral(first_bin, second_bin);
      
      signal_from_fit = dphi_signal_integral - flow_signal;
      if (signal_from_fit < 0) signal_from_fit = 0;

      nerr = sqrt(dphi_signal_integral);

      h_pt1_pt2_signal_inclusive->SetBinContent(gbin1, signal_from_fit);
      h_pt1_pt2_signal_inclusive->SetBinError(gbin1, nerr);
      h_pt1_pt2_signal_inclusive->SetBinContent(gbin2, signal_from_fit);
      h_pt1_pt2_signal_inclusive->SetBinError(gbin2, nerr);
      
      h_flat_data_pt1pt2_signal_inclusive->SetBinContent(ggbin1, signal_from_fit);
      h_flat_data_pt1pt2_signal_inclusive->SetBinContent(ggbin2, signal_from_fit);
      h_flat_data_pt1pt2_signal_inclusive->SetBinError(ggbin1, nerr);
      h_flat_data_pt1pt2_signal_inclusive->SetBinError(ggbin2, nerr);
	  
	  }
  } // end of leading/subleading bin loop
  
  std::cout << "subtracting_zyam" << std::endl;  
  // now the regular signal is done, not to just subtract the zyam
  h_pt1_pt2_signal_zyam_exclusive->Add(h_pt1_pt2_zyam_exclusive, -1);
  h_pt1_pt2_signal_zyam_inclusive->Add(h_pt1_pt2_zyam_inclusive, -1);
  h_flat_data_pt1pt2_signal_zyam_exclusive->Add(h_flat_data_pt1pt2_zyam_exclusive, -1);
  h_flat_data_pt1pt2_signal_zyam_inclusive->Add(h_flat_data_pt1pt2_zyam_inclusive, -1);

  int nbins_1 = h_flat_data_pt1pt2_signal_zyam_exclusive->GetXaxis()->GetNbins();
  int nbins_2 = h_pt1_pt2_signal_zyam_exclusive->GetXaxis()->GetNbins();
  std::cout << "checking_zyam" << std::endl;  
  for (int ix = 0; ix < nbins; ix++)
  {
    for (int iy = 0; iy < nbins; iy++)
	  {
      int gbin = 1 + ix + nbins*iy;

      if (h_flat_data_pt1pt2_signal_zyam_exclusive->GetBinContent(gbin) < 0)
      {
        h_flat_data_pt1pt2_signal_zyam_exclusive->SetBinContent(gbin, 0);
        h_flat_data_pt1pt2_signal_zyam_exclusive->SetBinError(gbin, 0);
      }
      if (h_flat_data_pt1pt2_signal_zyam_inclusive->GetBinContent(gbin) < 0)
      {
        h_flat_data_pt1pt2_signal_zyam_inclusive->SetBinContent(gbin, 0);
        h_flat_data_pt1pt2_signal_zyam_inclusive->SetBinError(gbin, 0);
      }
      int sbin = std::min(ix, iy);
      float nen = h_flat_data_pt1pt2_signal_zyam_exclusive->GetBinContent(gbin);
      float nerr = h_flat_data_pt1pt2_signal_zyam_exclusive->GetBinError(gbin);

      nen /= h_pt2_prob->GetBinContent(sbin + 1);
      nerr /= h_pt2_prob->GetBinContent(sbin + 1);
        
      h_flat_data_pt1pt2_signal_zyam_exclusive->SetBinContent(gbin, nen);
      h_flat_data_pt1pt2_signal_zyam_exclusive->SetBinError(gbin, nerr);
	  }
  } // end of loop over bins
 
  std::cout << "checking_zyam" << std::endl;  
  for (int ix = 0 ; ix < nbins_2; ix++)
  {
    for (int iy = 0 ; iy < nbins_2; iy++)
	  {

      int gbin = h_pt1_pt2_signal_zyam_exclusive->GetBin(ix+1, iy+1);
      if (h_pt1_pt2_signal_zyam_exclusive->GetBinContent(gbin) < 0)
        {
          h_pt1_pt2_signal_zyam_exclusive->SetBinContent(gbin, 0);
          h_pt1_pt2_signal_zyam_exclusive->SetBinError(gbin, 0);
        }
      if (h_pt1_pt2_signal_zyam_inclusive->GetBinContent(gbin) < 0)
        {
          h_pt1_pt2_signal_zyam_inclusive->SetBinContent(gbin, 0);
          h_pt1_pt2_signal_zyam_inclusive->SetBinError(gbin, 0);
        }
      int sbin = std::min(ix, iy);
      float nen = h_pt1_pt2_signal_zyam_exclusive->GetBinContent(gbin);
      float nerr = h_pt1_pt2_signal_zyam_exclusive->GetBinError(gbin);

      nen /= h_pt2_prob->GetBinContent(sbin + 1);
      nerr /= h_pt2_prob->GetBinContent(sbin + 1);

      h_pt1_pt2_signal_zyam_exclusive->SetBinContent(gbin, nen);
      h_pt1_pt2_signal_zyam_exclusive->SetBinError(gbin, nerr);

    }
  } // end of loop over bins

  std::cout << "writing" << std::endl;
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_preload_" + system_string + "_r0" + std::to_string(cone_size);

  unfoldpath += "_" + sys_name;
  unfoldpath += ".root";

  TFile *fout = new TFile(unfoldpath.Data(), "recreate");
  h_centrality->Write();
  h_mbd_vertex->Write();
  h_sumeT->Write();
  h_flat_data_pt1pt2_signal_zyam_exclusive->Write();
  h_flat_data_pt1pt2_signal_exclusive->Write();
  h_flat_data_pt1pt2_zyam_exclusive->Write();
  h_pt1_pt2_signal_zyam_exclusive->Write();
  h_pt1_pt2_signal_exclusive->Write();
  h_pt1_pt2_zyam_exclusive->Write();
  h_flat_data_pt1pt2_signal_zyam_inclusive->Write();
  h_flat_data_pt1pt2_signal_inclusive->Write();
  h_flat_data_pt1pt2_zyam_inclusive->Write();
  h_pt1_pt2_signal_zyam_inclusive->Write();
  h_pt1_pt2_signal_inclusive->Write();
  h_pt1_pt2_zyam_inclusive->Write();
  for (int i = 0; i < nbins; i++)
  {
    for (int j = 0; j < nbins; j++)
	  {

      h_dphi_exclusive[i][j]->Write();
      h_dphi_eta_exclusive[i][j]->Write();	  
      h_dphi_inclusive[i][j]->Write();
      h_dphi_eta_inclusive[i][j]->Write();
    }
  } // end of loop over bins
  
  fout->Close();
  
}
