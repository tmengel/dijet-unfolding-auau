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


void preprocess(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string & configfile = "binning_AA.config",  
  const std::string & infile = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root"
)
{
  read_binning rb(configfile.c_str());

  std::string system_string = rb.get_system_string(centrality_bin);
  std::cout << "System string: " << system_string << std::endl;

  const double dphicut = rb.get_dphicut();
  std::cout << "Delta phi cut: " << dphicut << std::endl;
  
  const float vertex_cut = static_cast<float>( rb.get_vtx_cut() );
  const float etacut = rb.get_abs_eta_acceptance( static_cast<float>(cone_size) * 0.1 );
  const float etacut_bkg = rb.get_eta_cut_bkg();
  const float reco_cut = static_cast<float>( rb.get_reco_pt_min_cut() );
  std::cout << "Vertex cut: " << vertex_cut << std::endl;
  std::cout << "Eta cut: " << etacut << std::endl;
  std::cout << "Eta cut bkg: " << etacut_bkg << std::endl;
  std::cout << "Reco cut: " << reco_cut << std::endl;
  
  double flow_v22_scale = rb.get_flow_sys();
  double flow_v33_scale = rb.get_flow_v33_sys();
  const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;

  std::string sys_name = "nominal";
  if (flow_sys){ sys_name = rb.get_flow_systematic_name(); }
    
  const int n_cent_bins = rb.get_number_centrality_bins();  
  const int n_pt_bins = rb.get_nbins();

  float pt_bins[ n_pt_bins + 1 ];
  float xj_bins[ n_pt_bins + 1 ];
  float cent_bins[ n_cent_bins + 1 ];

  rb.get_pt_bins( pt_bins );
  rb.get_xj_bins( xj_bins );
  rb.get_centrality_bins( cent_bins );
  for (int i = 0; i < n_cent_bins + 1; i++)
  {
    std::cout << "Centrality bin " << i << ": " << cent_bins[i] << std::endl;
  }
  for (int i = 0 ; i < n_pt_bins + 1; i++)
  {
    std::cout << pt_bins[i] << " -- " << xj_bins[i] << std::endl;
  }


  const float reco_leading_cut = rb.get_reco_leading_cut();
  const float reco_subleading_cut = rb.get_reco_subleading_cut();
  const int max_reco_bin = rb.get_maximum_reco_bin();
  std::cout << "Max reco: " << max_reco_bin << " - ( " << pt_bins[max_reco_bin] << ")" << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;

  float low_fit = rb.get_flow_fit_low();
  float high_fit = rb.get_flow_fit_high();
  float zyam_integral_low = rb.get_zyam_region_low();
  float zyam_integral_high = rb.get_zyam_region_high();
  float dphicut_z_low = rb.get_zyam_low();
  float dphicut_z_high = rb.get_zyam_high();
  float nbins_diff = (TMath::Pi() - dphicut)/(TMath::Pi()/32.);
  float nbins_zyam = (zyam_integral_high - zyam_integral_low)/(TMath::Pi()/32.);
  std::cout << "Flow fit range: " << low_fit << " -- " << high_fit << std::endl;
  std::cout << "ZYAM integral range: " << zyam_integral_low << " -- " << zyam_integral_high << std::endl;
  std::cout << "ZYAM dphi range: " << dphicut_z_low << " -- " << dphicut_z_high << std::endl;
  std::cout << "Number of bins in dphi difference: " << nbins_diff << std::endl;


  std::string prob_path = Form("%s/unfolding_hists/probability_hists_AA_r0%d.root", rb.get_code_location().c_str(), cone_size);
  auto * fprob = TFile::Open(prob_path.c_str(), "READ");
  if (!fprob || fprob->IsZombie()) 
  {
    std::cerr << "Missing probability-correction file: " << prob_path << std::endl;
    return;
  } 
  auto * h_pt2_prob = (TH1D*) fprob -> Get(Form("h_pt2_bin_log_correction_%d", centrality_bin));
  if (!h_pt2_prob)
  {
    std::cerr << "Missing correction histogram for centrality bin " << centrality_bin << " in " << prob_path << std::endl;
    return;
  }
  h_pt2_prob -> SetDirectory(0);
  fprob -> Close();

  auto * f = new TFile(infile.c_str(), "READ");
  auto * t = (TTree*) f->Get("T");
  if (!t)
  {
    std::cerr << "Could not find tree T in " << infile << std::endl;
    return;
  }
  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_unsub = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;
  int centrality;
  float mbd_vertex_z;
  float sumeT;
  if ( t )
  {
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("cent", 1);
    t->SetBranchStatus("zvrtx", 1);
    t->SetBranchStatus("sumeT", 1);
    t->SetBranchStatus("jet_pT", 1);
    t->SetBranchStatus("jet_unsub_pT", 1);
    t->SetBranchStatus("jet_E", 1);
    t->SetBranchStatus("jet_unsub_E", 1);
    t->SetBranchStatus("jet_eta", 1);
    t->SetBranchStatus("jet_phi", 1);
    t->SetCacheSize(256 * 1024 * 1024);
    t->AddBranchToCache("*", kTRUE);
    t->SetBranchAddress("cent", &centrality);
    t->SetBranchAddress("jet_pT", &reco_jet_pt);
    t->SetBranchAddress("jet_unsub_pT", &reco_jet_pt_unsub);
    t->SetBranchAddress("jet_E", &reco_jet_e);
    t->SetBranchAddress("jet_unsub_E", &reco_jet_e_unsub);
    t->SetBranchAddress("jet_eta", &reco_jet_eta);
    t->SetBranchAddress("jet_phi", &reco_jet_phi);
    t->SetBranchAddress("zvrtx", &mbd_vertex_z);
    t->SetBranchAddress("sumeT", &sumeT);
  }

  int entries = t -> GetEntries();
  std::cout << "Input tree has " << entries << " entries" << std::endl;
  std::cout << "Input tree has no minbias or GL1 trigger branch; assuming the input is preselected." << std::endl;

  TH1D *h_dphi_exclusive_all = new TH1D("h_dphi_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_exclusive_all = new TH1D("h_dphi_eta_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_inclusive_all = new TH1D("h_dphi_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_inclusive_all = new TH1D("h_dphi_eta_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_exclusive[n_pt_bins][n_pt_bins];
  TH1D *h_dphi_eta_exclusive[n_pt_bins][n_pt_bins];

  TH1D *h_dphi_inclusive[n_pt_bins][n_pt_bins];
  TH1D *h_dphi_eta_inclusive[n_pt_bins][n_pt_bins];
  
  for (int i = 0; i < n_pt_bins; i++)
  {
    for (int j = 0; j < n_pt_bins; j++)
    {
      h_dphi_exclusive[i][j] = new TH1D(Form("h_dphi_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      h_dphi_eta_exclusive[i][j] = new TH1D(Form("h_dphi_eta_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      
      h_dphi_inclusive[i][j] = new TH1D(Form("h_dphi_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
      h_dphi_eta_inclusive[i][j] = new TH1D(Form("h_dphi_eta_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
    }
  }

  TH2D *h_pt1_pt2_signal_zyam_exclusive = new TH2D("h_pt1_pt2_signal_zyam_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_signal_exclusive = new TH2D("h_pt1_pt2_signal_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_zyam_exclusive = new TH2D("h_pt1_pt2_zyam_exclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);

  TH2D *h_pt1_pt2_signal_zyam_inclusive = new TH2D("h_pt1_pt2_signal_zyam_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_signal_inclusive = new TH2D("h_pt1_pt2_signal_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_zyam_inclusive = new TH2D("h_pt1_pt2_zyam_inclusive",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);

  // background estimate (formerly getBackground.C), built from the same
  // exclusive leading/subleading pairing as the histograms above
  TH2D *h_pt1_pt2_ZYAM = new TH2D("h_pt1_pt2_ZYAM",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_sub = new TH2D("h_pt1_pt2_sub",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);
  TH2D *h_pt1_pt2_ratio = new TH2D("h_pt1_pt2_ratio",";#it{p}_{T,1}; #it{p}_{T,2};", n_pt_bins, pt_bins, n_pt_bins, pt_bins);

  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
  TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);

  TH1D *h_flat_data_pt1pt2_signal_zyam_exclusive = new TH1D("h_data_flat_pt1pt2_signal_zyam_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);
  TH1D *h_flat_data_pt1pt2_signal_exclusive = new TH1D("h_data_flat_pt1pt2_signal_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);
  TH1D *h_flat_data_pt1pt2_zyam_exclusive = new TH1D("h_data_flat_pt1pt2_zyam_exclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);

  TH1D *h_flat_data_pt1pt2_signal_zyam_inclusive = new TH1D("h_data_flat_pt1pt2_signal_zyam_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);
  TH1D *h_flat_data_pt1pt2_signal_inclusive = new TH1D("h_data_flat_pt1pt2_signal_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);
  TH1D *h_flat_data_pt1pt2_zyam_inclusive = new TH1D("h_data_flat_pt1pt2_zyam_inclusive",";#it{p}_{T,1, smear} + #it{p}_{T,2, smear}", n_pt_bins*n_pt_bins, 0, n_pt_bins*n_pt_bins);
  
  std::vector< std::pair< size_t, float > > accepted_jets;
  
  auto get_pt_bin = [&pt_bins, n_pt_bins](float pt) -> int
  {
    for (int ib = 0; ib < n_pt_bins; ib++)
    {
      if (pt >= pt_bins[ib] && pt < pt_bins[ib+1])
      {
        return ib;
      }
    }
    return n_pt_bins; // return an out-of-range value if not found
  };

  auto calc_dphi = [](float phi1, float phi2) -> float
  {
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > TMath::Pi())
    {
      dphi = 2*TMath::Pi() - dphi;
    }
    return dphi;
  };

  auto calc_abs_deta = [](float eta1, float eta2) -> float
  {
    return std::fabs(eta1 - eta2);
  };

  auto accept_jet_ue_sub = [](float pt_unsub, int centrality) -> bool
  {
    return pt_unsub >= 40.0*TMath::Exp(-0.038*centrality);
  };

  for (int i = 0; i < entries; i++)
  {

    t -> GetEntry( i );
    if (i % int(entries/10) == 0){ std::cout << "Event: " << i << " / " << entries << "\r" << std::flush; }
    if ( fabs(mbd_vertex_z) > vertex_cut ){ continue; }
    if ( centrality < cent_bins[ centrality_bin ] || centrality >= cent_bins[ centrality_bin+1 ] ){ continue; }
    
    // clear the accepted jets for this event
    accepted_jets.clear();

    int n_accepted_jets = 0;
    bool accept_event = true; 

    for ( size_t ijet = 0; ijet  < reco_jet_pt->size(); ijet++ )
    {

      if ( reco_jet_pt->at(ijet) <= reco_cut ){ continue;}
      if ( reco_jet_e->at(ijet) <= 0 ){ continue;}
      if ( reco_jet_e_unsub->at(ijet) <= 0 ){ continue;}
      if ( fabs( reco_jet_eta->at(ijet) ) > etacut ){ continue;}
      
      float pt_unsub = reco_jet_pt_unsub->at(ijet) - reco_jet_pt->at(ijet);
      if ( !accept_jet_ue_sub( pt_unsub, centrality ) )
      {
        accept_event = false;
        continue;
      }

      n_accepted_jets++;

      accepted_jets.push_back( std::make_pair( ijet, reco_jet_pt->at(ijet) ) );
      
      if ( reco_jet_pt->at(ijet) >= pt_bins[max_reco_bin] )
      {
        accept_event = false;
        continue;
      }

    } // end of jet loop

    if ( !accept_event || n_accepted_jets < 2 ){ continue; }

    // find leading 
    std::sort( accepted_jets.begin(), accepted_jets.end(), 
      [](const std::pair<size_t, float> & a, const std::pair<size_t, float> & b) 
    { return a.second > b.second; } );

    // check pt1
    if ( accepted_jets[0].second < reco_leading_cut ){ continue; }

    // fill event level histograms
    h_centrality -> Fill( centrality );
    h_mbd_vertex -> Fill( mbd_vertex_z );
    h_sumeT -> Fill( sumeT );

    int pt1_reco_bin = get_pt_bin( accepted_jets[0].second );

    // only iterate from 1 to n; the first candidate that clears the
    // subleading cut is the exclusive (leading + immediate next-highest)
    // partner, every later one only counts toward the inclusive fills
    bool filled_exclusive_pair = false;
    const float phi1 = reco_jet_phi->at( accepted_jets[0].first );
    const float eta1 = reco_jet_eta->at( accepted_jets[0].first );
    for ( int ijet = 1; ijet < (int)accepted_jets.size(); ijet++ )
    {

      if ( accepted_jets[ijet].second < reco_subleading_cut ){ continue; }

      const bool is_subleading = !filled_exclusive_pair;
      filled_exclusive_pair = true;

      auto dphi = calc_dphi( phi1, reco_jet_phi->at( accepted_jets[ijet].first ) );
      auto deta = calc_abs_deta( eta1, reco_jet_eta->at( accepted_jets[ijet].first ) );

      int pt2_bin = get_pt_bin( accepted_jets[ijet].second );

      bool signal_region = ( dphi >= dphicut );
      bool prep_region  = ( dphi >= dphicut_z_low && dphi < dphicut_z_high );

      if ( is_subleading ) 
      {
      
        h_dphi_exclusive_all -> Fill( dphi );
        h_dphi_exclusive[pt1_reco_bin][pt2_bin] -> Fill( dphi );
      
        if ( deta > etacut_bkg )
        {
          h_dphi_eta_exclusive_all -> Fill( dphi );
          h_dphi_eta_exclusive[pt1_reco_bin][pt2_bin] -> Fill( dphi );
        }

        if ( signal_region )
        {
          h_flat_data_pt1pt2_signal_zyam_exclusive -> Fill( pt1_reco_bin + n_pt_bins*pt2_bin );
          h_flat_data_pt1pt2_signal_zyam_exclusive -> Fill( pt2_bin + n_pt_bins*pt1_reco_bin );
          h_pt1_pt2_signal_zyam_exclusive -> Fill( accepted_jets[0].second, accepted_jets[ijet].second );
          h_pt1_pt2_signal_zyam_exclusive -> Fill( accepted_jets[ijet].second, accepted_jets[0].second );
        }
      
        if ( prep_region )
        {
          h_flat_data_pt1pt2_zyam_exclusive -> Fill( pt1_reco_bin + n_pt_bins*pt2_bin );
          h_flat_data_pt1pt2_zyam_exclusive -> Fill( pt2_bin + n_pt_bins*pt1_reco_bin );
          h_pt1_pt2_zyam_exclusive -> Fill( accepted_jets[0].second, accepted_jets[ijet].second );
          h_pt1_pt2_zyam_exclusive -> Fill( accepted_jets[ijet].second, accepted_jets[0].second );
        }
        
      }  // end of subleading check

      h_dphi_inclusive_all -> Fill( dphi );
      h_dphi_inclusive[pt1_reco_bin][pt2_bin] -> Fill( dphi );
      if ( deta > etacut_bkg )
      {
        h_dphi_eta_inclusive_all -> Fill( dphi );
        h_dphi_eta_inclusive[pt1_reco_bin][pt2_bin] -> Fill( dphi );
      }
      if ( signal_region )
      {
        h_flat_data_pt1pt2_signal_zyam_inclusive -> Fill( pt1_reco_bin + n_pt_bins*pt2_bin );
        h_flat_data_pt1pt2_signal_zyam_inclusive -> Fill( pt2_bin + n_pt_bins*pt1_reco_bin );
        h_pt1_pt2_signal_zyam_inclusive -> Fill( accepted_jets[0].second, accepted_jets[ijet].second );
        h_pt1_pt2_signal_zyam_inclusive -> Fill( accepted_jets[ijet].second, accepted_jets[0].second );
      }
      if ( prep_region )
      {
        h_flat_data_pt1pt2_zyam_inclusive -> Fill( pt1_reco_bin + n_pt_bins*pt2_bin );
        h_flat_data_pt1pt2_zyam_inclusive -> Fill( pt2_bin + n_pt_bins*pt1_reco_bin );
        h_pt1_pt2_zyam_inclusive -> Fill( accepted_jets[0].second, accepted_jets[ijet].second );
        h_pt1_pt2_zyam_inclusive -> Fill( accepted_jets[ijet].second, accepted_jets[0].second );
      }

    } // end of subleading loop
    
	  
  } // end of event loop

  // raw signal-region counts, snapshotted before the ZYAM subtraction below
  // mutates h_pt1_pt2_signal_zyam_exclusive in place
  TH2D *h_pt1_pt2_Signal = (TH2D*) h_pt1_pt2_signal_zyam_exclusive->Clone("h_pt1_pt2_Signal");

  TF1 *fits_o[n_pt_bins][n_pt_bins];
  TF1 *fits_bkg[n_pt_bins][n_pt_bins];
  for (int i = 0; i < n_pt_bins; i++)
  {
    for (int j = 0; j < n_pt_bins; j++)
    {
      fits_o[i][j] = nullptr;
      fits_bkg[i][j] = nullptr;
    }
  }

  for ( int i = 0; i < n_pt_bins; i++ )
  {
    for ( int j = 0; j <= i; j++ )
    {

      std::cout << " ---------- " << i << " - " << j << " -----------" << std::endl;

      const int nbins_over = h_dphi_eta_inclusive[i][j]->FindBin(high_fit) - h_dphi_eta_inclusive[i][j]->FindBin(low_fit);
      const int ncounts    = h_dphi_eta_inclusive[i][j]->Integral(h_dphi_eta_inclusive[i][j]->FindBin(low_fit), h_dphi_eta_inclusive[i][j]->FindBin(high_fit));
      std::cout << ncounts << std::endl;

      TF1 * ffit = new TF1(Form("ffit_o_%d_%d", i, j), "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", low_fit, high_fit);
      ffit -> SetParLimits(0, 0, 10000);
      ffit -> SetParLimits(1, 0, 0.5);
      ffit -> SetParLimits(2, 0, 0.5);
      ffit -> SetParameters(0, 0, 0);
      if (ncounts < nbins_over)
      {
        ffit -> FixParameter(1, 0);
        ffit -> FixParameter(2, 0);
      }

      h_dphi_eta_inclusive[i][j] -> Fit(ffit, "0RlQ", "", low_fit, high_fit);

      if (flow_sys)
      {
        std::cout << sys_name << " flow systematic: scaling fitted (v22,v33) from ("
          << ffit->GetParameter(1) << "," << ffit->GetParameter(2) << ") to ("
          << flow_v22_scale * ffit->GetParameter(1) << ","
          << flow_v33_scale * ffit->GetParameter(2) << ")" << " for bin " << i << "," << j << std::endl;
        ffit -> SetParameter(1, flow_v22_scale * ffit->GetParameter(1));
        ffit -> SetParameter(2, flow_v33_scale * ffit->GetParameter(2));
      }
      fits_o[i][j] = ffit;

      // Renormalizes the flow shape from ffit to a dphi histogram's own ZYAM
      // plateau, then returns the counts (and Poisson error) left over in the
      // signal region once that renormalized flow is subtracted off, plus the
      // pieces (raw signal-region integral, flow-under-signal, flow-under-zyam,
      // and the renormalized TF1 itself) that the background estimate needs.
      struct FlowSignal { float signal, err, dphi_signal_integral, flow_signal, flow_zyam; TF1 * fnorm; };
      auto flow_subtracted_signal = [&](TH1D * h_dphi, const char * fit_name) -> FlowSignal
      {
        const int first_bin  = h_dphi->FindBin(zyam_integral_low);
        const int second_bin = h_dphi->FindBin(zyam_integral_high);
        const float bin_diff = 1.0*(second_bin - first_bin + 1);

        const float dphi_zyam_integral = h_dphi->Integral(first_bin, second_bin)/bin_diff;
        const float dphi_flow_integral = ffit->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high - zyam_integral_low);
        const float p0 = ffit->GetParameter(0) + (dphi_zyam_integral - dphi_flow_integral);

        TF1 * fnorm = new TF1(fit_name, "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", 0, TMath::Pi());
        fnorm->SetParameters(p0, ffit->GetParameter(1), ffit->GetParameter(2));

        const int first_signal_bin  = h_dphi->FindBin(dphicut);
        const int second_signal_bin = -1;

        const float flow_signal          = nbins_diff*fnorm->Integral(dphicut, TMath::Pi())/(TMath::Pi() - dphicut);
        const float flow_zyam            = nbins_zyam*fnorm->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high - zyam_integral_low);
        const float dphi_signal_integral = h_dphi->Integral(first_signal_bin, second_signal_bin);

        const float signal = std::max(0.0f, dphi_signal_integral - flow_signal);
        const float err    = std::sqrt(dphi_signal_integral);
        return { signal, err, dphi_signal_integral, flow_signal, flow_zyam, fnorm };
      };

      const int sbin = std::min(i, j);
      const float prob = h_pt2_prob->GetBinContent(sbin + 1);

      const FlowSignal res_ex = flow_subtracted_signal(h_dphi_exclusive[i][j], Form("ffit_%d_%d", i, j));
      const float signal_ex = res_ex.signal / prob;
      const float err_ex    = res_ex.err / prob;

      const FlowSignal res_in = flow_subtracted_signal(h_dphi_inclusive[i][j], Form("ffit_in_%d_%d", i, j));
      const float signal_in = res_in.signal;
      const float err_in    = res_in.err;
      fits_bkg[i][j] = res_ex.fnorm;

      const int gbin1 = h_pt1_pt2_signal_exclusive->GetBin(i+1, j+1);
      const int gbin2 = h_pt1_pt2_signal_exclusive->GetBin(j+1, i+1);
      h_pt1_pt2_signal_exclusive->SetBinContent(gbin1, signal_ex);
      h_pt1_pt2_signal_exclusive->SetBinError(gbin1, err_ex);
      h_pt1_pt2_signal_exclusive->SetBinContent(gbin2, signal_ex);
      h_pt1_pt2_signal_exclusive->SetBinError(gbin2, err_ex);

      h_pt1_pt2_signal_inclusive->SetBinContent(gbin1, signal_in);
      h_pt1_pt2_signal_inclusive->SetBinError(gbin1, err_in);
      h_pt1_pt2_signal_inclusive->SetBinContent(gbin2, signal_in);
      h_pt1_pt2_signal_inclusive->SetBinError(gbin2, err_in);

      const int ggbin1 = 1 + i + n_pt_bins*j;
      const int ggbin2 = 1 + j + n_pt_bins*i;
      h_flat_data_pt1pt2_signal_exclusive->SetBinContent(ggbin1, signal_ex);
      h_flat_data_pt1pt2_signal_exclusive->SetBinContent(ggbin2, signal_ex);
      h_flat_data_pt1pt2_signal_exclusive->SetBinError(ggbin1, err_ex);
      h_flat_data_pt1pt2_signal_exclusive->SetBinError(ggbin2, err_ex);

      h_flat_data_pt1pt2_signal_inclusive->SetBinContent(ggbin1, signal_in);
      h_flat_data_pt1pt2_signal_inclusive->SetBinContent(ggbin2, signal_in);
      h_flat_data_pt1pt2_signal_inclusive->SetBinError(ggbin1, err_in);
      h_flat_data_pt1pt2_signal_inclusive->SetBinError(ggbin2, err_in);

      // background estimate: flow-subtracted (unfolding-prob-uncorrected)
      // exclusive signal, its background fraction, and the ZYAM-window counts
      // extrapolated into the signal region through the same renormalized fit
      h_pt1_pt2_sub->SetBinContent(gbin1, res_ex.signal);
      h_pt1_pt2_sub->SetBinContent(gbin2, res_ex.signal);

      const float ratio = res_ex.dphi_signal_integral > 0 ? res_ex.signal/res_ex.dphi_signal_integral : 0.0F;
      h_pt1_pt2_ratio->SetBinContent(gbin1, ratio);
      h_pt1_pt2_ratio->SetBinContent(gbin2, ratio);

      const float raw_zyam_counts = h_pt1_pt2_zyam_exclusive->GetBinContent(gbin1);
      const float zyam_estimate = (res_ex.flow_zyam > 0 && res_ex.flow_signal > 0)
        ? raw_zyam_counts * res_ex.flow_signal/res_ex.flow_zyam : 0.0F;
      h_pt1_pt2_ZYAM->SetBinContent(gbin1, zyam_estimate);
      h_pt1_pt2_ZYAM->SetBinContent(gbin2, zyam_estimate);

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
  for (int ix = 0; ix < n_pt_bins; ix++)
  {
    for (int iy = 0; iy < n_pt_bins; iy++)
	  {
      int gbin = 1 + ix + n_pt_bins*iy;

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
  for (int i = 0; i < n_pt_bins; i++)
  {
    for (int j = 0; j < n_pt_bins; j++)
	  {

      h_dphi_exclusive[i][j]->Write();
      h_dphi_eta_exclusive[i][j]->Write();	  
      h_dphi_inclusive[i][j]->Write();
      h_dphi_eta_inclusive[i][j]->Write();
    }
  } // end of loop over bins

  fout->Close();

  std::cout << "writing background" << std::endl;
  TString bkgpath = rb.get_code_location() + "/unfolding_hists/background_hists_AA_cent_" + std::to_string(centrality_bin) + "_r0" + std::to_string(cone_size) + ".root";

  TFile *fout_bkg = new TFile(bkgpath.Data(), "recreate");
  h_pt1_pt2_sub->Write();
  h_pt1_pt2_Signal->Write();
  h_pt1_pt2_ZYAM->Write();
  h_pt1_pt2_ratio->Write();
  for (int i = 0; i < n_pt_bins; i++)
  {
    for (int j = 0; j < n_pt_bins; j++)
    {
      if (fits_o[i][j])   { fits_o[i][j]->Write(); }
      if (fits_bkg[i][j]) { fits_bkg[i][j]->Write(); }

      h_dphi_exclusive[i][j]->Write();
      h_dphi_eta_exclusive[i][j]->Write();
      h_dphi_inclusive[i][j]->Write();
      h_dphi_eta_inclusive[i][j]->Write();
    }
  } // end of loop over bins

  fout_bkg->Close();

}
