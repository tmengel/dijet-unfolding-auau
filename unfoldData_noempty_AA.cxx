
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "TNtuple.h"
#include "TStyle.h"
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int unfoldData_noempty_AA(
  const std::string configfile = "binning_AA.config", 
  const int niterations = 10, 
  const int cone_size = 3, 
  const int centrality_bin = 0, 
  const int primer = 0
)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  int preloaded = 1;
  
  bool ispp = (centrality_bin < 0);

  read_binning rb(configfile.c_str());
  std::string system_string = rb.get_system_string(centrality_bin);
  std::string data_file = std::getenv("TNUPLE_DATA_FILE");
  std::cout << "Using data file: " << data_file << std::endl;

  float mbd_vertex;
  float pt1_reco;
  float pt2_reco;
  float dphi_reco;
  float nrecojets;
  float trigger;
  float centrality;
  TFile *fin = nullptr;
  TNtuple *tn = nullptr;
  if (!preloaded)
  {
    fin = new TFile(data_file.c_str(), "r");
    tn = (TNtuple*) fin->Get("tn_dijet");
    if (!tn)
    {
      std::cerr << "Missing tn_dijet in " << data_file << std::endl;
      return 1;
    }
    tn->SetBranchAddress("pt1_reco", &pt1_reco);
    tn->SetBranchAddress("pt2_reco", &pt2_reco);
    tn->SetBranchAddress("dphi_reco", &dphi_reco);
    tn->SetBranchAddress("njets", &nrecojets);
    tn->SetBranchAddress("centrality", &centrality);
    tn->SetBranchAddress("trigger", &trigger);
    tn->SetBranchAddress("mbd_vertex", &mbd_vertex);
  }

  const int minentries = rb.get_minentries();
	std::cout << "Minimum entries: " << minentries << std::endl;

	const int nbins = rb.get_nbins();
	float ipt_bins[nbins+1];
	float ixj_bins[nbins+1];
	rb.get_pt_bins(ipt_bins);
	rb.get_xj_bins(ixj_bins);
	std::cout << "Pt bins -- Xj bins: " << std::endl;
	for (int i = 0 ; i < nbins + 1; i++)
	{
		std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}

	const int centrality_bins = rb.get_number_centrality_bins();
	float icentrality_bins[centrality_bins + 1];
	rb.get_centrality_bins(icentrality_bins);
	std::cout << "Centrality bins: " << std::endl;
	for (int i = 0 ; i < centrality_bins + 1; i++)
	{
		std::cout << icentrality_bins[i] << std::endl;
	}

	const int max_reco_bin = rb.get_maximum_reco_bin();
	std::cout << "Max reco bin: " << max_reco_bin << std::endl;

	const double dphicut = rb.get_dphicut();
	const double dphicuttruth = dphicut;
	std::cout << "dphicut: " << dphicut << std::endl;
	std::cout << "dphicuttruth: " << dphicuttruth << std::endl;
  
  const float dphicut_z_low = rb.get_zyam_low();
  const float dphicut_z_high = rb.get_zyam_high();
  std::cout << "ZYAM low: " << dphicut_z_low << std::endl;
  std::cout << "ZYAM high: " << dphicut_z_high << std::endl;

	const float truth_leading_cut = rb.get_truth_leading_cut();
	const float truth_subleading_cut = rb.get_truth_subleading_cut();
	std::cout << "Truth leading cut: " << truth_leading_cut << std::endl;
	std::cout << "Truth subleading cut: " << truth_subleading_cut << std::endl;

	const float reco_leading_cut = rb.get_reco_leading_cut();
	const float reco_subleading_cut = rb.get_reco_subleading_cut();
	std::cout << "Reco leading cut: " << reco_leading_cut << std::endl;
	std::cout << "Reco subleading cut: " << reco_subleading_cut << std::endl;

	const float measure_leading_cut = rb.get_measure_leading_cut();
	const float measure_subleading_cut = rb.get_measure_subleading_cut();
	std::cout << "Measure leading cut: " << measure_leading_cut << std::endl;
	std::cout << "Measure subleading cut: " << measure_subleading_cut << std::endl;

	const int measure_leading_bin = rb.get_measure_leading_bin();
	const int measure_subleading_bin = rb.get_measure_subleading_bin();
	std::cout << "Measure leading bin: " << measure_leading_bin << std::endl;
	std::cout << "Measure subleading bin: " << measure_subleading_bin << std::endl;

	float sample_boundary[4] = {0};
	std::cout << "Sample boundaries: " << std::endl;
	for (int ib = 0; ib < 4; ib++)
	{
		sample_boundary[ib] = rb.get_sample_boundary(ib);
		std::cout <<  sample_boundary[ib] << std::endl;
	}


	const double flow_v22_scale = rb.get_flow_sys();
	const double flow_v33_scale = rb.get_flow_v33_sys();
	const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
	
	const int zyam_sys = rb.get_zyam_sys();
	const int inclusive_sys = rb.get_inclusive_sys();
	const int flavor_sys = rb.get_flavor_sys();
	const double JES_sys = rb.get_jes_sys();
	const double JER_sys = rb.get_jer_sys();
	const bool prior_sys = rb.get_prior_sys();

	std::string sys_name = "nominal";
	std::cout << "JES = " << JES_sys << std::endl;
	std::cout << "JER = " << JER_sys << std::endl;
	
	int using_sys = 0;

	if ( prior_sys )
	{
		using_sys = 1;
		sys_name = "PRIOR";
	}
	if ( flow_sys )
	{
		using_sys = 1;
		sys_name = rb.get_flow_systematic_name();
	}
	if (inclusive_sys)
	{
		using_sys = 1;
		sys_name = "INCLUSIVE";
	}
	if (flavor_sys == 1)
	{
		using_sys = 1;
		sys_name = "QQ";
	}
	else if (flavor_sys == 2)
	{
		using_sys = 1;
		sys_name = "QGGG";
	}
	if (JER_sys != 0)
	{
		using_sys = 1;
		if (JER_sys < 0) sys_name = "negJER";
		else if (JER_sys > 0) sys_name = "posJER";
		std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
	}
	if (JES_sys != 0)
	{
		using_sys = 1;
		if (JES_sys < 0) sys_name = "negJES";
		else if (JES_sys > 0) sys_name = "posJES";
		std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
	}

  bool draw_reference_comparisons = true;
  if (const char *draw_ref_env = std::getenv("DIJET_DRAW_REFERENCES"))
  {
    draw_reference_comparisons = std::atoi(draw_ref_env) != 0;
  }
    
  double scale_to_signal[nbins][nbins];
  for (int i = 0; i < nbins; i++)
  {
    for (int j = 0; j < nbins; j++)
	  {
	    scale_to_signal[i][j] = 1;
	  }
  }

  auto * fprob = new TFile(Form("%s/unfolding_hists/probability_hists_AA_r0%d.root", rb.get_code_location().c_str(), cone_size),"read");
  if ( !fprob || fprob->IsZombie() )
  {
    std::cerr << "Missing probability file: " << fprob->GetName() << std::endl;
    return 1;
  }
  auto * h_pt2_prob = (TH1D*) fprob->Get(Form("h_pt2_bin_log_correction_%d", centrality_bin)) -> Clone(Form("h_pt2_bin_log_correction_%d", centrality_bin));
  h_pt2_prob -> SetDirectory(0);
  fprob -> Close();

  auto * fbkg = new TFile(Form("%s/unfolding_hists/background_hists_%s_r0%d.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size),"read");   
  if ( !fbkg || fbkg->IsZombie() )
  {
    std::cerr << "Missing background file: " << fbkg->GetName() << std::endl;
    return 1;
  }

  for (int i = 0; i < nbins; i++)
  {
    for (int j = 0 ; j < nbins; j++)
    {
      if (j > i) continue;
      auto * background_fits = (TF1*) fbkg->Get(Form("ffit_%d_%d", i , j));
      double flow_signal = background_fits->Integral(dphicut, TMath::Pi())/(TMath::Pi()/32.);
      std::cout << "ij : " << i << " / " << j << " " << flow_signal << std::endl;
      if (i == j)
      {
        scale_to_signal[i][j] = 2*flow_signal;
      }
      else
      {
        scale_to_signal[j][i] = flow_signal;
        scale_to_signal[i][j] = flow_signal;
      }
      
    }
  }

  fbkg -> Close();
  
  
  
  // TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);

  // if (primer)
  // {
  //   responsepath += "_PRIMER" + std::to_string(primer);
  // }
  // responsepath += "_" + sys_name;
  
  // responsepath += ".root";
  std::string responsepath = rb.get_code_location() + "/response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
  if (primer)
  {
    responsepath += "_PRIMER" + std::to_string(primer);
  }
  responsepath += "_" + sys_name + ".root";

  auto * fresponse = new TFile(responsepath.c_str(),"READ");
  if ( !fresponse || fresponse->IsZombie() )
  {
    std::cerr << "Missing response file: " << fresponse->GetName() << std::endl;
    return 1;
  }

  auto * rooResponse = (RooUnfoldResponse*) fresponse -> Get("response_noempty");
  if (!rooResponse)
  {
    std::cout << "no repsonse" << std::endl;
    return 1;
  }

  auto * h_flat_truth_pt1pt2 = (TH1D*) fresponse->Get("h_truth_flat_pt1pt2");
  if (!h_flat_truth_pt1pt2)
  {
    std::cout << "no truth" << std::endl;
    return 1;
  }

  auto * h_flat_reco_pt1pt2 = (TH1D*) fresponse->Get("h_reco_flat_pt1pt2");
  if (!h_flat_reco_pt1pt2)
  {
    std::cout << "no reco" << std::endl;
    return 1;
  }

  auto * h_flat_truth_mapping = (TH1D*) fresponse->Get("h_flat_truth_mapping");
  if (!h_flat_truth_mapping)
  {
    std::cout << "no truth" << std::endl;
    return 1;
  }

  auto * h_flat_reco_mapping = (TH1D*) fresponse->Get("h_flat_reco_mapping");
  if (!h_flat_reco_mapping)
  {
    std::cout << "no reco" << std::endl;
    return 1;
  }

  auto * h_flat_truth_skim = (TH1D*) fresponse->Get("h_flat_truth_skim");
  if (!h_flat_truth_skim)
  {
    std::cout << "no truth" << std::endl;
    return 1;
  }

  auto * h_flat_reco_skim = (TH1D*) fresponse->Get("h_flat_reco_skim");
  if (!h_flat_reco_skim)
  {
    std::cout << "no reco" << std::endl;
    return 1;
  }



  TString unfolding_path = rb.get_code_location() + "/unfolding_hists/unfolding_hists_preload_" + system_string + "_r0" + std::to_string(cone_size);

  if (flow_sys) unfolding_path += "_" + sys_name + ".root";
  else unfolding_path += "_nominal.root";
  
  TFile *funin = new TFile(unfolding_path.Data(), "r");
  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";x_{J};1/N", nbins, ixj_bins);

  TH1D *h_dphi_reco = new TH1D("h_reco_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  TH1D *h_dphi_reco_ZYAM = new TH1D("h_reco_ZYAM_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());
  TH1D *h_dphi_reco_SIGNAL = new TH1D("h_reco_SIGNAL_dphi",";#Delta#phi;Counts", 32, 0, TMath::Pi());

  TH1D *h_dphi_reco_sub[4];
  TH1D *h_dphi_reco_ZYAM_sub[4];
  TH1D *h_dphi_reco_SIGNAL_sub[4];
  for (int i = 0; i < 4; i++)
  {
    h_dphi_reco_sub[i] = new TH1D(Form("h_reco_dphi_sub_%d", i),";#Delta#phi;Counts", 32, 0, TMath::Pi());
    h_dphi_reco_ZYAM_sub[i] = new TH1D(Form("h_reco_ZYAM_dphi_sub_%d", i),";#Delta#phi;Counts", 32, 0, TMath::Pi());
    h_dphi_reco_SIGNAL_sub[i] = new TH1D(Form("h_reco_SIGNAL_dphi_sub_%d", i),";#Delta#phi;Counts", 32, 0, TMath::Pi());
  }
  
  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T1, data};p_{T2, data}", nbins, ipt_bins, nbins, ipt_bins);
  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T1, smear} + p_{T2, smear}", nbins*nbins, 0, nbins*nbins);

  TH2D *h_pt1pt2_ZYAM = new TH2D("h_pt1pt2_ZYAM",";p_{T1, ZYAM};p_{T2, ZYAM}", nbins, ipt_bins, nbins, ipt_bins);
  TH1D *h_flat_data_pt1pt2_ZYAM = new TH1D("h_data_flat_pt1pt2_ZYAM",";p_{T1, smear} + p_{T2, smear}", nbins*nbins, 0, nbins*nbins);

  if (preloaded)
  {
    // The empty TH1Ds created above were registered in funin's in-memory directory
    // (funin was gDirectory) and shadow the on-disk histograms of the same name,
    // so funin->Get() would return the empty ones. Delete them first.
    delete h_mbd_vertex;
    delete h_centrality;
    h_mbd_vertex = (TH1D*) funin->Get("h_mbd_vertex");
    h_centrality = (TH1D*) funin->Get("h_centrality");

    string bkg_type = (zyam_sys?"zyam_":"");
    string sig_type = (inclusive_sys?"inclusive":"exclusive");
    h_pt1pt2 = (TH2D*) funin->Get(Form("h_pt1_pt2_signal_%s%s", bkg_type.c_str(), sig_type.c_str()));
    h_flat_data_pt1pt2 = (TH1D*) funin->Get(Form("h_data_flat_pt1pt2_signal_%s%s", bkg_type.c_str(), sig_type.c_str()));
    h_pt1pt2->SetName("h_pt1pt2");
    h_flat_data_pt1pt2->SetName("h_data_flat_pt1pt2");
  }

  int nbin_response = nbins*nbins;
  
  int entries = preloaded ? 0 : tn->GetEntries();
  if (!preloaded)
  {
    for (int i = 0; i < entries; i++)
    {
      tn->GetEntry(i);

      if (!ispp && centrality < icentrality_bins[centrality_bin] || centrality >= icentrality_bins[centrality_bin + 1]) continue;

      float maxi = std::max(pt1_reco, pt2_reco);
      float mini = std::min(pt1_reco, pt2_reco);

      float pt1_reco_bin = nbins;
      float pt2_reco_bin = nbins;

      float es1 = pt1_reco;
      float es2 = pt2_reco;
      if (es1 >= ipt_bins[max_reco_bin] ) continue;
      for (int ib = 0; ib < nbins; ib++)
        {

          if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
      {
        pt1_reco_bin = ib;
      }
          if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
      {
        pt2_reco_bin = ib;
      }
        }
      
    
      bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco >= dphicut);
      bool fill_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut);
      bool ZYAM_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco >= dphicut_z_low && dphi_reco < dphicut_z_high);

      if (fill_good)
        {
          h_dphi_reco->Fill(dphi_reco);
          if (maxi > 30)
      {
        if (pt2_reco_bin == 1)
          {
            h_dphi_reco_sub[0]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 2)
          {
            h_dphi_reco_sub[1]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 3 || pt2_reco_bin == 4)
          {
            h_dphi_reco_sub[2]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 5 || pt2_reco_bin == 6)
          {
            h_dphi_reco_sub[3]->Fill(dphi_reco);
          }
      }
        }
        
      if (reco_good)
        {
          h_pt1pt2->Fill(es1, es2);
          h_pt1pt2->Fill(es2, es1);
          h_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
          h_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
          h_reco_xj->Fill(mini/maxi);
          h_mbd_vertex->Fill(mbd_vertex);
          if (!ispp) h_centrality->Fill(centrality);
          h_dphi_reco_SIGNAL->Fill(dphi_reco);

          if (maxi >= 30)
      {
        if (pt2_reco_bin == 1)
          {
            h_dphi_reco_SIGNAL_sub[0]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 2)
          {
            h_dphi_reco_SIGNAL_sub[1]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 3 || pt2_reco_bin == 4)
          {
            h_dphi_reco_SIGNAL_sub[2]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 5 || pt2_reco_bin == 6)
          {
            h_dphi_reco_SIGNAL_sub[3]->Fill(dphi_reco);
          }
      }

        }
      else if (ZYAM_good)
        {

          h_pt1pt2_ZYAM->Fill(es1, es2);
          h_pt1pt2_ZYAM->Fill(es2, es1);
          h_flat_data_pt1pt2_ZYAM->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
          h_flat_data_pt1pt2_ZYAM->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
          h_dphi_reco_ZYAM->Fill(dphi_reco);
          if (maxi > 30)
      {
        if (pt2_reco_bin == 1)
          {
            h_dphi_reco_ZYAM_sub[0]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 2)
          {
            h_dphi_reco_ZYAM_sub[1]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 3 || pt2_reco_bin == 4)
          {
            h_dphi_reco_ZYAM_sub[2]->Fill(dphi_reco);
          }
        if (pt2_reco_bin == 5 || pt2_reco_bin == 6)
          {
            h_dphi_reco_ZYAM_sub[3]->Fill(dphi_reco);
          }
      }
        }
    }

    TH2D *h_pt1pt2_raw = (TH2D*) h_pt1pt2->Clone();
    h_pt1pt2_raw->SetName("h_pt1pt2_raw");

    TH1D *h_flat_data_pt1pt2_raw = (TH1D*) h_flat_data_pt1pt2->Clone();
    h_flat_data_pt1pt2_raw->SetName("h_flat_data_pt1pt2_raw");

    for (int ix = 0; ix < h_pt1pt2->GetXaxis()->GetNbins(); ix++)
    {
      for (int iy = 0; iy < h_pt1pt2->GetXaxis()->GetNbins(); iy++)
      {
        int xybin = h_pt1pt2->GetBin(ix+1, iy+1);

        int sbin = std::min(ix, iy);
          
        float n_bkg = h_pt1pt2->GetBinContent(xybin);
        float nerr = h_pt1pt2->GetBinError(xybin);

        if (n_bkg < scale_to_signal[ix][iy]) h_pt1pt2->SetBinContent(xybin, 0);
        else
        {
          float nen = n_bkg - scale_to_signal[ix][iy];
          nen /= h_pt2_prob->GetBinContent(sbin + 1);
          nerr /= h_pt2_prob->GetBinContent(sbin + 1);

          h_pt1pt2->SetBinContent(xybin, nen);
          h_pt1pt2->SetBinError(xybin, nerr);
        }
      }
    }
    
    for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
      {
        int xybin = 1 + ix + nbins*iy;
        int sbin = std::min(ix, iy);
          
        float nen = h_flat_data_pt1pt2->GetBinContent(xybin);
        float nerr = h_flat_data_pt1pt2->GetBinError(xybin);

        if (nen < scale_to_signal[ix][iy]) h_flat_data_pt1pt2->SetBinContent(xybin, 0);
        else
        {
          nen -= scale_to_signal[ix][iy];
          nen /= h_pt2_prob->GetBinContent(sbin + 1);
          nerr /= h_pt2_prob->GetBinContent(sbin + 1);
          h_flat_data_pt1pt2->SetBinContent(xybin, nen);
          h_flat_data_pt1pt2->SetBinError(xybin, nerr);
        }
      }
    }
	  
  }
 
  //  h_flat_data_pt1pt2->Add(h_flat_data_pt1pt2_ZYAM, -1);
  h_flat_data_pt1pt2->Scale(.5);

  TH1D *h_flat_data_skim = (TH1D*) h_flat_reco_skim->Clone();
  h_flat_data_skim->SetName("h_flat_data_skim");
  h_flat_data_skim->Reset();
  histo_opps::skim_down_histo(h_flat_data_skim, h_flat_data_pt1pt2, h_flat_reco_mapping); 
  
  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 1;
  for (int iter = 0; iter < niterations; iter++ )
  {
    
    // handleFakes=true -- see createResponse_exclusive_AA.cxx for why this
    // must be explicit: RooUnfoldBayes defaults it to false and silently
    // unfolds fake-jet content in the response as signal otherwise, instead
    // of excluding it as intended.
    RooUnfoldBayes   unfold (rooResponse, h_flat_data_skim, iter + 1, false, true);    // OR
    h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
    std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
    h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
    h_flat_unfold_pt1pt2[iter]->Reset();
    h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
    histo_opps::fill_up_histo(h_flat_unfold_skim[iter], h_flat_unfold_pt1pt2[iter], h_flat_truth_mapping);
  }

  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  for (int iter = 0; iter < niterations; iter++)
  {
    dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
    dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
  }
  
  dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
  dlutility::SetLineAtt(h_flat_data_pt1pt2, kBlue, 2, 1);

  h_flat_truth_pt1pt2->Draw("hist");
  h_flat_data_pt1pt2->Draw("hist same");
  h_flat_unfold_pt1pt2[niter]->Draw("same p");

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
  {
    h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
    h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
  }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
  {
    h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
  }

  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  h_pt1pt2_data->SetTitle(";Data p_{T, 1} [GeV]; Data p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");
      

  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
  {
    histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  }

  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_reco_xj, nbins);
  for (int iter = 0; iter < niterations; iter++)
  {
    histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
  }



  
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 1, 24);


  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  //h_xj_truth->GetXaxis()->SetRangeUser(0.3, 1.001);
  TH1D *h_xj_frame = draw_reference_comparisons ? h_xj_truth : h_xj_data;
  h_xj_frame->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_frame->SetMaximum(6);
  h_xj_frame->SetMinimum(0);
  h_xj_frame->Draw("p");
  h_xj_unfold[niter]->Draw("same p");
  h_xj_data->Draw("hist same");
  if (draw_reference_comparisons)
  {
    h_reco_xj->Draw("hist same");
    h_xj_truth->Draw("same hist");
  }
  h_xj_data->Draw("p same");
  h_xj_unfold[niter]->Draw("same hist");
  h_xj_unfold[niter]->Draw("same p");
  if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  else dlutility::DrawSPHENIXpp(0.22, 0.85);
  //dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 7/8 #pi", 0.22, 0.59);
  // dlutility::drawText("\\mathscr{L} = 25.7 pb^{-1}", 0.22, 0.54);

  TLegend *leg = new TLegend(0.22, 0.35, 0.4, 0.50);
  leg->SetLineWidth(0);
  leg->AddEntry(h_xj_data, "Data");
  if (draw_reference_comparisons) leg->AddEntry(h_xj_truth, "Pythia8");
  leg->AddEntry(h_xj_unfold[niter], "Unfolded");
  leg->Draw("same");
  dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
    
  cxj->cd(2);

  TH1D *h_data_compare = (TH1D*) h_xj_unfold[niter]->Clone();
  h_data_compare->Divide(draw_reference_comparisons ? h_xj_truth : h_xj_data);
  //h_data_compare->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_data_compare->SetTitle(draw_reference_comparisons ? ";x_{J}; Unfold / Truth" : ";x_{J}; Unfold / Data");
  dlutility::SetFont(h_data_compare, 42, 0.06);
  dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);

 
  h_data_compare->Draw("p");
  TLine *line = new TLine(0.3, 1, 1, 1);
  line->SetLineStyle(4);
  line->SetLineColor(kRed + 3);
  line->SetLineWidth(2);
  line->Draw("same");
  cxj->Print(
    Form("%s/unfolding_plots/unfolded_xj_%s_r%02d_%s%s.pdf",
    rb.get_code_location().c_str(), system_string.c_str(), cone_size,
    sys_name.c_str(), draw_reference_comparisons ? "" : "_no_ref")
  );

  TCanvas *c_iter = new TCanvas("c_iter", "c_iter");
  TH1D *h_closure[niterations];

  for (int iter = 0; iter < niterations; iter++)
  {
    h_closure[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
    h_closure[iter]->SetName(Form("h_closure_%d", iter));
  
    h_closure[iter]->Divide(h_xj_truth);
    h_closure[iter]->SetTitle(";x_{J}; Unfold / Truth");
    dlutility::SetFont(h_closure[iter], 42, 0.06);
  }
  int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 

  for (int iter = 0; iter < 5; iter++)
  {
    dlutility::SetLineAtt(h_closure[1 + 2*iter], colors[iter], 1,1);
    dlutility::SetMarkerAtt(h_closure[1 + 2*iter], colors[iter], 1,8);
  }
  h_closure[1]->Draw("p");
  h_closure[3]->Draw("p same");
  h_closure[5]->Draw("p same");
  h_closure[7]->Draw("p same");
  h_closure[9]->Draw("p same");
  TLine *line2 = new TLine(0.1, 1, 1, 1);
  line2->SetLineStyle(4);
  line2->SetLineColor(kRed + 3);
  line2->SetLineWidth(2);
  line2->Draw("same");



  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlack, 1, 24);
  
  TCanvas *cproj = new TCanvas("cproj","cproj", 500, 700);
  dlutility::ratioPanelCanvas(cproj);

  cproj->cd(1);
  h_xj_data->SetMaximum(3);
  h_xj_data->SetTitle(";x_{J}; #frac{1}{N_{pairs}}#frac{dN_{pair}}{dx_{J}}");
  dlutility::SetFont(h_xj_data, 42, 0.04);
  
  h_xj_data->Draw();
  h_reco_xj->Draw("same");
  if (!ispp) dlutility::DrawSPHENIX(0.22, 0.85);
  else dlutility::DrawSPHENIXpp(0.22, 0.85);
  //  dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T2}^{lead} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
  TLegend *legp = new TLegend(0.22, 0.45, 0.4, 0.55);
  legp->SetLineWidth(0);
  legp->SetTextFont(42);
  legp->SetTextSize(0.04);
  legp->AddEntry(h_reco_xj, "Data Filled");
  legp->AddEntry(h_xj_data, "Data Projected");
  legp->Draw("same");

  cproj->cd(2);
  TH1D *h_fillproj_compare = (TH1D*) h_xj_data->Clone();

  dlutility::SetFont(h_fillproj_compare, 42, 0.06);
  h_fillproj_compare->Divide(h_reco_xj);

  h_fillproj_compare->SetMaximum(1.2);
  
  h_fillproj_compare->SetMinimum(0.8);
  h_fillproj_compare->SetTitle(";x_{J};Projected / Filled");
  dlutility::SetLineAtt(h_fillproj_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_fillproj_compare, kBlack, 1, 8);

  h_fillproj_compare->Draw();
  TLine *line3 = new TLine(0.1, 1, 1, 1);
  line3->SetLineStyle(4);
  line3->SetLineColor(kRed + 3);
  line3->SetLineWidth(2);
  line3->Draw("same");
  cproj->Print(Form("%s/unfolding_plots/proj_compar_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
  
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/unfolding_hists_" + system_string + "_r0" + std::to_string(cone_size);

  if (primer)
  {
    unfoldpath += "_PRIMER" + std::to_string(primer);
  }
  unfoldpath += "_" + sys_name;
  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");
  
  h_flat_data_skim->Write();
  h_flat_data_pt1pt2->Write();
  //h_flat_data_pt1pt2_raw->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_mbd_vertex->Write();
  h_centrality->Write();
  
  for (int iter = 0; iter < niterations; ++iter)
  {
    h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
    h_flat_unfold_pt1pt2[iter]->Write();
  }
  fout->Close();
  

  return 0;
  
}
