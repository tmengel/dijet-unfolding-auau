#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>

#include "TNtuple.h"
#include "TRandom.h"
#include "TStyle.h"

using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"


// int createResponse_noempty_AA (
// 	const std::string configfile = "binning.config", 
// 	const int full_or_half = 0, 
// 	const int niterations = 10, 
// 	const int cone_size = 4, 
// 	const int centrality_bin = 0, 
// 	const int primer = 0,
// 	const bool OPTION_A_SETTINGS = false
// )
// {

// 	gStyle->SetOptStat(0);
// 	dlutility::SetyjPadStyle();
// 	if (full_or_half < 0 || full_or_half > 2)
// 	{
// 		std::cerr << "Closure mode must be 0 (analysis), 1 (half closure), or 2 (full closure)"	<< std::endl;
// 		return 1;
// 	}
	
// 	const bool half_closure = (full_or_half == 1);
// 	const bool full_closure = (full_or_half == 2);
// 	const bool full_or_half_closure = (full_or_half != 0);

// 	bool ispp = (centrality_bin < 0);
// 	read_binning rb(configfile.c_str());

// 	std::string system_string = (ispp ? "pp" : "AA_cent_" + std::to_string(centrality_bin));
// 	std::string j10_file = std::getenv("TNUPLE_SIM_FILE_JET10");
// 	std::string j20_file = std::getenv("TNUPLE_SIM_FILE_JET20");
// 	std::string j30_file = std::getenv("TNUPLE_SIM_FILE_JET30");
// 	std::cout << "Using matched simulation files:" << std::endl;
// 	std::cout << "  Jet10: " << j10_file << std::endl;
// 	std::cout << "  Jet20: " << j20_file << std::endl;
// 	std::cout << "  Jet30: " << j30_file << std::endl;

// 	float maxpttruth[3];
// 	float pt1_truth[3];
// 	float pt2_truth[3];
// 	float dphi_truth[3];
// 	float pt1_reco[3];
// 	float pt2_reco[3];
// 	float nrecojets[3];
// 	float dphi_reco[3];
// 	float match[3];
// 	float mbd_vertex[3];
// 	float centrality[3];
// 	float sumeT[3] = {0, 0, 0};

// 	float n_events[3];
// 	float b_n_events = 0;
// 	float event_weight[3] = {1, 1, 1};
// 	bool has_event_weight[3] = {false, false, false};
// 	bool has_sumeT[3] = {false, false, false};
	
// 	TFile *fin[3];
// 	fin[0] = new TFile(j10_file.c_str(), "r");
// 	fin[1] = new TFile(j20_file.c_str(), "r");
// 	fin[2] = new TFile(j30_file.c_str(), "r");
// 	TNtuple *tn[3];
// 	for (int i = 0 ; i < 3; i++)
// 	{
// 		if (!fin[i] || fin[i]->IsZombie())
// 		{
// 			std::cerr << "Cannot open matched simulation file " << fin[i]->GetName() << std::endl;
// 			return 1;
// 		}

// 		tn[i] = (TNtuple*) fin[i]->Get("tn_match");
// 		if (!tn[i])
// 		{
// 			std::cerr << "Missing tn_match in " << fin[i]->GetName() << std::endl;
// 			return 1;
// 		}

// 		tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
// 		tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
// 		tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
// 		tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
// 		tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
// 		tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
// 		tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
// 		tn[i]->SetBranchAddress("nrecojets", &nrecojets[i]);
// 		tn[i]->SetBranchAddress("matched", &match[i]);
// 		tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
// 		tn[i]->SetBranchAddress("centrality", &centrality[i]);
// 		has_sumeT[i] = tn[i]->GetBranch("sumeT") != nullptr;
// 		if (has_sumeT[i]) tn[i]->SetBranchAddress("sumeT", &sumeT[i]);
// 		else if (!ispp) std::cerr << "Warning: missing sumeT in " << fin[i]->GetName() << "; sumeT reweighting is disabled for this sample." << std::endl;
		
// 		has_event_weight[i] = tn[i]->GetBranch("weight") != nullptr;
// 		if (has_event_weight[i]) tn[i]->SetBranchAddress("weight", &event_weight[i]);
		

// 		TNtuple *tn_stats = (TNtuple*) fin[i]->Get("tn_stats");
// 		if (!tn_stats)
// 		{
// 			std::cerr << "Missing tn_stats in " << fin[i]->GetName() << std::endl;
// 			return 1;
// 		}
// 		tn_stats->SetBranchAddress("nevents", &b_n_events);
// 		tn_stats->GetEntry(0);
// 		n_events[i] = b_n_events;

// 	}
// 	std::cout << "n_events: " << n_events[0] << " " << n_events[1] << " " << n_events[2] << std::endl;
// 	std::cout << "event_weight: " << event_weight[0] << " " << event_weight[1] << " " << event_weight[2] << std::endl;
// 	std::cout << "has_event_weight: " << has_event_weight[0] << " " << has_event_weight[1] << " " << has_event_weight[2] << std::endl;
// 	const float cs_10 = 0.000003997;
//   	const float cs_20 = 6.218e-8;
//   	const float cs_30 = 2.505e-9;
  
// 	float scale_factor[3];
// 	scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
// 	scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30; 
// 	scale_factor[2] = 1;

// 	Int_t minentries = rb.get_minentries();
// 	Int_t read_nbins = rb.get_nbins();
// 	//Int_t primer = rb.get_primer();

// 	Double_t dphicut = rb.get_dphicut();
// 	Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

// 	TF1 *fgaus = new TF1("fgaus", "gaus");
// 	fgaus->SetRange(-0.5, 0.5);
// 	Int_t zyam_sys = rb.get_zyam_sys();
// 	Double_t flow_v22_scale = rb.get_flow_sys();
// 	Double_t flow_v33_scale = rb.get_flow_v33_sys();
// 	const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
// 	Int_t inclusive_sys = rb.get_inclusive_sys();
// 	Double_t JES_sys = rb.get_jes_sys();
// 	Double_t JER_sys = rb.get_jer_sys();
// 	Int_t prior_sys = rb.get_prior_sys();
// 	int using_sys = 0;

// 	TF1 * f_smear = nullptr;
// 	if (!ispp)
// 	{
// 		f_smear  = (TF1*) rb.get_smear_function(centrality_bin);
// 		if (!f_smear)
// 		{
// 			std::cout << " NO SMEAR " << std::endl;
// 			exit(-1);
// 		}
// 	}

// 	std::string sys_name = "nominal";
// 	std::cout << "JES = " << JES_sys << std::endl;
// 	std::cout << "JER = " << JER_sys << std::endl;
// 	float width = 0.8 + JER_sys;
// 	fgaus->SetParameters(1, 0, width);

// 	if (prior_sys)
// 	{
// 		using_sys = 1;
// 		sys_name = "PRIOR";
// 	}
// 	if (zyam_sys)
// 	{
// 		using_sys = 1;
// 		sys_name = "ZYAM";
// 	}
// 	if (flow_sys)
// 	{
// 		using_sys = 1;
// 		sys_name = rb.get_flow_systematic_name();
// 	}
// 	if (inclusive_sys)
// 	{
// 		using_sys = 1;
// 		sys_name = "INCLUSIVE";
// 	}
// 	if (JER_sys != 0)
// 	{
// 		using_sys = 1;
// 		if (JER_sys < 0) sys_name = "negJER";
// 		else if (JER_sys > 0) sys_name = "posJER";
// 		std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
// 	}
// 	if (JES_sys != 0)
// 	{
// 		using_sys = 1;
// 		if (JES_sys < 0) sys_name = "negJES";
// 		else if (JES_sys > 0) sys_name = "posJES";
// 		std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
// 	}

//   	const std::string diagnostic_name = half_closure ? "HALF_" + sys_name : (full_closure ? "FULL_" + sys_name : sys_name);

// 	const int nbins = read_nbins;
// 	float ipt_bins[nbins+1];
// 	float ixj_bins[nbins+1];

// 	rb.get_pt_bins(ipt_bins);
// 	rb.get_xj_bins(ixj_bins);
// 	for (int i = 0 ; i < nbins + 1; i++)
// 	{
// 		std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
// 	}

// 	Int_t max_reco_bin = rb.get_maximum_reco_bin();

// 	const int centrality_bins = rb.get_number_centrality_bins();
// 	float icentrality_bins[centrality_bins + 1];
// 	rb.get_centrality_bins(icentrality_bins);
	
// 	int prior_iteration = 1;
  
// 	TH1D * h_centrality_reweight = nullptr;
// 	TH1D * h_mbd_reweight = nullptr;
// 	TH1D * h_sumeT_reweight = nullptr;
	
// 	auto * fcent = new TFile(Form("centrality/centrality_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
// 	if ( !fcent || fcent->IsZombie() )
// 	{
// 		std::cerr << "Missing required centrality reweight for " << system_string << " " << sys_name << std::endl;
// 		return 1;
// 	}
// 	h_centrality_reweight = (TH1D*) fcent->Get("h_centrality_reweight") -> Clone(Form("h_centrality_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
// 	h_centrality_reweight -> SetDirectory(0);
// 	fcent -> Close();
	
// 	auto * fvtx = new TFile(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
// 	if ( !fvtx || fvtx->IsZombie() )
// 	{
// 		std::cerr << "Missing required vertex reweight for " << system_string << " "<< sys_name << std::endl;
// 		return 1;
// 	}
// 	h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight") -> Clone(Form("h_mbd_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
// 	h_mbd_reweight -> SetDirectory(0);
// 	fvtx -> Close();
	
// 	auto * fsumeT = new TFile(Form("sumeT/sumeT_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
// 	if ( !fsumeT || fsumeT->IsZombie() )
// 	{
// 		std::cerr << "Missing required sumeT reweight for " << system_string << " "<< sys_name << std::endl;
// 		return 1;
// 	}
// 	h_sumeT_reweight = (TH1D*) fsumeT->Get("h_sumeT_reweight") -> Clone(Form("h_sumeT_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
// 	h_sumeT_reweight -> SetDirectory(0);
// 	fsumeT -> Close();

// 	float truth_leading_cut = rb.get_truth_leading_cut();
// 	float truth_subleading_cut = rb.get_truth_subleading_cut();

// 	float reco_leading_cut = rb.get_reco_leading_cut();
// 	float reco_subleading_cut = rb.get_reco_subleading_cut();

// 	float measure_leading_cut = rb.get_measure_leading_cut();
// 	float measure_subleading_cut = rb.get_measure_subleading_cut();

// 	int measure_leading_bin = rb.get_measure_leading_bin();
// 	int measure_subleading_bin = rb.get_measure_subleading_bin();

// 	float sample_boundary[4] = {0};
// 	for (int ib = 0; ib < 4; ib++)
// 	{
// 		sample_boundary[ib] = rb.get_sample_boundary(ib);
// 		std::cout <<  sample_boundary[ib] << std::endl;
// 	}
// 	std::cout << "Max reco bin: " << max_reco_bin << std::endl;
// 	std::cout << "Truth1: " << truth_leading_cut << std::endl;
// 	std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
// 	std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
// 	std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
// 	std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
// 	std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

//   	TH1D * h_truth_lead_sample[3];
// 	for (int i = 0; i < 3; i++)
// 	{
// 		h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	}
  
// 	TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

// 	TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_match_reco_lead = new TH1D("h_match_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
// 	TH1D *h_match_reco_sublead = new TH1D("h_match_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	
// 	TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
// 	TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
// 	TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);
	
// 	// pure fills
// 	TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
// 	TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

// 	TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
// 	TH1D *h_linear_reco_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

// 	TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
// 	TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

// 	TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
// 	TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
// 	TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

// 	TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
// 	TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
// 	TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

// 	TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

// 	// For the prior sensitivity
// 	TH1D *h_flatreweight_pt1pt2 = new TH1D("h_unfold_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

// 	const double prior_fraction = prior_sys ? 0.5 : 1.0;

	
// 	if (!full_or_half && !prior_sys && !primer)
// 	{
// 		std::cout << "doing prior" << std::endl;
// 		TFile *fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
// 		TH1D *h_unfold_flat = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
// 		TFile *ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
// 		std::cout << "doing prior" << std::endl;
// 		TH1D *h_truth_flat = (TH1D*) ftr->Get("h_truth_flat_pt1pt2");
// 		TH2D *h2_flatreweight_pt1pt2 = new TH2D("h2_flatreweight_pt1pt2",";p_{T,1}^{truth} [GeV]; p_{T,2}^{truth} [GeV] ; Reweight Factor", nbins, ipt_bins, nbins, ipt_bins);
		
// 		if (h_unfold_flat->Integral() != 0) h_unfold_flat->Scale(1./h_unfold_flat->Integral());
// 		if (h_truth_flat->Integral() != 0)  h_truth_flat->Scale(1./h_truth_flat->Integral());
// 		for (int ibin = 0; ibin < nbins*nbins; ibin++)
// 		{
// 			float v = h_unfold_flat->GetBinContent(ibin+1);
// 			float b = h_truth_flat->GetBinContent(ibin+1);

// 			float pt1_bin = ibin/nbins;
// 			float pt2_bin = ibin%nbins;

// 			std::cout << "ibin : " << ibin << " : " << v << " / " << b << " = " << v/b << std::endl;
// 			int gbin = h2_flatreweight_pt1pt2->GetBin(pt1_bin+1, pt2_bin+1);
// 			if (b > 0)
// 			{
// 				h_flatreweight_pt1pt2->SetBinContent(ibin+1, v/b);
// 				h2_flatreweight_pt1pt2->SetBinContent(gbin, v/b);
// 			}
// 			else
// 			{
// 				h_flatreweight_pt1pt2->SetBinContent(ibin+1, 1);
// 				h2_flatreweight_pt1pt2->SetBinContent(gbin, 0);
// 			}
// 		}

// 		TCanvas *cre = new TCanvas("cre","cre", 500, 500);
// 		cre->SetLeftMargin(0.1);
// 		cre->SetRightMargin(0.19);
// 		// h2_flatreweight_pt1pt2->GetYaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
// 		// h2_flatreweight_pt1pt2->GetXaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
// 		h2_flatreweight_pt1pt2->Draw("colz");
// 		dlutility::DrawSPHENIX(0.2, 0.87);
// 		dlutility::drawText("Prior Reweighting Matrix", 0.2, 0.77);
// 		dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.2, 0.72);
// 		cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
// 	}
	

//   	const int nbin_response = nbins*nbins;
//   	RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

// 	TRandom * rng = new TRandom(0);
//   	for (int isample = 0; isample < 3; isample++)
//     {
// 		std::cout << "Sample " << isample << std::endl;
// 		int entries2 = tn[isample]->GetEntries();    
// 		for (int i = 0; i < entries2; i++)
// 		{

// 			tn[isample] -> GetEntry(i);
			
// 			int inrecojets = nrecojets[isample];
			
// 			double event_scale 			= scale_factor[isample];
// 			double mbd_vertex_scale 	= 1.0;
// 			double sumeT_scale 			= 1.0;
// 			double centrality_scale 	= 1.0;

// 			const bool in_maxpttruth_range = ( maxpttruth[isample] >= sample_boundary[isample] && maxpttruth[isample] < sample_boundary[isample+1] );
// 			if ( !in_maxpttruth_range ) 
// 			{ 
// 				continue; 
// 			}
			
// 			// const bool in_centrality_bin = ( centrality[isample] >= icentrality_bins[centrality_bin] && centrality[isample] < icentrality_bins[centrality_bin+1] );
// 			// if (!ispp && !is_in_centrality_bin) 
// 			// { 
// 			// 	continue; 
// 			// }

// 			bool fill_response = false;
// 			if ( full_or_half && rng->Uniform() >= 0.5 )
// 			{
// 				fill_response = true;
// 			}

// 			const int imbd_bin = h_mbd_reweight -> FindBin( mbd_vertex[isample] );
// 			if ( imbd_bin < 1 || imbd_bin > h_mbd_reweight->GetNbinsX() )
// 			{
// 				std::cerr << "Warning: mbd_vertex " << mbd_vertex[isample] << " is out of range for reweighting histogram" << std::endl;
// 				continue;
// 			}
// 			mbd_vertex_scale 	= h_mbd_reweight -> GetBinContent(imbd_bin);
			

// 			if ( !ispp )
// 			{
// 				const int isumeT_bin = h_sumeT_reweight -> FindBin( sumeT[isample] );
// 				const int icent_bin = h_centrality_reweight -> FindBin( centrality[isample] );
	
// 				if ( isumeT_bin < 1 || isumeT_bin > h_sumeT_reweight->GetNbinsX() )
// 				{
// 					// std::cerr << "Warning: sumeT " << sumeT[isample] << " is out of range for reweighting histogram" << std::endl;
// 					continue;
// 				}
// 				if ( icent_bin < 1 || icent_bin > h_centrality_reweight->GetNbinsX() )
// 				{
// 					// std::cerr << "Warning: centrality " << centrality[isample] << " is out of range for reweighting histogram" << std::endl;
// 					continue;
// 				}
				
// 				sumeT_scale 		= h_sumeT_reweight -> GetBinContent(isumeT_bin);
// 				centrality_scale 	= h_centrality_reweight -> GetBinContent(icent_bin);
// 			}


// 			const bool use_for_response = !half_closure || fill_response; // fill full closure, or fill r
// 			const bool use_for_test = full_closure || (half_closure && !fill_response);
// 			// const bool use_for_response = !(full_or_half && !fill_response); // fill full closure, or fill response half closure
// 			// const bool use_for_test = !(full_or_half && fill_response); // fill full closure

			
// 			// if ( !full_or_half && !prior_sys && !primer)
// 			if (primer != 1 && !full_or_half)
// 			{
// 				// event_scale *= mbd_vertex_scale * sumeT_scale * centrality_scale;
// 				event_scale *= mbd_vertex_scale * sumeT_scale;
// 			}

		
// 			float max_truth = 0;
// 			float max_reco = 0;
// 			float min_truth = 0;
// 			float min_reco = 0;

// 			if (pt1_truth[isample] >= pt2_truth[isample])
// 			{
// 				max_truth = pt1_truth[isample];
// 				max_reco = pt1_reco[isample];
// 				min_truth = pt2_truth[isample];
// 				min_reco = pt2_reco[isample];
// 			}
// 			else
// 			{
// 				max_truth = pt2_truth[isample];
// 				max_reco = pt2_reco[isample];
// 				min_truth = pt1_truth[isample];
// 				min_reco = pt1_reco[isample];
// 			}
			
// 			float pt1_truth_bin = nbins;
// 			float pt2_truth_bin = nbins;
// 			float pt1_reco_bin = nbins;
// 			float pt2_reco_bin = nbins;

// 			float e1 = max_truth;
// 			float e2 = min_truth;
// 			float es1 = max_reco;
// 			float es2 = min_reco;

// 			h_truth_lead_sample[isample]->Fill(e1, event_scale);

// 			double smear1 = 0;
// 			double smear2 = 0;

// 			if (ispp && !full_or_half)
// 			{
// 				smear1 = fgaus->GetRandom();
// 				smear2 = fgaus->GetRandom();
// 			}
// 			else if (!full_or_half)
// 			{
// 				fgaus->SetParameter(2, f_smear->Eval(e1));
// 				smear1 = fgaus->GetRandom();
// 				fgaus->SetParameter(2, f_smear->Eval(e2));
// 				smear2 = fgaus->GetRandom();
// 			}
		
// 			if (JES_sys != 0)
// 			{
// 				es1 = es1 + (JES_sys + smear1)*e1;
// 				es2 = es2 + (JES_sys + smear2)*e2; 
// 			}
// 			else
// 			{
// 				es1 = es1 + smear1*e1;
// 				es2 = es2 + smear2*e2; 	  
// 			}
		
// 			float maxi = std::max(es1, es2);
// 			float mini = std::min(es1, es2);
// 			float maxit = std::max(e1, e2);
// 			float minit = std::min(e1, e2);

// 			if (maxi >= ipt_bins[max_reco_bin]){ continue; }

// 			if (maxit > truth_leading_cut){ h_truth_lead->Fill(e1, event_scale);}
// 			if (minit >  truth_subleading_cut){ h_truth_sublead->Fill(e2, event_scale);}

// 			if (maxi >  reco_leading_cut){ h_reco_lead->Fill(maxi, event_scale);}
// 			if (mini >  reco_subleading_cut){ h_reco_sublead->Fill(mini, event_scale);}


// 			for (int ib = 0; ib < nbins; ib++)
// 			{
// 				if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
// 				{
// 					pt1_truth_bin = ib;
// 				}
// 				if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
// 				{
// 					pt2_truth_bin = ib;
// 				}
// 				if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
// 				{
// 					pt1_reco_bin = ib;
// 				}
// 				if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
// 				{
// 					pt2_reco_bin = ib;
// 				}
// 			}
		
// 			bool truth_good = (maxit >= sample_boundary[1]/* truth_leading_cut*/ && minit >= truth_subleading_cut && dphi_truth[isample] >= dphicuttruth);
// 			bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] >= dphicut);

// 			if (!truth_good && !reco_good) continue;

// 			if (!prior_sys && !primer && !full_or_half)
// 			{
// 				int recorrectbin = pt1_truth_bin*nbins + pt2_truth_bin;	      
// 				event_scale *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
// 			}
			
// 			if (truth_good && !reco_good)
// 			{
// 				if (use_for_response)
// 				{
// 					h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
// 					h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);

// 					rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 				}

// 				if (use_for_test)
// 				{
// 					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
// 					{
// 						h_truth_xj->Fill(minit/maxit, event_scale);
// 						h_linear_truth_xj->Fill(minit/maxit, event_scale);
// 					}
// 					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					h_e1e2->Fill(e1, e2, event_scale);
// 					h_e1e2->Fill(e2, e1, event_scale);
// 				}
// 			} 
			
// 			if (match[isample] && reco_good && truth_good)
// 			{
			
// 				if (maxit > truth_leading_cut) h_match_truth_lead->Fill(maxit, event_scale);
// 				if (minit >  truth_subleading_cut) h_match_truth_sublead->Fill(minit, event_scale);
// 				if (maxi >  reco_leading_cut) h_match_reco_lead->Fill(maxi, event_scale);
// 				if (mini >  reco_subleading_cut) h_match_reco_sublead->Fill(mini, event_scale);
				
			
// 				if (use_for_response)
// 				{
					
// 					h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
// 					h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
// 					h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
// 					h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
// 					h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
// 					h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);	      
				
// 				}
// 				// else if (use_for_test)
// 				if (use_for_test)
// 				{
// 					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
// 					{
// 						h_truth_xj->Fill(minit/maxit, event_scale);
// 						h_linear_truth_xj->Fill(minit/maxit, event_scale);
// 					}
// 					h_pt1pt2->Fill(es1, es2, event_scale);
// 					h_pt1pt2->Fill(es2, es1, event_scale);
// 					h_e1e2->Fill(e1, e2, event_scale);
// 					h_e1e2->Fill(e2, e1, event_scale);

// 					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
// 					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
// 					h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
// 					h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
// 				}
				
// 				if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
// 				{
// 					h_reco_xj->Fill(mini/maxi, event_scale);
// 					h_linear_reco_xj->Fill(mini/maxi, event_scale);
// 				}

// 				h_mbd_vertex->Fill(mbd_vertex[isample], event_scale);
// 				if (!ispp && has_sumeT[isample]){ h_sumeT->Fill(sumeT[isample], event_scale);}
// 				if (!ispp) {h_centrality->Fill(centrality[isample], event_scale);}
				
// 			} 
			
// 		}
//     }

// 	h_flat_reco_pt1pt2->Scale(.5);
// 	h_flat_truth_pt1pt2->Scale(.5);
// 	h_flat_reco_to_response_pt1pt2->Scale(.5);
// 	h_flat_truth_to_response_pt1pt2->Scale(.5);
// 	h_flat_response_pt1pt2->Scale(.5);

// 	int number_of_mins = 1;
// 	if (minentries == 0) number_of_mins = 20;
//   	for (int imin = minentries; imin < minentries + number_of_mins; imin++)       
//     {
// 		std::cout << "Minimum Entries == " << imin << std::endl;
// 		int nbinsx = h_flat_response_pt1pt2->GetXaxis()->GetNbins();
// 		int nbinsy = h_flat_response_pt1pt2->GetYaxis()->GetNbins();

// 		TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
// 		TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
// 		if (minentries)
// 		{
// 			h_flat_truth_mapping->SetName("h_flat_truth_mapping");
// 			h_flat_reco_mapping->SetName("h_flat_reco_mapping");
// 		}
// 		else
// 		{
// 			h_flat_truth_mapping->SetName(Form("h_flat_truth_mapping_min%d", imin));
// 			h_flat_reco_mapping->SetName(Form("h_flat_reco_mapping_min%d", imin));
// 		}
      
// 		h_flat_reco_mapping->Reset();
// 		h_flat_truth_mapping->Reset();
	
// 		int nempty_reco = 0;
// 		std::vector<int> binnumbers_reco{};
// 		int nempty_truth = 0;
// 		std::vector<int> binnumbers_truth{};
// 		int nrecobins = 0;
// 		int ntruthbins = 0;
// 		for (int ib = 0; ib < nbins*nbins; ib++)
// 		{
// 			if (h_count_flat_reco_pt1pt2->GetBinContent(ib+1) < imin)
// 			{
// 				nempty_reco++;
// 				binnumbers_reco.push_back(0);
// 				h_flat_reco_mapping->SetBinContent(ib+1, 0);
// 			}
// 	  		else
// 			{
// 				nrecobins++;
// 				binnumbers_reco.push_back(nrecobins);
// 				h_flat_reco_mapping->SetBinContent(ib+1, nrecobins);
// 			}
// 	  		if (h_count_flat_truth_pt1pt2->GetBinContent(ib+1) < imin)
// 	    	{
// 				nempty_truth++;
// 				binnumbers_truth.push_back(0);
// 				h_flat_truth_mapping->SetBinContent(ib+1, 0);
// 			}
// 			else
// 			{
// 				ntruthbins++;
// 				binnumbers_truth.push_back(ntruthbins);
// 				h_flat_truth_mapping->SetBinContent(ib+1, ntruthbins);
// 			}
// 		}
// 		TH1D *h_flat_truth_skim = new TH1D(Form("h_flat_truth_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
// 		TH1D *h_flat_reco_skim = new TH1D(Form("h_flat_reco_skim_min%d", imin),"", nrecobins, 0, nrecobins);
// 		TH1D *h_flat_truth_to_unfold_skim = new TH1D(Form("h_flat_truth_to_unfold_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
// 		TH1D *h_flat_reco_to_unfold_skim = new TH1D(Form("h_flat_reco_to_unfold_skim_min%d", imin),"", nrecobins, 0, nrecobins);

// 		TH2D *h_flat_response_skim = new TH2D(Form("h_flat_response_skim_min%d", imin),"", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

// 		if (minentries)
// 		{
// 			h_flat_response_skim->SetName("h_flat_response_skim");
// 			h_flat_truth_skim->SetName("h_flat_truth_skim");
// 			h_flat_reco_skim->SetName("h_flat_reco_skim");
// 			h_flat_truth_to_unfold_skim->SetName("h_flat_truth_to_unfold_skim");
// 			h_flat_reco_to_unfold_skim->SetName("h_flat_reco_to_unfold_skim");
// 		}

//       	for (int ib = 0; ib < nbins*nbins; ib++)
// 		{
// 			int bintruth = binnumbers_truth.at(ib);      
// 			int binreco = binnumbers_reco.at(ib);
// 			if (binreco)
// 			{
// 				h_flat_reco_skim->SetBinContent(binreco, h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1));
// 				h_flat_reco_skim->SetBinError(binreco, h_flat_reco_to_response_pt1pt2->GetBinError(ib+1));
// 				h_flat_reco_to_unfold_skim->SetBinContent(binreco, h_flat_reco_pt1pt2->GetBinContent(ib+1));
// 				h_flat_reco_to_unfold_skim->SetBinError(binreco, h_flat_reco_pt1pt2->GetBinError(ib+1));

// 			}

// 			if (bintruth)
// 			{
			
// 				h_flat_truth_skim->SetBinContent(bintruth, h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1));
// 				h_flat_truth_skim->SetBinError(bintruth, h_flat_truth_to_response_pt1pt2->GetBinError(ib+1));
// 				h_flat_truth_to_unfold_skim->SetBinContent(bintruth, h_flat_truth_pt1pt2->GetBinContent(ib+1));
// 				h_flat_truth_to_unfold_skim->SetBinError(bintruth, h_flat_truth_pt1pt2->GetBinError(ib+1));

// 				for (int ibr = 0; ibr < nbins*nbins; ibr++)
// 				{
// 					int binreco2 = binnumbers_reco.at(ibr);
// 					if (binreco2)
// 					{
// 						int rbin = h_flat_response_skim->GetBin(binreco2, bintruth);
// 						int tbin = h_flat_response_pt1pt2->GetBin(ibr+1, ib+1);
// 						h_flat_response_skim->SetBinContent(rbin, h_flat_response_pt1pt2->GetBinContent(tbin));
// 						h_flat_response_skim->SetBinError(rbin, h_flat_response_pt1pt2->GetBinError(tbin));
// 					}
// 				}
// 	    	}
// 		}

// 		std::cout <<" Nbins skim reco = "<<h_flat_reco_skim->GetNbinsX()<<std::endl;
// 		std::cout <<" Nbins skim reco = "<<h_flat_truth_skim->GetNbinsX()<<std::endl;
// 		std::cout <<" Nbins skim reco = "<<h_flat_response_skim->GetXaxis()->GetNbins()<<"  --  " << h_flat_response_skim->GetYaxis()->GetNbins()<<std::endl; 
// 		RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);
// 		if (minentries)
// 		{
// 			rooResponsehist.SetName("response_noempty");
// 		}
// 		else
// 		{
// 			rooResponsehist.SetName(Form("response_noempty_min%d", imin));
// 		}

// 		if (minentries)
// 		{
// 			TH1D* h_flat_unfold_skim[niterations];
// 			TH1D* h_flat_unfold_pt1pt2[niterations];
// 			TH1D* h_flat_refold_skim[niterations];
// 			TH1D* h_flat_refold_pt1pt2[niterations];
// 			TH1D* h_flat_truth_fold_skim;
// 			TH1D* h_flat_truth_fold_pt1pt2;
// 			int niter = 3;
			
// 			h_flat_truth_fold_skim = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_truth_skim, "hTruth_Folded");

// 			h_flat_truth_fold_pt1pt2 = (TH1D*) h_flat_truth_pt1pt2->Clone();
// 			h_flat_truth_fold_pt1pt2->Reset();
// 			h_flat_truth_fold_pt1pt2->SetName("h_flat_truth_fold_pt1pt2");
// 			for (int ib = 0; ib < nbins*nbins; ib++)
// 			{
// 				int bin = binnumbers_truth.at(ib);
// 				if (bin)
// 				{
// 					h_flat_truth_fold_pt1pt2->SetBinContent(ib+1, h_flat_truth_fold_skim->GetBinContent(bin));
// 					h_flat_truth_fold_pt1pt2->SetBinError(ib+1, h_flat_truth_fold_skim->GetBinError(bin));
// 				}
// 			}

// 			for (int iter = 0; iter < niterations; iter++ )
// 			{
		
// 				RooUnfoldBayes   unfold (&rooResponsehist, h_flat_reco_to_unfold_skim, iter + 1);    // OR

// 				h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
// 				std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
// 				h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
// 				h_flat_unfold_pt1pt2[iter]->Reset();
// 				h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
// 				for (int ib = 0; ib < nbins*nbins; ib++)
// 				{
// 					int bin = binnumbers_truth.at(ib);
// 					if (bin)
// 					{
// 						h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin));
// 						h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin));
// 					}
// 				}

// 				h_flat_refold_skim[iter] = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_unfold_skim[iter], "hRefolded");

// 				std::cout <<" Nbins skim reco = "<<h_flat_refold_skim[iter]->GetNbinsX()<<std::endl;
// 				h_flat_refold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
// 				h_flat_refold_pt1pt2[iter]->Reset();
// 				h_flat_refold_pt1pt2[iter]->SetName(Form("h_flat_refold_pt1pt2_%d",iter));
// 				for (int ib = 0; ib < nbins*nbins; ib++)
// 				{
// 					int bin = binnumbers_truth.at(ib);
// 					if (bin)
// 					{
// 						h_flat_refold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_refold_skim[iter]->GetBinContent(bin));
// 						h_flat_refold_pt1pt2[iter]->SetBinError(ib+1, h_flat_refold_skim[iter]->GetBinError(bin));
// 					}
// 				}

// 	    	}


// 			TCanvas *c = new TCanvas("c","c", 500, 500);

// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
// 				dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
// 			}
		
// 			dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
// 			dlutility::SetLineAtt(h_flat_reco_pt1pt2, kBlue, 2, 1);

// 			h_flat_truth_pt1pt2->Draw("hist");
// 			h_flat_reco_pt1pt2->Draw("hist same");
// 			h_flat_unfold_pt1pt2[niter]->Draw("same p");

// 			TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
// 			TH2D *h_pt1pt2_truth_fold = new TH2D("h_pt1pt2_truth_fold", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
// 			TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
// 			TH2D *h_pt1pt2_unfold[niterations];
// 			TH2D *h_pt1pt2_refold[niterations];
// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
// 				h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
// 				h_pt1pt2_refold[iter] = new TH2D("h_pt1pt2_refold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
// 				h_pt1pt2_refold[iter]->SetName(Form("h_pt1pt2_refold_iter%d", iter));
// 			}
// 			TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
// 			TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
// 			TH1D *h_xj_truth_fold = new TH1D("h_xj_truth_fold", ";x_{J};",nbins, ixj_bins);
// 			TH1D *h_xj_truth_direct = new TH1D("h_xj_truth_direct", ";x_{J};",nbins, ixj_bins);
// 			TH1D *h_xj_unfold[niterations];
// 			TH1D *h_xj_refold[niterations];
// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
// 				h_xj_refold[iter] = new TH1D(Form("h_xj_refold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
// 			}

// 			h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
// 			h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
// 			h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");

// 			histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
// 			histo_opps::make_sym_pt1pt2(h_flat_truth_fold_pt1pt2, h_pt1pt2_truth_fold, nbins);
// 			histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
// 				histo_opps::make_sym_pt1pt2(h_flat_refold_pt1pt2[iter], h_pt1pt2_refold[iter], nbins);
// 			}

// 			TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
// 			cpt1pt2->Divide(3, 1);
// 			cpt1pt2->cd(1);
// 			gPad->SetLogz();
// 			gPad->SetRightMargin(0.2);
// 			h_pt1pt2_reco->Draw("colz");
// 			dlutility::DrawSPHENIXpp(0.22, 0.85);
		
// 			cpt1pt2->cd(2);
// 			gPad->SetRightMargin(0.2);
// 			gPad->SetLogz();
// 			h_pt1pt2_truth->Draw("colz");
// 			cpt1pt2->cd(3);
// 			gPad->SetRightMargin(0.2);
// 			gPad->SetLogz();
// 			h_pt1pt2_unfold[niter]->Draw("colz");


// 			cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, diagnostic_name.c_str()));

// 			histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 			histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 			histo_opps::project_xj(h_pt1pt2_truth_fold, h_xj_truth_fold, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 			histo_opps::project_xj(h_e1e2, h_xj_truth_direct, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 				histo_opps::project_xj(h_pt1pt2_refold[iter], h_xj_refold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
// 			}
		
// 			TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
// 			dlutility::ratioPanelCanvas(cxj);
// 			cxj->cd(1);
// 			dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
// 			dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

// 			dlutility::SetLineAtt(h_truth_xj, kRed, 2, 1);
// 			dlutility::SetMarkerAtt(h_truth_xj, kRed, 2, 24);

// 			dlutility::SetLineAtt(h_linear_truth_xj, kRed, 2, 1);
// 			dlutility::SetMarkerAtt(h_linear_truth_xj, kRed, 2, 25);

// 			dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
// 			dlutility::SetMarkerAtt(h_reco_xj, kBlue, 2, 24);

// 			dlutility::SetLineAtt(h_linear_reco_xj, kBlue, 2, 1);
// 			dlutility::SetMarkerAtt(h_linear_reco_xj, kBlue, 2, 25);

		
// 			dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
// 			dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

// 			dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
// 			dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
		
// 			histo_opps::normalize_histo(h_xj_truth, nbins);
// 			histo_opps::normalize_histo(h_xj_reco, nbins);
// 			histo_opps::normalize_histo(h_truth_xj, nbins);
// 			histo_opps::normalize_histo(h_reco_xj, nbins);

// 			//histo_opps::normalize_histo(h_linear_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
// 			h_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
// 			//h_linear_reco_xj->Scale(1./h_linear_reco_xj->Integral(), "width");

// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
// 			}

// 			dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
// 			h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
// 			//h_xj_truth->SetMaximum(5);
// 			h_xj_truth->Draw("p");
// 			h_xj_unfold[niter]->Draw("same p");
// 			h_xj_reco->Draw("hist same");
// 			h_xj_truth->Draw("same hist");
// 			h_xj_reco->Draw("p same");
// 			h_xj_unfold[niter]->Draw("same hist");
// 			h_xj_unfold[niter]->Draw("same p");
// 			dlutility::DrawSPHENIXpp(0.22, 0.84);
// 			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
// 			dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
// 			dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
// 			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
// 			TLegend *leg = new TLegend(0.22, 0.3, 0.4, 0.55);
// 			leg->SetLineWidth(0);
// 			leg->AddEntry(h_xj_reco, "Pythia Reco");
// 			leg->AddEntry(h_xj_truth, "Pythia8 Truth");
// 			leg->AddEntry(h_xj_unfold[niter], "Unfolded");
// 			leg->Draw("same");
// 			dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
			
// 			cxj->cd(2);

// 			TH1D *h_reco_compare = (TH1D*) h_xj_unfold[niter]->Clone();
// 			h_reco_compare->Divide(h_xj_truth);
// 			h_reco_compare->SetTitle(";x_{J}; Unfold / Truth");
// 			dlutility::SetFont(h_reco_compare, 42, 0.06);
// 			dlutility::SetLineAtt(h_reco_compare, kBlack, 1,1);
// 			dlutility::SetMarkerAtt(h_reco_compare, kBlack, 1,8);

		
// 			h_reco_compare->Draw("p");
// 			TLine *line = new TLine(0.1, 1, 1, 1);
// 			line->SetLineStyle(4);
// 			line->SetLineColor(kRed + 3);
// 			line->SetLineWidth(2);
// 			line->Draw("same");

// 			TCanvas *c_iter = new TCanvas("c_iter", "c_iter");
// 			TH1D *h_closure[niterations];

// 			for (int iter = 0; iter < niterations; iter++)
// 			{
// 				h_closure[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
// 				h_closure[iter]->SetName(Form("h_closure_%d", iter));
				
// 				h_closure[iter]->Divide(h_xj_truth);
// 				h_closure[iter]->SetTitle(";x_{J}; Unfold / Truth");
// 				dlutility::SetFont(h_closure[iter], 42, 0.06);
// 			}
// 			int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 

// 			for (int iter = 0; iter < 5; iter++)
// 			{
// 				dlutility::SetLineAtt(h_closure[1 + 2*iter], colors[iter], 1,1);
// 				dlutility::SetMarkerAtt(h_closure[1 + 2*iter], colors[iter], 1,8);
// 			}
// 			h_closure[1]->Draw("p");
// 			h_closure[3]->Draw("p same");
// 			h_closure[5]->Draw("p same");
// 			h_closure[7]->Draw("p same");
// 			h_closure[9]->Draw("p same");
// 			TLine *line2 = new TLine(0.1, 1, 1, 1);
// 			line2->SetLineStyle(4);
// 			line2->SetLineColor(kRed + 3);
// 			line2->SetLineWidth(2);
// 			line2->Draw("same");


// 			TCanvas *cjet = new TCanvas("cjet","cjet", 700, 500);
// 			dlutility::createCutCanvas(cjet);
// 			cjet->cd(1);
// 			gPad->SetLogy();

// 			dlutility::SetLineAtt(h_truth_lead, kRed, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_lead, kRed, 1, 24);

// 			dlutility::SetLineAtt(h_truth_sublead, kBlue, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_sublead, kBlue, 1, 24);

// 			dlutility::SetLineAtt(h_reco_lead, kRed, 1, 1);
// 			dlutility::SetMarkerAtt(h_reco_lead, kRed, 1, 20);

// 			dlutility::SetLineAtt(h_reco_sublead, kBlue, 1, 1);
// 			dlutility::SetMarkerAtt(h_reco_sublead, kBlue, 1, 20);

// 			h_truth_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
// 			h_truth_lead->Draw();
// 			h_truth_sublead->Draw("same");
// 			h_reco_lead->Draw("same");
// 			h_reco_sublead->Draw("same");
// 			cjet->cd(2);
// 			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
// 			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.1, 0.74, 0, kBlack, 0.1);
// 			dlutility::drawText("All Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
// 			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);

// 			TLegend *leg33 = new TLegend(0.01, 0.13, 0.7, 0.62);
// 			leg33->SetLineWidth(0);
// 			leg33->AddEntry(h_reco_lead, "Leading Reco");
// 			leg33->AddEntry(h_reco_sublead, "Subleading Reco");
// 			leg33->AddEntry(h_truth_lead, "Leading Truth");
// 			leg33->AddEntry(h_truth_sublead, "Subleading Truth");
// 			leg33->SetTextSize(0.1);
// 			leg33->Draw("same");

// 			TCanvas *cjetdiv = new TCanvas("cjetdiv","cjetdiv", 700, 500);
// 			dlutility::createCutCanvas(cjetdiv);
// 			cjetdiv->cd(1);
// 			gPad->SetLogy();

// 			dlutility::SetLineAtt(h_truth_lead, kBlack, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_lead, kBlack, 0.8, 24);

// 			dlutility::SetLineAtt(h_truth_lead_sample[0], kRed, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_lead_sample[0], kRed, 1, 20);
// 			dlutility::SetLineAtt(h_truth_lead_sample[1], kGreen, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_lead_sample[1], kGreen, 1, 20);
// 			dlutility::SetLineAtt(h_truth_lead_sample[2], kBlue, 1, 1);
// 			dlutility::SetMarkerAtt(h_truth_lead_sample[2], kBlue, 1, 20);

// 			h_truth_lead->SetTitle(";Leading Jet p_{T} [GeV];counts * lumiscale");
// 			h_truth_lead_sample[0]->SetMinimum(1);
// 			h_truth_lead_sample[0]->Draw("p");
// 			h_truth_lead->Draw("p same");
// 			h_truth_lead_sample[1]->Draw("same p");
// 			h_truth_lead_sample[2]->Draw("same p");
// 			h_truth_lead->Draw(" same p");


// 			TLine  *line1 = new TLine(14, 1, 14, h_truth_lead->GetBinContent(15));
// 			dlutility::SetLineAtt(line1, kRed, 1.5, 1);
// 			line1->Draw("same");
// 			TLine  *lin2 = new TLine(20, 1, 20, h_truth_lead->GetBinContent(21));
// 			dlutility::SetLineAtt(lin2, kGreen, 1.5, 1);
// 			lin2->Draw("same");
// 			TLine  *lin3 = new TLine(30, 1, 30, h_truth_lead->GetBinContent(31));
// 			dlutility::SetLineAtt(lin3, kBlue, 1.5, 1);
// 			lin3->Draw("same");
// 			cjetdiv->cd(2);

// 			dlutility::DrawSPHENIXppsize(0.05, 0.84, 0.08, 1, 0, 1, "Pythia8");
// 			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.05, 0.74, 0, kBlack, 0.08);
// 			dlutility::drawText("All Dijet Pairs", 0.05, 0.69, 0, kBlack, 0.08);
// 			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.05, 0.64, 0, kBlack, 0.08);

// 			TLegend *leg3 = new TLegend(0.01, 0.13, 0.7, 0.62);
// 			leg3->SetLineWidth(0);
// 			leg3->AddEntry(h_truth_lead, "Combined Sample");
// 			leg3->AddEntry(h_truth_lead_sample[0], "10 GeV Sample");
// 			leg3->AddEntry(h_truth_lead_sample[1], "20 GeV Sample");
// 			leg3->AddEntry(h_truth_lead_sample[2], "30 GeV Sample");
// 			leg3->SetTextSize(0.08);
// 			leg3->Draw("same");

// 			cjetdiv->Print(Form("%s/unfolding_plots/combined_sample_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, diagnostic_name.c_str()));

// 			TCanvas *cmjet = new TCanvas("cmjet","cmjet", 700, 500);
// 			dlutility::createCutCanvas(cmjet);
// 			cmjet->cd(1);
// 			gPad->SetLogy();

// 			dlutility::SetLineAtt(h_match_truth_lead, kRed, 1, 1);
// 			dlutility::SetMarkerAtt(h_match_truth_lead, kRed, 1, 24);

// 			dlutility::SetLineAtt(h_match_truth_sublead, kBlue, 1, 1);
// 			dlutility::SetMarkerAtt(h_match_truth_sublead, kBlue, 1, 24);

// 			dlutility::SetLineAtt(h_match_reco_lead, kRed, 1, 1);
// 			dlutility::SetMarkerAtt(h_match_reco_lead, kRed, 1, 20);

// 			dlutility::SetLineAtt(h_match_reco_sublead, kBlue, 1, 1);
// 			dlutility::SetMarkerAtt(h_match_reco_sublead, kBlue, 1, 20);

// 			h_match_reco_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
// 			h_match_reco_lead->Draw();
// 			h_match_truth_sublead->Draw("same");
// 			h_match_truth_lead->Draw("same");
// 			h_match_reco_sublead->Draw("same");
// 			cmjet->cd(2);
// 			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
// 			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.1, 0.74, 0, kBlack, 0.1);
// 			dlutility::drawText("Matched Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
// 			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);
// 			TLegend *lef3 = new TLegend(0.01, 0.13, 0.7, 0.62);
// 			lef3->SetLineWidth(0);
// 			lef3->AddEntry(h_match_reco_lead, "Leading Reco");
// 			lef3->AddEntry(h_match_reco_sublead, "Subleading Reco");
// 			lef3->AddEntry(h_match_truth_lead, "Leading Truth");
// 			lef3->AddEntry(h_match_truth_sublead, "Subleading Truth");
// 			lef3->SetTextSize(0.1);
// 			lef3->Draw("same");


// 			TCanvas *cresponse = new TCanvas("fd","fd", 500, 500);
// 			h_flat_response_pt1pt2->Draw("colz");
// 			TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 500);
// 			cresponseskim->SetLogz();
// 			h_flat_response_skim->Draw("colz");
// 			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);


// 			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
			
// 			if (primer > 0)
// 			{
// 				responsepath += "_PRIMER" + std::to_string(primer);;
// 			}
// 			if (half_closure)
// 			{
// 				responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_HALF";
// 			}
// 			else if (full_closure)
// 			{
// 				responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_FULL";
// 			}
// 			responsepath += "_" + sys_name;

// 			responsepath += ".root";
// 			TFile *fr = new TFile(responsepath.Data(),"recreate");
// 			rooResponsehist.Write();
// 			h_flat_reco_pt1pt2->Write();
// 			h_flat_truth_pt1pt2->Write();
// 			h_flat_truth_to_response_pt1pt2->Write();
// 			h_flat_truth_fold_pt1pt2->Write();
// 			h_flat_response_pt1pt2->Write();
// 			h_flat_truth_mapping->Write();
// 			h_flat_reco_mapping->Write();
// 			h_flat_truth_skim->Write();
// 			h_flat_reco_skim->Write();
// 			h_flat_response_skim->Write();
// 			h_mbd_vertex->Write();
// 			h_centrality->Write();
// 			h_sumeT->Write();
// 			h_pt1pt2->Write();
// 			h_e1e2->Write();
// 			// Full closure drawing
// 			h_linear_truth_xj->Write();
// 			h_xj_truth->Write();
// 			h_xj_truth_fold->Write();
// 			h_xj_truth_direct->Write();
// 			h_truth_xj->Write();
// 			for (int iter = 0 ; iter < niterations; iter++)
// 			{
// 				h_flat_unfold_pt1pt2[iter]->Write();
// 				h_xj_unfold[iter]->Write();
// 				h_flat_refold_pt1pt2[iter]->Write();
// 				h_xj_refold[iter]->Write();
// 			}
// 			h_xj_reco->Write();
				
// 			fr->Close();

// 		}
// 		else
// 		{
// 			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_min_" + std::to_string(imin) + ".root";

// 			TFile *fr = new TFile(responsepath.Data(),"recreate");
// 			rooResponsehist.Write();
			
// 			TH1D *h1 = (TH1D*)h_flat_reco_pt1pt2->Clone();
// 			TH1D *h2 = (TH1D*)h_flat_truth_pt1pt2->Clone();
// 			TH2D *h3 = (TH2D*)h_flat_response_pt1pt2->Clone();
// 			TH1D *h4 = (TH1D*)h_flat_truth_mapping->Clone();
// 			TH1D *h5 = (TH1D*)h_flat_reco_mapping->Clone();
// 			TH1D *h6 = (TH1D*)h_flat_truth_skim->Clone();
// 			TH1D *h7 = (TH1D*)h_flat_reco_skim->Clone();
// 			TH2D *h8 = (TH2D*)h_flat_response_skim->Clone();
// 			h1->Write();
// 			h2->Write();
// 			h3->Write();
// 			h4->Write();
// 			h5->Write();
// 			h6->Write();
// 			h7->Write();
// 			h8->Write();
						
// 			fr->Close();
// 		}
// 		std::cout << "Reco empty: " << nempty_reco << std::endl;
// 		std::cout << "Truth empty: " << nempty_truth << std::endl;
//     }
  	
// 	return 0;
// }


// Variant of createResponse_noempty_AA.cxx that applies the PRIOR
// systematic as a deterministic partial (default 50%) blend toward the
// PRIMER1 data-driven prior, instead of skipping the correction entirely.
// Uses the same PRIMER1 source histograms, denominator, and bin indexing
// as the reference framework's nominal correction (see createResponse_noempty_AA.cxx);
// the only new ingredient is w_partial = 1 + prior_fraction*(w_full-1).
// Output is tagged "PRIOR50" (not "PRIOR") so it cannot collide with the
// existing full-vs-none PRIOR systematic files.
int createResponse_noempty_AA_prior50 (
	const std::string configfile = "binning.config",
	const int full_or_half = 0,
	const int niterations = 10,
	const int cone_size = 4,
	const int centrality_bin = 0,
	const int primer = 0,
	const double prior_fraction = 0.5)
{

	gStyle->SetOptStat(0);
	dlutility::SetyjPadStyle();
	if (full_or_half < 0 || full_or_half > 2)
	{
		std::cerr << "Closure mode must be 0 (analysis), 1 (half closure), or 2 (full closure)"	<< std::endl;
		return 1;
	}
	
	const bool half_closure = (full_or_half == 1);
	const bool full_closure = (full_or_half == 2);
	const bool full_or_half_closure = (full_or_half != 0);

	bool ispp = (centrality_bin < 0);
	read_binning rb(configfile.c_str());

	std::string system_string = (ispp ? "pp" : "AA_cent_" + std::to_string(centrality_bin));
	std::string j10_file = std::getenv("TNUPLE_SIM_FILE_JET10");
	std::string j20_file = std::getenv("TNUPLE_SIM_FILE_JET20");
	std::string j30_file = std::getenv("TNUPLE_SIM_FILE_JET30");
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
	float nrecojets[3];
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
	for (int i = 0 ; i < 3; i++)
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
		tn[i]->SetBranchAddress("nrecojets", &nrecojets[i]);
		tn[i]->SetBranchAddress("matched", &match[i]);
		tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
		tn[i]->SetBranchAddress("centrality", &centrality[i]);
		has_sumeT[i] = tn[i]->GetBranch("sumeT") != nullptr;
		if (has_sumeT[i]) tn[i]->SetBranchAddress("sumeT", &sumeT[i]);
		else if (!ispp) std::cerr << "Warning: missing sumeT in " << fin[i]->GetName() << "; sumeT reweighting is disabled for this sample." << std::endl;
		
		has_event_weight[i] = tn[i]->GetBranch("weight") != nullptr;
		if (has_event_weight[i]) tn[i]->SetBranchAddress("weight", &event_weight[i]);

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
		
	// float cs_10 = (2.889e-6);
	// float cs_20 = 5.4067742e-8;
	// float cs_30 = (2.505e-9);
	float cs_10 = 0.000003997;
  	float cs_20 = 6.218e-8;
  	float cs_30 = 2.505e-9;
  
	float scale_factor[3];
	scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
	scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30; 
	scale_factor[2] = 1;

	Int_t minentries = rb.get_minentries();
	Int_t read_nbins = rb.get_nbins();
	//Int_t primer = rb.get_primer();

	Double_t dphicut = rb.get_dphicut();
	Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

	TF1 *fgaus = new TF1("fgaus", "gaus");
	fgaus->SetRange(-0.5, 0.5);
	Int_t zyam_sys = rb.get_zyam_sys();
	Double_t flow_v22_scale = rb.get_flow_sys();
	Double_t flow_v33_scale = rb.get_flow_v33_sys();
	const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
	Int_t inclusive_sys = rb.get_inclusive_sys();
	Double_t JES_sys = rb.get_jes_sys();
	Double_t JER_sys = rb.get_jer_sys();
	Int_t prior_sys = rb.get_prior_sys();
	int using_sys = 0;

	TF1 *f_smear = nullptr;
	if (!ispp)
	{
		f_smear  = (TF1*) rb.get_smear_function(centrality_bin);
		if (!f_smear)
		{
			std::cout << " NO SMEAR " << std::endl;
			exit(-1);
		}
	}

	std::string sys_name = "nominal";
	std::cout << "JES = " << JES_sys << std::endl;
	std::cout << "JER = " << JER_sys << std::endl;
	float width = 0.8 + JER_sys;
	fgaus->SetParameters(1, 0, width);

	if (prior_sys)
	{
		using_sys = 1;
		sys_name = "PRIOR50";
	}
	if (zyam_sys)
	{
		using_sys = 1;
		sys_name = "ZYAM";
	}
	if (flow_sys)
	{
		using_sys = 1;
		sys_name = rb.get_flow_systematic_name();
	}
	if (inclusive_sys)
	{
		using_sys = 1;
		sys_name = "INCLUSIVE";
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

  	const std::string diagnostic_name = half_closure ? "HALF_" + sys_name : (full_closure ? "FULL_" + sys_name : sys_name);

	const int nbins = read_nbins;

	float ipt_bins[nbins+1];
	float ixj_bins[nbins+1];

	rb.get_pt_bins(ipt_bins);
	rb.get_xj_bins(ixj_bins);
	for (int i = 0 ; i < nbins + 1; i++)
	{
		std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}

	Int_t max_reco_bin = rb.get_maximum_reco_bin();

	const int centrality_bins = rb.get_number_centrality_bins();

	float icentrality_bins[centrality_bins + 1];
	rb.get_centrality_bins(icentrality_bins);
	

	int prior_iteration = 1;
  
 	// Vertex Reweight
	std::vector<std::pair<float, float>> centrality_scales;
	if (primer != 1 && !full_or_half && !ispp)
	{
		TFile *fcent = new TFile(Form("centrality/centrality_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()),"r");
		TH1D *h_centrality_reweight = (TH1D*) fcent->Get("h_centrality_reweight");
		if (!fcent || fcent->IsZombie() || !h_centrality_reweight)
		{
			std::cerr << "Missing required centrality reweight for " << system_string
						<< " " << sys_name << std::endl;
			return 1;
		}
		for (int ib = 0; ib < h_centrality_reweight->GetNbinsX(); ib++)
		{
			centrality_scales.push_back(std::make_pair(h_centrality_reweight->GetBinLowEdge(ib+1) + h_centrality_reweight->GetBinWidth(ib+1), h_centrality_reweight->GetBinContent(ib+1)));
		}
	} // cent scaled

	std::vector<std::pair<float, float>> vertex_scales;

	Int_t vtx_sys = rb.get_vtx_sys();
	if (primer != 1 && !full_or_half)
	{

		TFile *fvtx = new TFile(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()),"r");
		TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");
		if (!fvtx || fvtx->IsZombie() || !h_mbd_reweight)
		{
			std::cerr << "Missing required vertex reweight for " << system_string
						<< " " << sys_name << std::endl;
			return 1;
		}
		for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
		{
			vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
		}
	} // vtx scaled

	std::vector<std::pair<float, float>> sumeT_scales;
	if (primer != 1 && !full_or_half && !ispp)
	{
		TFile *fsumeT = new TFile(Form("sumeT/sumeT_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()),"r");
		TH1D *h_sumeT_reweight = fsumeT ? (TH1D*) fsumeT->Get("h_sumeT_reweight") : nullptr;
		if (h_sumeT_reweight)
		{
			for (int ib = 0; ib < h_sumeT_reweight->GetNbinsX(); ib++)
			{
				sumeT_scales.push_back(std::make_pair(h_sumeT_reweight->GetBinLowEdge(ib+1) + h_sumeT_reweight->GetBinWidth(ib+1), h_sumeT_reweight->GetBinContent(ib+1)));
			}
		}
		else
		{
			std::cerr << "Missing required sumeT reweight file/hist for " << system_string
						<< " r0" << cone_size << " " << sys_name
						<< std::endl;
			return 1;
		}
	} // sumeT scaled

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
	for (int ib = 0; ib < 4; ib++)
	{
		sample_boundary[ib] = rb.get_sample_boundary(ib);
		std::cout <<  sample_boundary[ib] << std::endl;
	}
	std::cout << "Max reco bin: " << max_reco_bin << std::endl;
	std::cout << "Truth1: " << truth_leading_cut << std::endl;
	std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
	std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
	std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
	std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
	std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  	TH1D *h_truth_lead_sample[3];
	for (int i = 0; i < 3; i++)
	{
		h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	}
  
	TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

	TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_reco_lead = new TH1D("h_match_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_reco_sublead = new TH1D("h_match_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
	TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
	TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);
	// pure fills
	TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
	TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

	TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
	TH1D *h_linear_reco_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

	TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
	TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

	TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

	TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

	TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

	// For the prior sensitivity
	TH1D *h_flatreweight_pt1pt2 = new TH1D("h_unfold_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

	// prior_sys no longer skips the correction (as it does in the reference
	// framework); it instead scales it down to prior_fraction of the full
	// nominal correction, deterministically, bin by bin. prior_fraction=1.0
	// would reproduce the nominal correction exactly; prior_fraction=0.0
	// would reproduce the reference's "skip entirely" PRIOR behavior.
	if (!full_or_half && !primer)
	{
		std::cout << "doing prior (fraction=" << (prior_sys ? prior_fraction : 1.0) << ")" << std::endl;
		TFile *fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
		TH1D *h_unfold_flat = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
		TFile *ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "r");
		std::cout << "doing prior" << std::endl;
		TH1D *h_truth_flat = (TH1D*) ftr->Get("h_truth_flat_pt1pt2");
		TH2D *h2_flatreweight_pt1pt2 = new TH2D("h2_flatreweight_pt1pt2",";p_{T,1}^{truth} [GeV]; p_{T,2}^{truth} [GeV] ; Reweight Factor", nbins, ipt_bins, nbins, ipt_bins);

		if (h_unfold_flat->Integral() != 0) h_unfold_flat->Scale(1./h_unfold_flat->Integral());
		if (h_truth_flat->Integral() != 0)  h_truth_flat->Scale(1./h_truth_flat->Integral());
		for (int ibin = 0; ibin < nbins*nbins; ibin++)
		{
			float v = h_unfold_flat->GetBinContent(ibin+1);
			float b = h_truth_flat->GetBinContent(ibin+1);

			float pt1_bin = ibin/nbins;
			float pt2_bin = ibin%nbins;

			int gbin = h2_flatreweight_pt1pt2->GetBin(pt1_bin+1, pt2_bin+1);
			if (b > 0)
			{
				const double w_full = v/b;
				const double w = prior_sys ? (1.0 + prior_fraction*(w_full - 1.0)) : w_full;
				std::cout << "ibin : " << ibin << " : " << v << " / " << b << " = " << w_full << " -> applied " << w << std::endl;
				h_flatreweight_pt1pt2->SetBinContent(ibin+1, w);
				h2_flatreweight_pt1pt2->SetBinContent(gbin, w);
			}
			else
			{
				h_flatreweight_pt1pt2->SetBinContent(ibin+1, 1);
				h2_flatreweight_pt1pt2->SetBinContent(gbin, 0);
			}
		}

		TCanvas *cre = new TCanvas("cre","cre", 500, 500);
		cre->SetLeftMargin(0.1);
		cre->SetRightMargin(0.19);
		h2_flatreweight_pt1pt2->GetYaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
		h2_flatreweight_pt1pt2->GetXaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
		h2_flatreweight_pt1pt2->Draw("colz");
		dlutility::DrawSPHENIX(0.2, 0.87);
		dlutility::drawText("Prior Reweighting Matrix", 0.2, 0.77);
		dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.2, 0.72);
		cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
	}
  	int nbin_response = nbins*nbins;
  
  	RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

	// A fixed seed makes the half-closure train/test assignment reproducible.
	TRandom rng(314159);
	int max_matched_entries = 0;
	if (const char *max_matched_env = std::getenv("DIJET_MAX_MATCHED_EVENTS"))
	{
		max_matched_entries = std::atoi(max_matched_env);
		if (max_matched_entries > 0)
		{
			std::cout << "Proof/sample mode: limiting each matched simulation sample to "
			<< max_matched_entries << " entries via DIJET_MAX_MATCHED_EVENTS." << std::endl;
		}
	}
	
  	for (int isample = 0; isample < 3; isample++)
    {
		std::cout << "Sample " << isample << std::endl;
		int entries0 = 0;
		int entries1 = tn[isample]->GetEntries()/2;    
		int entries2 = tn[isample]->GetEntries();
		if (max_matched_entries > 0) entries2 = std::min(entries2, max_matched_entries);
		//if (full_or_half == 1) entries2 = entries1;

		for (int i = entries0; i < entries2; i++)
		{
			tn[isample]->GetEntry(i);
			int inrecojets = nrecojets[isample];
			double event_scale = has_event_weight[isample] ? event_weight[isample] : scale_factor[isample];

			if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;	  

			if (!(ispp) &&  centrality[isample] < icentrality_bins[centrality_bin] || centrality[isample] >= icentrality_bins[centrality_bin+1]) continue;
			// Vertex Rewieghting
			bool fill_response = false;
			if (full_or_half_closure)
			{
				if (rng.Uniform() >= 0.5)
				{
					fill_response = true;
				}
			}

			// Matches the reference framework's original condition
			// (!(full_or_half && !fill_response) / !(full_or_half && fill_response)):
			// both are unconditionally true when full_or_half==0, so normal
			// (non-closure) running fills both the response and the
			// truth/reco "test" histograms from every event.
			const bool use_for_response = !(full_or_half_closure && !fill_response);
			const bool use_for_test = !(full_or_half_closure && fill_response);

			if (primer != 1 && !full_or_half)
			{
				for (int ib = 0; ib < (int)vertex_scales.size(); ib++)
				{
					if (mbd_vertex[isample] < vertex_scales.at(ib).first)
					{
						event_scale *= vertex_scales.at(ib).second;
						//if (i < 100) std:cout << "found z = " << vertex_scales.at(ib).first << " " <<vertex_scales.at(ib).second<<std::endl;
						break;
					}
				} //mbd vertex reweighting

				if (!ispp && has_sumeT[isample] && !sumeT_scales.empty())
				{
					for (int ib = 0; ib < (int)sumeT_scales.size(); ib++)
					{
						if (sumeT[isample] < sumeT_scales.at(ib).first)
						{
							event_scale *= sumeT_scales.at(ib).second;
							break;
						}
					}
				} // sumeT reweighting

				if (!ispp && !full_or_half)
				{
					for (int ib = 0; ib < (int)centrality_scales.size(); ib++)
					{
						if (centrality[isample] < centrality_scales.at(ib).first)
						{
							event_scale *= centrality_scales.at(ib).second;
							//if (i < 100) std:cout << "found z = " << vertex_scales.at(ib).first << " " <<vertex_scales.at(ib).second<<std::endl;
							break;
						}
					}
				} // centrality reweighting
			} // end of primer != 1 && !full_or_half

			float max_truth = 0;
			float max_reco = 0;
			float min_truth = 0;
			float min_reco = 0;

			if (pt1_truth[isample] >= pt2_truth[isample])
			{
				max_truth = pt1_truth[isample];
				max_reco = pt1_reco[isample];
				min_truth = pt2_truth[isample];
				min_reco = pt2_reco[isample];
			}
			else
			{
				max_truth = pt2_truth[isample];
				max_reco = pt2_reco[isample];
				min_truth = pt1_truth[isample];
				min_reco = pt1_reco[isample];
			}
			
			float pt1_truth_bin = nbins;
			float pt2_truth_bin = nbins;
			float pt1_reco_bin = nbins;
			float pt2_reco_bin = nbins;

			float e1 = max_truth;
			float e2 = min_truth;
			float es1 = max_reco;
			float es2 = min_reco;

			h_truth_lead_sample[isample]->Fill(e1, event_scale);

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
				es1 = es1 + (JES_sys + smear1)*e1;
				es2 = es2 + (JES_sys + smear2)*e2; 
			}
			else
			{
				es1 = es1 + smear1*e1;
				es2 = es2 + smear2*e2; 	  
			}
		
			float maxi = std::max(es1, es2);
			float mini = std::min(es1, es2);

			float maxit = std::max(e1, e2);
			float minit = std::min(e1, e2);

			if (maxi >= ipt_bins[max_reco_bin]) continue;

			if (maxit > truth_leading_cut) h_truth_lead->Fill(e1, event_scale);
			if (minit >  truth_subleading_cut) h_truth_sublead->Fill(e2, event_scale);

			if (maxi >  reco_leading_cut) h_reco_lead->Fill(maxi, event_scale);
			if (mini >  reco_subleading_cut) h_reco_sublead->Fill(mini, event_scale);


			for (int ib = 0; ib < nbins; ib++)
			{
				if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
				{
					pt1_truth_bin = ib;
				}
				if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
				{
					pt2_truth_bin = ib;
				}
				if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
				{
					pt1_reco_bin = ib;
				}
				if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
				{
					pt2_reco_bin = ib;
				}
			}
		
			bool truth_good = (maxit >= sample_boundary[1]/* truth_leading_cut*/ && minit >= truth_subleading_cut && dphi_truth[isample] >= dphicuttruth);

			bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] >= dphicut);

			if (!truth_good && !reco_good) continue;

			if (!primer && !full_or_half)
			{
				int recorrectbin = pt1_truth_bin*nbins + pt2_truth_bin;
				event_scale *= h_flatreweight_pt1pt2->GetBinContent(recorrectbin);
			}
			//RooUnfoldResponse rooResponse_histo(*h_flat_reco_to_response_pt1pt2, *h_flat_truth_to_response, *flat_response_pt1pt2);

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

				// if (!(full_or_half && !fill_response))
				// {
				// 	h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
				// 	h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);

				// 	rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);

				// }

				if (use_for_test)
				{
					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
					{
						h_truth_xj->Fill(minit/maxit, event_scale);
						h_linear_truth_xj->Fill(minit/maxit, event_scale);
					}
					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_e1e2->Fill(e1, e2, event_scale);
					h_e1e2->Fill(e2, e1, event_scale);
				}

				// if (!(full_or_half && fill_response))
				// {
				// 	if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
				// 	{
				// 		h_truth_xj->Fill(minit/maxit, event_scale);
				// 		h_linear_truth_xj->Fill(minit/maxit, event_scale);
				// 	}
				// 	h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	h_e1e2->Fill(e1, e2, event_scale);
				// 	h_e1e2->Fill(e2, e1, event_scale);

				// }
				continue;
			} // end of truth_good && !reco_good
		
			// if (!truth_good && reco_good)
			//   {
			//     h_pt1pt2->Fill(es1, es2, event_scale);
			//     h_pt1pt2->Fill(es2, es1, event_scale);
			
			//     h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
			//     if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
			// 	{
			// 	  h_reco_xj->Fill(mini/maxi, event_scale);
			// 	  h_linear_reco_xj->Fill(mini/maxi, event_scale);
			// 	}
			//     h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
			//     rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
			//     continue;
			//   }
			
			if (match[isample] && reco_good && truth_good)
			{
			
				if (maxit > truth_leading_cut) h_match_truth_lead->Fill(maxit, event_scale);
				if (minit >  truth_subleading_cut) h_match_truth_sublead->Fill(minit, event_scale);
				if (maxi >  reco_leading_cut) h_match_reco_lead->Fill(maxi, event_scale);
				if (mini >  reco_subleading_cut) h_match_reco_sublead->Fill(mini, event_scale);
				
			
				if (use_for_response)
				{
					
					h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
					h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
					h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
					h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
					h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);	      
				
				}
				else if (use_for_test)
				{
					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
					{
						h_truth_xj->Fill(minit/maxit, event_scale);
						h_linear_truth_xj->Fill(minit/maxit, event_scale);
					}
					h_pt1pt2->Fill(es1, es2, event_scale);
					h_pt1pt2->Fill(es2, es1, event_scale);
					h_e1e2->Fill(e1, e2, event_scale);
					h_e1e2->Fill(e2, e1, event_scale);

					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
				}

				// if (!(full_or_half && !fill_response))
				// {
				
				// 	h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
				// 	h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
				// 	h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
				// 	h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
				// 	h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
				// 	h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);	      
					
				// }
				// else if (!(full_or_half && fill_response))
				// {
				// 	if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
				// 	{
				// 		h_truth_xj->Fill(minit/maxit, event_scale);
				// 		h_linear_truth_xj->Fill(minit/maxit, event_scale);
				// 	}
				// 	h_pt1pt2->Fill(es1, es2, event_scale);
				// 	h_pt1pt2->Fill(es2, es1, event_scale);
				// 	h_e1e2->Fill(e1, e2, event_scale);
				// 	h_e1e2->Fill(e2, e1, event_scale);

				// 	h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
				// 	h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				// 	h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
				// 	h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
				// }
			
			
				
				if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
				{
					h_reco_xj->Fill(mini/maxi, event_scale);
					h_linear_reco_xj->Fill(mini/maxi, event_scale);
				}
				h_mbd_vertex->Fill(mbd_vertex[isample], event_scale);
				if (!ispp && has_sumeT[isample]) h_sumeT->Fill(sumeT[isample], event_scale);
				if (!ispp) h_centrality->Fill(centrality[isample], event_scale);
				
				continue;
			} 
			// else
			//   {
			//     h_truth_xj->Fill(e2/e1, event_scale);
			//     h_linear_truth_xj->Fill(e2/e1, event_scale);
			//     h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
			//     h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);

			//     h_pt1pt2->Fill(es1, es2, event_scale);
			//     h_pt1pt2->Fill(es2, es1, event_scale);

			//     h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
			//     h_reco_xj->Fill(mini/maxi, event_scale);
			//     h_linear_reco_xj->Fill(mini/maxi, event_scale);
			//     rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
			//     rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
			//     rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);

			//     h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
			//     h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
			//     h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
			//     h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
			//     continue;
			//   }

			//std::cout << "something not right" << std::endl;
		}
    }

	h_flat_reco_pt1pt2->Scale(.5);
	h_flat_truth_pt1pt2->Scale(.5);

	h_flat_reco_to_response_pt1pt2->Scale(.5);
	h_flat_truth_to_response_pt1pt2->Scale(.5);
	h_flat_response_pt1pt2->Scale(.5);

	int number_of_mins = 1;
	if (minentries == 0) number_of_mins = 20;
  	for (int imin = minentries; imin < minentries + number_of_mins; imin++)       
    {
		std::cout << "Minimum Entries == " << imin << std::endl;
		int nbinsx = h_flat_response_pt1pt2->GetXaxis()->GetNbins();
		int nbinsy = h_flat_response_pt1pt2->GetYaxis()->GetNbins();

		TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
		TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
		if (minentries)
		{
			h_flat_truth_mapping->SetName("h_flat_truth_mapping");
			h_flat_reco_mapping->SetName("h_flat_reco_mapping");
		}
		else
		{
			h_flat_truth_mapping->SetName(Form("h_flat_truth_mapping_min%d", imin));
			h_flat_reco_mapping->SetName(Form("h_flat_reco_mapping_min%d", imin));
		}
      
		h_flat_reco_mapping->Reset();
		h_flat_truth_mapping->Reset();
	
		int nempty_reco = 0;
		std::vector<int> binnumbers_reco{};
		int nempty_truth = 0;
		std::vector<int> binnumbers_truth{};
		int nrecobins = 0;
		int ntruthbins = 0;
		for (int ib = 0; ib < nbins*nbins; ib++)
		{
			if (h_count_flat_reco_pt1pt2->GetBinContent(ib+1) < imin)
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
	  		if (h_count_flat_truth_pt1pt2->GetBinContent(ib+1) < imin)
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
		TH1D *h_flat_truth_skim = new TH1D(Form("h_flat_truth_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
		TH1D *h_flat_reco_skim = new TH1D(Form("h_flat_reco_skim_min%d", imin),"", nrecobins, 0, nrecobins);
		TH1D *h_flat_truth_to_unfold_skim = new TH1D(Form("h_flat_truth_to_unfold_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
		TH1D *h_flat_reco_to_unfold_skim = new TH1D(Form("h_flat_reco_to_unfold_skim_min%d", imin),"", nrecobins, 0, nrecobins);

		TH2D *h_flat_response_skim = new TH2D(Form("h_flat_response_skim_min%d", imin),"", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

		if (minentries)
		{
			h_flat_response_skim->SetName("h_flat_response_skim");
			h_flat_truth_skim->SetName("h_flat_truth_skim");
			h_flat_reco_skim->SetName("h_flat_reco_skim");
			h_flat_truth_to_unfold_skim->SetName("h_flat_truth_to_unfold_skim");
			h_flat_reco_to_unfold_skim->SetName("h_flat_reco_to_unfold_skim");
		}

      	for (int ib = 0; ib < nbins*nbins; ib++)
		{
			int bintruth = binnumbers_truth.at(ib);      
			int binreco = binnumbers_reco.at(ib);
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
					int binreco2 = binnumbers_reco.at(ibr);
					if (binreco2)
					{
						int rbin = h_flat_response_skim->GetBin(binreco2, bintruth);
						int tbin = h_flat_response_pt1pt2->GetBin(ibr+1, ib+1);
						h_flat_response_skim->SetBinContent(rbin, h_flat_response_pt1pt2->GetBinContent(tbin));
						h_flat_response_skim->SetBinError(rbin, h_flat_response_pt1pt2->GetBinError(tbin));
					}
				}
	    	}
		}

		std::cout <<" Nbins skim reco = "<<h_flat_reco_skim->GetNbinsX()<<std::endl;
		std::cout <<" Nbins skim reco = "<<h_flat_truth_skim->GetNbinsX()<<std::endl;
		std::cout <<" Nbins skim reco = "<<h_flat_response_skim->GetXaxis()->GetNbins()<<"  --  " << h_flat_response_skim->GetYaxis()->GetNbins()<<std::endl; 
		RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);
		if (minentries)
		{
			rooResponsehist.SetName("response_noempty");
		}
		else
		{
			rooResponsehist.SetName(Form("response_noempty_min%d", imin));
		}

		if (minentries)
		{
			TH1D* h_flat_unfold_skim[niterations];
			TH1D* h_flat_unfold_pt1pt2[niterations];
			TH1D* h_flat_refold_skim[niterations];
			TH1D* h_flat_refold_pt1pt2[niterations];
			TH1D* h_flat_truth_fold_skim;
			TH1D* h_flat_truth_fold_pt1pt2;
			int niter = 3;
			
			h_flat_truth_fold_skim = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_truth_skim, "hTruth_Folded");

			h_flat_truth_fold_pt1pt2 = (TH1D*) h_flat_truth_pt1pt2->Clone();
			h_flat_truth_fold_pt1pt2->Reset();
			h_flat_truth_fold_pt1pt2->SetName("h_flat_truth_fold_pt1pt2");
			for (int ib = 0; ib < nbins*nbins; ib++)
			{
				int bin = binnumbers_truth.at(ib);
				if (bin)
				{
					h_flat_truth_fold_pt1pt2->SetBinContent(ib+1, h_flat_truth_fold_skim->GetBinContent(bin));
					h_flat_truth_fold_pt1pt2->SetBinError(ib+1, h_flat_truth_fold_skim->GetBinError(bin));
				}
			}

			for (int iter = 0; iter < niterations; iter++ )
			{
		
				RooUnfoldBayes   unfold (&rooResponsehist, h_flat_reco_to_unfold_skim, iter + 1);    // OR

				h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
				std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
				h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
				h_flat_unfold_pt1pt2[iter]->Reset();
				h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
				for (int ib = 0; ib < nbins*nbins; ib++)
				{
					int bin = binnumbers_truth.at(ib);
					if (bin)
					{
						h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin));
						h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin));
					}
				}

				h_flat_refold_skim[iter] = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_unfold_skim[iter], "hRefolded");

				std::cout <<" Nbins skim reco = "<<h_flat_refold_skim[iter]->GetNbinsX()<<std::endl;
				h_flat_refold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
				h_flat_refold_pt1pt2[iter]->Reset();
				h_flat_refold_pt1pt2[iter]->SetName(Form("h_flat_refold_pt1pt2_%d",iter));
				for (int ib = 0; ib < nbins*nbins; ib++)
				{
					int bin = binnumbers_truth.at(ib);
					if (bin)
					{
						h_flat_refold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_refold_skim[iter]->GetBinContent(bin));
						h_flat_refold_pt1pt2[iter]->SetBinError(ib+1, h_flat_refold_skim[iter]->GetBinError(bin));
					}
				}

	    	}


			TCanvas *c = new TCanvas("c","c", 500, 500);

			for (int iter = 0; iter < niterations; iter++)
			{
				dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
				dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
			}
		
			dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
			dlutility::SetLineAtt(h_flat_reco_pt1pt2, kBlue, 2, 1);

			h_flat_truth_pt1pt2->Draw("hist");
			h_flat_reco_pt1pt2->Draw("hist same");
			h_flat_unfold_pt1pt2[niter]->Draw("same p");

			TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_truth_fold = new TH2D("h_pt1pt2_truth_fold", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_unfold[niterations];
			TH2D *h_pt1pt2_refold[niterations];
			for (int iter = 0; iter < niterations; iter++)
			{
				h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
				h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
				h_pt1pt2_refold[iter] = new TH2D("h_pt1pt2_refold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
				h_pt1pt2_refold[iter]->SetName(Form("h_pt1pt2_refold_iter%d", iter));
			}
			TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
			TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_truth_fold = new TH1D("h_xj_truth_fold", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_truth_direct = new TH1D("h_xj_truth_direct", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_unfold[niterations];
			TH1D *h_xj_refold[niterations];
			for (int iter = 0; iter < niterations; iter++)
			{
				h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
				h_xj_refold[iter] = new TH1D(Form("h_xj_refold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
			}

			h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
			h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
			h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");

			histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
			histo_opps::make_sym_pt1pt2(h_flat_truth_fold_pt1pt2, h_pt1pt2_truth_fold, nbins);
			histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
				histo_opps::make_sym_pt1pt2(h_flat_refold_pt1pt2[iter], h_pt1pt2_refold[iter], nbins);
			}

			TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
			cpt1pt2->Divide(3, 1);
			cpt1pt2->cd(1);
			gPad->SetLogz();
			gPad->SetRightMargin(0.2);
			h_pt1pt2_reco->Draw("colz");
			dlutility::DrawSPHENIXpp(0.22, 0.85);
		
			cpt1pt2->cd(2);
			gPad->SetRightMargin(0.2);
			gPad->SetLogz();
			h_pt1pt2_truth->Draw("colz");
			cpt1pt2->cd(3);
			gPad->SetRightMargin(0.2);
			gPad->SetLogz();
			h_pt1pt2_unfold[niter]->Draw("colz");


			cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, diagnostic_name.c_str()));

			histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_pt1pt2_truth_fold, h_xj_truth_fold, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_e1e2, h_xj_truth_direct, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
				histo_opps::project_xj(h_pt1pt2_refold[iter], h_xj_refold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			}
		
			TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
			dlutility::ratioPanelCanvas(cxj);
			cxj->cd(1);
			dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
			dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

			dlutility::SetLineAtt(h_truth_xj, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_truth_xj, kRed, 2, 24);

			dlutility::SetLineAtt(h_linear_truth_xj, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_linear_truth_xj, kRed, 2, 25);

			dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_reco_xj, kBlue, 2, 24);

			dlutility::SetLineAtt(h_linear_reco_xj, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_linear_reco_xj, kBlue, 2, 25);

		
			dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

			dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
		
			histo_opps::normalize_histo(h_xj_truth, nbins);
			histo_opps::normalize_histo(h_xj_reco, nbins);
			histo_opps::normalize_histo(h_truth_xj, nbins);
			histo_opps::normalize_histo(h_reco_xj, nbins);

			//histo_opps::normalize_histo(h_linear_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
			h_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
			//h_linear_reco_xj->Scale(1./h_linear_reco_xj->Integral(), "width");

			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
			}

			dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
			h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
			//h_xj_truth->SetMaximum(5);
			h_xj_truth->Draw("p");
			h_xj_unfold[niter]->Draw("same p");
			h_xj_reco->Draw("hist same");
			h_xj_truth->Draw("same hist");
			h_xj_reco->Draw("p same");
			h_xj_unfold[niter]->Draw("same hist");
			h_xj_unfold[niter]->Draw("same p");
			dlutility::DrawSPHENIXpp(0.22, 0.84);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
			dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
			dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
			TLegend *leg = new TLegend(0.22, 0.3, 0.4, 0.55);
			leg->SetLineWidth(0);
			leg->AddEntry(h_xj_reco, "Pythia Reco");
			leg->AddEntry(h_xj_truth, "Pythia8 Truth");
			leg->AddEntry(h_xj_unfold[niter], "Unfolded");
			leg->Draw("same");
			dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
			
			cxj->cd(2);

			TH1D *h_reco_compare = (TH1D*) h_xj_unfold[niter]->Clone();
			h_reco_compare->Divide(h_xj_truth);
			h_reco_compare->SetTitle(";x_{J}; Unfold / Truth");
			dlutility::SetFont(h_reco_compare, 42, 0.06);
			dlutility::SetLineAtt(h_reco_compare, kBlack, 1,1);
			dlutility::SetMarkerAtt(h_reco_compare, kBlack, 1,8);

		
			h_reco_compare->Draw("p");
			TLine *line = new TLine(0.1, 1, 1, 1);
			line->SetLineStyle(4);
			line->SetLineColor(kRed + 3);
			line->SetLineWidth(2);
			line->Draw("same");

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


			TCanvas *cjet = new TCanvas("cjet","cjet", 700, 500);
			dlutility::createCutCanvas(cjet);
			cjet->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_truth_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead, kRed, 1, 24);

			dlutility::SetLineAtt(h_truth_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_truth_sublead, kBlue, 1, 24);

			dlutility::SetLineAtt(h_reco_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_reco_lead, kRed, 1, 20);

			dlutility::SetLineAtt(h_reco_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_reco_sublead, kBlue, 1, 20);

			h_truth_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
			h_truth_lead->Draw();
			h_truth_sublead->Draw("same");
			h_reco_lead->Draw("same");
			h_reco_sublead->Draw("same");
			cjet->cd(2);
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.1, 0.74, 0, kBlack, 0.1);
			dlutility::drawText("All Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);

			TLegend *leg33 = new TLegend(0.01, 0.13, 0.7, 0.62);
			leg33->SetLineWidth(0);
			leg33->AddEntry(h_reco_lead, "Leading Reco");
			leg33->AddEntry(h_reco_sublead, "Subleading Reco");
			leg33->AddEntry(h_truth_lead, "Leading Truth");
			leg33->AddEntry(h_truth_sublead, "Subleading Truth");
			leg33->SetTextSize(0.1);
			leg33->Draw("same");

			TCanvas *cjetdiv = new TCanvas("cjetdiv","cjetdiv", 700, 500);
			dlutility::createCutCanvas(cjetdiv);
			cjetdiv->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_truth_lead, kBlack, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead, kBlack, 0.8, 24);

			dlutility::SetLineAtt(h_truth_lead_sample[0], kRed, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[0], kRed, 1, 20);
			dlutility::SetLineAtt(h_truth_lead_sample[1], kGreen, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[1], kGreen, 1, 20);
			dlutility::SetLineAtt(h_truth_lead_sample[2], kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[2], kBlue, 1, 20);

			h_truth_lead->SetTitle(";Leading Jet p_{T} [GeV];counts * lumiscale");
			h_truth_lead_sample[0]->SetMinimum(1);
			h_truth_lead_sample[0]->Draw("p");
			h_truth_lead->Draw("p same");
			h_truth_lead_sample[1]->Draw("same p");
			h_truth_lead_sample[2]->Draw("same p");
			h_truth_lead->Draw(" same p");


			TLine  *line1 = new TLine(14, 1, 14, h_truth_lead->GetBinContent(15));
			dlutility::SetLineAtt(line1, kRed, 1.5, 1);
			line1->Draw("same");
			TLine  *lin2 = new TLine(20, 1, 20, h_truth_lead->GetBinContent(21));
			dlutility::SetLineAtt(lin2, kGreen, 1.5, 1);
			lin2->Draw("same");
			TLine  *lin3 = new TLine(30, 1, 30, h_truth_lead->GetBinContent(31));
			dlutility::SetLineAtt(lin3, kBlue, 1.5, 1);
			lin3->Draw("same");
			cjetdiv->cd(2);

			dlutility::DrawSPHENIXppsize(0.05, 0.84, 0.08, 1, 0, 1, "Pythia8");
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.05, 0.74, 0, kBlack, 0.08);
			dlutility::drawText("All Dijet Pairs", 0.05, 0.69, 0, kBlack, 0.08);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.05, 0.64, 0, kBlack, 0.08);

			TLegend *leg3 = new TLegend(0.01, 0.13, 0.7, 0.62);
			leg3->SetLineWidth(0);
			leg3->AddEntry(h_truth_lead, "Combined Sample");
			leg3->AddEntry(h_truth_lead_sample[0], "10 GeV Sample");
			leg3->AddEntry(h_truth_lead_sample[1], "20 GeV Sample");
			leg3->AddEntry(h_truth_lead_sample[2], "30 GeV Sample");
			leg3->SetTextSize(0.08);
			leg3->Draw("same");

			cjetdiv->Print(Form("%s/unfolding_plots/combined_sample_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, diagnostic_name.c_str()));

			TCanvas *cmjet = new TCanvas("cmjet","cmjet", 700, 500);
			dlutility::createCutCanvas(cmjet);
			cmjet->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_match_truth_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_match_truth_lead, kRed, 1, 24);

			dlutility::SetLineAtt(h_match_truth_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_match_truth_sublead, kBlue, 1, 24);

			dlutility::SetLineAtt(h_match_reco_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_match_reco_lead, kRed, 1, 20);

			dlutility::SetLineAtt(h_match_reco_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_match_reco_sublead, kBlue, 1, 20);

			h_match_reco_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
			h_match_reco_lead->Draw();
			h_match_truth_sublead->Draw("same");
			h_match_truth_lead->Draw("same");
			h_match_reco_sublead->Draw("same");
			cmjet->cd(2);
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.1, 0.74, 0, kBlack, 0.1);
			dlutility::drawText("Matched Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);
			TLegend *lef3 = new TLegend(0.01, 0.13, 0.7, 0.62);
			lef3->SetLineWidth(0);
			lef3->AddEntry(h_match_reco_lead, "Leading Reco");
			lef3->AddEntry(h_match_reco_sublead, "Subleading Reco");
			lef3->AddEntry(h_match_truth_lead, "Leading Truth");
			lef3->AddEntry(h_match_truth_sublead, "Subleading Truth");
			lef3->SetTextSize(0.1);
			lef3->Draw("same");


			TCanvas *cresponse = new TCanvas("fd","fd", 500, 500);
			h_flat_response_pt1pt2->Draw("colz");
			TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 500);
			cresponseskim->SetLogz();
			h_flat_response_skim->Draw("colz");
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);


			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
			
			if (primer > 0)
			{
				responsepath += "_PRIMER" + std::to_string(primer);;
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
			TFile *fr = new TFile(responsepath.Data(),"recreate");
			rooResponsehist.Write();
			h_flat_reco_pt1pt2->Write();
			h_flat_truth_pt1pt2->Write();
			h_flat_truth_to_response_pt1pt2->Write();
			h_flat_truth_fold_pt1pt2->Write();
			h_flat_response_pt1pt2->Write();
			h_flat_truth_mapping->Write();
			h_flat_reco_mapping->Write();
			h_flat_truth_skim->Write();
			h_flat_reco_skim->Write();
			h_flat_response_skim->Write();
			h_mbd_vertex->Write();
			h_centrality->Write();
			h_sumeT->Write();
			h_pt1pt2->Write();
			h_e1e2->Write();
			// Full closure drawing
			h_linear_truth_xj->Write();
			h_xj_truth->Write();
			h_xj_truth_fold->Write();
			h_xj_truth_direct->Write();
			h_truth_xj->Write();
			for (int iter = 0 ; iter < niterations; iter++)
			{
				h_flat_unfold_pt1pt2[iter]->Write();
				h_xj_unfold[iter]->Write();
				h_flat_refold_pt1pt2[iter]->Write();
				h_xj_refold[iter]->Write();
			}
			h_xj_reco->Write();
				
			fr->Close();

		}
		else
		{
			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_min_" + std::to_string(imin) + ".root";

			TFile *fr = new TFile(responsepath.Data(),"recreate");
			rooResponsehist.Write();
			
			TH1D *h1 = (TH1D*)h_flat_reco_pt1pt2->Clone();
			TH1D *h2 = (TH1D*)h_flat_truth_pt1pt2->Clone();
			TH2D *h3 = (TH2D*)h_flat_response_pt1pt2->Clone();
			TH1D *h4 = (TH1D*)h_flat_truth_mapping->Clone();
			TH1D *h5 = (TH1D*)h_flat_reco_mapping->Clone();
			TH1D *h6 = (TH1D*)h_flat_truth_skim->Clone();
			TH1D *h7 = (TH1D*)h_flat_reco_skim->Clone();
			TH2D *h8 = (TH2D*)h_flat_response_skim->Clone();
			h1->Write();
			h2->Write();
			h3->Write();
			h4->Write();
			h5->Write();
			h6->Write();
			h7->Write();
			h8->Write();
						
			fr->Close();
		}
		std::cout << "Reco empty: " << nempty_reco << std::endl;
		std::cout << "Truth empty: " << nempty_truth << std::endl;
    }
  	
	return 0;
}
