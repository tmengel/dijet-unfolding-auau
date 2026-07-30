
#include <cstdlib>
#include <iostream>
#include "TTree.h"
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

using std::cout;
using std::endl;

struct jet
{
  int id;
  float pt = 0;
  int ptbin = 0;
  float eta = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
};

const float etacut = 0.8;
const float eta_cut_dijet = 0.8;

const float vertex_cut = 60;

void getBackground(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string configfile = "binning_AA.config",
  const std::string infile =
    "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root")
{

  std::cout << "Staring" << std::endl;
  read_binning rb(configfile.c_str());
  std::cout << "Read" << std::endl;
  Int_t read_nbins = rb.get_nbins();
  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;

  float dphicut_z_low = rb.get_zyam_low();
  float dphicut_z_high = rb.get_zyam_high();

  float low_fit = 0.0;  
  float high_fit = 2.5;  

  float zyam_integral_low = dphicut_z_low;
  float zyam_integral_high = high_fit;//dphicut_z_high;;
  const int n_centrality_bins = rb.get_number_centrality_bins();  
  float icentrality_bins[n_centrality_bins+1];

  rb.get_centrality_bins(icentrality_bins);
  
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


  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

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
  t->SetBranchStatus("jet_pT", 1);
  t->SetBranchStatus("jet_unsub_pT", 1);
  t->SetBranchStatus("jet_E", 1);
  t->SetBranchStatus("jet_unsub_E", 1);
  t->SetBranchStatus("jet_eta", 1);
  t->SetBranchStatus("jet_phi", 1);
  t->SetCacheSize(256 * 1024 * 1024);
  t->AddBranchToCache("*", kTRUE);

  std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_unsub = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;

  int centrality;
  float mbd_vertex_z;

  t->SetBranchAddress("cent", &centrality);
  t->SetBranchAddress("jet_pT", &reco_jet_pt);
  t->SetBranchAddress("jet_unsub_pT", &reco_jet_pt_unsub);
  t->SetBranchAddress("jet_E", &reco_jet_e);
  t->SetBranchAddress("jet_unsub_E", &reco_jet_e_unsub);
  t->SetBranchAddress("jet_eta", &reco_jet_eta);
  t->SetBranchAddress("jet_phi", &reco_jet_phi);
  t->SetBranchAddress("zvrtx", &mbd_vertex_z);



  std::vector<struct jet> myrecojets;
  std::vector<struct jet> myrecojets_all;

  int total_entries = t->GetEntries();

  std::pair<int, float> id_leaders[2];

  TH1D *h_dphi_exclusive_dphi_first = new TH1D("h_dphi_exclusive_dphi_first",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, TMath::Pi()/2., TMath::Pi());
  TH1D *h_dphi_exclusive_dphi_second = new TH1D("h_dphi_exclusive_dphi_second",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, TMath::Pi()/2., TMath::Pi());
  
  TH1D *h_dphi_exclusive_all = new TH1D("h_dphi_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_exclusive_all = new TH1D("h_dphi_eta_exclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_inclusive_all = new TH1D("h_dphi_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_inclusive_all = new TH1D("h_dphi_eta_inclusive_all",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());


  TH1D *h_dphi_exclusive[nbins][nbins]; //new TH1D("h_dphi_exclusive",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_exclusive[nbins][nbins]; //new TH1D("h_dphi_eta_exclusive",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

  TH1D *h_dphi_inclusive[nbins][nbins]; //new TH1D("h_dphi_inclusive",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
  TH1D *h_dphi_eta_inclusive[nbins][nbins]; //new TH1D("h_dphi_eta_inclusive",";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());

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

  TH2D *h_pt1_pt2_Signal = new TH2D("h_pt1_pt2_Signal",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1_pt2_ZYAM = new TH2D("h_pt1_pt2_ZYAM",";#it{p}_{T,1}; #it{p}_{T,2};", nbins, ipt_bins, nbins, ipt_bins);
  
  float njet_lead[nbins];
  for (int i = 0; i < nbins; i++)
    {
      njet_lead[i] = 0;
    }
  float njet_lead_all[nbins];
  for (int i = 0; i < nbins; i++)
    {
      njet_lead_all[i] = 0;
    }

  TF1 *fcut = new TF1("fcut","[0]+[1]*TMath::Exp(-[2]*x)",0.0,100.0);
  fcut->SetParameters(2.5,36.2,0.035);

  for (int i = 0; i < total_entries; i++)
    {
		t->GetEntry(i);

		if (i % 100000 == 0) std::cout << "Event: " << i << " / " << total_entries << "\r" << std::flush;
		if (fabs(mbd_vertex_z) > vertex_cut) continue;

		if (centrality < icentrality_bins[centrality_bin] || centrality >= icentrality_bins[centrality_bin+1]) continue;

		int nrecojets = reco_jet_pt->size();
		// this is the reco index for the leading and subleading jet
		myrecojets.clear();
		myrecojets_all.clear();

		int bad_event = 0;
		bool found_reco_dijet = false;
		int njet_good = 0;
		int njet_good_all = 0;
		float cut_value = fcut->Eval(centrality);
      	for (int j = 0; j < nrecojets;j++)
		{
			if (reco_jet_pt->at(j) < reco_subleading_cut) continue;
			if (reco_jet_e->at(j) < 0) continue;
			if (reco_jet_e_unsub->at(j) < 0) continue;
			if (fabs(reco_jet_eta->at(j)) > etacut) continue;
			float pt_unsub = reco_jet_pt_unsub->at(j) - reco_jet_pt->at(j);

	  		if (pt_unsub > cut_value)
			{
			bad_event = 1;
			continue;
			}

	  njet_good++;

	  float es1 = reco_jet_pt->at(j);
	  int bin = nbins;
	  if (es1 > ipt_bins[nbins] ) continue;
	  for (int ib = 0; ib < nbins; ib++)
	    {

	      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
		{
		  bin = ib;
		}
	    }

	  struct jet tempjet1;

	  tempjet1.pt = es1;
	  tempjet1.ptbin = bin;
	  tempjet1.eta = reco_jet_eta->at(j);
	  tempjet1.phi = reco_jet_phi->at(j);
	  tempjet1.id = j;
	  myrecojets.push_back(tempjet1);

	}
      
      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt->at(j) < reco_subleading_cut) continue;
	  if (reco_jet_e->at(j) < 0) continue;
	  if (reco_jet_e_unsub->at(j) < 0) continue;
	  if (fabs(reco_jet_eta->at(j)) > etacut) continue;
	  float pt_unsub = reco_jet_pt_unsub->at(j) - reco_jet_pt->at(j);

	  if (pt_unsub > cut_value)
	    {
	      bad_event = 1;
	      continue;
	    }

	  njet_good_all++;

	  float es1 = reco_jet_pt->at(j);
	  int bin = nbins;
	  if (es1 > ipt_bins[nbins] ) continue;
	  for (int ib = 0; ib < nbins; ib++)
	    {

	      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
		{
		  bin = ib;
		}
	    }

	  struct jet tempjet1;

	  tempjet1.pt = es1;
	  tempjet1.ptbin = bin;
	  tempjet1.eta = reco_jet_eta->at(j);
	  tempjet1.phi = reco_jet_phi->at(j);
	  tempjet1.id = j;
	  myrecojets_all.push_back(tempjet1);

	}

      if (bad_event) continue;
      
      if (njet_good >= 1)
	{

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) {return a.pt > b.pt; });


	  double dphir = 0;
	  double detar = 0;
	  auto leading_iter = myrecojets.begin();

	
	  // Check if the reco jet satisfies the cuts
	  // check dphi
	  auto subleading_iter = myrecojets.begin() + 1;

	  if (leading_iter->pt < reco_leading_cut) continue;

	  njet_lead[leading_iter->ptbin] += 1.0;
      
	  if (subleading_iter != myrecojets.end())
	    {
	      //detar = fabs(leading_iter->eta - subleading_iter->eta);
	      dphir = fabs(leading_iter->phi - subleading_iter->phi);

	      if (dphir > TMath::Pi())
		{
		  dphir = 2*TMath::Pi() - dphir;
		}
	  
	      if (dphir >= dphicut)
		{
		  h_pt1_pt2_Signal->Fill(leading_iter->pt, subleading_iter->pt);
		}

	      if (dphir >= dphicut_z_low && dphir < dphicut_z_high)
		{
		  h_pt1_pt2_ZYAM->Fill(leading_iter->pt, subleading_iter->pt);
		}

	      h_dphi_exclusive_all->Fill(dphir);
	      
	      h_dphi_exclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);
	      
	    }

	  while (subleading_iter != myrecojets.end())
	    {
	      //detar = fabs(leading_iter->eta - subleading_iter->eta);
	      dphir = fabs(leading_iter->phi - subleading_iter->phi);
	      if (dphir > TMath::Pi())
		{
		  dphir = 2*TMath::Pi() - dphir;
		}


	      h_dphi_inclusive_all->Fill(dphir);
	      h_dphi_inclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);

	      subleading_iter++;
	    }
	}
      if (njet_good_all >= 1)
	{

	  std::sort(myrecojets_all.begin(), myrecojets_all.end(), [] (auto a, auto b) {return a.pt > b.pt; });

	  double dphir = 0;
	  double detar = 0;
	  
	  auto leading_iter = myrecojets_all.begin();

	
	  // Check if the reco jet satisfies the cuts
	  // check dphi
	  auto subleading_iter = myrecojets_all.begin() + 1;

	  if (leading_iter->pt < reco_leading_cut) continue;

	  njet_lead_all[leading_iter->ptbin] += 1.0;
      
	  if (subleading_iter != myrecojets_all.end())
	    {
	      detar = fabs(leading_iter->eta - subleading_iter->eta);
	      dphir = fabs(leading_iter->phi - subleading_iter->phi);
	      if (dphir > TMath::Pi())
		{
		  dphir = 2*TMath::Pi() - dphir;
		}

	  
	      if (fabs(detar) > eta_cut_dijet)
		{
		  h_dphi_eta_exclusive_all->Fill(dphir);	  
		  h_dphi_eta_exclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);	  
		}
	    }

	  while (subleading_iter != myrecojets_all.end())
	    {
	      detar = fabs(leading_iter->eta - subleading_iter->eta);
	      dphir = fabs(leading_iter->phi - subleading_iter->phi);
	      if (dphir > TMath::Pi())
		{
		  dphir = 2*TMath::Pi() - dphir;
		}


	      if (fabs(detar) > eta_cut_dijet)
		{
		  h_dphi_eta_inclusive_all->Fill(dphir);	  
		  h_dphi_eta_inclusive[leading_iter->ptbin][subleading_iter->ptbin]->Fill(dphir);	  
		}
	      subleading_iter++;
	    }
	}
    }

  // fit to get the v2,v3 modulation of background
  // then calculate the ZYAM that falls in the signal range
  // Bakcground to signal ratio
  
  TF1 *ffit = new TF1("ffit", "[0] * ( 1 + 2 * [1] * cos(2.0 * x) + 2 * [2] * cos(3.0 * x))", low_fit, high_fit);
  TF1 *ffitline = new TF1("ffitline", "[0]", low_fit, high_fit);

  TF1 *fits[nbins][nbins];
  TF1 *fits_o[nbins][nbins];

  ffit->SetParLimits(0, 0, 10000);
  ffit->SetParLimits(1, 0, 0.5);
  ffit->SetParLimits(2, 0, 0.5);

  TH2D *h_pt1_pt2_sub = (TH2D*) h_pt1_pt2_Signal->Clone();
  h_pt1_pt2_sub->SetName("h_pt1_pt2_sub");

  TH2D *h_pt1_pt2_ratio = (TH2D*) h_pt1_pt2_Signal->Clone();
  h_pt1_pt2_ratio->SetName("h_pt1_pt2_ratio");
  

  for (int i = 0; i < nbins; i++)
    {
      for (int j = 0; j < nbins; j++)
	{
	  std::cout << " ---------- " << i << " - " << j << " -----------" << std::endl;
	  
	  int nbins_over = h_dphi_eta_inclusive[i][j]->FindBin(high_fit) - h_dphi_eta_inclusive[i][j]->FindBin(low_fit);
	  int ncounts = h_dphi_eta_inclusive[i][j]->Integral(h_dphi_eta_inclusive[i][j]->FindBin(low_fit), h_dphi_eta_inclusive[i][j]->FindBin(high_fit));
	  if (ncounts < nbins_over)
	    {
	      fits[i][j] = (TF1*) ffitline->Clone();
	      fits[i][j]->SetName(Form("ffit_%d_%d", i, j));

	      fits[i][j]->SetParameter(0, 0);
	      ffitline->SetParameter(0, 0);

	      h_dphi_eta_inclusive[i][j]->Fit("ffitline", "0RlQ");

	      fits_o[i][j] = (TF1*) ffitline->Clone();
	      fits_o[i][j]->SetName(Form("ffit_o_%d_%d", i, j));

	      fits_o[i][j]->SetParameter(0, ffitline->GetParameter(0));
	      
	      if (njet_lead[i] > 0)
		{
		  //h_dphi_inclusive[i][j]->Scale(1./njet_lead[i]);
		  // h_dphi_eta_inclusive[i][j]->Scale(1./njet_lead[i]);
		}

	      int first_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_low);
	      int second_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_high);

	      float bin_diff = second_bin - first_bin;

	      float dphi_zyam_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin)/bin_diff;
	      float dphi_flow_integral = ffitline->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high -  zyam_integral_low);


	      float par1 = ffitline->GetParameter(0);

	      float integral_diff = dphi_zyam_integral - dphi_flow_integral;
	      par1 += integral_diff;

	      fits[i][j]->SetParameter(0,par1);
	      fits[i][j]->SetRange(0, TMath::Pi());

	      first_bin = h_dphi_exclusive[i][j]->FindBin(dphicut);
	      second_bin = -1;

	      float nbins_diff = (TMath::Pi() - dphicut)/(TMath::Pi()/32.);
	      float nbins_zyam = (zyam_integral_high - zyam_integral_low)/(TMath::Pi()/32.);
	      float flow_signal = nbins_diff*fits[i][j]->Integral(dphicut, TMath::Pi())/(TMath::Pi() - dphicut);
	      float flow_zyam = nbins_zyam * ffitline->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high - zyam_integral_low);;

	      float dphi_signal_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin);

	      int gbin = h_pt1_pt2_ZYAM->GetBin(i+1, j+1);

	      if (flow_zyam > 0 && flow_signal > 0)
		{
		  h_pt1_pt2_ZYAM->SetBinContent(gbin, h_pt1_pt2_ZYAM->GetBinContent(gbin) * flow_signal/flow_zyam);
		}
	      else
		{
		  h_pt1_pt2_ZYAM->SetBinContent(gbin, 0);
		}


	      if (j > i) continue;
	      std::cout << "Fitting Summary: " << i << "-" << j << std::endl;
	      std::cout << "    Parameters fit: " << ffitline->GetParameter(0) <<  std::endl;
	      std::cout << "    Parameters exc: " << fits[i][j]->GetParameter(0) << std::endl;
	      std::cout << "    Integrals: Flow in signal : " << flow_signal << std::endl;
	      std::cout << "               Flow in zyam   : " << flow_zyam << std::endl;
	      std::cout << "               Hist in exc    : " << dphi_signal_integral << std::endl;
	      std::cout << "               Hist in inc    : " << dphi_zyam_integral << " (per_bin) "  << std::endl;
	      std::cout << "    Total signal : " << dphi_signal_integral - flow_signal << std::endl;
	      float sig = dphi_signal_integral - flow_signal;
	      if (sig > 0)
		{
		  h_pt1_pt2_sub->SetBinContent(gbin, sig);	      
		}
	      else
		{
		  h_pt1_pt2_sub->SetBinContent(gbin, 0);	      
		}
	      h_pt1_pt2_ratio->SetBinContent(gbin, sig/dphi_signal_integral);

	    }
	  else
	    {
	      fits[i][j] = (TF1*) ffit->Clone();
	      fits[i][j]->SetName(Form("ffit_%d_%d", i, j));

	      fits[i][j]->SetParameters(0,0,0);
	      ffit->SetParameters(0,0,0);

	      h_dphi_eta_inclusive[i][j]->Fit("ffit", "0RlQ");

	      fits_o[i][j] = (TF1*) ffit->Clone();
	      fits_o[i][j]->SetName(Form("ffit_o_%d_%d", i, j));
	      fits_o[i][j]->SetParameters(ffit->GetParameter(0), ffit->GetParameter(1), ffit->GetParameter(2));
	      // normalize
	      if (njet_lead[i] > 0)
		{
		  //h_dphi_inclusive[i][j]->Scale(1./njet_lead[i]);
		  //h_dphi_eta_inclusive[i][j]->Scale(1./njet_lead[i]);
		}

	      int first_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_low);
	      int second_bin = h_dphi_exclusive[i][j]->FindBin(zyam_integral_high);

	      float bin_diff = second_bin - first_bin;

	      float dphi_zyam_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin)/bin_diff;
	      float dphi_flow_integral = ffit->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high -  zyam_integral_low);


	      float par1 = ffit->GetParameter(0);

	      if (std::isfinite(dphi_flow_integral) && dphi_flow_integral > 0)
	        {
	          par1 *= dphi_zyam_integral / dphi_flow_integral;
	        }
	      else
	        {
	          // Sparse pT bins can produce a zero constant fit.  Preserve the
	          // intended normalization to the observed sideband without 0/0.
	          par1 = std::isfinite(dphi_zyam_integral)
	            ? dphi_zyam_integral : 0.0F;
	        }

	      fits[i][j]->SetParameters(par1, ffit->GetParameter(1), ffit->GetParameter(2));
	      fits[i][j]->SetRange(0, TMath::Pi());

	      first_bin = h_dphi_exclusive[i][j]->FindBin(dphicut);
	      second_bin = -1;

	      float nbins_diff = (TMath::Pi() - dphicut)/(TMath::Pi()/32.);
	      float nbins_zyam = (zyam_integral_high - zyam_integral_low)/(TMath::Pi()/32.);
	      float flow_signal = nbins_diff*fits[i][j]->Integral(dphicut, TMath::Pi())/(TMath::Pi() - dphicut);
	      float flow_zyam = nbins_zyam * ffit->Integral(zyam_integral_low, zyam_integral_high)/(zyam_integral_high - zyam_integral_low);;

	      float dphi_signal_integral = h_dphi_exclusive[i][j]->Integral(first_bin, second_bin);

	      int gbin = h_pt1_pt2_ZYAM->GetBin(i+1, j+1);

	      if (flow_zyam > 0 && flow_signal > 0)
		{
		  h_pt1_pt2_ZYAM->SetBinContent(gbin, h_pt1_pt2_ZYAM->GetBinContent(gbin) * flow_signal/flow_zyam);
		}
	      else
		{
		  h_pt1_pt2_ZYAM->SetBinContent(gbin, 0);
		}

	      if (j > i) continue;
	      std::cout << "Fitting Summary: " << i << "-" << j << std::endl;
	      std::cout << "    Parameters fit: " << ffit->GetParameter(0) <<  " ,  " << ffit->GetParameter(1) <<  " ,  " << ffit->GetParameter(2) <<  std::endl;
	      std::cout << "    Parameters exc: " << fits[i][j]->GetParameter(0) <<  " ,  " << fits[i][j]->GetParameter(1) <<  " ,  " << fits[i][j]->GetParameter(2) <<  std::endl;
	      std::cout << "    Integrals: Flow in signal : " << flow_signal << std::endl;
	      std::cout << "               Flow in zyam   : " << flow_zyam << std::endl;
	      std::cout << "               Hist in exc    : " << dphi_signal_integral << std::endl;
	      std::cout << "               Hist in inc    : " << dphi_zyam_integral << " (per_bin) "  << std::endl;
	      std::cout << "    Total signal : " << dphi_signal_integral - flow_signal << std::endl;
	      float sig = dphi_signal_integral - flow_signal;
	      if (sig > 0)
		{
		  h_pt1_pt2_sub->SetBinContent(gbin, sig);	      
		}
	      else
		{
		  h_pt1_pt2_sub->SetBinContent(gbin, 0);	      
		}

	      h_pt1_pt2_ratio->SetBinContent(
	        gbin, dphi_signal_integral > 0 ? sig/dphi_signal_integral : 0.0F);
	    }
	}
    }
  

  
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/background_hists_AA_cent_" + std::to_string(centrality_bin) + "_r0" + std::to_string(cone_size);
  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");

  h_pt1_pt2_sub->Write();
  h_pt1_pt2_Signal->Write();
  h_pt1_pt2_ZYAM->Write();
  h_pt1_pt2_ratio->Write();
  for (int i = 0; i < nbins; i++)
    {
      for (int j = 0; j < nbins; j++)
	{

	  fits[i][j]->Write();
	  fits_o[i][j]->Write();
	  h_dphi_exclusive[i][j]->Write();//new TH1D(Form("_dphi_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
	  h_dphi_eta_exclusive[i][j]->Write();//new TH1D(Form("_dphi_eta_exclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
	  
	  h_dphi_inclusive[i][j]->Write();//new TH1D(Form("_dphi_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
	  h_dphi_eta_inclusive[i][j]->Write();//new TH1D(Form("_dphi_eta_inclusive_%d_%d", i, j),";#Delta#phi; #frac{1}{N_{lead}}#frac{dN_{pair}}{d#Delta#phi}", 32, 0, TMath::Pi());
	}
    }
  fout->Close();

}
