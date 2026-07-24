#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include <algorithm>
#include <cmath>

const int color_unfold = kAzure - 6;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

const int color_ZYAM = kOrange + 7;
const float marker_ZYAM = 20;
const float msize_ZYAM = 0.9;
const float lsize_ZYAM = 1.1;

const int color_INCLUSIVE = kRed + 1;
const float marker_INCLUSIVE = 20;
const float msize_INCLUSIVE = 0.9;
const float lsize_INCLUSIVE = 1.1;

const int color_prior = kRed - 2;
const float marker_prior = 20;
const float msize_prior = 0.9;
const float lsize_prior = 1.1;

const int color_negjer = kPink + 2;
const float marker_negjer = 20;
const float msize_negjer = 0.9;
const float lsize_negjer = 1.1;

const int color_posjer = kMagenta + 2;
const float marker_posjer = 20;
const float msize_posjer = 0.9;
const float lsize_posjer = 1.1;

const int color_pjes = kCyan + 1;
const float marker_pjes = 20;
const float msize_pjes = 0.9;
const float lsize_pjes = 1.1;

const int color_njes = kGreen  - 2;
const float marker_njes = 20;
const float msize_njes = 0.9;
const float lsize_njes = 1.1;

const int color_pythia = kBlack;
const float msize_pythia = 0.7;
const float marker_pythia = 21;
const float lsize_pythia = 1.1;
const int color_reco = kRed - 2;
const float marker_reco = 21;
const float msize_reco = 0.7;
const float lsize_reco = 1.1;
const int color_data = kBlue - 9;
const float marker_data = 20;
const float msize_data = 0.9;
const float lsize_data = 1.1;
void drawSys_AA(const int cone_size = 3, const int centrality_bin = 0)
{

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(std::getenv("AUAU_CONFIG"));

  bool ispp = (centrality_bin < 0);
  std::string system_string = (ispp?"pp":"AA_cent_" + std::to_string(centrality_bin));
  
  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins+1];
  rb.get_centrality_bins(icentrality_bins);

  Double_t first_xj = rb.get_first_xj();
  int first_bin = 0;
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];
  double  dxj_bins[nbins+1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] >= first_xj && first_bin == 0) first_bin = i;
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

  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int measure_bins[10] = {0};
  int subleading_measure_bins[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
      subleading_measure_bins[ir] = rb.get_subleading_measure_region(ir);
      std::cout << measure_bins[ir] << " -- " <<  subleading_measure_bins[ir] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  
  const int niterations = 10;
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_%s_r%02d_nominal.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TProfile *h_xj_rms[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_xj_rms[iter] = (TProfile*) finu->Get(Form("hp_xj_%d", iter));
    }
  TFile *fin = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_nominal.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }

  // JES +
  TFile *finpjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_posJES.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_pjes_flat_data_pt1pt2 = (TH1D*) finpjes->Get("h_data_flat_pt1pt2");
  h_pjes_flat_data_pt1pt2->SetName("h_pjes_flat_data_pt1pt2");
  TH1D *h_pjes_flat_reco_pt1pt2 = (TH1D*) finpjes->Get("h_reco_flat_pt1pt2");
  h_pjes_flat_reco_pt1pt2->SetName("h_pjes_flat_reco_pt1pt2");
  TH1D *h_pjes_flat_truth_pt1pt2 = (TH1D*) finpjes->Get("h_truth_flat_pt1pt2");
  h_pjes_flat_truth_pt1pt2->SetName("h_pjes_flat_truth_pt1pt2");
  TH1D *h_pjes_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_flat_unfold_pt1pt2[iter] = (TH1D*) finpjes->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_pjes_flat_unfold_pt1pt2[iter]->SetName(Form("h_pjes_flat_unfold_pt1pt2_%d", iter));
    }
  // JES +
  TFile *finnjes = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_negJES.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_njes_flat_data_pt1pt2 = (TH1D*) finnjes->Get("h_data_flat_pt1pt2");
  h_njes_flat_data_pt1pt2->SetName("h_njes_flat_data_pt1pt2");
  TH1D *h_njes_flat_reco_pt1pt2 = (TH1D*) finnjes->Get("h_reco_flat_pt1pt2");
  h_njes_flat_reco_pt1pt2->SetName("h_njes_flat_reco_pt1pt2");
  TH1D *h_njes_flat_truth_pt1pt2 = (TH1D*) finnjes->Get("h_truth_flat_pt1pt2");
  h_njes_flat_truth_pt1pt2->SetName("h_njes_flat_truth_pt1pt2");
  TH1D *h_njes_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_flat_unfold_pt1pt2[iter] = (TH1D*) finnjes->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_njes_flat_unfold_pt1pt2[iter]->SetName(Form("h_njes_flat_unfold_pt1pt2_%d", iter));
    }
  // posJER
  TFile *finposjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_posJER.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_posjer_flat_data_pt1pt2 = (TH1D*) finposjer->Get("h_data_flat_pt1pt2");
  h_posjer_flat_data_pt1pt2->SetName("h_posjer_flat_data_pt1pt2");
  TH1D *h_posjer_flat_reco_pt1pt2 = (TH1D*) finposjer->Get("h_reco_flat_pt1pt2");
  h_posjer_flat_reco_pt1pt2->SetName("h_posjer_flat_reco_pt1pt2");
  TH1D *h_posjer_flat_truth_pt1pt2 = (TH1D*) finposjer->Get("h_truth_flat_pt1pt2");
  h_posjer_flat_truth_pt1pt2->SetName("h_posjer_flat_truth_pt1pt2");
  TH1D *h_posjer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_flat_unfold_pt1pt2[iter] = (TH1D*) finposjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_posjer_flat_unfold_pt1pt2[iter]->SetName(Form("h_posjer_flat_unfold_pt1pt2_%d", iter));
    }
  // negJER
  TFile *finnegjer = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_negJER.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_negjer_flat_data_pt1pt2 = (TH1D*) finnegjer->Get("h_data_flat_pt1pt2");
  h_negjer_flat_data_pt1pt2->SetName("h_negjer_flat_data_pt1pt2");
  TH1D *h_negjer_flat_reco_pt1pt2 = (TH1D*) finnegjer->Get("h_reco_flat_pt1pt2");
  h_negjer_flat_reco_pt1pt2->SetName("h_negjer_flat_reco_pt1pt2");
  TH1D *h_negjer_flat_truth_pt1pt2 = (TH1D*) finnegjer->Get("h_truth_flat_pt1pt2");
  h_negjer_flat_truth_pt1pt2->SetName("h_negjer_flat_truth_pt1pt2");
  TH1D *h_negjer_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_flat_unfold_pt1pt2[iter] = (TH1D*) finnegjer->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_negjer_flat_unfold_pt1pt2[iter]->SetName(Form("h_negjer_flat_unfold_pt1pt2_%d", iter));
    }

  /* // VTX */
  /* TFile *finvtx = new TFile("%s/unfolding_hists_VTX.root","r"); */

  /* TH1D *h_vtx_flat_data_pt1pt2 = (TH1D*) finvtx->Get("h_data_flat_pt1pt2"); */
  /* h_vtx_flat_data_pt1pt2->SetName("h_vtx_flat_data_pt1pt2"); */
  /* TH1D *h_vtx_flat_reco_pt1pt2 = (TH1D*) finvtx->Get("h_reco_flat_pt1pt2"); */
  /* h_vtx_flat_reco_pt1pt2->SetName("h_vtx_flat_reco_pt1pt2"); */
  /* TH1D *h_vtx_flat_truth_pt1pt2 = (TH1D*) finvtx->Get("h_truth_flat_pt1pt2"); */
  /* h_vtx_flat_truth_pt1pt2->SetName("h_vtx_flat_truth_pt1pt2"); */
  /* TH1D *h_vtx_flat_unfold_pt1pt2[niterations]; */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_flat_unfold_pt1pt2[iter] = (TH1D*) finvtx->Get(Form("h_flat_unfold_pt1pt2_%d", iter)); */
  /*     h_vtx_flat_unfold_pt1pt2[iter]->SetName(Form("h_vtx_flat_unfold_pt1pt2_%d", iter)); */
  /*   } */

  // Coherent v22/v33 COMB envelope

  TFile *finZYAM = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_COMBDown.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");
  TFile *finCOMBUp = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_COMBUp.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_ZYAM_flat_data_pt1pt2 = (TH1D*) finZYAM->Get("h_data_flat_pt1pt2");
  
  TH1D *h_ZYAM_flat_reco_pt1pt2 = (TH1D*) finZYAM->Get("h_reco_flat_pt1pt2");
  TH1D *h_ZYAM_flat_truth_pt1pt2 = (TH1D*) finZYAM->Get("h_truth_flat_pt1pt2");

  TH1D *h_ZYAM_flat_unfold_pt1pt2[niterations];
  TH1D *h_COMBUp_flat_unfold_pt1pt2[niterations];
  if (!ispp)
    {
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_ZYAM_flat_unfold_pt1pt2[iter] = (TH1D*) finZYAM->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_ZYAM_flat_unfold_pt1pt2[iter]->SetName(Form("h_ZYAM_flat_unfold_pt1pt2_%d", iter));
	  h_COMBUp_flat_unfold_pt1pt2[iter] = (TH1D*) finCOMBUp->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_COMBUp_flat_unfold_pt1pt2[iter]->SetName(Form("h_COMBUp_flat_unfold_pt1pt2_%d", iter));
	  
	}
      h_ZYAM_flat_truth_pt1pt2->SetName("h_ZYAM_flat_truth_pt1pt2");
      h_ZYAM_flat_data_pt1pt2->SetName("h_ZYAM_flat_data_pt1pt2");
      h_ZYAM_flat_reco_pt1pt2->SetName("h_ZYAM_flat_reco_pt1pt2");
    }
  // INCLUSIVE 

  TFile *finINCLUSIVE = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_INCLUSIVE.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_INCLUSIVE_flat_data_pt1pt2 = (TH1D*) finINCLUSIVE->Get("h_data_flat_pt1pt2");
  
  TH1D *h_INCLUSIVE_flat_reco_pt1pt2 = (TH1D*) finINCLUSIVE->Get("h_reco_flat_pt1pt2");
  TH1D *h_INCLUSIVE_flat_truth_pt1pt2 = (TH1D*) finINCLUSIVE->Get("h_truth_flat_pt1pt2");

  TH1D *h_INCLUSIVE_flat_unfold_pt1pt2[niterations];
  if (!ispp)
    {
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_INCLUSIVE_flat_unfold_pt1pt2[iter] = (TH1D*) finINCLUSIVE->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_INCLUSIVE_flat_unfold_pt1pt2[iter]->SetName(Form("h_INCLUSIVE_flat_unfold_pt1pt2_%d", iter));
	  
	}
      h_INCLUSIVE_flat_truth_pt1pt2->SetName("h_INCLUSIVE_flat_truth_pt1pt2");
      h_INCLUSIVE_flat_data_pt1pt2->SetName("h_INCLUSIVE_flat_data_pt1pt2");
      h_INCLUSIVE_flat_reco_pt1pt2->SetName("h_INCLUSIVE_flat_reco_pt1pt2");
    }
  // PRIOR 
  TFile *finprior = new TFile(Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_PRIOR.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size),"r");

  TH1D *h_prior_flat_data_pt1pt2 = (TH1D*) finprior->Get("h_data_flat_pt1pt2");
  h_prior_flat_data_pt1pt2->SetName("h_prior_flat_data_pt1pt2");
  TH1D *h_prior_flat_reco_pt1pt2 = (TH1D*) finprior->Get("h_reco_flat_pt1pt2");
  h_prior_flat_reco_pt1pt2->SetName("h_prior_flat_reco_pt1pt2");
  TH1D *h_prior_flat_truth_pt1pt2 = (TH1D*) finprior->Get("h_truth_flat_pt1pt2");
  h_prior_flat_truth_pt1pt2->SetName("h_prior_flat_truth_pt1pt2");
  TH1D *h_prior_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_flat_unfold_pt1pt2[iter] = (TH1D*) finprior->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_prior_flat_unfold_pt1pt2[iter]->SetName(Form("h_prior_flat_unfold_pt1pt2_%d", iter));
    }

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_xj_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_data_range[irange] = new TH1D(Form("h_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  // pJES
  TH2D *h_pjes_pt1pt2_reco = new TH2D("h_pjes_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pjes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_pt1pt2_unfold[iter] = new TH2D("h_pjes_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pjes_pt1pt2_unfold[iter]->SetName(Form("h_pjes_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_pjes_xj_reco = new TH1D("h_pjes_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_pjes_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_pjes_xj_unfold[iter] = new TH1D(Form("h_pjes_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  
  TH1D *h_pjes_xj_reco_range[mbins];
  TH1D *h_pjes_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjes_xj_reco_range[irange] = new TH1D(Form("h_pjes_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjes_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjes_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  // nJES

  TH2D *h_njes_pt1pt2_reco = new TH2D("h_njes_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_njes_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_pt1pt2_unfold[iter] = new TH2D("h_njes_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_njes_pt1pt2_unfold[iter]->SetName(Form("h_njes_pt1pt2_unfold_iter%d", iter));
    }
  TH1D *h_njes_xj_reco = new TH1D("h_njes_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_njes_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_njes_xj_unfold[iter] = new TH1D(Form("h_njes_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  
  TH1D *h_njes_xj_reco_range[mbins];
  TH1D *h_njes_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njes_xj_reco_range[irange] = new TH1D(Form("h_njes_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njes_xj_unfold_range[irange][iter] = new TH1D(Form("h_njes_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  //POSJER
  TH2D *h_posjer_pt1pt2_reco = new TH2D("h_posjer_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_posjer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_pt1pt2_unfold[iter] = new TH2D("h_posjer_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_posjer_pt1pt2_unfold[iter]->SetName(Form("h_posjer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_posjer_xj_reco = new TH1D("h_posjer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_posjer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_posjer_xj_unfold[iter] = new TH1D(Form("h_posjer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_posjer_xj_reco_range[mbins];
  TH1D *h_posjer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_posjer_xj_reco_range[irange] = new TH1D(Form("h_posjer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_posjer_xj_unfold_range[irange][iter] = new TH1D(Form("h_posjer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  //NEGJER
  TH2D *h_negjer_pt1pt2_reco = new TH2D("h_negjer_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_negjer_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_pt1pt2_unfold[iter] = new TH2D("h_negjer_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_negjer_pt1pt2_unfold[iter]->SetName(Form("h_negjer_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_negjer_xj_reco = new TH1D("h_negjer_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_negjer_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_negjer_xj_unfold[iter] = new TH1D(Form("h_negjer_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_negjer_xj_reco_range[mbins];
  TH1D *h_negjer_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_negjer_xj_reco_range[irange] = new TH1D(Form("h_negjer_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_negjer_xj_unfold_range[irange][iter] = new TH1D(Form("h_negjer_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  /* //VTX */
  /* TH2D *h_vtx_pt1pt2_reco = new TH2D("h_vtx_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins); */
  /* TH2D *h_vtx_pt1pt2_unfold[niterations]; */

  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_pt1pt2_unfold[iter] = new TH2D("h_vtx_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins); */
  /*     h_vtx_pt1pt2_unfold[iter]->SetName(Form("h_vtx_pt1pt2_unfold_iter%d", iter)); */
  /*   } */

  /* TH1D *h_vtx_xj_reco = new TH1D("h_vtx_xj_reco", ";x_{J};", nbins, ixj_bins); */
  /* TH1D *h_vtx_xj_unfold[niterations]; */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_vtx_xj_unfold[iter] = new TH1D(Form("h_vtx_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins); */
  /*   } */

  /* TH1D *h_vtx_xj_reco_range[mbins]; */
  /* TH1D *h_vtx_xj_unfold_range[mbins][niterations]; */
  /* for (int irange = 0; irange < mbins; irange++) */
  /*   { */
  /*     h_vtx_xj_reco_range[irange] = new TH1D(Form("h_vtx_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins); */
  /*     for (int iter = 0; iter < niterations; iter++) */
  /* 	{ */
  /* 	  h_vtx_xj_unfold_range[irange][iter] = new TH1D(Form("h_vtx_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins); */
  /* 	} */
  /*   } */
  //ZYAM
  TH2D *h_ZYAM_pt1pt2_reco = new TH2D("h_ZYAM_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_ZYAM_pt1pt2_unfold[niterations];
  TH2D *h_COMBUp_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_ZYAM_pt1pt2_unfold[iter] = new TH2D("h_ZYAM_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_ZYAM_pt1pt2_unfold[iter]->SetName(Form("h_ZYAM_pt1pt2_unfold_iter%d", iter));
      h_COMBUp_pt1pt2_unfold[iter] = new TH2D("h_COMBUp_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_COMBUp_pt1pt2_unfold[iter]->SetName(Form("h_COMBUp_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_ZYAM_xj_reco = new TH1D("h_ZYAM_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_ZYAM_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_ZYAM_xj_unfold[iter] = new TH1D(Form("h_ZYAM_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_ZYAM_xj_reco_range[mbins];
  TH1D *h_ZYAM_xj_unfold_range[mbins][niterations];
  TH1D *h_COMBUp_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_ZYAM_xj_reco_range[irange] = new TH1D(Form("h_ZYAM_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_ZYAM_xj_unfold_range[irange][iter] = new TH1D(Form("h_ZYAM_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_COMBUp_xj_unfold_range[irange][iter] = new TH1D(Form("h_COMBUp_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //INCLUSIVE
  TH2D *h_INCLUSIVE_pt1pt2_reco = new TH2D("h_INCLUSIVE_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_INCLUSIVE_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_INCLUSIVE_pt1pt2_unfold[iter] = new TH2D("h_INCLUSIVE_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_INCLUSIVE_pt1pt2_unfold[iter]->SetName(Form("h_INCLUSIVE_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_INCLUSIVE_xj_reco = new TH1D("h_INCLUSIVE_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_INCLUSIVE_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_INCLUSIVE_xj_unfold[iter] = new TH1D(Form("h_INCLUSIVE_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_INCLUSIVE_xj_reco_range[mbins];
  TH1D *h_INCLUSIVE_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_INCLUSIVE_xj_reco_range[irange] = new TH1D(Form("h_INCLUSIVE_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_INCLUSIVE_xj_unfold_range[irange][iter] = new TH1D(Form("h_INCLUSIVE_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  //PRIOR
  TH2D *h_prior_pt1pt2_reco = new TH2D("h_prior_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_prior_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_pt1pt2_unfold[iter] = new TH2D("h_prior_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_prior_pt1pt2_unfold[iter]->SetName(Form("h_prior_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_prior_xj_reco = new TH1D("h_prior_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_prior_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_prior_xj_unfold[iter] = new TH1D(Form("h_prior_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_prior_xj_reco_range[mbins];
  TH1D *h_prior_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_prior_xj_reco_range[irange] = new TH1D(Form("h_prior_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_prior_xj_unfold_range[irange][iter] = new TH1D(Form("h_prior_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }

  //pjes
  histo_opps::make_sym_pt1pt2(h_pjes_flat_reco_pt1pt2, h_pjes_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_pjes_flat_unfold_pt1pt2[iter], h_pjes_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pjes_pt1pt2_reco, h_pjes_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pjes_pt1pt2_unfold[iter], h_pjes_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_pjes_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_pjes_xj_unfold[iter], nbins);
    }
  //njes
  histo_opps::make_sym_pt1pt2(h_njes_flat_reco_pt1pt2, h_njes_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_njes_flat_unfold_pt1pt2[iter], h_njes_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_njes_pt1pt2_reco, h_njes_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_njes_pt1pt2_unfold[iter], h_njes_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_njes_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_njes_xj_unfold[iter], nbins);
    }
  //posjer
  histo_opps::make_sym_pt1pt2(h_posjer_flat_reco_pt1pt2, h_posjer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_posjer_flat_unfold_pt1pt2[iter], h_posjer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_posjer_pt1pt2_reco, h_posjer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_posjer_pt1pt2_unfold[iter], h_posjer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_posjer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_posjer_xj_unfold[iter], nbins);
    }
  //negjer
  histo_opps::make_sym_pt1pt2(h_negjer_flat_reco_pt1pt2, h_negjer_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_negjer_flat_unfold_pt1pt2[iter], h_negjer_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_negjer_pt1pt2_reco, h_negjer_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_negjer_pt1pt2_unfold[iter], h_negjer_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_negjer_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_negjer_xj_unfold[iter], nbins);
    }
  /* //vtx */
  /* histo_opps::make_sym_pt1pt2(h_vtx_flat_reco_pt1pt2, h_vtx_pt1pt2_reco, nbins); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::make_sym_pt1pt2(h_vtx_flat_unfold_pt1pt2[iter], h_vtx_pt1pt2_unfold[iter], nbins); */
  /*   } */

  /* histo_opps::project_xj(h_vtx_pt1pt2_reco, h_vtx_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::project_xj(h_vtx_pt1pt2_unfold[iter], h_vtx_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2); */
  /*   } */
  
  /* histo_opps::normalize_histo(h_vtx_xj_reco, nbins); */
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     histo_opps::normalize_histo(h_vtx_xj_unfold[iter], nbins); */
  /*   } */

  //ZYAM
  if (!ispp)
    {
      histo_opps::make_sym_pt1pt2(h_ZYAM_flat_reco_pt1pt2, h_ZYAM_pt1pt2_reco, nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::make_sym_pt1pt2(h_ZYAM_flat_unfold_pt1pt2[iter], h_ZYAM_pt1pt2_unfold[iter], nbins);
	  histo_opps::make_sym_pt1pt2(h_COMBUp_flat_unfold_pt1pt2[iter], h_COMBUp_pt1pt2_unfold[iter], nbins);
	}

      histo_opps::project_xj(h_ZYAM_pt1pt2_reco, h_ZYAM_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_ZYAM_pt1pt2_unfold[iter], h_ZYAM_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	}
  
      histo_opps::normalize_histo(h_ZYAM_xj_reco, nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_ZYAM_xj_unfold[iter], nbins);
	}
    }
  //INCLUSIVE
  if (!ispp)
    {
      histo_opps::make_sym_pt1pt2(h_INCLUSIVE_flat_reco_pt1pt2, h_INCLUSIVE_pt1pt2_reco, nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::make_sym_pt1pt2(h_INCLUSIVE_flat_unfold_pt1pt2[iter], h_INCLUSIVE_pt1pt2_unfold[iter], nbins);
	}

      histo_opps::project_xj(h_INCLUSIVE_pt1pt2_reco, h_INCLUSIVE_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_INCLUSIVE_pt1pt2_unfold[iter], h_INCLUSIVE_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	}
  
      histo_opps::normalize_histo(h_INCLUSIVE_xj_reco, nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_INCLUSIVE_xj_unfold[iter], nbins);
	}
    }
  //prior
  histo_opps::make_sym_pt1pt2(h_prior_flat_reco_pt1pt2, h_prior_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_prior_flat_unfold_pt1pt2[iter], h_prior_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_prior_pt1pt2_reco, h_prior_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_prior_pt1pt2_unfold[iter], h_prior_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_prior_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_prior_xj_unfold[iter], nbins);
    }
  
  TH1D *h_final_xj_data_range[mbins];
  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_data_range[irange] = new TH1D(Form("h_final_xj_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //pjes
  TH1D *h_pjes_final_xj_reco_range[mbins];
  TH1D *h_pjes_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_pjes_final_xj_reco_range[irange] = new TH1D(Form("h_pjes_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_pjes_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_pjes_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //njes
  TH1D *h_njes_final_xj_reco_range[mbins];
  TH1D *h_njes_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_njes_final_xj_reco_range[irange] = new TH1D(Form("h_njes_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_njes_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_njes_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //posjer
  TH1D *h_posjer_final_xj_reco_range[mbins];
  TH1D *h_posjer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_posjer_final_xj_reco_range[irange] = new TH1D(Form("h_posjer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_posjer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_posjer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //negjer
  TH1D *h_negjer_final_xj_reco_range[mbins];
  TH1D *h_negjer_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_negjer_final_xj_reco_range[irange] = new TH1D(Form("h_negjer_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_negjer_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_negjer_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  /* //vtx */
  /* TH1D *h_vtx_final_xj_reco_range[mbins]; */
  /* TH1D *h_vtx_final_xj_unfold_range[mbins][niterations]; */
  /* for (int irange = 0; irange < mbins; irange++) */
  /*   { */
  /*     h_vtx_final_xj_reco_range[irange] = new TH1D(Form("h_vtx_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins); */
  /*     for (int iter = 0; iter < niterations; iter++) */
  /* 	{ */
  /* 	  h_vtx_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_vtx_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins); */
  /* 	} */
  /*   } */

  //ZYAM
  TH1D *h_ZYAM_final_xj_reco_range[mbins];
  TH1D *h_ZYAM_final_xj_unfold_range[mbins][niterations];
  TH1D *h_COMBUp_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_ZYAM_final_xj_reco_range[irange] = new TH1D(Form("h_ZYAM_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_ZYAM_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_ZYAM_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_COMBUp_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_COMBUp_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  //INCLUSIVE
  TH1D *h_INCLUSIVE_final_xj_reco_range[mbins];
  TH1D *h_INCLUSIVE_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_INCLUSIVE_final_xj_reco_range[irange] = new TH1D(Form("h_INCLUSIVE_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_INCLUSIVE_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_INCLUSIVE_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  //prior
  TH1D *h_prior_final_xj_reco_range[mbins];
  TH1D *h_prior_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_prior_final_xj_reco_range[irange] = new TH1D(Form("h_prior_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_prior_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_prior_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }


  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_data, h_xj_data_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_truth_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_data_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_data_range[irange], h_final_xj_data_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

    }

  for (int irange = 0; irange < mbins; irange++)
    {
      //pjes
      histo_opps::project_xj(h_pjes_pt1pt2_reco, h_pjes_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pjes_pt1pt2_unfold[iter], h_pjes_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_pjes_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_pjes_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_pjes_xj_reco_range[irange], h_pjes_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_pjes_xj_unfold_range[irange][iter], h_pjes_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}
      //njes
      histo_opps::project_xj(h_njes_pt1pt2_reco, h_njes_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_njes_pt1pt2_unfold[iter], h_njes_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_njes_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_njes_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_njes_xj_reco_range[irange], h_njes_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_njes_xj_unfold_range[irange][iter], h_njes_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}
      //posjer
      histo_opps::project_xj(h_posjer_pt1pt2_reco, h_posjer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_posjer_pt1pt2_unfold[iter], h_posjer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_posjer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_posjer_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_posjer_xj_reco_range[irange], h_posjer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_posjer_xj_unfold_range[irange][iter], h_posjer_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}
      //negjer
      histo_opps::project_xj(h_negjer_pt1pt2_reco, h_negjer_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_negjer_pt1pt2_unfold[iter], h_negjer_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_negjer_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_negjer_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_negjer_xj_reco_range[irange], h_negjer_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_negjer_xj_unfold_range[irange][iter], h_negjer_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}
      /* //vtx */
      /* histo_opps::project_xj(h_vtx_pt1pt2_reco, h_vtx_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::project_xj(h_vtx_pt1pt2_unfold[iter], h_vtx_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2); */
      /* 	} */


      /* histo_opps::normalize_histo(h_vtx_xj_reco_range[irange], nbins); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::normalize_histo(h_vtx_xj_unfold_range[irange][iter], nbins); */
      /* 	} */

      /* histo_opps::finalize_xj(h_vtx_xj_reco_range[irange], h_vtx_final_xj_reco_range[irange], nbins, first_xj); */
      /* for (int iter = 0; iter < niterations; iter++) */
      /* 	{ */
      /* 	  histo_opps::finalize_xj(h_vtx_xj_unfold_range[irange][iter], h_vtx_final_xj_unfold_range[irange][iter], nbins, first_xj); */
      /* 	} */
      //ZYAM
      if (!ispp)
	{
	  histo_opps::project_xj(h_ZYAM_pt1pt2_reco, h_ZYAM_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::project_xj(h_ZYAM_pt1pt2_unfold[iter], h_ZYAM_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	      histo_opps::project_xj(h_COMBUp_pt1pt2_unfold[iter], h_COMBUp_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	    }
	  
	  
	  histo_opps::normalize_histo(h_ZYAM_xj_reco_range[irange], nbins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::normalize_histo(h_ZYAM_xj_unfold_range[irange][iter], nbins);
	      histo_opps::normalize_histo(h_COMBUp_xj_unfold_range[irange][iter], nbins);
	    }
	  
	  histo_opps::finalize_xj(h_ZYAM_xj_reco_range[irange], h_ZYAM_final_xj_reco_range[irange], nbins, first_xj);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::finalize_xj(h_ZYAM_xj_unfold_range[irange][iter], h_ZYAM_final_xj_unfold_range[irange][iter], nbins, first_xj);
	      histo_opps::finalize_xj(h_COMBUp_xj_unfold_range[irange][iter], h_COMBUp_final_xj_unfold_range[irange][iter], nbins, first_xj);
	    }
	}
            //INCLUSIVE
      if (!ispp)
	{
	  histo_opps::project_xj(h_INCLUSIVE_pt1pt2_reco, h_INCLUSIVE_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::project_xj(h_INCLUSIVE_pt1pt2_unfold[iter], h_INCLUSIVE_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	    }
	  
	  
	  histo_opps::normalize_histo(h_INCLUSIVE_xj_reco_range[irange], nbins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::normalize_histo(h_INCLUSIVE_xj_unfold_range[irange][iter], nbins);
	    }
	  
	  histo_opps::finalize_xj(h_INCLUSIVE_xj_reco_range[irange], h_INCLUSIVE_final_xj_reco_range[irange], nbins, first_xj);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::finalize_xj(h_INCLUSIVE_xj_unfold_range[irange][iter], h_INCLUSIVE_final_xj_unfold_range[irange][iter], nbins, first_xj);
	    }
	}
      //prior
      histo_opps::project_xj(h_prior_pt1pt2_reco, h_prior_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_prior_pt1pt2_unfold[iter], h_prior_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_prior_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_prior_xj_unfold_range[irange][iter], nbins);
	}

      histo_opps::finalize_xj(h_prior_xj_reco_range[irange], h_prior_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_prior_xj_unfold_range[irange][iter], h_prior_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

    }
    

  // JES/JER systematics
  TH1D *h_sys_posjer[mbins][niterations];
  TH1D *h_sys_negjer[mbins][niterations];
  TH1D *h_sys_pjes[mbins][niterations];
  TH1D *h_sys_njes[mbins][niterations];

  // ZYAM systematics
  TH1D *h_sys_ZYAM[mbins][niterations];
  // INCLUSIVE systematics
  TH1D *h_sys_INCLUSIVE[mbins][niterations];
  // PRIOR systematics
  TH1D *h_sys_prior[mbins][niterations];
  // VTX reweight systeamtic
  //  TH1D *h_sys_vtx[mbins][niterations];
  // total systematics
  TH1D *h_total_sys_range[mbins][niterations];
  TH1D *h_total_sys_neg_range[mbins][niterations];

  
  TCanvas *cxjjes = new TCanvas("cxjjes","cxjjes", 500, 700);
  dlutility::ratioPanelCanvas(cxjjes);

  for (int niter = 0; niter < niterations; niter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjjes->cd(1);
	  dlutility::SetLineAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_pjes_final_xj_unfold_range[irange][niter], color_pjes, msize_pjes, marker_pjes);
	  dlutility::SetLineAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, lsize_njes, 1);
	  dlutility::SetMarkerAtt(h_njes_final_xj_unfold_range[irange][niter], color_njes, msize_njes, marker_njes);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_pjes_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_njes_final_xj_unfold_range[irange][niter]->Draw("same p");

	  if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);}
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59);
	  if (!ispp) dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.55, 0.7, 0.87, 0.84);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_pjes_final_xj_unfold_range[irange][niter], "Unfolded +3% JES");
	  leg->AddEntry(h_njes_final_xj_unfold_range[irange][niter], "Unfolded -3% JES");
	  leg->Draw("same");

    
	  cxjjes->cd(2);

	  TH1D *h_pjes_compare = (TH1D*) h_pjes_final_xj_unfold_range[irange][niter]->Clone();
	  h_pjes_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  h_pjes_compare->SetTitle(";x_{J}; Unfold #pm3% JES / Unfold");
	  dlutility::SetFont(h_pjes_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_pjes_compare, color_pjes, 1,1);
	  dlutility::SetMarkerAtt(h_pjes_compare, color_pjes, 1,8);

	  h_pjes_compare->SetMaximum(1.5);
	  h_pjes_compare->SetMinimum(0.5);
	  h_pjes_compare->Draw("p");

	  TH1D *h_njes_compare = (TH1D*) h_njes_final_xj_unfold_range[irange][niter]->Clone();
	  h_njes_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  dlutility::SetLineAtt(h_njes_compare, color_njes, 1,1);
	  dlutility::SetMarkerAtt(h_njes_compare, color_njes, 1,8);
      
	  h_njes_compare->Draw("p same");

	  h_sys_njes[irange][niter] = (TH1D*) h_njes_compare->Clone();
	  h_sys_pjes[irange][niter] = (TH1D*) h_pjes_compare->Clone();
      
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjjes->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	  cxjjes->SaveAs(Form("%s/systematic_plots/h_JES_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	}

  
      TCanvas *cxjjer = new TCanvas("cxjjer","cxjjer", 500, 700);
      dlutility::ratioPanelCanvas(cxjjer);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjjer->cd(1);
	  dlutility::SetLineAtt(h_posjer_final_xj_unfold_range[irange][niter], color_posjer, lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_posjer_final_xj_unfold_range[irange][niter], color_posjer, msize_posjer, marker_posjer);

	  dlutility::SetLineAtt(h_negjer_final_xj_unfold_range[irange][niter], color_negjer, lsize_negjer, 1);
	  dlutility::SetMarkerAtt(h_negjer_final_xj_unfold_range[irange][niter], color_negjer, msize_negjer, marker_negjer);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  
	  h_posjer_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_negjer_final_xj_unfold_range[irange][niter]->Draw("same p");

	  if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);}
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59);
	  if (!ispp) dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);//	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.55, 0.64, 0.87, 0.84);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_posjer_final_xj_unfold_range[irange][niter], "Unfolded +2% JER");
	  leg->AddEntry(h_negjer_final_xj_unfold_range[irange][niter], "Unfolded -2% JER");
	  leg->Draw("same");

    
	  cxjjer->cd(2);

	  TH1D *h_posjer_compare = (TH1D*) h_posjer_final_xj_unfold_range[irange][niter]->Clone();
	  h_posjer_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  h_posjer_compare->SetTitle(";x_{J}; Unfold + Sys. JER / Unfold");
	  dlutility::SetFont(h_posjer_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_posjer_compare, color_posjer, 1,1);
	  dlutility::SetMarkerAtt(h_posjer_compare, color_posjer, 1,8);
	  TH1D *h_negjer_compare = (TH1D*) h_negjer_final_xj_unfold_range[irange][niter]->Clone();
	  h_negjer_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  h_negjer_compare->SetTitle(";x_{J}; Unfold - Sys. JER / Unfold");
	  dlutility::SetFont(h_negjer_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_negjer_compare, color_negjer, 1,1);
	  dlutility::SetMarkerAtt(h_negjer_compare, color_negjer, 1,8);

	  h_posjer_compare->SetMaximum(1.5);
	  h_posjer_compare->SetMinimum(0.5);
	  h_posjer_compare->Draw("p");
	  h_negjer_compare->Draw("p same");
	  h_sys_posjer[irange][niter] = (TH1D*) h_posjer_compare->Clone();
	  h_sys_negjer[irange][niter] = (TH1D*) h_negjer_compare->Clone();
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjjer->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	  cxjjer->SaveAs(Form("%s/systematic_plots/h_JER_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));

	}
      /* TCanvas *cxjvtx = new TCanvas("cxjvtx","cxjvtx", 500, 700); */
      /* dlutility::ratioPanelCanvas(cxjvtx); */

      /* for (int irange = 0; irange < mbins; irange++) */
      /* 	{ */
      /* 	  cxjvtx->cd(1); */
      /* 	  dlutility::SetLineAtt(h_vtx_final_xj_unfold_range[irange][niter], color_vtx, lsize_vtx, 1); */
      /* 	  dlutility::SetMarkerAtt(h_vtx_final_xj_unfold_range[irange][niter], color_vtx, msize_vtx, marker_vtx); */

      /* 	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1); */
      /* 	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold); */

      /* 	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1); */
      /* 	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia); */

      /* 	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04); */
      /* 	  h_final_xj_truth_range[irange]->SetMaximum(4.5); */
      /* 	  h_final_xj_truth_range[irange]->SetMinimum(0); */
      /* 	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}"); */

      /* 	  h_final_xj_truth_range[irange]->Draw("p"); */
      /* 	  h_final_xj_unfold_range[irange][niter]->Draw("same p"); */
      /* 	  h_vtx_final_xj_unfold_range[irange][niter]->Draw("same p"); */

      /* 	  if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);} */
      /* 	  dlutility::drawText("anti-k_{t} R = 0.4", 0.22, 0.74); */
      /* 	  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69); */
      /* 	  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64); */
      /* 	  dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59); */
      /* 	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54); */
      /* 	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49); */
      /* 	  TLegend *leg = new TLegend(0.55, 0.7, 0.87, 0.84); */
      /* 	  leg->SetLineWidth(0); */
      /* 	  leg->SetTextSize(0.04); */
      /* 	  leg->SetTextFont(42); */
      /* 	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight"); */
      /* 	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded"); */
      /* 	  leg->AddEntry(h_vtx_final_xj_unfold_range[irange][niter], "Unfolded + z_{vtx} reweight"); */
      /* 	  leg->Draw("same"); */

    
      /* 	  cxjvtx->cd(2); */

      /* 	  TH1D *h_vtx_compare = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone(); */
      /* 	  h_vtx_compare->Divide(h_vtx_final_xj_unfold_range[irange][niter]); */
      /* 	  h_vtx_compare->SetTitle(";x_{J}; Unfold / Unfold + z_{vtx} reweight"); */
      /* 	  dlutility::SetFont(h_vtx_compare, 42, 0.06); */
      /* 	  dlutility::SetLineAtt(h_vtx_compare, color_vtx, 1,1); */
      /* 	  dlutility::SetMarkerAtt(h_vtx_compare, color_vtx, 1,8); */

      /* 	  h_vtx_compare->SetMaximum(2.0); */
      /* 	  h_vtx_compare->SetMinimum(0.0); */
      /* 	  h_vtx_compare->Draw("p"); */
      /* 	  h_sys_vtx[irange][niter] = (TH1D*) h_vtx_compare->Clone(); */
      /* 	  TLine *line = new TLine(0.1, 1, 1, 1); */
      /* 	  line->SetLineStyle(4); */
      /* 	  line->SetLineColor(kRed + 3); */
      /* 	  line->SetLineWidth(2); */
      /* 	  line->Draw("same"); */
      /* 	  cxjvtx->SaveAs(Form("systematic_plots/h_VTX_xj_unfolded_range_%d_iter_%d.png", irange, niter)); */
      /* 	  cxjvtx->SaveAs(Form("systematic_plots/h_VTX_xj_unfolded_range_%d_iter_%d.pdf", irange, niter)); */
      /* 	} */

      if (!ispp)
	{
	  TCanvas *cxjZYAM = new TCanvas("cxjZYAM","cxjZYAM", 500, 700);
	  dlutility::ratioPanelCanvas(cxjZYAM);

	  for (int irange = 0; irange < mbins; irange++)
	    {
	      cxjZYAM->cd(1);
	      dlutility::SetLineAtt(h_ZYAM_final_xj_unfold_range[irange][niter], color_ZYAM, lsize_ZYAM, 1);
	      dlutility::SetMarkerAtt(h_ZYAM_final_xj_unfold_range[irange][niter], color_ZYAM, msize_ZYAM, marker_ZYAM);
	      dlutility::SetLineAtt(h_COMBUp_final_xj_unfold_range[irange][niter], kRed + 1, lsize_ZYAM, 2);
	      dlutility::SetMarkerAtt(h_COMBUp_final_xj_unfold_range[irange][niter], kRed + 1, msize_ZYAM, 24);

	      dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	      dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	      dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	      dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	      dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	      h_final_xj_truth_range[irange]->SetMaximum(4.5);
	      h_final_xj_truth_range[irange]->SetMinimum(0);
	      h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	      h_final_xj_truth_range[irange]->Draw("p");
	      h_final_xj_unfold_range[irange][niter]->Draw("same p");
	      h_ZYAM_final_xj_unfold_range[irange][niter]->Draw("same p");
	      h_COMBUp_final_xj_unfold_range[irange][niter]->Draw("same p");

	      if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);}
	      dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	      dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59);
	      if (!ispp) dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);//dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	      TLegend *leg = new TLegend(0.55, 0.7, 0.87, 0.84);
	      leg->SetLineWidth(0);
	      leg->SetTextSize(0.04);
	      leg->SetTextFont(42);
	      leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight");
	      leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded w/ Flow Bkg.");
	      leg->AddEntry(h_ZYAM_final_xj_unfold_range[irange][niter], "COMBDown: 0.7(v_{2,2},v_{3,3})");
	      leg->AddEntry(h_COMBUp_final_xj_unfold_range[irange][niter], "COMBUp: 1.3(v_{2,2},v_{3,3})");
	      leg->Draw("same");

    
	      cxjZYAM->cd(2);

	      TH1D *h_ZYAM_compare = (TH1D*) h_ZYAM_final_xj_unfold_range[irange][niter]->Clone();
	      h_ZYAM_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	      TH1D *h_COMBUp_compare = (TH1D*) h_COMBUp_final_xj_unfold_range[irange][niter]->Clone();
	      h_COMBUp_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	      h_ZYAM_compare->SetTitle(";x_{J};COMB variation / nominal");
	      dlutility::SetFont(h_ZYAM_compare, 42, 0.06);
	      dlutility::SetLineAtt(h_ZYAM_compare, color_ZYAM, 1,1);
	      dlutility::SetMarkerAtt(h_ZYAM_compare, color_ZYAM, 1,8);

	      h_ZYAM_compare->SetMaximum(1.5);
	      h_ZYAM_compare->SetMinimum(0.5);
	      h_ZYAM_compare->Draw("p");
	      dlutility::SetLineAtt(h_COMBUp_compare, kRed + 1, 1,2);
	      dlutility::SetMarkerAtt(h_COMBUp_compare, kRed + 1, 1,24);
	      h_COMBUp_compare->Draw("p same");

	      h_sys_ZYAM[irange][niter] = (TH1D*) h_ZYAM_compare->Clone();
	      h_sys_ZYAM[irange][niter]->SetName(Form("h_COMB_envelope_range_%d_iter_%d", irange, niter));
	      for (int ibin = 1; ibin <= h_sys_ZYAM[irange][niter]->GetNbinsX(); ++ibin)
	        if (std::fabs(h_COMBUp_compare->GetBinContent(ibin) - 1.0) >
	            std::fabs(h_ZYAM_compare->GetBinContent(ibin) - 1.0))
	          h_sys_ZYAM[irange][niter]->SetBinContent(
	            ibin, h_COMBUp_compare->GetBinContent(ibin));
	      TLine *line = new TLine(0.1, 1, 1, 1);
	      line->SetLineStyle(4);
	      line->SetLineColor(kRed + 3);
	      line->SetLineWidth(2);
	      line->Draw("same");
	      cxjZYAM->SaveAs(Form("%s/systematic_plots/h_COMB_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	      cxjZYAM->SaveAs(Form("%s/systematic_plots/h_COMB_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	    }
	  TCanvas *cxjINCLUSIVE = new TCanvas("cxjINCLUSIVE","cxjINCLUSIVE", 500, 700);
	  dlutility::ratioPanelCanvas(cxjINCLUSIVE);

	  for (int irange = 0; irange < mbins; irange++)
	    {
	      cxjINCLUSIVE->cd(1);
	      dlutility::SetLineAtt(h_INCLUSIVE_final_xj_unfold_range[irange][niter], color_INCLUSIVE, lsize_INCLUSIVE, 1);
	      dlutility::SetMarkerAtt(h_INCLUSIVE_final_xj_unfold_range[irange][niter], color_INCLUSIVE, msize_INCLUSIVE, marker_INCLUSIVE);

	      dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	      dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	      dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	      dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	      dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	      h_final_xj_truth_range[irange]->SetMaximum(4.5);
	      h_final_xj_truth_range[irange]->SetMinimum(0);
	      h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	      h_final_xj_truth_range[irange]->Draw("p");
	      h_final_xj_unfold_range[irange][niter]->Draw("same p");
	      h_INCLUSIVE_final_xj_unfold_range[irange][niter]->Draw("same p");

	      if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);}
	      dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	      dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59);
	      if (!ispp) dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);//dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	      dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	      TLegend *leg = new TLegend(0.55, 0.7, 0.87, 0.84);
	      leg->SetLineWidth(0);
	      leg->SetTextSize(0.04);
	      leg->SetTextFont(42);
	      leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight");
	      leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded w/ Exclusive");
	      leg->AddEntry(h_INCLUSIVE_final_xj_unfold_range[irange][niter], "Unfolded w/ Inclusive ");
	      leg->Draw("same");

    
	      cxjINCLUSIVE->cd(2);

	      TH1D *h_INCLUSIVE_compare = (TH1D*) h_INCLUSIVE_final_xj_unfold_range[irange][niter]->Clone();
	      h_INCLUSIVE_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	      h_INCLUSIVE_compare->SetTitle(";x_{J}; Unfold w/ Inclusive / Unfold");
	      dlutility::SetFont(h_INCLUSIVE_compare, 42, 0.06);
	      dlutility::SetLineAtt(h_INCLUSIVE_compare, color_INCLUSIVE, 1,1);
	      dlutility::SetMarkerAtt(h_INCLUSIVE_compare, color_INCLUSIVE, 1,8);

	      h_INCLUSIVE_compare->SetMaximum(1.5);
	      h_INCLUSIVE_compare->SetMinimum(0.5);
	      h_INCLUSIVE_compare->Draw("p");

	      h_sys_INCLUSIVE[irange][niter] = (TH1D*) h_INCLUSIVE_compare->Clone();
	      TLine *line = new TLine(0.1, 1, 1, 1);
	      line->SetLineStyle(4);
	      line->SetLineColor(kRed + 3);
	      line->SetLineWidth(2);
	      line->Draw("same");
	      cxjINCLUSIVE->SaveAs(Form("%s/systematic_plots/h_INCLUSIVE_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	      cxjINCLUSIVE->SaveAs(Form("%s/systematic_plots/h_INCLUSIVE_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	    }

	}

      TCanvas *cxjprior = new TCanvas("cxjprior","cxjprior", 500, 700);
      dlutility::ratioPanelCanvas(cxjprior);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjprior->cd(1);
	  dlutility::SetLineAtt(h_prior_final_xj_unfold_range[irange][niter], color_prior, lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_prior_final_xj_unfold_range[irange][niter], color_prior, msize_prior, marker_prior);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.04);
	  h_final_xj_truth_range[irange]->SetMaximum(4.5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");

	  h_final_xj_truth_range[irange]->Draw("p");
	  h_final_xj_unfold_range[irange][niter]->Draw("same p");
	  h_prior_final_xj_unfold_range[irange][niter]->Draw("same p");

	  if (!ispp) {dlutility::DrawSPHENIX(0.22, 0.84);} else {dlutility::DrawSPHENIXpp(0.22, 0.84);}
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.22, 0.59);
	  if (!ispp) 
	    dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);//	  dlutility::drawText("L_{int} = 25.7 pb^{-1}", 0.22, 0.54);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
	  TLegend *leg = new TLegend(0.55, 0.7, 0.87, 0.84);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Pythia8 - Reweight");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Unfolded");
	  leg->AddEntry(h_prior_final_xj_unfold_range[irange][niter], "Unfolded w/o Prior");
	  leg->Draw("same");
    
	  cxjprior->cd(2);

	  TH1D *h_prior_compare = (TH1D*) h_prior_final_xj_unfold_range[irange][niter]->Clone();
	  h_prior_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  h_prior_compare->SetTitle(";x_{J}; Unfold w/o Prior / Unfold");
	  dlutility::SetFont(h_prior_compare, 42, 0.06);
	  dlutility::SetLineAtt(h_prior_compare, color_prior, 1,1);
	  dlutility::SetMarkerAtt(h_prior_compare, color_prior, 1,8);

	  h_prior_compare->SetMaximum(1.5);
	  h_prior_compare->SetMinimum(0.5);
	  h_prior_compare->Draw("p");

	  h_sys_prior[irange][niter] = (TH1D*) h_prior_compare->Clone();
	  TLine *line = new TLine(0.1, 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxjprior->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	  cxjprior->SaveAs(Form("%s/systematic_plots/h_PRIOR_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	}

      // all systematics on one plot
      TCanvas *cxjsys = new TCanvas("cxjsys","cxjsys", 750, 500);

      for (int irange = 0; irange < mbins; irange++)
	{
	  cxjsys->Clear();
	  dlutility::createCutCanvas(cxjsys, 0.65);
	  cxjsys->cd(1);
	  
	  TH1D *h_sys_jer = (TH1D*) h_sys_posjer[irange][niter]->Clone();
	  TH1D *h_sys_jer_flip = (TH1D*) h_sys_posjer[irange][niter]->Clone();

	  //TH1D *h_sys_vtx_flip = (TH1D*) h_sys_vtx[irange][niter]->Clone();
	  TH1D *h_sys_ZYAM_flip;
	  if (!ispp)
	    {
	      h_sys_ZYAM_flip = (TH1D*) h_sys_ZYAM[irange][niter]->Clone();
	    }
	  TH1D *h_sys_INCLUSIVE_flip;
	  if (!ispp)
	    {
	      h_sys_INCLUSIVE_flip = (TH1D*) h_sys_INCLUSIVE[irange][niter]->Clone();
	    }

	  TH1D *h_sys_prior_flip = (TH1D*) h_sys_prior[irange][niter]->Clone();

	  TH1D *h_sys_jes = (TH1D*) h_sys_njes[irange][niter]->Clone();
	  TH1D *h_sys_jes_flip = (TH1D*) h_sys_njes[irange][niter]->Clone();	  

	  TH1D *h_total_sys = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys->Reset();
	  h_total_sys->SetName("h_total_sys");
	  h_total_sys->SetTitle(";x_{J}; Variation / Nominal");
	  TH1D *h_total_sys_flip = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_flip->Reset();
	  h_total_sys_flip->SetName("h_total_sys_flip");
	  h_total_sys_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_range[irange][niter]->Reset();
	  h_total_sys_range[irange][niter]->SetName(Form("h_total_sys_range_%d_iter_%d", irange, niter));
	  h_total_sys_neg_range[irange][niter] = (TH1D*) h_sys_pjes[irange][niter]->Clone();
	  h_total_sys_neg_range[irange][niter]->Reset();
	  h_total_sys_neg_range[irange][niter]->SetName(Form("h_total_sys_neg_range_%d_iter_%d", irange, niter));

      
	  for (int i = 0; i < h_sys_prior_flip->GetNbinsX(); i++)
	    {
	  
	      if (!ispp) h_sys_ZYAM_flip->SetBinContent(i+1, 2 - h_sys_ZYAM[irange][niter]->GetBinContent(i+1));
	      if (!ispp) h_sys_INCLUSIVE_flip->SetBinContent(i+1, 2 - h_sys_INCLUSIVE[irange][niter]->GetBinContent(i+1));
	      h_sys_prior_flip->SetBinContent(i+1, 2 - h_sys_prior[irange][niter]->GetBinContent(i+1));
	      float maxjer = h_sys_posjer[irange][niter]->GetBinContent(i+1) - 1;
	      float minjer = h_sys_negjer[irange][niter]->GetBinContent(i+1) - 1;

	      if (maxjer < 0 && minjer > 0)
		{
		  maxjer = h_sys_negjer[irange][niter]->GetBinContent(i+1) - 1;
		  minjer = h_sys_posjer[irange][niter]->GetBinContent(i+1) - 1;
		}
	      else if(maxjer > 0 && minjer > 0)
		{
		  maxjer = std::max(maxjer, minjer);
		  minjer = 0;
		}
	      else if(maxjer < 0 && minjer < 0)
		{
		  minjer = std::min(maxjer, minjer);
		  maxjer = 0;
		}
	      
	      float maxjes = h_sys_pjes[irange][niter]->GetBinContent(i+1) - 1;
	      float minjes = h_sys_njes[irange][niter]->GetBinContent(i+1) - 1;
	      if (maxjes < 0 && minjes > 0)
		{
		  maxjes = h_sys_njes[irange][niter]->GetBinContent(i+1) - 1;
		  minjes = h_sys_pjes[irange][niter]->GetBinContent(i+1) - 1;
		}
	      else if(maxjes > 0 && minjes > 0)
		{
		  maxjes = std::max(maxjes, minjes);
		  minjes = 0;
		}
	      else if(maxjes < 0 && minjes < 0)
		{
		  minjes = std::min(maxjes, minjes);
		  maxjes = 0;
		}

	      std::cout << "iter: " << niter << " / range "<< irange << " / JES " << maxjes << " /  JER " << maxjer << std::endl; 

	      h_sys_jer->SetBinContent(i+1, maxjer + 1);
	      h_sys_jer_flip->SetBinContent(i+1, 1-fabs(minjer));

	      h_sys_jes->SetBinContent(i+1, maxjes + 1);
	      h_sys_jes_flip->SetBinContent(i+1, 1 - fabs(minjes));

	      float toal_neg = TMath::Power(minjer, 2) +
		TMath::Power(minjes, 2) +
		TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1) - 1, 2);
	      if (!ispp)
		{
		  toal_neg += TMath::Power(h_sys_ZYAM[irange][niter]->GetBinContent(i+1) - 1, 2);
		  toal_neg += TMath::Power(h_sys_INCLUSIVE[irange][niter]->GetBinContent(i+1) - 1, 2);
		}
	      
	      toal_neg = sqrt(toal_neg);
	      
	      float toal_pos = TMath::Power(maxjer, 2) +
		TMath::Power(maxjes, 2) +		
		TMath::Power(h_sys_prior[irange][niter]->GetBinContent(i+1) - 1, 2);
	      if (!ispp)
		{
		  toal_pos += TMath::Power(h_sys_ZYAM[irange][niter]->GetBinContent(i+1) - 1, 2);
		  toal_pos += TMath::Power(h_sys_INCLUSIVE[irange][niter]->GetBinContent(i+1) - 1, 2);
		}
	      toal_pos = sqrt(toal_pos);
	      h_total_sys_range[irange][niter]->SetBinContent(i+1, toal_pos);
	      h_total_sys_neg_range[irange][niter]->SetBinContent(i+1, toal_neg);
	      h_total_sys->SetBinContent(i+1, 1 + toal_pos);
	      h_total_sys_flip->SetBinContent(i+1, 1 - toal_neg);

	    }

	  dlutility::SetLineAtt(h_sys_jer, color_posjer, 1 + lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer, color_posjer, msize_posjer, 1);

	  dlutility::SetLineAtt(h_sys_jer_flip, color_posjer, 1 + lsize_posjer, 1);
	  dlutility::SetMarkerAtt(h_sys_jer_flip, color_posjer, msize_posjer, 1);

	  dlutility::SetLineAtt(h_sys_jes, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes, color_pjes, msize_pjes, 1);
	  dlutility::SetLineAtt(h_sys_jes_flip, color_pjes, 1 + lsize_pjes, 1);
	  dlutility::SetMarkerAtt(h_sys_jes_flip, color_pjes, msize_pjes, 1);

	  if (!ispp)
	    {
	      dlutility::SetLineAtt(h_sys_ZYAM[irange][niter], color_ZYAM, 1 + lsize_ZYAM, 1);
	      dlutility::SetMarkerAtt(h_sys_ZYAM[irange][niter], color_ZYAM, msize_ZYAM, 1);
	      dlutility::SetLineAtt(h_sys_ZYAM_flip, color_ZYAM, 1 + lsize_ZYAM, 1);
	      dlutility::SetMarkerAtt(h_sys_ZYAM_flip, color_ZYAM, msize_ZYAM, 1);
	      dlutility::SetLineAtt(h_sys_INCLUSIVE[irange][niter], color_INCLUSIVE, 1 + lsize_INCLUSIVE, 1);
	      dlutility::SetMarkerAtt(h_sys_INCLUSIVE[irange][niter], color_INCLUSIVE, msize_INCLUSIVE, 1);
	      dlutility::SetLineAtt(h_sys_INCLUSIVE_flip, color_INCLUSIVE, 1 + lsize_INCLUSIVE, 1);
	      dlutility::SetMarkerAtt(h_sys_INCLUSIVE_flip, color_INCLUSIVE, msize_INCLUSIVE, 1);
	    }
	  dlutility::SetLineAtt(h_sys_prior[irange][niter], color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior[irange][niter], color_prior, msize_prior, 1);
	  dlutility::SetLineAtt(h_sys_prior_flip, color_prior, 1 + lsize_prior, 1);
	  dlutility::SetMarkerAtt(h_sys_prior_flip, color_prior, msize_prior, 1);
	  
	  dlutility::SetLineAtt(h_total_sys, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys, kBlack, 1, 8);
	  dlutility::SetLineAtt(h_total_sys_flip, kBlack, 1, 1);
	  dlutility::SetMarkerAtt(h_total_sys_flip, kBlack, 1, 8);


	  h_total_sys->SetMinimum(0.0);
	  h_total_sys->SetMaximum(2.0);
	  dlutility::SetFont(h_total_sys, 42, 0.06, 0.05, 0.05, 0.05);
	  h_total_sys->Rebin(nbins - first_bin, Form("h_rebin_sys_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_total_sys_flip->Rebin(nbins - first_bin, Form("h_rebin_sys_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jer->Rebin(nbins - first_bin, Form("h_rebin_jer_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jer_flip->Rebin(nbins - first_bin, Form("h_rebin_jer_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jes->Rebin(nbins - first_bin, Form("h_rebin_jes_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_jes_flip->Rebin(nbins - first_bin, Form("h_rebin_jes_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);

	  h_sys_prior[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_prior_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  h_sys_prior_flip->Rebin(nbins - first_bin, Form("h_rebin_prior_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	  if (!ispp)
	    {
	      h_sys_ZYAM[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_ZYAM_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	      h_sys_ZYAM_flip->Rebin(nbins - first_bin, Form("h_rebin_ZYAM_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	      h_sys_INCLUSIVE[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_INCLUSIVE_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);
	      h_sys_INCLUSIVE_flip->Rebin(nbins - first_bin, Form("h_rebin_INCLUSIVE_flip_%d_%d", irange, niter), &dxj_bins[first_bin]);

	    }
	  h_total_sys->GetXaxis()->SetRangeUser(dxj_bins[first_bin], dxj_bins[nbins]);
	  h_total_sys->Draw("p");
	  h_total_sys_flip->Draw("p same");

	  /* h_sys_ZYAM[irange][niter]->Draw("hist same"); */
	  /* h_sys_ZYAM_flip->Draw("hist same"); */

	  h_sys_jer->Draw("hist same");
	  h_sys_jer_flip->Draw("hist same");
	  h_sys_jes->Draw("hist same");
	  h_sys_jes_flip->Draw("hist same");
	  if (!ispp)
	    {
	      h_sys_ZYAM[irange][niter]->Draw("hist same");
	      h_sys_ZYAM_flip->Draw("hist same");
	      h_sys_INCLUSIVE[irange][niter]->Draw("hist same");
	      h_sys_INCLUSIVE_flip->Draw("hist same");

	    }
	  h_sys_prior[irange][niter]->Draw("hist same");
	  h_sys_prior_flip->Draw("hist same");



	  cxjsys->cd(2);
	  if (!ispp) dlutility::DrawSPHENIXcut(0.02, 0.9, 1, 0.09);
	  else dlutility::DrawSPHENIXpp(0.02, 0.9, 0.09);
	  dlutility::drawText(Form("anti-#it{k}_{t} #kern[-0.3]{#it{R}} = %0.1f", 0.1*cone_size), 0.02, 0.78, 0, kBlack, 0.09);
	  dlutility::drawText(Form("%2.1f #leq #kern[-0.15]{#it{p}_{T,1}} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.02, 0.72, 0, kBlack, 0.09);
	  dlutility::drawText(Form("#it{p}_{T,2} #kern[-0.15]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.02,0.66, 0, kBlack, 0.09);
	  dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 7#pi/8}", 0.02,0.60, 0, kBlack, 0.09);
	  if (!ispp) dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.02, 0.54, 0, kBlack, 0.09);
	  TLegend *leg4 = new TLegend(0.02, 0.1, 0.9, 0.53);
	  leg4->SetTextSize(0.09);
	  leg4->SetTextFont(42);
	  leg4->SetLineWidth(0);
	  leg4->AddEntry(h_total_sys, "Total Systematics","p");
	  leg4->AddEntry(h_sys_jes_flip, "JES Systematics");
	  leg4->AddEntry(h_sys_jer_flip, "JER Systematics");
	  if (!ispp) leg4->AddEntry(h_sys_ZYAM_flip, "COMB flow");
	  if (!ispp) leg4->AddEntry(h_sys_INCLUSIVE_flip, "Inclusive");
	  leg4->AddEntry(h_sys_prior_flip, "Prior");
	  leg4->Draw();
	  cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_%s_r%02d_range_%d_iter_%d.png",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));
	  cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_%s_r%02d_range_%d_iter_%d.pdf",  rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange, niter));

	}
    }
  TFile *fsys = new TFile(Form("%s/uncertainties/systematics_%s_r%02d.root",  rb.get_code_location().c_str(), system_string.c_str(), cone_size), "recreate");
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  h_total_sys_range[irange][iter]->Write();
	  h_total_sys_neg_range[irange][iter]->Write();
	  h_sys_pjes[irange][iter]->Write(Form("h_sys_posJES_range_%d_iter_%d", irange, iter));
	  h_sys_njes[irange][iter]->Write(Form("h_sys_negJES_range_%d_iter_%d", irange, iter));
	  h_sys_posjer[irange][iter]->Write(Form("h_sys_posJER_range_%d_iter_%d", irange, iter));
	  h_sys_negjer[irange][iter]->Write(Form("h_sys_negJER_range_%d_iter_%d", irange, iter));
	  h_sys_prior[irange][iter]->Write(Form("h_sys_PRIOR_range_%d_iter_%d", irange, iter));
	  if (!ispp)
	    {
	      h_sys_ZYAM[irange][iter]->Write(Form("h_sys_COMB_range_%d_iter_%d", irange, iter));
	      h_sys_INCLUSIVE[irange][iter]->Write(Form("h_sys_INCLUSIVE_range_%d_iter_%d", irange, iter));
	    }
	}
    }
  fsys->Close();
  return;
}
