#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

void drawPrior_AA(const int cone_size = 3, const int centrality_bin = 0)
{
  int color_nom = kBlack;
  int color_p1 = kBlue;
  int color_p2 = kRed;
  read_binning rb("binning_AA.config");

  Int_t read_nbins = rb.get_nbins();
  std::string dphi_string = rb.get_dphi_string();
  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins+1];
  rb.get_centrality_bins(icentrality_bins);

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);
  fgaus->SetParameters(1, 0, 0.1);

  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  if (JER_sys > 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
      fgaus->SetParameters(1, 0, JER_sys);
    }
  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
    }
  const int nbins = read_nbins;

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
  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int measure_bins[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;


  TString responsepath = Form("response_matrices/response_matrix_AA_cent_%d_r%02d_nominal.root", centrality_bin, cone_size);
  TFile *fr = new TFile(responsepath.Data(),"r");
  TH1D *h_nom_flat_reco_pt1pt2 = (TH1D*) fr->Get("h_reco_flat_pt1pt2");
  h_nom_flat_reco_pt1pt2->SetName("h_nom_flat_reco_pt1pt2");
  TH1D *h_nom_flat_truth_pt1pt2 = (TH1D*) fr->Get("h_truth_flat_pt1pt2");
  h_nom_flat_truth_pt1pt2->SetName("h_nom_flat_truth_pt1pt2");
  
  TString responsepath_p1 = Form("response_matrices/response_matrix_AA_cent_%d_r%02d_PRIMER1.root", centrality_bin, cone_size);
  TFile *fr_p1 = new TFile(responsepath_p1.Data(),"r");
  TH1D *h_p1_flat_reco_pt1pt2= (TH1D*) fr_p1->Get("h_reco_flat_pt1pt2");
  h_p1_flat_reco_pt1pt2->SetName("h_p1_flat_reco_pt1pt2");
  TH1D *h_p1_flat_truth_pt1pt2 = (TH1D*) fr_p1->Get("h_truth_flat_pt1pt2");
  h_p1_flat_truth_pt1pt2->SetName("h_p1_flat_truth_pt1pt2");

  TString responsepath_p2 = Form("response_matrices/response_matrix_AA_cent_%d_r%02d_PRIMER2.root", centrality_bin, cone_size);
  TFile *fr_p2 = new TFile(responsepath_p2.Data(),"r");
  TH1D *h_p2_flat_reco_pt1pt2= (TH1D*) fr_p2->Get("h_reco_flat_pt1pt2");
  h_p2_flat_reco_pt1pt2->SetName("h_p2_reco_flat_pt1pt2");
  TH1D *h_p2_flat_truth_pt1pt2 = (TH1D*) fr_p2->Get("h_truth_flat_pt1pt2");
  h_p2_flat_truth_pt1pt2->SetName("h_p2_flat_truth_pt1pt2");

  // make TH2s and TH1s for projection prep:

  TH2D *h_nom_pt1pt2_reco = new TH2D("h_nom_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_nom_pt1pt2_truth = new TH2D("h_nom_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_nom_xj_reco = new TH1D("h_nom_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_nom_xj_truth = new TH1D("h_nom_xj_truth", ";x_{J};",nbins, ixj_bins);

  TH2D *h_p1_pt1pt2_reco = new TH2D("h_p1_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_p1_pt1pt2_truth = new TH2D("h_p1_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_p1_xj_reco = new TH1D("h_p1_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_p1_xj_truth = new TH1D("h_p1_xj_truth", ";x_{J};",nbins, ixj_bins);

  TH2D *h_p2_pt1pt2_reco = new TH2D("h_p2_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_p2_pt1pt2_truth = new TH2D("h_p2_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_p2_xj_reco = new TH1D("h_p2_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_p2_xj_truth = new TH1D("h_p2_xj_truth", ";x_{J};",nbins, ixj_bins);
    
  histo_opps::make_sym_pt1pt2(h_nom_flat_truth_pt1pt2, h_nom_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_nom_flat_reco_pt1pt2, h_nom_pt1pt2_reco, nbins);

  histo_opps::make_sym_pt1pt2(h_p1_flat_truth_pt1pt2, h_p1_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_p1_flat_reco_pt1pt2, h_p1_pt1pt2_reco, nbins);

  histo_opps::make_sym_pt1pt2(h_p2_flat_truth_pt1pt2, h_p2_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_p2_flat_reco_pt1pt2, h_p2_pt1pt2_reco, nbins);

  histo_opps::project_xj(h_nom_pt1pt2_reco, h_nom_xj_reco, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_nom_pt1pt2_truth, h_nom_xj_truth, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);

  histo_opps::project_xj(h_p1_pt1pt2_reco, h_p1_xj_reco, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_p1_pt1pt2_truth, h_p1_xj_truth, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);

  histo_opps::project_xj(h_p2_pt1pt2_reco, h_p2_xj_reco, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_p2_pt1pt2_truth, h_p2_xj_truth, nbins, measure_bins[1], measure_bins[2], measure_subleading_bin, nbins - 2);


  histo_opps::normalize_histo(h_nom_xj_reco, nbins);
  histo_opps::normalize_histo(h_nom_xj_truth, nbins);

  histo_opps::normalize_histo(h_p1_xj_reco, nbins);
  histo_opps::normalize_histo(h_p1_xj_truth, nbins);

  histo_opps::normalize_histo(h_p2_xj_reco, nbins);
  histo_opps::normalize_histo(h_p2_xj_truth, nbins);


  dlutility::SetMarkerAtt(h_nom_xj_reco, color_nom, 1.2, 20);
  dlutility::SetLineAtt(h_nom_xj_reco, color_nom, 1, 1);
  dlutility::SetMarkerAtt(h_nom_xj_truth, color_nom, 1.2, 21);
  dlutility::SetLineAtt(h_nom_xj_truth, color_nom, 1, 1);

  dlutility::SetMarkerAtt(h_p1_xj_reco, color_p1, 1.2, 20);
  dlutility::SetLineAtt(h_p1_xj_reco, color_p1, 1, 1);
  dlutility::SetMarkerAtt(h_p1_xj_truth, color_p1, 1.2, 20);
  dlutility::SetLineAtt(h_p1_xj_truth, color_p1, 1, 1);

  dlutility::SetMarkerAtt(h_p2_xj_reco, color_p2, 1.2, 20);
  dlutility::SetLineAtt(h_p2_xj_reco, color_p2, 1, 1);
  dlutility::SetMarkerAtt(h_p2_xj_truth, color_p2, 1.2, 21);
  dlutility::SetLineAtt(h_p2_xj_truth, color_p2, 1, 1);

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  h_nom_xj_truth->SetMaximum(4.5);
  dlutility::SetFont(h_nom_xj_truth, 42, 0.05);
  //h_nom_xj_reco->Draw("p");
  h_nom_xj_truth->Draw("p");
  //h_p2_xj_reco->Draw("p same");
  h_p1_xj_truth->Draw("p same");

  int irange = 1;
  dlutility::DrawSPHENIX(0.22, 0.84);
  dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, 0.74);
  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText(Form("#Delta#phi #geq %s", dphi_string.c_str()), 0.22, 0.59);
  dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);

  TLegend *leg = new TLegend(0.65, 0.78, 0.85, 0.86);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.03);
  //leg->AddEntry(h_nom_xj_reco,"Nominal Reco");
  leg->AddEntry(h_p1_xj_truth,"Truth");
  leg->AddEntry(h_nom_xj_truth," Truth w/ Reweight");
  //leg->AddEntry(h_p2_xj_reco,"Prior Reco");

  leg->Draw("same");
  c->Print(Form("%s/unfolding_plots/prior_compare_AA_cent_%d_r%02d.pdf", rb.get_code_location().c_str(), centrality_bin, cone_size));
  return;
}
