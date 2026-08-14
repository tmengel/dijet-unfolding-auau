#include <iostream>
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

const float etacut = 1.1;
const float eta_cut_dijet = 0.8;

const float vertex_cut = 60;

void getProbability(const int cone_size = 3,  const std::string configfile = "binning_AA.config")
{
  std::cout << "Staring" << std::endl;
  read_binning rb(configfile.c_str());
  std::cout << "Read" << std::endl;
  Int_t read_nbins = rb.get_nbins();


  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();

  float jet_area =TMath::Pi()*TMath::Power(cone_size*0.1, 2);
  float dphi_area =2*(TMath::Pi() - dphicut)*1.6;
  float total_area = (TMath::Pi()*2*1.6);;
  float scale_factor = (total_area - dphi_area - jet_area) / total_area;
  float scale_factor_above = (total_area - jet_area) / total_area;

  
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


  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  // get centrality
  const int centrality_bins = rb.get_number_centrality_bins();

  float icentrality_bins[centrality_bins + 1];
  rb.get_centrality_bins(icentrality_bins);


  TFile *fun[centrality_bins];
  TH1D *h_centrality_data[centrality_bins];
  for (int i = 0;  i < centrality_bins; i++)
    {      
      fun[i] = new TFile(Form("%s/unfolding_hists/unfolding_hists_preload_AA_cent_%d_r%02d_nominal.root", rb.get_code_location().c_str(), i , cone_size), "r");
      if (!fun[i]) return;
      h_centrality_data[i] = (TH1D*) fun[i]->Get("h_centrality");
      if (!h_centrality_data[i]) return;
      h_centrality_data[i]->SetName(Form("h_centrality_data_%d", i));
    }
  
  std::string infile = "../../trees/TREE_DIJET_v10_2_502_2024p022_v001_gl10-all.root";

  TFile *f = new TFile(infile.c_str(), "r");

  TH1D *h_pt2_correction[100];
  TH1D *h_pt2_bin_correction[centrality_bins];
  TH1D *h_pt2_bin_log_correction[centrality_bins];

  TH1D *h_centrality = (TH1D*) f->Get("h_centrality");

  TH1D *h_jet_spectra[100];

  TH1D *h_jet_spectra_meas[centrality_bins];
  TH1D *h_scaled_events = new TH1D("h_scaled_events","", 4, -0.5, 3.5);
 
  int no_hist[4] = {0};
  for ( int i = 0; i < 100; i++)
    {

      h_pt2_correction[i] = new TH1D(Form("h_pt2_correction_%d", i),";#it{p}_{T};Subleading Efficiency", 45, 5, 50);

      h_jet_spectra[i] = (TH1D*) f->Get(Form("h_jet_spectra_etacut_%d", i));
      if (h_jet_spectra[i]->Integral() == 0) continue;
      for (int j = 0; j < 4; j++)
	{
	  std::cout << j << std::endl;      
	  if (icentrality_bins[j] <= i && icentrality_bins[j+1] > i)
	    {
	      int cbin = 0;
	      if (!no_hist[j])
		{
		  no_hist[j] = 1;
		  std::cout << "Making the histos " << j << std::endl;      
		  h_jet_spectra_meas[j] = (TH1D*) h_jet_spectra[i]->Clone();
		  h_jet_spectra_meas[j]->SetName(Form("h_jet_spectra_meas_%d", j));
		  h_jet_spectra_meas[j]->Reset();
		  cbin = h_centrality_data[j]->FindBin(i);
		  
		  float n_events_MB = h_centrality->GetBinContent(i+1);
		  float n_events_jet = h_centrality_data[j]->GetBinContent(cbin);
		  h_centrality_data[j]->Scale(n_events_MB/n_events_jet);
		  std::cout << n_events_MB << " / " << n_events_jet << std::endl;
		}

	      cbin = h_centrality_data[j]->FindBin(i);
	      float scale = h_centrality_data[j]->GetBinContent(cbin)/h_centrality->GetBinContent(i+1);
	      h_jet_spectra_meas[j]->Add(h_jet_spectra[i], scale);
	      h_scaled_events->Fill(j, h_centrality->GetBinContent(i+1)*scale);
	    }
	};
    }

  for (int i = 0; i < centrality_bins; i++)
    {
      h_jet_spectra_meas[i]->Scale(1./h_scaled_events->GetBinContent(i+1),"width");
      h_pt2_bin_correction[i] = new TH1D(Form("h_pt2_bin_correction_%d", i),";#it{p}_{T} [GeV];Subleading Efficiency", 45, 5, 50);
      h_pt2_bin_log_correction[i] = new TH1D(Form("h_pt2_bin_log_correction_%d", i),";#it{p}_{T} [GeV] ;Subleading Efficiency", nbins, ipt_bins);
      for (int j = 0 ; j < h_pt2_bin_correction[i]->GetNbinsX(); j++)
	{
	  double bin_low_edge = h_pt2_bin_correction[i]->GetBinLowEdge(j+1);
	  double bin_high_edge = h_pt2_bin_correction[i]->GetBinWidth(j+1) + bin_low_edge;
	  double bin_center = h_pt2_bin_correction[i]->GetBinCenter(j+1);
	  int bin_number1 = h_jet_spectra_meas[i]->FindBin(bin_center);
	  int bin_number2 = h_jet_spectra_meas[i]->FindBin(bin_high_edge);
	  double bin1_integral = h_jet_spectra_meas[i]->Integral(bin_number1, bin_number2, "width")*scale_factor;
	  double bin2_integral = h_jet_spectra_meas[i]->Integral(bin_number2, -1, "width")*scale_factor_above;
	  double prob = TMath::Exp(-1*(bin1_integral + bin2_integral));
	  h_pt2_bin_correction[i]->SetBinContent(j+1, prob);
	}
      for (int j = 0 ; j < h_pt2_bin_log_correction[i]->GetNbinsX(); j++)
	{
	  double bin_low_edge = h_pt2_bin_log_correction[i]->GetBinLowEdge(j+1);
	  double bin_high_edge = h_pt2_bin_log_correction[i]->GetBinWidth(j+1) + bin_low_edge;
	  double bin_center = h_pt2_bin_log_correction[i]->GetBinCenter(j+1);
	  int bin_number1 = h_jet_spectra_meas[i]->FindBin(bin_center);
	  int bin_number2 = h_jet_spectra_meas[i]->FindBin(bin_high_edge);
	  double bin1_integral = h_jet_spectra_meas[i]->Integral(bin_number1, bin_number2, "width")*scale_factor;
	  double bin2_integral = h_jet_spectra_meas[i]->Integral(bin_number2, -1, "width")*scale_factor_above;
	  double prob = TMath::Exp(-1*(bin1_integral + bin2_integral));
	  h_pt2_bin_log_correction[i]->SetBinContent(j+1, prob);
	}
    }

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  TCanvas *c = new TCanvas("c","c", 1000, 450);
  c->Divide(2, 1);

  c->cd(1);
  gPad->SetLogy();
  
  dlutility::SetLineAtt(h_jet_spectra_meas[0], kRed, 2, 1);
  dlutility::SetLineAtt(h_jet_spectra_meas[1], kBlue, 2, 1);
  dlutility::SetLineAtt(h_jet_spectra_meas[2], kGreen, 2, 1);
  dlutility::SetLineAtt(h_jet_spectra_meas[3], kOrange, 2, 1);
  dlutility::SetFont(h_jet_spectra_meas[0], 42, 0.05);
  gPad->SetLeftMargin(0.22);
  //  h_pt2_bin_correction[0]->SetMaximum(1.3);
  h_jet_spectra_meas[0]->GetXaxis()->SetRangeUser(5, 30);
  h_jet_spectra_meas[0]->SetTitle("; #it{p}_{T} [GeV] ; #frac{1}{N_{evt}}#frac{dN_{jet}}{d#it{p}_{T}}");
  h_jet_spectra_meas[0]->GetYaxis()->SetTitleOffset(1.6);

  h_jet_spectra_meas[0]->Draw("hist");
  h_jet_spectra_meas[1]->Draw("hist same");
  h_jet_spectra_meas[2]->Draw("hist same");
  h_jet_spectra_meas[3]->Draw("hist same");

  dlutility::DrawSPHENIX(0.5, 0.85);

  TLegend *leg = new TLegend(0.65, 0.6, 0.85, 0.75);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_jet_spectra_meas[0], "0-10%","l");
  leg->AddEntry(h_jet_spectra_meas[1], "10-30%","l");
  leg->AddEntry(h_jet_spectra_meas[2], "30-50%","l");
  leg->AddEntry(h_jet_spectra_meas[3], "50-90%","l");
  leg->Draw("same");
  
  c->cd(2);
  gPad->SetLeftMargin(0.22);
  dlutility::SetLineAtt(h_pt2_bin_correction[0], kRed, 2, 1);
  dlutility::SetLineAtt(h_pt2_bin_correction[1], kBlue, 2, 1);
  dlutility::SetLineAtt(h_pt2_bin_correction[2], kGreen, 2, 1);
  dlutility::SetLineAtt(h_pt2_bin_correction[3], kOrange, 2, 1);
  dlutility::SetFont(h_pt2_bin_correction[0], 42, 0.05);

  h_pt2_bin_correction[0]->SetMaximum(1.5);
  h_pt2_bin_correction[0]->GetXaxis()->SetRangeUser(5, 30);

  h_pt2_bin_correction[0]->Draw("hist");
  h_pt2_bin_correction[1]->Draw("hist same");
  h_pt2_bin_correction[2]->Draw("hist same");
  h_pt2_bin_correction[3]->Draw("hist same");

  TLine *l1 = new TLine(5, 1, 30, 1);
  l1->SetLineWidth(2);
  l1->SetLineStyle(4);
  l1->SetLineColor(kBlack);
  l1->Draw("same");
  dlutility::DrawSPHENIX(0.27, 0.85);

  leg = new TLegend(0.65, 0.73, 0.85, 0.9);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_pt2_bin_correction[0], "0-10%","l");
  leg->AddEntry(h_pt2_bin_correction[1], "10-30%","l");
  leg->AddEntry(h_pt2_bin_correction[2], "30-50%","l");
  leg->AddEntry(h_pt2_bin_correction[3], "50-90%","l");
  leg->Draw("same");
  
  c->Print(Form("%s/unfolding_plots/probabilities_AA_r%02d.pdf", rb.get_code_location().c_str(), cone_size));
	   
  TString unfoldpath = rb.get_code_location() + "/unfolding_hists/probability_hists_AA_r0" + std::to_string(cone_size);
  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");

  for (int i = 0 ; i < 100; i++)
    {
      h_pt2_correction[i]->Write();
    }
  for (int i = 0 ; i < centrality_bins; i++)
    {
      h_pt2_bin_correction[i]->Write();
      h_pt2_bin_log_correction[i]->Write();
    }
  fout->Close();

}
