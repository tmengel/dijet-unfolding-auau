#include "dlUtility.h"
#include "read_binning.h"
void drawResponse_AA(const int cone_size = 3, const int centrality_bin = 0)
{

  read_binning rb(std::getenv("AUAU_CONFIG"));

  Int_t read_nbins = rb.get_nbins();

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

  float sample_boundary[4] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;


  TString responsepath = Form("response_matrices/response_matrix_AA_cent_%d_r%02d_nominal.root", centrality_bin, cone_size);
  TFile *fr = new TFile(responsepath.Data(),"r");

  TH2D *h_flat_response_skim = (TH2D*) fr->Get("h_flat_response_skim");

  TH1D *h_xj_reco = (TH1D*) fr->Get("h_xj_reco");
  TH1D *h_xj_truth = (TH1D*) fr->Get("h_xj_truth");
  TH1D *h_xj_unfold[10];
  for (int iter =0; iter < 10; iter++)
    {
      h_xj_unfold[iter] = (TH1D*) fr->Get(Form("h_xj_unfold_iter%d", iter));
    }

  TString response_halfpath = Form("response_matrices/response_matrix_AA_cent_%d_r%02d_HALF.root", centrality_bin, cone_size);
  TFile *frh = new TFile(response_halfpath.Data(),"r");

  TH1D *h_xj_half_reco = (TH1D*) frh->Get("h_xj_reco");
  TH1D *h_xj_half_truth = (TH1D*) frh->Get("h_xj_truth");
  TH1D *h_xj_half_unfold[10];
  for (int iter =0; iter < 10; iter++)
    {
      h_xj_half_unfold[iter] = (TH1D*) frh->Get(Form("h_xj_unfold_iter%d", iter));
    }

  
  dlutility::SetyjPadStyle();
  TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 500);
  cresponseskim->SetLogz();

  gPad->SetRightMargin(0.05);
  h_flat_response_skim->SetTitle(";#it{p}_{T,1}^{reco, bin} #times nbins + #it{p}_{T,2}^{reco, bin};#it{p}_{T,1}^{truth, bin} #times nbins + #it{p}_{T,2}^{truth, bin}");
  h_flat_response_skim->Draw("col");
  //dlutility::DrawSPHENIX(0.1, 0.95, 0.04, 0, 1, 1);
  string sPHENIX_MARK = "#bf{#it{sPHENIX}} #it{Simulation}";
  string species = "Au+Au #kern[-0.2]{#sqrt{s_{NN}}} = 200 GeV";
  //string extratext = "#it{Internal}";
  

  dlutility::drawText(sPHENIX_MARK.c_str(), 0.15 ,0.95, 0, kBlack, 0.035);
  dlutility::drawText(species.c_str(), 0.85 ,0.95, 1, kBlack, 0.035);
  dlutility::drawText("HIJING+PYTHIA8", 0.19, 0.88);
  dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.19, 0.83);

  cresponseskim->Print(Form("%s/unfolding_plots/response_matrix_AA_cent_%d.pdf", rb.get_code_location().c_str(), centrality_bin));

  std::cout << "bins: " << h_flat_response_skim->GetXaxis()->GetNbins() << " by " << h_flat_response_skim->GetYaxis()->GetNbins() << std::endl;
  TCanvas *cresponseskimzoom = new TCanvas("fdsz","fdsz", 500, 500);
  cresponseskimzoom->SetLogz();

  gPad->SetRightMargin(0.05);
  h_flat_response_skim->GetXaxis()->SetRangeUser(46, 60);
  h_flat_response_skim->GetYaxis()->SetRangeUser(121, 140);
  h_flat_response_skim->Draw("col");
  dlutility::DrawSPHENIX_Prelim(0.93, 0.4);
  dlutility::drawText("Response matrix", 0.93, 0.3, 1);

  cresponseskimzoom->Print("response_matrix_zoom.pdf");

  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);

  dlutility::ratioPanelCanvas(cxj);
  for (int iter =0 ; iter < 10; iter++)
    {
      cxj->cd(1);
      dlutility::SetLineAtt(h_xj_unfold[iter], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_xj_unfold[iter], kBlack, 1, 8);
  
      dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
      dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);
      
      dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
      dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
  
      dlutility::SetFont(h_xj_truth, 42, 0.04);
      h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
      h_xj_truth->SetMaximum(5);
      h_xj_truth->Draw("p");
      h_xj_unfold[iter]->Draw("same p");
      h_xj_reco->Draw("p same");

      dlutility::DrawSPHENIX(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      dlutility::drawText(Form("Full Closure - N_{iter} = %d", iter + 1), 0.22, 0.54);
      TLegend *leg = new TLegend(0.22, 0.25, 0.4, 0.5);
      leg->SetLineWidth(0);
      leg->AddEntry(h_xj_reco, "PYTHIA-8 Reco");
      leg->AddEntry(h_xj_truth, "PYTHIA-8 Truth");
      leg->AddEntry(h_xj_unfold[iter], "Unfolded");
      leg->Draw("same");

    
      cxj->cd(2);

      TH1D *h_reco_compare = (TH1D*) h_xj_unfold[iter]->Clone();
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
      cxj->Print(Form("full_closure_iter_%d.pdf", iter));
    }

  TCanvas *cxjh = new TCanvas("cxjh","cxjh", 500, 700);

  dlutility::ratioPanelCanvas(cxjh);
  for (int iter =0 ; iter < 10; iter++)
    {
      cxjh->cd(1);
      dlutility::SetLineAtt(h_xj_half_unfold[iter], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_xj_half_unfold[iter], kBlack, 1, 8);
  
      dlutility::SetLineAtt(h_xj_half_truth, kRed, 2, 1);
      dlutility::SetMarkerAtt(h_xj_half_truth, kRed, 1, 8);
      
      dlutility::SetLineAtt(h_xj_half_reco, kBlue, 2, 1);
      dlutility::SetMarkerAtt(h_xj_half_reco, kBlue, 1, 8);
  
      dlutility::SetFont(h_xj_half_truth, 42, 0.04);
      h_xj_half_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
      h_xj_half_truth->SetMaximum(5);
      h_xj_half_truth->Draw("p");
      h_xj_half_unfold[iter]->Draw("same p");
      h_xj_half_reco->Draw("p same");

      dlutility::DrawSPHENIX(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      dlutility::drawText(Form("Half Closure - N_{iter} = %d", iter + 1), 0.22, 0.54);
      TLegend *leg = new TLegend(0.22, 0.25, 0.4, 0.5);
      leg->SetLineWidth(0);
      leg->AddEntry(h_xj_reco, "PYTHIA-8 Reco");
      leg->AddEntry(h_xj_truth, "PYTHIA-8 Truth");
      leg->AddEntry(h_xj_unfold[iter], "Unfolded");
      leg->Draw("same");
    
      cxjh->cd(2);

      TH1D *h_reco_compare = (TH1D*) h_xj_half_unfold[iter]->Clone();
      h_reco_compare->Divide(h_xj_half_truth);
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
      cxjh->Print(Form("half_closure_iter_%d.pdf", iter));
    }

  return;
}
