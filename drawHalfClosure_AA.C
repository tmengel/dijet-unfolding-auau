#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
const bool NUCLEAR = true;
const int color_unfold_fill = kAzure - 4;
const int color_unfold = kAzure - 5;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

const int color_pythia = kRed;
const float msize_pythia = 0.9;
const float marker_pythia = 20;
const float lsize_pythia = 1.1;

const int color_herwig = kViolet;
const float msize_herwig = 0.9;
const float marker_herwig = 20;
const float lsize_herwig = 1.1;

const int color_reco = kRed;
const float marker_reco = 24;
const float msize_reco = 0.9;
const float lsize_reco = 1.1;
const int color_data = kAzure - 6;
const float marker_data = 24;
const float msize_data = 0.9;
const float lsize_data = 1.1;
void drawHalfClosure_pp(const int cone_size = 3)
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  const int niterations = 30;
  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();


  const int nbins = read_nbins;
  int first_bin = 0;
  float ipt_bins[nbins+1];
  double  dxj_bins[nbins+1];

  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] > 0.3 && first_bin == 0) first_bin = i;

    }

  
  Int_t max_reco_bin = rb.get_maximum_reco_bin();
  
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
      std::cout << ipt_bins[measure_bins[ir]] << " -- " <<  ipt_bins[subleading_measure_bins[ir]] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  

  TFile *fin = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_HALF_nominal.root ", rb.get_code_location().c_str(), cone_size),"r");
  if (!fin)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  TH1D *h_xj_truth = (TH1D*) fin->Get("h_xj_truth");
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = (TH1D*) fin->Get(Form("h_xj_unfold_iter%d", iter));
    }

  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);

  TH1D *h_closure_test[niterations];
  
  for (int niter = 0; niter < niterations; niter++)
    {
      //for (int irange = 0; irange < mbins; irange++)
      //{
      cxj->cd(1);
      dlutility::SetLineAtt(h_xj_unfold[niter], color_unfold, lsize_unfold, 1);
      dlutility::SetMarkerAtt(h_xj_unfold[niter], color_unfold, msize_unfold, marker_unfold);

      dlutility::SetLineAtt(h_xj_truth, color_pythia, lsize_pythia, 1);
      dlutility::SetMarkerAtt(h_xj_truth, color_pythia, msize_pythia, marker_pythia);

      dlutility::SetFont(h_xj_truth, 42, 0.05);

      h_xj_truth->SetMaximum(5);
      h_xj_truth->SetMinimum(0);
      h_xj_truth->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      TH1D *ht = (TH1D*) h_xj_truth->Clone();
      TH1D *hu = (TH1D*) h_xj_unfold[niter]->Clone();

      ht->Draw("p E1");
      hu->Draw("same p E1");


      dlutility::DrawSPHENIX(0.22, 0.84);
      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[1]], ipt_bins[measure_bins[2]]), 0.22, 0.69);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);

      TLegend *leg = new TLegend(0.2, 0.4, 0.4, 0.56);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_xj_truth, "Training Truth","p");
      leg->AddEntry(h_xj_unfold[niter], "Testing Unfold","p");
      leg->Draw("same");

      cxj->cd(2);

      TH1D *h_data_compare = (TH1D*) h_xj_truth->Clone();

      h_data_compare->Reset();
      for (int i = 0; i < h_data_compare->GetNbinsX(); i++)
	{
	  if (h_xj_truth->GetBinContent(i+1) == 0) continue;
	  h_data_compare->SetBinContent(i+1, (h_xj_unfold[niter]->GetBinContent(i+1) - h_xj_truth->GetBinContent(i+1))/(h_xj_truth->GetBinContent(i+1)));
	}
      h_data_compare->SetTitle(";x_{J}; Truth - Unfold / Truth");
      dlutility::SetFont(h_data_compare, 42, 0.1, 0.07, 0.07, 0.07);
      dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
      dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);
      h_data_compare->SetMaximum(0.2);
      h_data_compare->SetMinimum(-0.2);

      TH1D *hd = (TH1D*) h_data_compare->Clone();

      hd->Draw("p");

      h_closure_test[niter] = (TH1D*) hd->Clone();
      h_closure_test[niter]->SetName(Form("h_closure_test_%d", niter));
	  
      TLine *line = new TLine(hd->GetBinLowEdge(1), 1, 1, 1);
      line->SetLineStyle(4);
      line->SetLineColor(kRed + 3);
      line->SetLineWidth(2);
      line->Draw("same");
      cxj->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter %d.png", rb.get_code_location().c_str(), cone_size, 1, niter));
      cxj->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_%d.pdf", rb.get_code_location().c_str(), cone_size, 1, niter));
    }

  TCanvas *cxjclos = new TCanvas("cxjclos","cxjclos", 700, 500);
  dlutility::createCutCanvas(cxjclos);
  cxjclos->cd(1);

  dlutility::SetLineAtt(h_closure_test[0], kBlue, 2, 1);
  dlutility::SetLineAtt(h_closure_test[1], kCyan, 2, 1);
  dlutility::SetLineAtt(h_closure_test[3], kGreen, 2, 1);
  dlutility::SetLineAtt(h_closure_test[5], kViolet, 2, 1);
  dlutility::SetLineAtt(h_closure_test[7], kRed, 2, 1);
  dlutility::SetFont(h_closure_test[0], 42, 0.04);
  h_closure_test[0]->SetMinimum(-0.2);
  h_closure_test[0]->SetMaximum(0.2);
  h_closure_test[0]->Draw("hist");
  //h_closure_test[0]->SetTitleOffset(0.9);
  h_closure_test[0]->SetTitle(";x_{J}; Half Truth - Half Unfold / Half Truth");
  h_closure_test[1]->Draw("hist same");
  h_closure_test[3]->Draw("hist same");
  h_closure_test[5]->Draw("hist same");
  h_closure_test[7]->Draw("hist same");    

  cxjclos->cd(2);
  int mrange = 1;
  dlutility::DrawSPHENIXcut(0.02, 0.9, 0, 0.08);
  dlutility::drawText(Form("anti-#it{k}_{t} #kern[-0.3]{#it{R}} = %0.1f", 0.1*cone_size), 0.02, 0.78, 0, kBlack, 0.08);
  dlutility::drawText(Form("%2.1f #kern[-0.1]{#leq #it{p}_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[mrange]], ipt_bins[measure_bins[mrange+1]]), 0.02, 0.73, 0, kBlack, 0.08);
  dlutility::drawText(Form("#it{p}_{T,2} #kern[-0.15]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.02,0.68, 0, kBlack, 0.08);
  dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 7#pi/8}", 0.02,0.63, 0, kBlack, 0.08);

  dlutility::drawText("Closure Test", 0.02,0.53, 0, kBlack, 0.08);

  TLegend *leg4 = new TLegend(0.02, 0.25, 0.9, 0.48);
  leg4->SetTextSize(0.08);
  leg4->SetTextFont(42);
  leg4->SetLineWidth(0);
  leg4->AddEntry(h_closure_test[0],"1 iteration","l");
  leg4->AddEntry(h_closure_test[1],"2 iteration","l");
  leg4->AddEntry(h_closure_test[3],"4 iteration","l");
  leg4->AddEntry(h_closure_test[5],"6 iteration","l");
  leg4->AddEntry(h_closure_test[7],"8 iteration","l");
  leg4->Draw();
  cxjclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_1_iter_all.pdf", rb.get_code_location().c_str(), cone_size));
  cxjclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_1_iter_all.png", rb.get_code_location().c_str(), cone_size));

  TCanvas *cclos = new TCanvas("cxjclos","cxjclos", 1000, 800);
  cclos->Divide(3,2);

  int nbins1 = h_closure_test[0]->GetNbinsX();
  
  for (int i = 0; i < 6; i++)
    {
      cclos->cd(i+1);
     
      TH1D *h_niter_closure = new TH1D("h_niter_closure",";N_{iterations}; Rel. Error", niterations, 0, niterations);
      TH1D *h_niter_stat = new TH1D("h_niter_stat",";N_{iterations}; Rel. Error", niterations, 0, niterations);
      TH1D *h_niter_quad = new TH1D("h_niter_quad",";N_{iterations}; Rel. Error", niterations, 0, niterations);
		  
      for (int iter = 0; iter < niterations; iter++)
	{
	  std::cout << fabs(h_closure_test[iter]->GetBinCenter(nbins1 - i)) << std::endl;
	  float closure = fabs(h_closure_test[iter]->GetBinContent(nbins1 - i));
	  h_niter_closure->SetBinContent(iter + 1, fabs(h_closure_test[iter]->GetBinContent(nbins1 - i)));
	  float stat=fabs(h_xj_unfold[iter]->GetBinError(nbins1 - i)/h_xj_truth->GetBinContent(nbins1 - i));
	  h_niter_stat->SetBinContent(iter + 1, stat);
	  h_niter_quad->SetBinContent(iter + 1, sqrt(stat*stat + closure*closure));
	}

      h_niter_closure->SetLineColor(kBlack);
      //h_niter_stat->SetLineColor(kBlue);
      //h_niter_closure->SetLineColor(kBlack);
      h_niter_closure->SetMaximum(0.15);
      h_niter_closure->Draw("hist");

      dlutility::drawText(Form("x_{J} = %0.2f", h_closure_test[0]->GetBinCenter(nbins1 - i)), 0.22, 0.8);


    }

  cclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_1_iter_xjbin.pdf", rb.get_code_location().c_str(), cone_size));
      
  return;
}
