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
void drawFullClosure_AA(const int cone_size = 3, const int centrality_bin = 0)
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  const int niterations = 10;
  read_binning rb("binning_AA.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins+1];
  rb.get_centrality_bins(icentrality_bins);

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
  

  TFile *fin = new TFile(Form("%s/response_matrices/response_matrix_AA_cent_%d_r%02d.root", rb.get_code_location().c_str(), centrality_bin, cone_size),"r");
  if (!fin)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }

  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_reco, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }

  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_truth_range[irange], nbins);
      histo_opps::normalize_histo(h_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold_range[irange][iter], nbins);
	}


      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

    }

  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  for (int niter = 0; niter < niterations; niter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  cxj->cd(1);
	  dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	  dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	  dlutility::SetLineAtt(h_final_xj_reco_range[irange], color_reco, lsize_reco, 1);
	  dlutility::SetMarkerAtt(h_final_xj_reco_range[irange], color_reco, msize_reco, marker_reco);

	  dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.05);

	  h_final_xj_truth_range[irange]->SetMaximum(5);
	  h_final_xj_truth_range[irange]->SetMinimum(0);
	  h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	  TH1D *ht = (TH1D*) h_final_xj_truth_range[irange]->Rebin(nbins - first_bin, Form("h_rebin_truth_%d", irange), &dxj_bins[first_bin]);
	  TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);

	  ht->Draw("p E1");
	  hu->Draw("same p E1");

	  h_final_xj_reco_range[irange]->Draw("same p");

      
	  dlutility::DrawSPHENIX(0.22, 0.84);
	  dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, 0.74);
	  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, 0.69);
	  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
	  dlutility::drawText(Form("%d - %d %%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]), 0.22, 0.54);
	  TLegend *leg = new TLegend(0.2, 0.4, 0.4, 0.56);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.04);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final_xj_truth_range[irange], "Training Truth","p");
	  leg->AddEntry(h_final_xj_reco_range[irange], "Training Reco ","p");
	  leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Training Unfold","p");
	  leg->Draw("same");

	  cxj->cd(2);

	  TH1D *h_data_compare = (TH1D*) h_final_xj_truth_range[irange]->Clone();
	  h_data_compare->Divide(h_final_xj_unfold_range[irange][niter]);
	  h_data_compare->SetTitle(";x_{J}; Testing / Training");
	  dlutility::SetFont(h_data_compare, 42, 0.1, 0.07, 0.07, 0.07);
	  dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
	  dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);
	  h_data_compare->SetMaximum(2.3);
	  h_data_compare->SetMinimum(0.0);
	  TH1D *hd = (TH1D*) h_data_compare->Rebin(nbins - first_bin, "h_rebin_compare", &dxj_bins[first_bin]);

	  hd->Draw("p");

	  TLine *line = new TLine(hd->GetBinLowEdge(1), 1, 1, 1);
	  line->SetLineStyle(4);
	  line->SetLineColor(kRed + 3);
	  line->SetLineWidth(2);
	  line->Draw("same");
	  cxj->Print(Form("%s/unfolding_plots/h_xj_full_closure_AA_cent_%d_r%02d_range_%d_iter %d.png", rb.get_code_location().c_str(), centrality_bin, cone_size, irange, niter));
	  cxj->Print(Form("%s/unfolding_plots/h_xj_full_closure_AA_cent_%d_r%02d_range_%d_iter_%d.pdf", rb.get_code_location().c_str(), centrality_bin, cone_size, irange, niter));
	}
    }
  return;
}
