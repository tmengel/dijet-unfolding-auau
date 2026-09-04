#include "dlUtility.h"
#include "PlotUtils.h"
#include "read_binning.h"
#include "histo_opps.h"

const bool NUCLEAR = true;

const int color_pp_unfold_fill = kBlack;
const int color_pp_unfold = kBlack;
const float marker_pp_unfold = 24;
const float msize_pp_unfold = 0.9;
const float lsize_pp_unfold = 1.1;

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
const int color_pp_data = kAzure - 6;
const float marker_pp_data = 24;
const float msize_pp_data = 0.9;
const float lsize_pp_data = 1.1;

// Lays four square pads with no internal gap onto canvas c, in reading order
// (top-left, top-right, bottom-left, bottom-right), so the resulting 2x2
// grid shares one x-axis (bottom row only) and one y-axis (left column only).
// void make2x2SharedAxisPads(TCanvas *c, TPad *pads[4], const char *prefix,
// 			    double padL = 0.14, double padB = 0.14,
// 			    double padT = 0.06, double padR = 0.035)
// {
//   const double plotW = 1.0 - padL - padR;
//   const double plotH = 1.0 - padB - padT;
//   const double xmid = padL + plotW / 2.0;
//   const double ymid = padB + plotH / 2.0;

//   c->cd();
//   pads[0] = new TPad(Form("%s_tl", prefix), "", 0.0, ymid, xmid, 1.0);
//   pads[1] = new TPad(Form("%s_tr", prefix), "", xmid, ymid, 1.0, 1.0);
//   pads[2] = new TPad(Form("%s_bl", prefix), "", 0.0, 0.0, xmid, ymid);
//   pads[3] = new TPad(Form("%s_br", prefix), "", xmid, 0.0, 1.0, ymid);

//   // SetMargin(left, right, bottom, top); inner edges (facing the other
//   // panels) get zero margin so the four plot areas butt up against
//   // each other with no gap.
//   pads[0]->SetMargin(padL / xmid, 0.0, 0.0, padT / (1.0 - ymid));
//   pads[1]->SetMargin(0.0, padR / (1.0 - xmid), 0.0, padT / (1.0 - ymid));
//   pads[2]->SetMargin(padL / xmid, 0.0, padB / ymid, 0.0);
//   pads[3]->SetMargin(0.0, padR / (1.0 - xmid), padB / ymid, 0.0);

//   for (int ip = 0; ip < 4; ip++)
//     {
//       pads[ip]->SetTickx(1);
//       pads[ip]->SetTicky(1);
//       pads[ip]->Draw();
//     }
// }

// Make a clean 2x2 grid with touching panel boundaries.
// Only the outside of the grid has margins; internal boundaries touch.
void make2x2SharedAxisPads(TCanvas *c, TPad *pads[4],
                           const char *prefix,
                           double left   = 0.13,
                           double right  = 0.035,
                           double bottom = 0.12,
                           double top    = 0.045)
{
    const double xmid = 0.5;
    const double ymid = 0.5;

    c->cd();

    // top-left, top-right, bottom-left, bottom-right
    pads[0] = new TPad(Form("%s_tl", prefix), "", 0.0, ymid, xmid, 1.0);
    pads[1] = new TPad(Form("%s_tr", prefix), "", xmid, ymid, 1.0, 1.0);
    pads[2] = new TPad(Form("%s_bl", prefix), "", 0.0, 0.0, xmid, ymid);
    pads[3] = new TPad(Form("%s_br", prefix), "", xmid, 0.0, 1.0, ymid);

    // Outer margins only.
    // Internal boundaries have zero margin so the panels touch.
    pads[0]->SetMargin(left / xmid, 0.0, 0.0, top / (1.0 - ymid));
    pads[1]->SetMargin(0.0, right / (1.0 - xmid), 0.0, top / (1.0 - ymid));

    pads[2]->SetMargin(left / xmid, 0.0, bottom / ymid, 0.0);
    pads[3]->SetMargin(0.0, right / (1.0 - xmid), bottom / ymid, 0.0);

    for (int i = 0; i < 4; ++i)
    {
        pads[i]->SetBorderMode(0);
        pads[i]->SetFrameBorderMode(0);
        pads[i]->SetFillStyle(0);
        pads[i]->SetTickx(1);
        pads[i]->SetTicky(1);
        pads[i]->Draw();
    }
}
void drawFinalUnfold_AA_cent_grid(const int cone_size = 3, const std::string configfile = "configs/binning_AA.config")
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  // dlutility::SetyjPadStyle();
  
  PlotUtils::set_sphenix_style();


  int color_unfold_fill[5] = {kRed + 1, kAzure + 2, kGreen + 3, kViolet - 1, kBlack};
  int color_unfold[5] = {kRed + 1, kAzure + 2, kGreen+3, kViolet -1, kBlack};
  int marker_unfold[5] = {20, 21, 33, 34, 24};
  float msize_unfold[5] = {1.1, 1, 1.5, 1.2, 1.1};
  float lsize_unfold = 1.8;


  read_binning rb(configfile.c_str());

  Double_t first_xj = rb.get_first_xj();

  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  std::string dphi_string = rb.get_dphi_string();

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins+1];
  rb.get_centrality_bins(icentrality_bins);

  if (cent_bins != 4)
    {
      std::cout << "drawFinalUnfold_AA_cent_grid expects exactly 4 centrality bins to fill a 2x2 grid, got " << cent_bins << std::endl;
      return;
    }

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
      if (dxj_bins[i] > 0.2 && first_bin == 0) first_bin = i;

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


  TFile *finjw = new TFile(Form("%s/truth_hists/JEWEL_for_Dan.root",  rb.get_code_location().c_str()),"r");
  if (!finjw)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }

  TH1D *h_jewel_xj_med = (TH1D*) finjw->Get("h1_medium");
  TH1D *h_jewel_xj_vac = (TH1D*) finjw->Get("h1_vacuum");
  dlutility::SetLineAtt(h_jewel_xj_med, kGreen + 2, 2, 1);
  dlutility::SetLineAtt(h_jewel_xj_vac, kBlack, 2, 1);


  TFile *fintr = new TFile(Form("%s/truth_hists/truth_hist_r%02d.root",  rb.get_code_location().c_str(),  cone_size),"r");
  if (!fintr)
    {
      std::cout << "no truth hists" << std::endl;
      return;
    }

  TH1D *h_linear_truth_xj[3];
  for (int i = 0; i < 3; i++)
    {
      h_linear_truth_xj[i] = (TH1D*) fintr->Get(Form("h_linear_truth_xj_%d", i));
      h_linear_truth_xj[i]->SetName(Form("h_linear_pythia_xj_%d", i));
      if (!h_linear_truth_xj[i])
	{
	  std::cout << "AHH"<< std::endl;
	  return;
	}
    }


  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  const int niterations = 10;
  // pp nominal iteration (no centrality dependence, matches createResponse.C).
  const int niter = 1;
  // AA nominal iteration (0-indexed) -> N_iter = 2; matches the prior_iteration
  // constant in createResponse_noempty_AA.cxx. Kept as a per-centrality helper
  // so the AA and pp choices stay independent at each call site.
  auto niter_for_cent = [](int /*cent*/) { return 1; };
  //pp results
  TProfile *hp_xj_rms_pp[mbins][niterations];
  TFile *finupp = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finupp)
    {
      std::cout << " no unc " << std::endl;
      return;
    }

  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  hp_xj_rms_pp[irange][iter] = (TProfile*) finupp->Get(Form("hp_xj_range_%d_%d", irange, iter));
	}
    }

  TProfile *hp_xj_rms[cent_bins][mbins][niterations];

  TFile *finu[cent_bins];
  for (int i = 0; i < cent_bins; i++)
    {
      finu[i] = new TFile(Form("%s/uncertainties/uncertainties_AA_cent_%d_r%02d_nominal.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!finu[i])
	{
	  std::cout << " no unc " << std::endl;
	  return;
	}

      for (int irange = 0; irange < mbins; irange++)
	{
	  for (int iter = 0; iter < niterations; ++iter)
	    {
	      hp_xj_rms[i][irange][iter] = (TProfile*) finu[i]->Get(Form("hp_xj_range_%d_%d", irange, iter));
	      hp_xj_rms[i][irange][iter]->SetName(Form("hp_xj_rms_%d_%d_%d", i, irange, iter));
	    }
	}
    }
  TH1D *h_total_sys_range[cent_bins][mbins][niterations];
  TH1D *h_total_sys_neg_range[cent_bins][mbins][niterations];

  TFile *fins[cent_bins];
  for (int i = 0; i < cent_bins; i++)
    {
      fins[i] = new TFile(Form("%s/uncertainties/systematics_AA_cent_%d_r%02d.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!fins[i])
	{
	  std::cout << " no sys fin " << i << std::endl;
	  return;
	}
      for (int irange = 0; irange < mbins; irange++)
	{
	  for (int iter = 0; iter < niterations; ++iter)
	    {
	      h_total_sys_range[i][irange][iter] = (TH1D*) fins[i]->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	      h_total_sys_neg_range[i][irange][iter] = (TH1D*) fins[i]->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));

	    }
	}
    }
  TH1D *h_total_pp_sys_range[mbins][niterations];
  TH1D *h_total_pp_sys_neg_range[mbins][niterations];

  TFile *fins_pp = new TFile(Form("%s/uncertainties/systematics_pp_r%02d.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!fins_pp)
    {
      std::cout << " no sys pp"  << std::endl;
      return;
    }
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0; iter < niterations; ++iter)
	{
	  h_total_pp_sys_range[irange][iter] = (TH1D*) fins_pp->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	  h_total_pp_sys_neg_range[irange][iter] = (TH1D*) fins_pp->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));
	  h_total_pp_sys_range[irange][iter]->SetName(Form("h_total_pp_sys_range_%d_%d", irange, iter));
	  h_total_pp_sys_neg_range[irange][iter]->SetName(Form("h_total_pp_sys_neg_range_%d_%d", irange, iter));
	}
    }


  TFile *fin[cent_bins];
  TH1D *h_flat_data_pt1pt2[cent_bins];
  TH1D *h_flat_unfold_pt1pt2[cent_bins][niterations];

  for (int i = 0; i <  cent_bins; i++)
    {
      fin[i]= new TFile(Form("%s/unfolding_hists/unfolding_hists_AA_cent_%d_r%02d_nominal.root",  rb.get_code_location().c_str(), i, cone_size),"r");
      if (!fin[i])
	{
	  std::cout << " no file " << std::endl;
	  return;
	}
      h_flat_data_pt1pt2[i] = (TH1D*) fin[i]->Get("h_data_flat_pt1pt2");
      h_flat_data_pt1pt2[i]->SetName(Form("h_flat_data_pt1pt2_%d", i));
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_flat_unfold_pt1pt2[i][iter] = (TH1D*) fin[i]->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	  h_flat_unfold_pt1pt2[i][iter]->SetName(Form("h_flat_unfold_pt1pt2_%d_%d", i, iter));
	}
    }


  TFile *finpp = new TFile(Form("%s/unfolding_hists/unfolding_hists_pp_r%02d_nominal.root",  rb.get_code_location().c_str(), cone_size),"r");
  if (!finpp)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  TH1D *h_flat_pp_data_pt1pt2 = (TH1D*) finpp->Get("h_data_flat_pt1pt2");
  h_flat_pp_data_pt1pt2->SetName("h_flat_pp_data_pt1pt2");

  TH1D *h_flat_pp_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_pp_unfold_pt1pt2[iter] = (TH1D*) finpp->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_pp_unfold_pt1pt2[iter]->SetName(Form("h_flat_pp_unfold_pt1pt2_%d", iter));
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  TH2D *h_pt1pt2_data[cent_bins];
  TH2D *h_pt1pt2_pp_data = new TH2D("h_pt1pt2_pp_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[cent_bins][niterations];
  TH2D *h_pt1pt2_pp_unfold[niterations];
  for (int i = 0; i < cent_bins; i++)
    {
      h_pt1pt2_data[i] = new TH2D(Form("h_pt1pt2_data_%d", i), ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int i = 0; i < cent_bins; i++)
	{
	  h_pt1pt2_unfold[i][iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
	  h_pt1pt2_unfold[i][iter]->SetName(Form("h_pt1pt2_unfold_%d_iter%d", i, iter));
	}
      h_pt1pt2_pp_unfold[iter] = new TH2D("h_pt1pt2_pp_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_pp_unfold[iter]->SetName(Form("h_pt1pt2_pp_unfold_iter%d", iter));

    }

  TH1D *h_xj_data[cent_bins];
  TH1D *h_xj_pp_data = new TH1D("h_xj_pp_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[cent_bins][niterations];
  TH1D *h_xj_pp_unfold[niterations];
  for (int i = 0; i < cent_bins; i++)
    {
      h_xj_data[i]= new TH1D(Form("h_xj_data_%d", i), ";x_{J};", nbins, ixj_bins);
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      for (int i = 0; i < cent_bins; i++)
	{
	  h_xj_unfold[i][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", i, iter), ";x_{J};",nbins, ixj_bins);
	}
      h_xj_pp_unfold[iter] = new TH1D(Form("h_xj_pp_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data_range[cent_bins][mbins];
  TH1D *h_xj_pp_data_range[mbins];
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_xj_pp_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {

      h_xj_pp_data_range[irange] = new TH1D(Form("h_xj_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0 ; i < cent_bins; i++)
	{
	  h_xj_data_range[i][irange] = new TH1D(Form("h_xj_data_range_%d_%d",i, irange), ";x_{J};", nbins, ixj_bins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_xj_unfold_range[i][irange][iter] = new TH1D(Form("h_xj_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_pp_unfold_range[irange][iter] = new TH1D(Form("h_xj_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  TH1D *h_xjunc_data[cent_bins];
  TH1D *h_xjunc_pp_data = new TH1D("h_xjunc_pp_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_truth = new TH1D("h_xjunc_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xjunc_unfold[cent_bins][niterations];
  TH1D *h_xjunc_pp_unfold[niterations];
  for (int i = 0 ; i < cent_bins; i++)
    {
      h_xjunc_data[i]= new TH1D(Form("h_xjunc_data_%d", i), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xjunc_unfold[i][iter] = new TH1D(Form("h_xjunc_unfold_%d_iter%d", i, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xjunc_pp_unfold[iter] = new TH1D(Form("h_xjunc_pp_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TH1D *h_xjunc_data_range[cent_bins][mbins];
  TH1D *h_xjunc_pp_data_range[mbins];
  TH1D *h_xjunc_reco_range[mbins];
  TH1D *h_xjunc_truth_range[mbins];
  TH1D *h_xjunc_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_xjunc_pp_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {

      h_xjunc_pp_data_range[irange] = new TH1D(Form("h_xjunc_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_reco_range[irange] = new TH1D(Form("h_xjunc_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xjunc_truth_range[irange] = new TH1D(Form("h_xjunc_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0; i < cent_bins; i++)
	{
	  h_xjunc_data_range[i][irange] = new TH1D(Form("h_xjunc_data_range_%d_%d", i, irange), ";x_{J};", nbins, ixj_bins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_xjunc_unfold_range[i][irange][iter] = new TH1D(Form("h_xjunc_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	  for (int iter = 0; iter < niterations; iter++)
	    h_xjunc_pp_unfold_range[irange][iter] = new TH1D(Form("h_xjunc_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  histo_opps::make_sym_pt1pt2(h_flat_pp_data_pt1pt2, h_pt1pt2_pp_data, nbins);
  for (int i = 0; i < cent_bins; i++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2[i], h_pt1pt2_data[i], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[i][iter], h_pt1pt2_unfold[i][iter], nbins);
	}

    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_pp_unfold_pt1pt2[iter], h_pt1pt2_pp_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  histo_opps::project_xj(h_pt1pt2_pp_data, h_xj_pp_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);

  for (int i = 0; i < cent_bins; i++)
    {
      histo_opps::project_xj(h_pt1pt2_data[i], h_xj_data[i], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[i][iter], h_xj_unfold[i][iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_pp_unfold[iter], h_xj_pp_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }

  histo_opps::normalize_histo(h_xj_pp_data, nbins);
  for (int i = 0 ; i < cent_bins; i++)
    {
      histo_opps::normalize_histo(h_xj_data[i], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_unfold[i][iter], nbins);
	}
    }
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_pp_unfold[iter], nbins);
    }
  if (NUCLEAR) std::cout << __LINE__ << std::endl;

  TH1D *h_final_xj_pp_data_range[mbins];
  TH1D *h_final_xj_pp_unfold_range[mbins][niterations];

  TH1D *h_final_xj_pp_systematics[mbins][niterations];
  TGraphAsymmErrors *g_final_xj_pp_systematics[mbins][niterations];

  TH1D *h_final_xj_data_range[cent_bins][mbins];
  TH1D *h_final_xj_unfold_range[cent_bins][mbins][niterations];
  TH1D *h_final_xj_systematics[cent_bins][mbins][niterations];
  TGraphAsymmErrors *g_final_xj_systematics[cent_bins][mbins][niterations];

  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_pp_data_range[irange] = new TH1D(Form("h_final_xj_pp_data_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int i = 0; i < cent_bins; i++)
	{
	  h_final_xj_data_range[i][irange] = new TH1D(Form("h_final_xj_data_range_%d_%d", i, irange), ";x_{J};", nbins, ixj_bins);

	  for (int iter = 0; iter < niterations; iter++)
	    {
	      h_final_xj_unfold_range[i][irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_%d_iter%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	      h_final_xj_systematics[i][irange][iter] = new TH1D(Form("h_final_xj_systematics_%d_%d_%d", i, irange, iter), ";x_{J};",nbins, ixj_bins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{

	  h_final_xj_pp_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_pp_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	  h_final_xj_pp_systematics[irange][iter] = new TH1D(Form("h_final_xj_pp_systematics_%d_%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int irange = 0; irange < mbins; irange++)
    {

      histo_opps::project_xj(h_pt1pt2_pp_data, h_xj_pp_data_range[irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);

      for (int i = 0; i < cent_bins; i++)
	{
	  histo_opps::project_xj(h_pt1pt2_data[i], h_xj_data_range[i][irange], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::project_xj(h_pt1pt2_unfold[i][iter], h_xj_unfold_range[i][irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_pp_unfold[iter], h_xj_pp_unfold_range[irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}


      histo_opps::normalize_histo(h_xj_pp_data_range[irange], nbins);
      for (int i = 0 ; i < cent_bins; i++)
	{
	  histo_opps::normalize_histo(h_xj_data_range[i][irange], nbins);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::normalize_histo(h_xj_unfold_range[i][irange][iter], nbins);
	      histo_opps::normalize_histo(hp_xj_rms[i][irange][iter], nbins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_xj_pp_unfold_range[irange][iter], nbins);
	  histo_opps::normalize_histo(hp_xj_rms_pp[irange][iter], nbins);
	}


      histo_opps::finalize_xj(h_xj_pp_data_range[irange], h_final_xj_pp_data_range[irange], nbins, first_xj);
      for (int i = 0; i < cent_bins; i++)
	{
	  histo_opps::finalize_xj(h_xj_data_range[i][irange], h_final_xj_data_range[i][irange], nbins, first_xj);
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      histo_opps::finalize_xj(h_xj_unfold_range[i][irange][iter], h_final_xj_unfold_range[i][irange][iter], nbins, first_xj);
	      histo_opps::set_xj_errors(h_final_xj_unfold_range[i][irange][iter], hp_xj_rms[i][irange][iter], nbins);


	      h_final_xj_systematics[i][irange][iter] = (TH1D*) h_final_xj_unfold_range[i][irange][iter]->Clone();
	      h_final_xj_systematics[i][irange][iter]->SetName(Form("h_final_xj_systematics_%d_%d_%d", i, irange, iter));
	      g_final_xj_systematics[i][irange][iter] = new TGraphAsymmErrors(h_final_xj_systematics[i][irange][iter]);
	      g_final_xj_systematics[i][irange][iter]->SetName(Form("g_final_xj_systematics_%d_%d_%d", i, irange, iter));
	      histo_opps::get_xj_systematics(g_final_xj_systematics[i][irange][iter], h_total_sys_neg_range[i][irange][iter], h_total_sys_range[i][irange][iter], nbins);
	    }
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_pp_unfold_range[irange][iter], h_final_xj_pp_unfold_range[irange][iter], nbins, first_xj);
	  histo_opps::set_xj_errors(h_final_xj_pp_unfold_range[irange][iter], hp_xj_rms_pp[irange][iter], nbins);

	  h_final_xj_pp_systematics[irange][iter] = (TH1D*) h_final_xj_pp_unfold_range[irange][iter]->Clone();
	  h_final_xj_pp_systematics[irange][iter]->SetName(Form("h_final_xj_pp_systematics_%d_%d", irange, iter));
	  g_final_xj_pp_systematics[irange][iter] = new TGraphAsymmErrors(h_final_xj_pp_systematics[irange][iter]);
	  g_final_xj_pp_systematics[irange][iter]->SetName(Form("g_final_xj_pp_systematics_%d_%d", irange, iter));
	  histo_opps::get_xj_systematics(g_final_xj_pp_systematics[irange][iter], h_total_pp_sys_neg_range[irange][iter], h_total_pp_sys_range[irange][iter], nbins);

	}

    }

  // Override the pp central values / systematics computed above with the pre-built
  // final pp results in final_plots_pp_r03.root, instead of the values just derived
  // from uncertainties_pp/systematics_pp/unfolding_hists_pp.
  {
    TFile *fpp_final = new TFile(Form("%s/final_plots_pp_r03.root", rb.get_code_location().c_str()), "r");
    if (!fpp_final || fpp_final->IsZombie())
      {
	std::cout << " no pp final results file " << std::endl;
	return;
      }
    // Each pre-built pp point is placed in the bin that contains its own x,
    // rather than counting up from a fixed starting bin. final_plots_pp_r03.root
    // is produced separately, so the bin its graphs start at moves whenever that
    // file is regenerated; assuming a fixed offset silently shifts the markers
    // relative to their systematic boxes.
    for (int irange = 0; irange < mbins; irange++)
      {
	for (int iter = 0; iter < niterations; iter++)
	  {
	    TGraphAsymmErrors *g_stat = (TGraphAsymmErrors*) fpp_final->Get(Form("g_final_xj_statistics_%d_%d", irange, iter));
	    TGraphAsymmErrors *g_sys  = (TGraphAsymmErrors*) fpp_final->Get(Form("g_final_xj_systematics_%d_%d", irange, iter));

	    TH1D *h_pp = h_final_xj_pp_unfold_range[irange][iter];
	    h_pp->Reset();
	    g_final_xj_pp_systematics[irange][iter] = new TGraphAsymmErrors(h_pp);
	    for (int ib = 1; ib <= nbins; ib++)
	      {
		g_final_xj_pp_systematics[irange][iter]->SetPoint(ib - 1, h_pp->GetBinCenter(ib), 0);
		g_final_xj_pp_systematics[irange][iter]->SetPointError(ib - 1, h_pp->GetBinWidth(ib)/2., h_pp->GetBinWidth(ib)/2., 0, 0);
	      }

	    int npts = g_stat ? g_stat->GetN() : 0;
	    for (int j = 0; j < npts; j++)
	      {
		double x, y, sx, sy;
		g_stat->GetPoint(j, x, y);
		const int bin = h_pp->FindBin(x);
		if (bin < 1 || bin > nbins) continue;
		h_pp->SetBinContent(bin, y);
		h_pp->SetBinError(bin, g_stat->GetErrorYhigh(j));

		g_sys->GetPoint(j, sx, sy);
		// Anchor the systematic box on the same bin centre as the marker.
		g_final_xj_pp_systematics[irange][iter]->SetPoint(bin - 1, h_pp->GetBinCenter(bin), sy);
		double exlow = g_final_xj_pp_systematics[irange][iter]->GetErrorXlow(bin - 1);
		double exhigh = g_final_xj_pp_systematics[irange][iter]->GetErrorXhigh(bin - 1);
		g_final_xj_pp_systematics[irange][iter]->SetPointError(bin - 1, exlow, exhigh, g_sys->GetErrorYlow(j), g_sys->GetErrorYhigh(j));
	      }
	  }
      }

    for (int irange = 0; irange < mbins; irange++)
      {
	TH1D *h_data_file = (TH1D*) fpp_final->Get(Form("h_final_xj_data_range_%d", irange));
	if (h_data_file) h_final_xj_pp_data_range[irange] = (TH1D*) h_data_file->Clone(Form("h_final_xj_pp_data_range_%d", irange));
      }
  }

  // Draw the four centrality bins as a 2x2 grid of pads sharing one x-axis
  // (bottom row) and one y-axis (left column), reading order top-left ->
  // bottom-right = most central -> most peripheral. One canvas per pT1 range.
  const char *cent_label[4] = {"0 - 10%", "10 - 30%", "30 - 50%", "50 - 90%"};
  const int quad_to_cent[4] = {0, 1, 2, 3};

  for (int irange = 0; irange < mbins; irange++)
    {
      TCanvas *cgrid = new TCanvas(Form("cxj_grid_range_%d", irange), "cxj_grid", 1000, 850 );
      cgrid->Divide(2, 2);
      // TPad *pads[4];
      // make2x2SharedAxisPads(cgrid, pads, Form("pad_r%d", irange));

      for (int i = 1; i <= 4; ++i) {
          cgrid->cd(i);
          gPad->SetLeftMargin(i == 1 || i == 3 ? 0.14 : 0.0);
          gPad->SetRightMargin(i == 2 || i == 4 ? 0.035 : 0.0);
          gPad->SetBottomMargin(i == 3 || i == 4 ? 0.14 : 0.0);
          gPad->SetTopMargin(i == 1 || i == 2 ? 0.05 : 0.0);
      }

      for (int iq = 0; iq < 4; iq++)
	{
    cgrid->cd(iq + 1);
	  // const int ic = quad_to_cent[iq];
	  // const bool bottomRow = (iq == 2 || iq == 3);
	  // const bool leftCol   = (iq == 0 || iq == 2);

     const int ic = quad_to_cent[iq];
    const bool bottomRow = (iq == 2 || iq == 3);
    const bool leftCol   = (iq == 0 || iq == 2);
	  // pads[iq]->cd();

	  TH1D *hh = new TH1D(
        Form("hh_grid_%d_%d", irange, ic),
        "",
        nbins,
        ixj_bins
    );

      hh->SetMinimum(0);
      hh->SetMaximum(3);

      hh->SetTitle(Form(";%s;%s",
                        bottomRow ? "x_{J}" : "",
                        leftCol ? "#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}" : ""));

      dlutility::SetFont(hh, 42, 0.075, 0.075, 0.065, 0.065);

      // Better axis spacing
      hh->GetXaxis()->SetTitleOffset(1.0);
      hh->GetYaxis()->SetTitleOffset(1.0);

      // Labels only on shared axes
      if (!bottomRow)
      {
          hh->GetXaxis()->SetLabelSize(0);
          hh->GetXaxis()->SetTitleSize(0);
      }

      if (!leftCol)
      {
          hh->GetYaxis()->SetLabelSize(0);
          hh->GetYaxis()->SetTitleSize(0);
      }

      // Make ticks consistent
      hh->GetXaxis()->SetTickLength(0.025);
      hh->GetYaxis()->SetTickLength(0.025);

      hh->Draw("hist");

	  dlutility::SetLineAtt(h_final_xj_pp_unfold_range[irange][niter], color_pp_unfold, lsize_pp_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_pp_unfold_range[irange][niter], color_pp_unfold, msize_pp_unfold, marker_pp_unfold);
	  dlutility::SetLineAtt(g_final_xj_pp_systematics[irange][niter], color_pp_unfold, lsize_pp_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_xj_pp_systematics[irange][niter], color_pp_unfold, msize_pp_unfold, marker_pp_unfold);
	  g_final_xj_pp_systematics[irange][niter]->SetFillColorAlpha(color_pp_unfold_fill, 0.25);

	  dlutility::SetLineAtt(h_final_xj_unfold_range[ic][irange][niter_for_cent(ic)], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[ic][irange][niter_for_cent(ic)], color_unfold[ic], msize_unfold[ic], marker_unfold[ic]);
	  dlutility::SetLineAtt(g_final_xj_systematics[ic][irange][niter_for_cent(ic)], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_xj_systematics[ic][irange][niter_for_cent(ic)], color_unfold[ic], msize_unfold[ic], marker_unfold[ic]);
	  g_final_xj_systematics[ic][irange][niter_for_cent(ic)]->SetFillColorAlpha(color_unfold_fill[ic], 0.3);

	  g_final_xj_pp_systematics[irange][niter]->Draw("same p E2");
	  h_final_xj_pp_unfold_range[irange][niter]->Draw("same p E1");
	  g_final_xj_systematics[ic][irange][niter_for_cent(ic)]->Draw("same p E2");
	  h_final_xj_unfold_range[ic][irange][niter_for_cent(ic)]->Draw("same p E1");

	  dlutility::drawText(cent_label[ic], 0.25, 0.88, 0, kBlack, 0.09);

	  if (iq == 0)
	    {
	      dlutility::DrawSPHENIXboth(0.25, 0.78, 1, 0.075);
	      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size * 0.1), 0.25, 0.63, 0, kBlack, 0.065);
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange + 1]]), 0.25, 0.55, 0, kBlack, 0.065);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.25, 0.47, 0, kBlack, 0.065);
	      dlutility::drawText(Form("#Delta#phi #geq %s", dphi_string.c_str()), 0.25, 0.39, 0, kBlack, 0.065);

	      TLegend *leg = new TLegend(0.25, 0.14, 0.85, 0.34);
	      leg->SetLineWidth(0);
	      leg->SetTextSize(0.065);
	      leg->SetTextFont(42);
	      leg->AddEntry(h_final_xj_unfold_range[ic][irange][niter_for_cent(ic)], "Au+Au data", "pf");
	      leg->AddEntry(h_final_xj_pp_unfold_range[irange][niter], "#it{p}+#it{p} data", "pf");
	      leg->Draw("same");
	    }
	}

      cgrid->cd();
      cgrid->Print(Form("%s/final_plots/h_final_xj_unfolded_AA_cent_grid_r%02d_range_%d.png", rb.get_code_location().c_str(), cone_size, irange));
      cgrid->Print(Form("%s/final_plots/h_final_xj_unfolded_AA_cent_grid_r%02d_range_%d.pdf", rb.get_code_location().c_str(), cone_size, irange));
    }

  return;
}
