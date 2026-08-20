#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TStyle.h"

#include <iostream>
#include <string>

namespace {

TH1D *require_hist(TFile *file, const char *name, const char *path)
{
  TH1D *hist = file ? static_cast<TH1D *>(file->Get(name)) : nullptr;
  if (!hist)
    std::cerr << "Missing histogram " << name << " in " << path << std::endl;
  return hist;
}

TProfile *require_profile(TFile *file, const char *name, const char *path)
{
  TProfile *hist = file ? static_cast<TProfile *>(file->Get(name)) : nullptr;
  if (!hist)
    std::cerr << "Missing profile " << name << " in " << path << std::endl;
  return hist;
}

TFile *open_required(const TString &path)
{
  TFile *file = TFile::Open(path, "READ");
  if (!file || file->IsZombie())
    {
      std::cerr << "Could not open " << path << std::endl;
      return nullptr;
    }
  return file;
}

}

void drawFinalUnfold_AA_only(const int cone_size = 3,
                             const int centrality_bin = 0,
                             const std::string configfile = "binning_AA.config",
                             const int niter_arg = -1)
{
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.0001);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());
  const std::string system_string = "AA_cent_" + std::to_string(centrality_bin);
  const int nbins = rb.get_nbins();
  const int mbins = rb.get_measure_bins();
  const double first_xj = rb.get_first_xj();
  const std::string dphi_string = rb.get_dphi_string();

  float ipt_bins[nbins + 1];
  float ixj_bins[nbins + 1];
  double dxj_bins[nbins + 1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  int first_bin = 0;
  for (int i = 0; i < nbins + 1; ++i)
    {
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] >= first_xj && first_bin == 0)
        first_bin = i;
    }

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins + 1];
  rb.get_centrality_bins(icentrality_bins);

  int measure_bins[10] = {0};
  for (int irange = 0; irange < mbins + 1; ++irange)
    measure_bins[irange] = rb.get_measure_region(irange);

  const int measure_subleading_bin = rb.get_measure_subleading_bin();
  const int niterations = 10;
  // Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2; matches the
  // prior_iteration constant in createResponse_noempty_AA.cxx.
  const int niter = (niter_arg >= 0) ? niter_arg : 1;
  if (niter < 0 || niter >= niterations)
    {
      std::cerr << "Requested iteration " << niter << " is outside [0, "
                << niterations - 1 << "]" << std::endl;
      return;
    }

  TString unfold_path = Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_nominal.root",
                             rb.get_code_location().c_str(), system_string.c_str(), cone_size);
  TString unc_path = Form("%s/uncertainties/uncertainties_%s_r%02d_nominal.root",
                          rb.get_code_location().c_str(), system_string.c_str(), cone_size);
  TString sys_path = Form("%s/uncertainties/systematics_%s_r%02d.root",
                          rb.get_code_location().c_str(), system_string.c_str(), cone_size);

  TFile *funfold = open_required(unfold_path);
  TFile *func = open_required(unc_path);
  TFile *fsys = open_required(sys_path);
  if (!funfold || !func || !fsys)
    return;

  TH1D *h_flat_data_pt1pt2 = require_hist(funfold, "h_data_flat_pt1pt2", unfold_path);
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    h_flat_unfold_pt1pt2[iter] = require_hist(funfold, Form("h_flat_unfold_pt1pt2_%d", iter), unfold_path);

  if (!h_flat_data_pt1pt2 || !h_flat_unfold_pt1pt2[niter])
    return;

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data_AA_only", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    h_pt1pt2_unfold[iter] = new TH2D(Form("h_pt1pt2_unfold_AA_only_%d", iter), ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);

  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; ++iter)
    if (h_flat_unfold_pt1pt2[iter])
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);

  TCanvas *canvas = new TCanvas("c_final_AA_only", "c_final_AA_only", 500, 700);
  dlutility::ratioPanelCanvas(canvas);
  TFile *final_output = TFile::Open(Form(
    "%s/final_plots/final_xj_%s_r%02d_AAonly.root",
    rb.get_code_location().c_str(), system_string.c_str(), cone_size), "RECREATE");

  for (int irange = 0; irange < mbins; ++irange)
    {
      TProfile *hp_xj_rms = require_profile(func, Form("hp_xj_range_%d_%d", irange, niter), unc_path);
      TH1D *h_total_sys = require_hist(fsys, Form("h_total_sys_range_%d_iter_%d", irange, niter), sys_path);
      TH1D *h_total_sys_neg = require_hist(fsys, Form("h_total_sys_neg_range_%d_iter_%d", irange, niter), sys_path);
      if (!hp_xj_rms || !h_total_sys || !h_total_sys_neg)
        continue;

      TH1D *h_xj_data = new TH1D(Form("h_xj_data_AA_only_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      TH1D *h_xj_unfold = new TH1D(Form("h_xj_unfold_AA_only_range_%d", irange), ";x_{J};", nbins, ixj_bins);

      histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_bins[irange], measure_bins[irange + 1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_unfold[niter], h_xj_unfold, nbins, measure_bins[irange], measure_bins[irange + 1], measure_subleading_bin, nbins - 2);

      histo_opps::normalize_histo(h_xj_data, nbins);
      histo_opps::normalize_histo(h_xj_unfold, nbins);
      histo_opps::normalize_histo(hp_xj_rms, nbins);

      TH1D *h_final_data = new TH1D(Form("h_final_data_AA_only_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      TH1D *h_final_unfold = new TH1D(Form("h_final_unfold_AA_only_range_%d", irange), ";x_{J};", nbins, ixj_bins);

      histo_opps::finalize_xj(h_xj_data, h_final_data, nbins, first_xj);
      histo_opps::finalize_xj(h_xj_unfold, h_final_unfold, nbins, first_xj);
      histo_opps::set_xj_errors(h_final_unfold, hp_xj_rms, nbins);

      TGraphAsymmErrors *g_systematics = new TGraphAsymmErrors(h_final_unfold);
      g_systematics->SetName(Form("g_final_xj_systematics_range_%d", irange));
      histo_opps::get_xj_systematics(g_systematics, h_total_sys_neg, h_total_sys, nbins);

      TH1D *frame = static_cast<TH1D *>(h_final_unfold->Rebin(nbins - first_bin, Form("h_final_frame_AA_only_%d", irange), &dxj_bins[first_bin]));
      TH1D *data = static_cast<TH1D *>(h_final_data->Rebin(nbins - first_bin, Form("h_final_data_rebin_AA_only_%d", irange), &dxj_bins[first_bin]));

      dlutility::SetLineAtt(frame, kAzure - 6, 1.1, 1);
      dlutility::SetMarkerAtt(frame, kAzure - 6, 1.0, 20);
      dlutility::SetLineAtt(data, kGray + 2, 1.0, 1);
      dlutility::SetMarkerAtt(data, kGray + 2, 0.9, 24);
      dlutility::SetLineAtt(g_systematics, kAzure - 6, 1.0, 1);
      dlutility::SetMarkerAtt(g_systematics, kAzure - 6, 1.0, 20);
      g_systematics->SetFillColorAlpha(kAzure - 6, 0.3);

      canvas->cd(1);
      frame->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
      frame->SetMinimum(0);
      frame->SetMaximum(4.5);
      dlutility::SetFont(frame, 42, 0.05);
      frame->Draw("p E1");
      g_systematics->Draw("same p E2");
      data->Draw("same p E1");
      frame->Draw("same p E1");

      dlutility::DrawSPHENIXboth(0.22, 0.84, 1);
      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size * 0.1), 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange + 1]]), 0.22, 0.69);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText(Form("#Delta#phi #geq %s", dphi_string.c_str()), 0.22, 0.59);
      dlutility::drawText(Form("%d - %d %%", static_cast<int>(icentrality_bins[centrality_bin]), static_cast<int>(icentrality_bins[centrality_bin + 1])), 0.22, 0.54);

      TLegend *legend = new TLegend(0.58, 0.69, 0.86, 0.89);
      legend->SetLineWidth(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(frame, "AA unfolded", "p");
      legend->AddEntry(data, "AA reco data", "p");
      legend->Draw("same");

      canvas->cd(2);
      TH1D *compare = static_cast<TH1D *>(frame->Clone(Form("h_compare_AA_only_%d", irange)));
      compare->Divide(data);
      compare->SetTitle(";x_{J};Unfold / Reco Data");
      compare->SetMinimum(0);
      compare->SetMaximum(2.3);
      dlutility::SetFont(compare, 42, 0.1, 0.07, 0.07, 0.07);
      dlutility::SetLineAtt(compare, kBlack, 1, 1);
      dlutility::SetMarkerAtt(compare, kBlack, 1, 8);
      compare->Draw("p");
      TLine *line = new TLine(compare->GetBinLowEdge(1), 1, 1, 1);
      line->SetLineStyle(4);
      line->SetLineColor(kRed + 3);
      line->SetLineWidth(2);
      line->Draw("same");

      canvas->Print(Form("%s/final_plots/h_final_xj_unfolded_%s_r%02d_range_%d_AAonly.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange));
      canvas->Print(Form("%s/final_plots/h_final_xj_unfolded_%s_r%02d_range_%d_AAonly.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, irange));
      if (final_output && !final_output->IsZombie())
        {
          final_output->cd();
          frame->Write(Form("h_final_xj_unfold_range_%d", irange));
          data->Write(Form("h_final_xj_data_range_%d", irange));
          g_systematics->Write();
        }
    }
  if (final_output)
    {
      final_output->Write();
      final_output->Close();
    }
}
