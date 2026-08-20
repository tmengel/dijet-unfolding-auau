#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

namespace {

// Draw one iteration-scan panel (sigma_sim, sigma_data, sigma_bin-by-bin and
// their quadrature sum vs N_iter) and write it to <stem>.pdf and <stem>.png.
// yMax > 0 fixes the vertical range; otherwise it is auto-scaled with enough
// headroom that the N_iter = 1 point clears the sPHENIX/centrality labels.
void drawIterationScan(TH1D *h_stat, TH1D *h_unfold, TH1D *h_binbybin,
                       TH1D *h_total, const char *canvas_name,
                       const std::string &cent_label, const std::string &stem,
                       double yMax = -1)
{
  TCanvas *c = new TCanvas(canvas_name, canvas_name, 500, 500);

  h_total->SetMinimum(0.);
  h_total->SetMaximum(yMax > 0 ? yMax : 1.6*h_total->GetMaximum());

  dlutility::SetMarkerAtt(h_total, kBlack, 1, 8);
  dlutility::SetLineAtt(h_total, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_binbybin, kRed, 1, 8);
  dlutility::SetLineAtt(h_binbybin, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_stat, kBlue, 1, 8);
  dlutility::SetLineAtt(h_stat, kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_unfold, kGreen, 1, 8);
  dlutility::SetLineAtt(h_unfold, kGreen, 1, 1);

  h_total->Draw("p hist");
  h_stat->Draw("same p hist");
  h_unfold->Draw("same p hist");
  h_binbybin->Draw("same p hist");
  h_total->Draw("p hist same");

  dlutility::DrawSPHENIX(0.22, 0.87);
  dlutility::drawText(cent_label.c_str(), 0.22, 0.77);

  TLegend *leg = new TLegend(0.218, 0.566, 0.397, 0.726);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_stat, "#sigma_{sim}","p");
  leg->AddEntry(h_unfold, "#sigma_{data}","p");
  leg->AddEntry(h_binbybin, "#sigma_{bin-by-bin}","p");
  leg->AddEntry(h_total, "#sigma_{conv} = #sqrt{#sigma^{2}_{sim} + #sigma^{2}_{data} + #sigma^{2}_{bin-by-bin}}","p");
  leg->Draw("same");

  c->SaveAs((stem + ".pdf").c_str());
  c->SaveAs((stem + ".png").c_str());
}

}

void makeIterationPlot_AA(const int cone_size = 3, const int centrality_bin = 0, const int prior = 0)
{
  const int niterations = 10;

  std::string sysname = "nominal";
  if (prior == 1) sysname = "PRIMER1";
  if (prior == 2) sysname = "PRIMER2";
  
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  // read_binning rb("binning_AA.config");
  read_binning rb(std::getenv("AUAU_CONFIG"));

  Int_t read_nbins = rb.get_nbins();
  
  Double_t dphicut = rb.get_dphicut();

  const int n_centrality_bins = rb.get_number_centrality_bins();  
  float icentrality_bins[n_centrality_bins+1];

  rb.get_centrality_bins(icentrality_bins);

  const int nbins = read_nbins;
  const int nbins2 = nbins*nbins;

  Double_t ipt_bins[nbins+1];
  Double_t ixj_bins[nbins+1];

  float fipt_bins[nbins+1];
  float fixj_bins[nbins+1];

  rb.get_pt_bins(fipt_bins);
  rb.get_xj_bins(fixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      ipt_bins[i] = fipt_bins[i];
      ixj_bins[i] = fixj_bins[i];
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
  float low_trigger[3] = {0};

  const int mbins = rb.get_measure_bins();
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

  
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ib = 0; ib < 3; ib++)
    {
      low_trigger[ib] = rb.get_low_trigger(ib);
      std::cout <<  low_trigger[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  
  std::cout << __LINE__ << std::endl;

  TFile *f_uncertainties = new TFile(Form("%s/uncertainties/uncertainties_AA_cent_%d_r%02d_%s.root", rb.get_code_location().c_str(), centrality_bin, cone_size, sysname.c_str()),"r");

  if (!f_uncertainties)
    {
      std::cout << "no uncertainities" << std::endl;
      return;
    }
  TProfile *hp_xj[niterations];
  TProfile *hp_pt1pt2[niterations];
    
  TFile *fin = new TFile(Form("%s/unfolding_hists/unfolding_hists_AA_cent_%d_r%02d_%s.root", rb.get_code_location().c_str(), centrality_bin, cone_size, sysname.c_str()),"r");
  if (!fin)
    {
      std::cout << "no file" << std::endl;
      return;
    }

  TH1D *h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      hp_xj[iter] = (TProfile*) f_uncertainties->Get(Form("hp_xj_%d", iter));
      hp_pt1pt2[iter] = (TProfile*) f_uncertainties->Get(Form("hp_pt1pt2_%d", iter));
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }
  int nbins_pt1pt2 = hp_pt1pt2[0]->GetNbinsX();
  
  std::cout << __LINE__ << std::endl;
  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";#it{p}_{T,1};#it{p}_{T,2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }
  std::cout << __LINE__ << std::endl;
  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  TH1D *h_xj_profile_unfold[niterations];
  TH1D *h_pt1pt2_profile_unfold[niterations];  
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
      h_xj_profile_unfold[iter] = new TH1D(Form("h_xj_profile_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
      h_pt1pt2_profile_unfold[iter] = new TH1D(Form("h_pt1pt2_profile_unfold_iter%d", iter), ";#it{p}_{T,1} + #it{p}_{T,2} bin;",nbins_pt1pt2, 0, nbins_pt1pt2);
    }

  std::cout << __LINE__ << std::endl;
  TH1D *h_xjunc_data = new TH1D("h_xjunc_data", ";x_{J};", nbins, ixj_bins);
  std::cout << __LINE__ << std::endl;
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      std::cout << __LINE__ << std::endl;
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  std::cout << __LINE__ << std::endl;
  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin + 1, measure_leading_bin + 3, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin + 1, measure_leading_bin + 3, measure_subleading_bin, nbins - 2);//measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
      histo_opps::tprofile_to_histo(hp_xj[iter], h_xj_profile_unfold[iter], nbins);
      histo_opps::tprofile_to_histo(hp_pt1pt2[iter], h_pt1pt2_profile_unfold[iter], nbins_pt1pt2);
    }
  histo_opps::normalize_histo(h_xj_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
      histo_opps::normalize_histo(h_xj_profile_unfold[iter], nbins);
    }

  TH1D *h_final_xj_data = (TH1D*) h_xj_data->Clone();
  h_final_xj_data->SetName("h_final_xj_data");
  h_final_xj_data->Reset();
  histo_opps::finalize_xj(h_xj_data,h_final_xj_data, nbins, 0.4);
  TH1D * h_final_xj_unfold[niterations];
  TH1D * h_final_xj_profile_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_final_xj_unfold[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
      h_final_xj_unfold[iter]->SetName(Form("h_final_xj_unfold_%d", iter));
      h_final_xj_unfold[iter]->Reset();
      h_final_xj_profile_unfold[iter] = (TH1D*) h_xj_profile_unfold[iter]->Clone();
      h_final_xj_profile_unfold[iter]->SetName(Form("h_final_xj_profile_unfold_%d", iter));
      h_final_xj_profile_unfold[iter]->Reset();

      histo_opps::finalize_xj(h_xj_unfold[iter], h_final_xj_unfold[iter], nbins, 0.4);
      histo_opps::finalize_xj(h_xj_profile_unfold[iter], h_final_xj_profile_unfold[iter], nbins, 0.4);
    }  


  int nbins_flat = h_flat_unfold_pt1pt2[0]->GetNbinsX();
  
  // Three variants of the same scan:
  //   _abs  -- absolute changes/errors, in raw counts.
  //   _rel  -- fractional (per-bin) changes/errors, kept for comparison. On a
  //            steeply falling p_{T,1},p_{T,2} spectrum this is dominated by the
  //            sparse high-p_{T} bins, which is why it is not used for the choice.
  //   _norm -- each distribution first normalised to unit integral, then
  //            absolute changes/errors. This is the AA note's definition
  //            (absolute change of the normalised distribution) and reproduces
  //            the note's iteration figure; the nominal N_iter comes from it.
  const char *abs_axes = "; N_{iter}; #sigma_{conv} [counts]";
  const char *rel_axes = "; N_{iter}; #sigma_{conv} (fractional)";
  const char *norm_axes = "; N_{iter}; #sigma_{conv}";

  TH1D *h_statistical_uncertainties = new TH1D("h_statistical_uncertainties_abs", abs_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_unfold_uncertainties = new TH1D("h_unfold_uncertainties_abs", abs_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_binbybin_uncertainties = new TH1D("h_binbybin_uncertainties_abs", abs_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_total_uncertainties = new TH1D("h_total_uncertainties_abs", abs_axes, niterations + 1, -0.5, niterations + 0.5);

  TH1D *h_statistical_uncertainties_rel = new TH1D("h_statistical_uncertainties_rel", rel_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_unfold_uncertainties_rel = new TH1D("h_unfold_uncertainties_rel", rel_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_binbybin_uncertainties_rel = new TH1D("h_binbybin_uncertainties_rel", rel_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_total_uncertainties_rel = new TH1D("h_total_uncertainties_rel", rel_axes, niterations + 1, -0.5, niterations + 0.5);

  TH1D *h_statistical_uncertainties_norm = new TH1D("h_statistical_uncertainties_norm", norm_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_unfold_uncertainties_norm = new TH1D("h_unfold_uncertainties_norm", norm_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_binbybin_uncertainties_norm = new TH1D("h_binbybin_uncertainties_norm", norm_axes, niterations + 1, -0.5, niterations + 0.5);
  TH1D *h_total_uncertainties_norm = new TH1D("h_total_uncertainties_norm", norm_axes, niterations + 1, -0.5, niterations + 0.5);
  std::cout << __LINE__ << std::endl;

  // Normalisation factors for the _norm variant. Held separately rather than
  // scaling the histograms in place, so the _abs and _rel variants keep raw counts.
  auto safe_integral = [](double v) { return v != 0 ? v : 1.0; };
  const double norm_data = safe_integral(h_flat_data_pt1pt2->Integral());
  double norm_unfold[niterations];
  double norm_sim[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    {
      norm_unfold[iter] = safe_integral(h_flat_unfold_pt1pt2[iter]->Integral());
      norm_sim[iter] = safe_integral(h_pt1pt2_profile_unfold[iter]->Integral());
    }

  for (int iter = 0; iter < niterations; ++iter)
    {
      Double_t stat_unc = 0;
      Double_t unfold_unc = 0;
      Double_t binbybin_unc = 0;

      Double_t stat_unc_rel = 0;
      Double_t unfold_unc_rel = 0;
      Double_t binbybin_unc_rel = 0;

      Double_t stat_unc_norm = 0;
      Double_t unfold_unc_norm = 0;
      Double_t binbybin_unc_norm = 0;

      for (int ibin = 0; ibin < nbins_flat; ibin++)
	{

	  float bin_cont = h_flat_unfold_pt1pt2[iter]->GetBinContent(ibin+1);
	  float prev_cont = h_flat_data_pt1pt2->GetBinContent(ibin+1);

	  if (iter > 0)
	    {
	      prev_cont = h_flat_unfold_pt1pt2[iter - 1]->GetBinContent(ibin+1);
	    }

	  if (bin_cont == 0) continue;

	  const float unfold_err = h_flat_unfold_pt1pt2[iter]->GetBinError(ibin + 1);
	  const float sim_err = h_pt1pt2_profile_unfold[iter]->GetBinError(ibin + 1);
	  const float sim_cont = h_pt1pt2_profile_unfold[iter]->GetBinContent(ibin + 1);

	  // Absolute (not fractional) change: the p_{T,1},p_{T,2} spectrum is
	  // steeply falling, so relative differences are dominated by the sparse
	  // high-p_{T} bins. Matches the AA note's unfolding description.
	  float err1 = fabs(prev_cont - bin_cont);

	  binbybin_unc += TMath::Power(err1,2);
	  stat_unc+= TMath::Power(sim_err, 2);
	  unfold_unc+= TMath::Power(unfold_err, 2);

	  // Fractional counterpart, for the comparison panel only.
	  binbybin_unc_rel += TMath::Power(err1/bin_cont, 2);
	  unfold_unc_rel += TMath::Power(unfold_err/bin_cont, 2);
	  if (sim_cont != 0) stat_unc_rel += TMath::Power(sim_err/sim_cont, 2);

	  // Normalised counterpart: absolute change of the unit-integral
	  // distributions, i.e. the note's definition.
	  const float bin_cont_norm = bin_cont/norm_unfold[iter];
	  const float prev_cont_norm = (iter > 0) ? prev_cont/norm_unfold[iter - 1]
	                                          : prev_cont/norm_data;
	  binbybin_unc_norm += TMath::Power(prev_cont_norm - bin_cont_norm, 2);
	  unfold_unc_norm += TMath::Power(unfold_err/norm_unfold[iter], 2);
	  stat_unc_norm += TMath::Power(sim_err/norm_sim[iter], 2);
	}

      std::cout << stat_unc << " + " << unfold_unc <<" + " << binbybin_unc << " = " <<  sqrt(stat_unc + unfold_unc + binbybin_unc) << std::endl;


      Double_t total_unc = sqrt(binbybin_unc + unfold_unc + stat_unc);
      Double_t total_unc_rel = sqrt(binbybin_unc_rel + unfold_unc_rel + stat_unc_rel);
      Double_t total_unc_norm = sqrt(binbybin_unc_norm + unfold_unc_norm + stat_unc_norm);

      stat_unc = sqrt(stat_unc);
      unfold_unc = sqrt(unfold_unc);
      binbybin_unc = sqrt(binbybin_unc);

      stat_unc_rel = sqrt(stat_unc_rel);
      unfold_unc_rel = sqrt(unfold_unc_rel);
      binbybin_unc_rel = sqrt(binbybin_unc_rel);

      stat_unc_norm = sqrt(stat_unc_norm);
      unfold_unc_norm = sqrt(unfold_unc_norm);
      binbybin_unc_norm = sqrt(binbybin_unc_norm);

      h_statistical_uncertainties->Fill(iter + 1, stat_unc);
      h_unfold_uncertainties->Fill(iter + 1, unfold_unc);
      h_binbybin_uncertainties->Fill(iter + 1, binbybin_unc);
      h_total_uncertainties->Fill(iter + 1, total_unc);

      h_statistical_uncertainties_rel->Fill(iter + 1, stat_unc_rel);
      h_unfold_uncertainties_rel->Fill(iter + 1, unfold_unc_rel);
      h_binbybin_uncertainties_rel->Fill(iter + 1, binbybin_unc_rel);
      h_total_uncertainties_rel->Fill(iter + 1, total_unc_rel);

      h_statistical_uncertainties_norm->Fill(iter + 1, stat_unc_norm);
      h_unfold_uncertainties_norm->Fill(iter + 1, unfold_unc_norm);
      h_binbybin_uncertainties_norm->Fill(iter + 1, binbybin_unc_norm);
      h_total_uncertainties_norm->Fill(iter + 1, total_unc_norm);
    }


  const std::string cent_label = Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]);

  // Absolute panel keeps the historical file name -- downstream scripts
  // (get_plots.sh, get_note_plots.sh) copy it by that name.
  drawIterationScan(h_statistical_uncertainties, h_unfold_uncertainties,
                    h_binbybin_uncertainties, h_total_uncertainties,
                    "c_unc", cent_label,
                    Form("%s/unfolding_plots/iteration_tune_AA_cent_%d_r%02d_%s",
                         rb.get_code_location().c_str(), centrality_bin, cone_size, sysname.c_str()));

  drawIterationScan(h_statistical_uncertainties_rel, h_unfold_uncertainties_rel,
                    h_binbybin_uncertainties_rel, h_total_uncertainties_rel,
                    "c_unc_rel", cent_label,
                    Form("%s/unfolding_plots/iteration_tune_relative_AA_cent_%d_r%02d_%s",
                         rb.get_code_location().c_str(), centrality_bin, cone_size, sysname.c_str()));

  drawIterationScan(h_statistical_uncertainties_norm, h_unfold_uncertainties_norm,
                    h_binbybin_uncertainties_norm, h_total_uncertainties_norm,
                    "c_unc_norm", cent_label,
                    Form("%s/unfolding_plots/iteration_tune_normalized_AA_cent_%d_r%02d_%s",
                         rb.get_code_location().c_str(), centrality_bin, cone_size, sysname.c_str()),
                    1.2); // fixed range, matching the AA note's iteration figure
  return;
}
