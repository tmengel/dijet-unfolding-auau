#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include <string>
#include <vector>

// x_J: truth (no event reweighting) vs truth (reweighted) vs unfolded data, for
// the exclusive and the inclusive measurement overlaid on one canvas (and, with
// draw_separate, repeated one measurement per canvas).
//
// Exclusive and inclusive are separate measurements with different pair
// definitions, so each carries its own truth: the inclusive series is never
// compared against the exclusive truth.
//
// The overlay's ratio panel is inclusive over exclusive, taken series by
// series -- truth over truth, reweighted truth over reweighted truth, unfolded
// over unfolded -- so the two are only ever divided by their own counterpart.
// Both truths are drawn as histogram lines (no points) in both panels; only the
// unfolded data, which is the thing carrying an uncertainty worth showing, is
// drawn as points.
//
// "truth (no reweight)"  = h_truth_flat_pt1pt2 from the PRIMER1 pass, which runs
//                          before centrality/vertex/sumeT reweight histograms
//                          exist and so carries no event reweighting at all
//                          (see the `if ( !PRIMER1 || OVERRIDE_EVENT_WEIGHT )`
//                          block in createResponse_exclusive_v2_AA.cxx).
// "truth (reweight)"     = the same histogram from the final pass, with the
//                          centrality/vertex/sumeT event weights applied.
// "unfolded"             = h_flat_unfold_pt1pt2_<niter> from the final pass.
//
// The exclusive series comes from the "nominal" files (createResponse_exclusive_v2_AA.cxx)
// and the inclusive series from the "INCLUSIVE" files
// (createResponse_exclusive_v2_AA_inclusive.cxx, run_inclusive_sys_AA_exclusive.sh).

namespace
{
  // Colour separates the two measurements, marker separates the three series
  // within each measurement.
  const int color_exclusive = kAzure - 6;
  const int color_inclusive = kRed + 1;

  const float marker_truth_raw = 25;  // open square:   truth, no event reweighting
  const float marker_truth_rw = 21;   // filled square: truth, event reweighted
  const float marker_unfold = 20;     // filled circle: unfolded data

  const float msize = 0.9;
  const float lwidth = 2;  // the truths are lines now, so give them some weight

  // The overlay's ratio panel (inclusive/exclusive) drops colour-by-measurement
  // -- there's only one line there per series -- and uses colour-by-series
  // instead so the three ratios can be told apart.
  const int color_cross_truth_raw = kGray + 2;
  const int color_cross_truth_rw = kBlack;
  const int color_cross_unfold = kBlue + 2;

  // One measurement definition's worth of finalised, per-leading-pT-range x_J.
  struct xj_set
  {
    bool ok = false;
    std::vector<TH1D *> truth_raw;
    std::vector<TH1D *> truth_rw;
    std::vector<TH1D *> unfold;
  };

  TH1D *fetch(TFile *f, const char *hname, const char *newname)
  {
    if (!f) return nullptr;
    TH1D *h = (TH1D *) f->Get(hname);
    if (!h)
      {
        std::cerr << "Missing " << hname << " in " << f->GetName() << std::endl;
        return nullptr;
      }
    h = (TH1D *) h->Clone(newname);
    h->SetDirectory(0);
    return h;
  }

  // flat pt1pt2 -> symmetric 2D -> per-range x_J -> normalised -> trimmed.
  std::vector<TH1D *> flat_to_xj_ranges(TH1D *h_flat, const std::string &name,
                                        const int nbins, const float *ipt_bins, const float *ixj_bins,
                                        const int mbins, const int *measure_bins,
                                        const int measure_subleading_bin, const float first_xj)
  {
    std::vector<TH1D *> out;
    TH2D *h_pt1pt2 = new TH2D(Form("h2_%s", name.c_str()), ";#it{p}_{T,1};#it{p}_{T,2}",
                              nbins, ipt_bins, nbins, ipt_bins);
    histo_opps::make_sym_pt1pt2(h_flat, h_pt1pt2, nbins);

    for (int irange = 0; irange < mbins; irange++)
      {
        TH1D *h_xj = new TH1D(Form("h_xj_%s_range_%d", name.c_str(), irange), ";x_{J};", nbins, ixj_bins);
        histo_opps::project_xj(h_pt1pt2, h_xj, nbins,
                               measure_bins[irange], measure_bins[irange + 1],
                               measure_subleading_bin, nbins - 2);
        histo_opps::normalize_histo(h_xj, nbins);

        TH1D *h_final = new TH1D(Form("h_final_xj_%s_range_%d", name.c_str(), irange), ";x_{J};", nbins, ixj_bins);
        histo_opps::finalize_xj(h_xj, h_final, nbins, first_xj);
        out.push_back(h_final);
      }
    return out;
  }

  xj_set build_xj_set(const std::string &tag, const std::string &sys_name,
                      const std::string &system_string, const std::string &code_location,
                      const int cone_size, const int niter,
                      const int nbins, const float *ipt_bins, const float *ixj_bins,
                      const int mbins, const int *measure_bins,
                      const int measure_subleading_bin, const float first_xj)
  {
    xj_set out;

    const TString finalpath = Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_%s.root",
                                   code_location.c_str(), system_string.c_str(), cone_size, sys_name.c_str());
    const TString primerpath = Form("%s/unfolding_hists/unfolding_hists_%s_r%02d_PRIMER1_%s.root",
                                    code_location.c_str(), system_string.c_str(), cone_size, sys_name.c_str());

    TFile *ffinal = TFile::Open(finalpath.Data(), "READ");
    if (!ffinal || ffinal->IsZombie())
      {
        std::cerr << "Missing " << finalpath << std::endl;
        return out;
      }
    TFile *fprimer = TFile::Open(primerpath.Data(), "READ");
    if (!fprimer || fprimer->IsZombie())
      {
        std::cerr << "Missing " << primerpath << std::endl;
        return out;
      }

    TH1D *h_flat_truth_rw = fetch(ffinal, "h_truth_flat_pt1pt2", Form("h_%s_truth_rw_flat", tag.c_str()));
    TH1D *h_flat_unfold = fetch(ffinal, Form("h_flat_unfold_pt1pt2_%d", niter), Form("h_%s_unfold_flat", tag.c_str()));
    TH1D *h_flat_truth_raw = fetch(fprimer, "h_truth_flat_pt1pt2", Form("h_%s_truth_raw_flat", tag.c_str()));

    ffinal->Close();
    fprimer->Close();

    if (!h_flat_truth_rw || !h_flat_unfold || !h_flat_truth_raw) return out;

    out.truth_raw = flat_to_xj_ranges(h_flat_truth_raw, tag + "_truth_raw", nbins, ipt_bins, ixj_bins,
                                      mbins, measure_bins, measure_subleading_bin, first_xj);
    out.truth_rw = flat_to_xj_ranges(h_flat_truth_rw, tag + "_truth_rw", nbins, ipt_bins, ixj_bins,
                                     mbins, measure_bins, measure_subleading_bin, first_xj);
    out.unfold = flat_to_xj_ranges(h_flat_unfold, tag + "_unfold", nbins, ipt_bins, ixj_bins,
                                   mbins, measure_bins, measure_subleading_bin, first_xj);
    out.ok = true;
    return out;
  }
}

void drawTruthUnfoldCompare_AA(
  const int cone_size = 3, 
  const int centrality_bin = 0, 
  const int niter = 1,                             
  const bool draw_separate = false,
  const std::string configfile = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());

  const bool ispp = (centrality_bin < 0);
  const std::string system_string = (ispp ? "pp" : "AA_cent_" + std::to_string(centrality_bin));
  const std::string code_location = rb.get_code_location();

  const int cent_bins = rb.get_number_centrality_bins();
  float icentrality_bins[cent_bins + 1];
  rb.get_centrality_bins(icentrality_bins);

  const float first_xj = rb.get_first_xj();
  const int nbins = rb.get_nbins();
  const std::string dphi_string = rb.get_dphi_string();

  float ipt_bins[nbins + 1];
  float ixj_bins[nbins + 1];
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  // Low edge of the first bin finalize_xj actually keeps. Drawing from here
  // rather than from first_xj stops the histogram outlines running along zero
  // across the empty bins below the measured range.
  float xj_plot_min = ixj_bins[0];
  for (int i = 0; i < nbins; i++)
    {
      if (ixj_bins[i] >= first_xj) { xj_plot_min = ixj_bins[i]; break; }
    }

  const int measure_subleading_bin = rb.get_measure_subleading_bin();
  const int mbins = rb.get_measure_bins();
  int measure_bins[10] = {0};
  for (int ir = 0; ir < mbins + 1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
    }

  // Exclusive: the nominal (createResponse_exclusive_v2_AA.cxx) chain.
  // Inclusive: the INCLUSIVE (createResponse_exclusive_v2_AA_inclusive.cxx) chain.
  // Each carries its own truth pair -- the inclusive panel is never drawn
  // against the exclusive truth.
  const std::string tags[2] = {"exclusive", "inclusive"};
  const std::string sys_names[2] = {"nominal", "INCLUSIVE"};
  const std::string labels[2] = {"Exclusive", "Inclusive"};
  const int set_colors[2] = {color_exclusive, color_inclusive};

  xj_set sets[2];
  for (int iset = 0; iset < 2; iset++)
    {
      sets[iset] = build_xj_set(tags[iset], sys_names[iset], system_string, code_location,
                                cone_size, niter, nbins, ipt_bins, ixj_bins,
                                mbins, measure_bins, measure_subleading_bin, first_xj);
      if (!sets[iset].ok)
        std::cerr << "Missing inputs for " << tags[iset] << " -- it will be left off the plots." << std::endl;
    }

  // Style: colour separates the two measurements, line/marker separates the
  // three series within each of them. Both truths are lines only -- marker
  // style 1 so nothing stray shows up if they are ever drawn with points.
  for (int iset = 0; iset < 2; iset++)
    {
      if (!sets[iset].ok) continue;
      for (int irange = 0; irange < mbins; irange++)
        {
          dlutility::SetLineAtt(sets[iset].truth_raw[irange], set_colors[iset], lwidth, 2);
          sets[iset].truth_raw[irange]->SetMarkerStyle(1);
          dlutility::SetLineAtt(sets[iset].truth_rw[irange], set_colors[iset], lwidth, 1);
          sets[iset].truth_rw[irange]->SetMarkerStyle(1);
          dlutility::SetLineAtt(sets[iset].unfold[irange], set_colors[iset], lwidth, 1);
          dlutility::SetMarkerAtt(sets[iset].unfold[irange], set_colors[iset], msize, marker_unfold);
        }
    }

  // Ratios used by the single-measurement canvases: against that set's OWN
  // reweighted truth, the prior its own unfolding used.
  TH1D *h_raw_ratio[2][10] = {{nullptr}};
  TH1D *h_unfold_ratio[2][10] = {{nullptr}};
  for (int iset = 0; iset < 2; iset++)
    {
      if (!sets[iset].ok) continue;
      for (int irange = 0; irange < mbins; irange++)
        {
          h_raw_ratio[iset][irange] = (TH1D *) sets[iset].truth_raw[irange]->Clone(Form("h_%s_truth_raw_ratio_range_%d", tags[iset].c_str(), irange));
          h_raw_ratio[iset][irange]->Divide(sets[iset].truth_rw[irange]);
          h_unfold_ratio[iset][irange] = (TH1D *) sets[iset].unfold[irange]->Clone(Form("h_%s_unfold_ratio_range_%d", tags[iset].c_str(), irange));
          h_unfold_ratio[iset][irange]->Divide(sets[iset].truth_rw[irange]);
        }
    }

  // Ratios used by the overlay canvas: inclusive over exclusive, series by
  // series. Truth is compared with truth and unfolded with unfolded -- the two
  // measurements are only ever divided by their own counterpart, never mixed.
  const bool have_both = (sets[0].ok && sets[1].ok);
  TH1D *h_cross_truth_raw[10] = {nullptr};
  TH1D *h_cross_truth_rw[10] = {nullptr};
  TH1D *h_cross_unfold[10] = {nullptr};
  if (have_both)
    {
      for (int irange = 0; irange < mbins; irange++)
        {
          h_cross_truth_raw[irange] = (TH1D *) sets[1].truth_raw[irange]->Clone(Form("h_cross_truth_raw_range_%d", irange));
          h_cross_truth_raw[irange]->Divide(sets[0].truth_raw[irange]);
          h_cross_truth_rw[irange] = (TH1D *) sets[1].truth_rw[irange]->Clone(Form("h_cross_truth_rw_range_%d", irange));
          h_cross_truth_rw[irange]->Divide(sets[0].truth_rw[irange]);
          h_cross_unfold[irange] = (TH1D *) sets[1].unfold[irange]->Clone(Form("h_cross_unfold_range_%d", irange));
          h_cross_unfold[irange]->Divide(sets[0].unfold[irange]);

          // Same vocabulary as the top pad -- dashed line = truth, solid line =
          // truth w/ reweight, points = unfolded -- but colour now separates
          // the three series instead of the two measurements, since a ratio
          // of the two measurements is neither of them.
          dlutility::SetLineAtt(h_cross_truth_raw[irange], color_cross_truth_raw, lwidth, 2);
          h_cross_truth_raw[irange]->SetMarkerStyle(1);
          dlutility::SetLineAtt(h_cross_truth_rw[irange], color_cross_truth_rw, lwidth, 1);
          h_cross_truth_rw[irange]->SetMarkerStyle(1);
          dlutility::SetLineAtt(h_cross_unfold[irange], color_cross_unfold, lwidth, 1);
          dlutility::SetMarkerAtt(h_cross_unfold[irange], color_cross_unfold, msize, marker_unfold);
        }
    }
  else
    {
      std::cerr << "Only one measurement available -- the overlay ratio panel will be empty." << std::endl;
    }

  // Draws the standard sPHENIX/kinematics caption block for one range.
  auto draw_captions = [&](const int irange)
  {
    if (!ispp) { dlutility::DrawSPHENIX(0.22, 0.84); }
    else { dlutility::DrawSPHENIXpp(0.22, 0.84); }
    dlutility::drawText(Form("anti-k_{t} R = %0.1f", 0.1 * cone_size), 0.22, 0.74);
    dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange + 1]]), 0.22, 0.69);
    dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
    dlutility::drawText(Form("#Delta#phi #geq %s", dphi_string.c_str()), 0.22, 0.59);
    if (!ispp) dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin + 1]), 0.22, 0.54);
    dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.22, 0.49);
  };

  auto draw_unity_line = [&]()
  {
    TLine *line = new TLine(xj_plot_min, 1, 1, 1);
    line->SetLineStyle(4);
    line->SetLineColor(kRed + 3);
    line->SetLineWidth(2);
    line->Draw("same");
  };

  // ---------------------------------------------------------------------
  // Both measurements on one canvas.
  // ---------------------------------------------------------------------
  TCanvas *cboth = new TCanvas("cboth", "cboth", 700, 700);
  dlutility::ratioPanelCanvas(cboth);

  for (int irange = 1; irange < 2; irange++)
    {
      cboth->cd(1);

      // Frame off whichever set is present, then overlay the rest.
      TH1D *h_frame = nullptr;
      for (int iset = 0; iset < 2 && !h_frame; iset++)
        if (sets[iset].ok) h_frame = sets[iset].truth_rw[irange];
      if (!h_frame) continue;

      dlutility::SetFont(h_frame, 42, 0.04);
      h_frame->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
      h_frame->SetMaximum(4.5);
      h_frame->SetMinimum(0);
      h_frame->GetXaxis()->SetRangeUser(xj_plot_min, 1.0);
      h_frame->Draw("hist");

      for (int iset = 0; iset < 2; iset++)
        {
          if (!sets[iset].ok) continue;
          sets[iset].truth_raw[irange]->Draw("same hist");
          sets[iset].truth_rw[irange]->Draw("same hist");
          sets[iset].unfold[irange]->Draw("same p");
        }

      draw_captions(irange);

      TLegend *leg = new TLegend(0.55, 0.56, 0.94, 0.86);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.032);
      leg->SetTextFont(42);
      for (int iset = 0; iset < 2; iset++)
        {
          if (!sets[iset].ok) continue;
          leg->AddEntry(sets[iset].truth_raw[irange], Form("%s Truth", labels[iset].c_str()), "l");
          leg->AddEntry(sets[iset].truth_rw[irange], Form("%s Truth w/ Reweight", labels[iset].c_str()), "l");
          leg->AddEntry(sets[iset].unfold[irange], Form("%s Unfolded Data", labels[iset].c_str()), "pl");
        }
      leg->Draw("same");

      cboth->cd(2);

      if (have_both)
        {
          TH1D *h_ratio_frame = h_cross_truth_rw[irange];
          h_ratio_frame->SetTitle(";x_{J}; Inclusive / Exclusive");
          dlutility::SetFont(h_ratio_frame, 42, 0.06);
          h_ratio_frame->SetMaximum(1.5);
          h_ratio_frame->SetMinimum(0.5);
          h_ratio_frame->GetXaxis()->SetRangeUser(xj_plot_min, 1.0);

          h_ratio_frame->Draw("hist");
          h_cross_truth_raw[irange]->Draw("hist same");
          h_cross_unfold[irange]->Draw("p same");
          draw_unity_line();
        }

      cboth->SaveAs(Form("%s/unfolding_plots/truth_unfold_compare_%s_r%02d_range_%d_iter_%d.png",
                         code_location.c_str(), system_string.c_str(), cone_size, irange, niter));
      cboth->SaveAs(Form("%s/unfolding_plots/truth_unfold_compare_%s_r%02d_range_%d_iter_%d.pdf",
                         code_location.c_str(), system_string.c_str(), cone_size, irange, niter));
    }

  // ---------------------------------------------------------------------
  // The same content split one measurement per canvas, for when the overlay
  // is too busy.
  // ---------------------------------------------------------------------
  if (!draw_separate) return;

  for (int iset = 0; iset < 2; iset++)
    {
      if (!sets[iset].ok) continue;

      TCanvas *c = new TCanvas(Form("c_%s", tags[iset].c_str()), tags[iset].c_str(), 600, 700);
      dlutility::ratioPanelCanvas(c);

      for (int irange = 0; irange < mbins; irange++)
        {
          c->cd(1);

          TH1D *h_rw = sets[iset].truth_rw[irange];
          dlutility::SetFont(h_rw, 42, 0.04);
          h_rw->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
          h_rw->SetMaximum(4.5);
          h_rw->SetMinimum(0);

          h_rw->GetXaxis()->SetRangeUser(xj_plot_min, 1.0);
          h_rw->Draw("hist");
          sets[iset].truth_raw[irange]->Draw("same hist");
          sets[iset].unfold[irange]->Draw("same p");

          draw_captions(irange);
          dlutility::drawText(Form("%s pairs", labels[iset].c_str()), 0.22, 0.44);

          TLegend *leg = new TLegend(0.58, 0.66, 0.93, 0.84);
          leg->SetLineWidth(0);
          leg->SetTextSize(0.038);
          leg->SetTextFont(42);
          leg->AddEntry(sets[iset].truth_raw[irange], "Truth", "l");
          leg->AddEntry(h_rw, "Truth w/ Reweight", "l");
          leg->AddEntry(sets[iset].unfold[irange], "Unfolded Data", "pl");
          leg->Draw("same");

          c->cd(2);

          h_raw_ratio[iset][irange]->SetTitle(";x_{J}; Ratio to Truth w/ Reweight");
          dlutility::SetFont(h_raw_ratio[iset][irange], 42, 0.06);
          h_raw_ratio[iset][irange]->SetMaximum(1.5);
          h_raw_ratio[iset][irange]->SetMinimum(0.5);
          h_raw_ratio[iset][irange]->GetXaxis()->SetRangeUser(xj_plot_min, 1.0);
          h_raw_ratio[iset][irange]->Draw("hist");
          h_unfold_ratio[iset][irange]->Draw("p same");
          draw_unity_line();

          c->SaveAs(Form("%s/unfolding_plots/truth_unfold_compare_%s_%s_r%02d_range_%d_iter_%d.png",
                         code_location.c_str(), tags[iset].c_str(), system_string.c_str(), cone_size, irange, niter));
          c->SaveAs(Form("%s/unfolding_plots/truth_unfold_compare_%s_%s_r%02d_range_%d_iter_%d.pdf",
                         code_location.c_str(), tags[iset].c_str(), system_string.c_str(), cone_size, irange, niter));
        }
    }

  return;
}
