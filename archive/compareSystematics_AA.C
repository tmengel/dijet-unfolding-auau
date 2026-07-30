// Compare the AA systematic-variation plot between the current and reference
// versions.  The reference file contains only total systematics, so its two
// total-uncertainty curves are drawn as dotted black lines.
//
// Example: root -l -q 'compareSystematics_AA.C(1,4)'

#include <algorithm>
#include <cmath>

#include "TCanvas.h"
#include "TError.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

#include "histo_opps.h"

namespace {
TH1D *getHist(TFile &file, const TString &name)
{
  TH1D *hist = dynamic_cast<TH1D *>(file.Get(name));
  if (!hist) Error("compareSystematics_AA", "Missing histogram %s in %s", name.Data(), file.GetName());
  return hist;
}

TH1D *copyHist(const TH1D *source, const TString &name)
{
  TH1D *copy = static_cast<TH1D *>(source->Clone(name));
  copy->SetDirectory(nullptr);
  return copy;
}

// The stored total uncertainties are absolute deviations.  Convert them to
// the Variation/Nominal representation used in the systematic plot.
TH1D *totalVariation(const TH1D *source, bool positive, const TString &name)
{
  TH1D *result = copyHist(source, name);
  for (int bin = 1; bin <= result->GetNbinsX(); ++bin)
    result->SetBinContent(bin, 1. + (positive ? 1. : -1.) * source->GetBinContent(bin));
  return result;
}

void makeEnvelope(const TH1D *up, const TH1D *down, TH1D *positive, TH1D *negative)
{
  for (int bin = 1; bin <= up->GetNbinsX(); ++bin) {
    double maximum = up->GetBinContent(bin) - 1.;
    double minimum = down->GetBinContent(bin) - 1.;
    if (maximum < 0. && minimum > 0.) std::swap(maximum, minimum);
    else if (maximum > 0. && minimum > 0.) { maximum = std::max(maximum, minimum); minimum = 0.; }
    else if (maximum < 0. && minimum < 0.) { minimum = std::min(maximum, minimum); maximum = 0.; }
    positive->SetBinContent(bin, 1. + maximum);
    negative->SetBinContent(bin, 1. - std::abs(minimum));
  }
}

void styleLine(TH1D *hist, int color, int style = 1)
{
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->SetLineStyle(style);
  hist->SetMarkerColor(color);
}

void drawPair(TH1D *positive, TH1D *negative, int color)
{
  styleLine(positive, color);
  styleLine(negative, color);
  positive->Draw("HIST SAME");
  negative->Draw("HIST SAME");
}

TH1D *graphAsHistogram(const TGraph *graph, const TH1D *templateHist, const TString &name)
{
  TH1D *hist = copyHist(templateHist, name);
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
    hist->SetBinContent(bin, graph->Eval(hist->GetXaxis()->GetBinCenter(bin)));
  return hist;
}

void saveIndividualSystematic(TH1D *currentUpper, TH1D *currentLower,
                              TH1D *oldUpper, TH1D *oldLower,
                              const char *label, const char *fileLabel,
                              int color, int range, int iter)
{
  TCanvas *canvas = new TCanvas(Form("c_%s", fileLabel), label, 700, 600);
  canvas->SetLeftMargin(0.15); canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.14); canvas->SetTopMargin(0.05); canvas->SetTicks(1, 1);
  currentUpper->SetTitle(Form("%s; x_{J}; Variation / Nominal", label));
  currentUpper->SetMinimum(0.); currentUpper->SetMaximum(2.);
  currentUpper->GetXaxis()->SetRangeUser(0.3, 1.0);
  currentUpper->GetXaxis()->SetTitleSize(0.055); currentUpper->GetYaxis()->SetTitleSize(0.055);
  currentUpper->GetXaxis()->SetLabelSize(0.045); currentUpper->GetYaxis()->SetLabelSize(0.045);
  currentUpper->SetLineColor(color); currentUpper->SetMarkerColor(color); currentUpper->SetMarkerStyle(20);
  currentLower->SetLineColor(color); currentLower->SetMarkerColor(color); currentLower->SetMarkerStyle(24);
  oldUpper->SetLineColor(color); oldUpper->SetLineStyle(3); oldUpper->SetLineWidth(2);
  oldLower->SetLineColor(color); oldLower->SetLineStyle(3); oldLower->SetLineWidth(2);
  currentUpper->Draw("HIST"); currentLower->Draw("HIST SAME");
  oldUpper->Draw("HIST SAME"); oldLower->Draw("HIST SAME");
  TLine unity(0.3, 1., 1., 1.); unity.SetLineStyle(2); unity.SetLineColor(kGray + 2); unity.Draw();
  TLegend *legend = new TLegend(0.47, 0.67, 0.94, 0.89);
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.035);
  legend->AddEntry(currentUpper, "New (solid)", "l");
  legend->AddEntry(oldUpper, "Old (dotted)", "l");
  legend->AddEntry(currentUpper, "Upper variation", "p");
  legend->AddEntry(currentLower, "Lower variation", "p");
  legend->Draw();
  canvas->SaveAs(Form("compare_systematic_%s_AA_cent_0_r03_range_%d_iter_%d.pdf", fileLabel, range, iter));
}

// Recreate the old x_J distribution from its saved flattened unfolded
// spectrum.  The old configuration is embedded as TEnv in its uncertainty
// file, so this does not depend on the current binning configuration.
TH1D *oldFinalXj(TFile &file, const TEnv &config, int range, int iter, const TString &name)
{
  const int nbins = config.GetValue("nbins", 0);
  const double minimum = config.GetValue("minimum", 0.);
  const double fixed = config.GetValue("fixed", 0.);
  const int bbins = config.GetValue("bbins", 0);
  if (nbins <= 0 || minimum <= 0. || fixed <= 0. || bbins <= 0) return nullptr;
  TH1D *flat = getHist(file, Form("h_flat_unfold_pt1pt2_%d", iter));
  if (!flat) return nullptr;
  const double alpha = std::pow(fixed / minimum, 1. / bbins);
  double pt[30], xj[30];
  for (int bin = 0; bin <= nbins; ++bin) pt[bin] = minimum * std::pow(alpha, bin);
  const double maximum = pt[nbins];
  for (int bin = 0; bin <= nbins; ++bin) xj[bin] = minimum / maximum * std::pow(alpha, bin);
  int leadingStart[4] = {0};
  for (int region = 0; region < 4; ++region) {
    const double threshold = config.GetValue(Form("pt1_%d", region), 1.);
    while (leadingStart[region] < nbins && pt[leadingStart[region]] < threshold) ++leadingStart[region];
  }
  const int subleadingStart = config.GetValue("measure_subleading_bin_buffer", 0);
  if (range < 0 || range >= config.GetValue("measure_bins", 0)) return nullptr;
  TH2D pt1pt2(Form("%s_pt1pt2", name.Data()), "", nbins, pt, nbins, pt);
  TH1D projected(Form("%s_projected", name.Data()), "", nbins, xj);
  TH1D *final = new TH1D(name, "", nbins, xj);
  histo_opps::make_sym_pt1pt2(flat, &pt1pt2, nbins);
  histo_opps::project_xj(&pt1pt2, &projected, nbins, leadingStart[range], leadingStart[range + 1],
                         subleadingStart, nbins - 2);
  histo_opps::normalize_histo(&projected, nbins);
  histo_opps::finalize_xj(&projected, final, nbins, config.GetValue("first_xj", 0.));
  final->SetDirectory(nullptr);
  return final;
}

TH1D *oldRatio(TFile &variation, TFile &nominal, const TEnv &config,
               int range, int iter, const TString &name)
{
  TH1D *varied = oldFinalXj(variation, config, range, iter, name + "_varied");
  TH1D *baseline = oldFinalXj(nominal, config, range, iter, name + "_baseline");
  if (!varied || !baseline) return nullptr;
  varied->Divide(baseline);
  varied->SetName(name);
  return varied;
}
}

void compareSystematics_AA(
    int range = 0,
    int iter = 0,
    const char *currentFile = "/home/tmengel/PPG14/dijet-unfolding-final/uncertainties/systematics_AA_cent_0_r03.root",
    const char *referenceFile = "/home/tmengel/PPG14/sphenix_digitized_root_bundle/"
                                "sphenix_xj_systematics_digitized.root")
{
  gStyle->SetOptStat(0);
  TFile current(currentFile, "READ");
  TFile reference(referenceFile, "READ");
  if (current.IsZombie() || reference.IsZombie()) {
    Error("compareSystematics_AA", "Could not open one or both input files.");
    return;
  }

  const TString totalName = Form("h_total_sys_range_%d_iter_%d", range, iter);
  const TString totalNegName = Form("h_total_sys_neg_range_%d_iter_%d", range, iter);
  TH1D *curTotalUp = getHist(current, totalName);
  TH1D *curTotalDown = getHist(current, totalNegName);
  TGraph *oldUp = dynamic_cast<TGraph *>(reference.Get("total_systematics_upper"));
  TGraph *oldDown = dynamic_cast<TGraph *>(reference.Get("total_systematics_lower"));
  TH1D *jesUp = getHist(current, Form("h_sys_posJES_range_%d_iter_%d", range, iter));
  TH1D *jesDown = getHist(current, Form("h_sys_negJES_range_%d_iter_%d", range, iter));
  TH1D *jerUp = getHist(current, Form("h_sys_posJER_range_%d_iter_%d", range, iter));
  TH1D *jerDown = getHist(current, Form("h_sys_negJER_range_%d_iter_%d", range, iter));
  TH1D *comb = getHist(current, Form("h_sys_COMB_range_%d_iter_%d", range, iter));
  TH1D *inclusive = getHist(current, Form("h_sys_INCLUSIVE_range_%d_iter_%d", range, iter));
  TH1D *prior = getHist(current, Form("h_sys_PRIOR_range_%d_iter_%d", range, iter));
  if (!curTotalUp || !curTotalDown || !oldUp || !oldDown || !jesUp || !jesDown ||
      !jerUp || !jerDown || !comb || !inclusive || !prior) return;

  TH1D *totalUp = totalVariation(curTotalUp, true, "totalUp");
  TH1D *totalDown = totalVariation(curTotalDown, false, "totalDown");
  TH1D *jesPositive = copyHist(jesUp, "jesPositive");
  TH1D *jesNegative = copyHist(jesUp, "jesNegative");
  TH1D *jerPositive = copyHist(jerUp, "jerPositive");
  TH1D *jerNegative = copyHist(jerUp, "jerNegative");
  makeEnvelope(jesUp, jesDown, jesPositive, jesNegative);
  makeEnvelope(jerUp, jerDown, jerPositive, jerNegative);
  TH1D *combNegative = copyHist(comb, "combNegative");
  TH1D *inclusiveNegative = copyHist(inclusive, "inclusiveNegative");
  TH1D *priorNegative = copyHist(prior, "priorNegative");
  for (int bin = 1; bin <= comb->GetNbinsX(); ++bin) {
    combNegative->SetBinContent(bin, 2. - comb->GetBinContent(bin));
    inclusiveNegative->SetBinContent(bin, 2. - inclusive->GetBinContent(bin));
    priorNegative->SetBinContent(bin, 2. - prior->GetBinContent(bin));
  }

  // Digitized curves provide every old systematic directly.  Old ZYAM is
  // represented by the digitized UE-background-subtraction pair.
  TH1D *oldJESPositive = copyHist(getHist(reference, "JES_upper"), "oldJESPositive");
  TH1D *oldJESNegative = copyHist(getHist(reference, "JES_lower"), "oldJESNegative");
  TH1D *oldJERPositive = copyHist(getHist(reference, "JER_upper"), "oldJERPositive");
  TH1D *oldJERNegative = copyHist(getHist(reference, "JER_lower"), "oldJERNegative");
  TH1D *oldComb = copyHist(getHist(reference, "UE_background_subtraction_upper"), "oldComb");
  TH1D *oldCombNegative = copyHist(getHist(reference, "UE_background_subtraction_lower"), "oldCombNegative");
  TH1D *oldInclusive = copyHist(getHist(reference, "inclusive_upper"), "oldInclusive");
  TH1D *oldInclusiveNegative = copyHist(getHist(reference, "inclusive_lower"), "oldInclusiveNegative");
  TH1D *oldPriorRatio = copyHist(getHist(reference, "unfolding_upper"), "oldPrior");
  TH1D *oldPriorNegative = copyHist(getHist(reference, "unfolding_lower"), "oldPriorNegative");
  if (!oldJESPositive || !oldJESNegative || !oldJERPositive || !oldJERNegative || !oldComb ||
      !oldCombNegative || !oldInclusive || !oldInclusiveNegative || !oldPriorRatio || !oldPriorNegative) return;

  TCanvas *canvas = new TCanvas("c_compareSystematics", "Systematics comparison", 1120, 650);
  TPad *plot = new TPad("plot", "", 0., 0., 0.62, 1.);
  TPad *info = new TPad("info", "", 0.62, 0., 1., 1.);
  plot->SetLeftMargin(0.15); plot->SetRightMargin(0.02); plot->SetBottomMargin(0.14); plot->SetTopMargin(0.04);
  plot->SetTicks(1, 1); plot->Draw(); info->Draw();

  plot->cd();
  totalUp->SetTitle(";x_{J};Variation / Nominal");
  totalUp->SetMinimum(0.0); totalUp->SetMaximum(2.0);
  // The old analysis begins at x_J = 0.3; use the common plotted domain.
  totalUp->GetXaxis()->SetRangeUser(0.3, 1.0);
  totalUp->GetXaxis()->SetTitleSize(0.06); totalUp->GetYaxis()->SetTitleSize(0.06);
  totalUp->GetXaxis()->SetLabelSize(0.05); totalUp->GetYaxis()->SetLabelSize(0.05);
  totalUp->GetYaxis()->SetTitleOffset(1.05);
  totalUp->SetMarkerStyle(20); totalUp->SetMarkerSize(0.8); totalUp->SetMarkerColor(kBlack);
  totalDown->SetMarkerStyle(20); totalDown->SetMarkerSize(0.8); totalDown->SetMarkerColor(kBlack);
  totalUp->Draw("P"); totalDown->Draw("P SAME");
  drawPair(jesPositive, jesNegative, kCyan + 1);
  drawPair(jerPositive, jerNegative, kMagenta + 1);
  drawPair(comb, combNegative, kOrange + 7);
  drawPair(inclusive, inclusiveNegative, kRed + 1);
  drawPair(prior, priorNegative, kRed - 2);
  // Dotted old components.  Old ZYAM is the predecessor of current COMB.
  drawPair(oldJESPositive, oldJESNegative, kCyan + 1);
  drawPair(oldJERPositive, oldJERNegative, kMagenta + 1);
  drawPair(oldComb, oldCombNegative, kOrange + 7);
  drawPair(oldInclusive, oldInclusiveNegative, kRed + 1);
  drawPair(oldPriorRatio, oldPriorNegative, kRed - 2);
  for (TH1D *hist : {oldJESPositive, oldJESNegative, oldJERPositive, oldJERNegative,
                      oldComb, oldCombNegative, oldInclusive, oldInclusiveNegative,
                      oldPriorRatio, oldPriorNegative}) hist->SetLineStyle(3);
  oldJESPositive->Draw("HIST SAME"); oldJESNegative->Draw("HIST SAME");
  oldJERPositive->Draw("HIST SAME"); oldJERNegative->Draw("HIST SAME");
  oldComb->Draw("HIST SAME"); oldCombNegative->Draw("HIST SAME");
  oldInclusive->Draw("HIST SAME"); oldInclusiveNegative->Draw("HIST SAME");
  oldPriorRatio->Draw("HIST SAME"); oldPriorNegative->Draw("HIST SAME");
  oldUp->SetLineColor(kBlack); oldUp->SetLineWidth(2); oldUp->SetLineStyle(3);
  oldDown->SetLineColor(kBlack); oldDown->SetLineWidth(2); oldDown->SetLineStyle(3);
  oldUp->Draw("L SAME"); oldDown->Draw("L SAME");
  totalUp->Draw("P SAME"); totalDown->Draw("P SAME");

  info->cd();
  TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.065);
  text.DrawLatex(0.08, 0.93, "#bf{sPHENIX}  #it{Preliminary}");
  text.DrawLatex(0.08, 0.86, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
  text.DrawLatex(0.08, 0.77, "anti-k_{t}  R = 0.3");
  text.DrawLatex(0.08, 0.70, Form("p_{T,1} range index = %d", range));
  text.DrawLatex(0.08, 0.63, "p_{T,2} > 10.1 GeV");
  text.DrawLatex(0.08, 0.56, "#Delta#phi #geq 7#pi/8");
  text.DrawLatex(0.08, 0.49, "0 - 10 %");
  TLegend *legend = new TLegend(0.08, 0.08, 0.96, 0.46);
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextFont(42); legend->SetTextSize(0.055);
  legend->AddEntry(totalUp, "Total current", "p");
  legend->AddEntry(oldUp, "Total reference (dotted)", "l");
  legend->AddEntry(jesPositive, "JES systematics", "l");
  legend->AddEntry(jerPositive, "JER systematics", "l");
  legend->AddEntry(comb, "UE bkg. sub.", "l");
  legend->AddEntry(inclusive, "Inclusive", "l");
  legend->AddEntry(prior, "Unfolding", "l");
  legend->Draw();
  canvas->SaveAs(Form("compare_systematics_AA_cent_0_r03_range_%d_iter_%d.pdf", range, iter));

  TH1D *oldTotalUp = graphAsHistogram(oldUp, totalUp, "oldTotalUp");
  TH1D *oldTotalDown = graphAsHistogram(oldDown, totalDown, "oldTotalDown");
  saveIndividualSystematic(totalUp, totalDown, oldTotalUp, oldTotalDown, "Total systematics", "total", kBlack, range, iter);
  saveIndividualSystematic(jesPositive, jesNegative, oldJESPositive, oldJESNegative, "JES systematics", "jes", kCyan + 1, range, iter);
  saveIndividualSystematic(jerPositive, jerNegative, oldJERPositive, oldJERNegative, "JER systematics", "jer", kMagenta + 1, range, iter);
  saveIndividualSystematic(comb, combNegative, oldComb, oldCombNegative, "UE background subtraction", "ue_background", kOrange + 7, range, iter);
  saveIndividualSystematic(inclusive, inclusiveNegative, oldInclusive, oldInclusiveNegative, "Inclusive systematic", "inclusive", kRed + 1, range, iter);
  saveIndividualSystematic(prior, priorNegative, oldPriorRatio, oldPriorNegative, "Unfolding systematic", "unfolding", kRed - 2, range, iter);
}
