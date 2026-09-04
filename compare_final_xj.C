#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
#include "PlotUtils.h"

namespace {

bool compatible_binning(const TH1 *a, const TH1 *b)
{
  if (!a || !b || a->GetNbinsX() != b->GetNbinsX()) return false;

  const int n = a->GetNbinsX();
  for (int i = 1; i <= n + 1; ++i) {
    const double edge_a = a->GetXaxis()->GetBinLowEdge(i);
    const double edge_b = b->GetXaxis()->GetBinLowEdge(i);
    if (std::abs(edge_a - edge_b) > 1.e-9) return false;
  }
  return true;
}

void style_hist(TH1 *h, int color, int marker)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  // h->SetMarkerSize(1.1);
  // h->SetLineWidth(2);
}

void style_graph(TGraphAsymmErrors *g, int color)
{
  if (!g) return;
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetFillColorAlpha(color, 0.25);
  g->SetLineWidth(1);
}

} // namespace


void compare_final_xj(const char *file_v1,
                      const char *file_v2,
                      int irange = 0,
                      int iter = 1,
                      int centrality_bin = 0,
                      const char *label_v1 = "version 1",
                      const char *label_v2 = "version 2",
                      const char *output = "compare_final_xj.pdf",
                      bool compare_pp = false)
{


    gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // #include "PlotUtils.h"
  PlotUtils::set_sphenix_style();
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  // dlutility::SetyjPadStyle();

  PlotUtils::set_sphenix_style();
  int cone_size = 3;

  read_binning rb("/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config");

  int arr_color_unfold_fill[5] = {kRed + 1, kAzure + 2, kGreen + 3, kViolet - 1, kBlack};
  int arr_color_unfold[5] = {kRed + 1, kAzure + 2, kGreen+3, kViolet -1, kBlack};
  int arr_marker_unfold[5] = {20, 21, 33, 34, 24};
  float arr_msize_unfold[5] = {1.1, 1, 1.5, 1.2, 1.1};
  std::cout << "centrality_bin: " << centrality_bin << std::endl;
  std::cout << "arr_color_unfold_fill[centrality_bin]: " << arr_color_unfold_fill[centrality_bin] << std::endl;
  const int color_unfold_fill = arr_color_unfold_fill[centrality_bin];
  const int color_unfold= arr_color_unfold[centrality_bin];
  const float marker_unfold= arr_marker_unfold[centrality_bin];

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  std::string dphi_string = rb.get_dphi_string();
  
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



  const int color1 = kGray + 2;
  int color2 = arr_color_unfold[centrality_bin];

  TFile *f1 = TFile::Open(file_v1, "READ");
  TFile *f2 = TFile::Open(file_v2, "READ");
  if (!f1 || f1->IsZombie()) {
    std::cerr << "ERROR: cannot open " << file_v1 << std::endl;
    return;
  }
  if (!f2 || f2->IsZombie()) {
    std::cerr << "ERROR: cannot open " << file_v2 << std::endl;
    f1->Close();
    return;
  }

  const TString hist_name = compare_pp
      ? Form("h_final_xj_pp_unfold_range_%d_iter_%d", irange, iter)
      : Form("h_final_xj_unfold_range_%d_iter_%d", irange, iter);

  const TString syst_name = compare_pp
      ? Form("g_final_xj_pp_systematics_range_%d_iter_%d", irange, iter)
      : Form("g_final_xj_systematics_range_%d_iter_%d", irange, iter);

  TH1D *h1_in = dynamic_cast<TH1D *>(f1->Get(hist_name));
  TH1D *h2_in = dynamic_cast<TH1D *>(f2->Get(hist_name));
  TGraphAsymmErrors *g1_in = dynamic_cast<TGraphAsymmErrors *>(f1->Get(syst_name));
  TGraphAsymmErrors *g2_in = dynamic_cast<TGraphAsymmErrors *>(f2->Get(syst_name));

  if (!h1_in || !h2_in) {
    std::cerr << "ERROR: could not find " << hist_name << " in both files.\n"
              << "  file 1: " << (h1_in ? "found" : "missing") << "\n"
              << "  file 2: " << (h2_in ? "found" : "missing") << std::endl;
    f1->Close();
    f2->Close();
    return;
  }

  if (!compatible_binning(h1_in, h2_in)) {
    std::cerr << "ERROR: histogram binning differs between the two files." << std::endl;
    f1->Close();
    f2->Close();
    return;
  }

  // Detach clones from their input files.
  TH1D *h1 = dynamic_cast<TH1D *>(h1_in->Clone("h_xj_version1"));
  TH1D *h2 = dynamic_cast<TH1D *>(h2_in->Clone("h_xj_version2"));
  h1->SetDirectory(nullptr);
  h2->SetDirectory(nullptr);

  TGraphAsymmErrors *g1 = g1_in
      ? dynamic_cast<TGraphAsymmErrors *>(g1_in->Clone("g_xj_syst_version1"))
      : nullptr;
  TGraphAsymmErrors *g2 = g2_in
      ? dynamic_cast<TGraphAsymmErrors *>(g2_in->Clone("g_xj_syst_version2"))
      : nullptr;

  f1->Close();
  f2->Close();
  
  style_hist(h1, color1, 20);
  style_hist(h2, color2, 24);
  style_graph(g1, color1);
  style_graph(g2, color2);
  // shift all the points in the second graph slightly to the right for visibility
  // if (g2) {
  //   const int n2 = g2->GetN();
  //   for (int i = 0; i < n2; ++i) {
  //     double x, y;
  //     g2->GetPoint(i, x, y);
  //     g2->SetPoint(i, x + 0.005, y);
  //   }
  // }
  
  TH1D *ratio = dynamic_cast<TH1D *>(h2->Clone("h_xj_ratio_v2_over_v1"));
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  ratio->SetTitle(";x_{J};Ratio");
  style_hist(ratio, color2, 20);

  for (int ibin = 1; ibin <= h1->GetNbinsX(); ++ibin) {
    const double a = h1->GetBinContent(ibin);
    const double b = h2->GetBinContent(ibin);
    const double ea = h1->GetBinError(ibin);
    const double eb = h2->GetBinError(ibin);

    if (a == 0.0) {
      ratio->SetBinContent(ibin, 0.0);
      ratio->SetBinError(ibin, 0.0);
      continue;
    }

    const double r = b / a;
    double er = 0.0;
    if (b != 0.0) {
      er = std::abs(r) * std::sqrt((ea / a) * (ea / a) +
                                   (eb / b) * (eb / b));
    } else {
      er = std::abs(eb / a);
    }
    ratio->SetBinContent(ibin, r);
    ratio->SetBinError(ibin, er);
  }

  const double ymax = 1.30 * std::max(h1->GetMaximum(), h2->GetMaximum());

  TCanvas *c = new TCanvas("c_compare_final_xj", "compare final xj", 650, 750);

  TPad *p_top = new TPad("p_top", "p_top", 0.0, 0.30, 1.0, 1.0);
  TPad *p_bot = new TPad("p_bot", "p_bot", 0.0, 0.00, 1.0, 0.30);
  p_top->SetLeftMargin(0.14);
  p_top->SetRightMargin(0.04);
  p_top->SetTopMargin(0.06);
  p_top->SetBottomMargin(0.02);
  p_bot->SetLeftMargin(0.14);
  p_bot->SetRightMargin(0.04);
  p_bot->SetTopMargin(0.02);
  p_bot->SetBottomMargin(0.34);
  p_top->Draw();
  p_bot->Draw();

  p_top->cd();
  h1->SetTitle(Form("%s, range %d, iteration %d;x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}",
                    compare_pp ? "p+p" : "Au+Au", irange, iter));
  h1->SetMinimum(0.0);
  h1->SetMaximum(ymax > 0.0 ? ymax : 1.0);
  h1->GetYaxis()->SetRangeUser(0.0, h1->GetMaximum() * 1.5);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetLabelSize(0.044);
  h1->GetXaxis()->SetLabelSize(0.0);
  h1 -> GetXaxis()->SetRangeUser(0.2, 1.0);
  h1->Draw("E");

  g1 -> SetFillColorAlpha(color1, 0.25);
  g2 -> SetFillColorAlpha(color2, 0.1);
  if (g1) g1->Draw("E2 SAME");
  if (g2) g2->Draw("E2 SAME");
  h1->Draw("E P SAME");
  h2->Draw("E P SAME");

  std::vector< std::string > my_tags;
  my_tags.push_back("#it{#bf{sPHENIX} Internal}");
  my_tags.push_back(Form("Au+Au#kern[0.02]{#sqrt{#it{s_{_{NN}}}} = 200 GeV}"));
  my_tags.push_back(Form("anti-#it{k}_{t} #kern[-0.1]{#it{R} = %0.1f}", cone_size*0.1));
  my_tags.push_back(Form("%2.1f #kern[-0.07]{#leq #it{p}_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]));
  my_tags.push_back(Form("#it{p}_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]));
  my_tags.push_back(Form("#Delta#phi #geq %s", dphi_string.c_str()));
  my_tags.push_back(Form("%d - %d%%", (int)icentrality_bins[centrality_bin], (int)icentrality_bins[centrality_bin+1]));
  double x1=0.18;
  double y1=0.88;
  double dy=0.05;
  for (size_t i = 0; i < my_tags.size(); ++i) {
    PlotUtils::draw_text(x1, y1 - i*dy, kBlack, my_tags[i].c_str(), 0.04);
  }
  TLegend *leg = new TLegend(0.6, 0.88 - 2*dy, 0.85, 0.88);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->Draw("same");

  leg->AddEntry(h1, label_v1, "lep");
  leg->AddEntry(h2, label_v2, "lep");
  leg->Draw();

  p_bot->cd();
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(1.05);
  ratio->GetXaxis()->SetLabelSize(0.10);
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.58);
  ratio->GetYaxis()->SetLabelSize(0.095);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetXaxis()->SetRangeUser(0.2, 1.0);
  ratio->Draw("E1");

  

  // TLine *unity = new TLine(ratio->GetXaxis()->GetXmin(), 1.0,
  //                          ratio->GetXaxis()->GetXmax(), 1.0);
  TLine *unity = new TLine(0.2, 1.0, 1.0, 1.0);
  unity->SetLineColor(kGray + 2);
  unity->SetLineStyle(2);
  unity->SetLineWidth(2);
  unity->Draw("SAME");

  c->SaveAs(output);

  std::cout << "Compared object: " << hist_name << "\n"
            << "Systematics object: " << syst_name << "\n"
            << "Saved: " << output << std::endl;
}