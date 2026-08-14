#include "PlotUtils.h"
void quick_comp()
{
    auto *f11 = TFile::Open("/home/tmengel/PPG14/version0/v001_20260715/response_matrices/response_matrix_AA_cent_0_r03_PRIMER1_PRIOR.root", "READ");
    auto *f12 = TFile::Open("/home/tmengel/PPG14/version0/v001_20260715/unfolding_hists/unfolding_hists_AA_cent_0_r03_PRIMER1_PRIOR.root", "READ");
    auto *f21 = TFile::Open("/home/tmengel/PPG14/version1/dijet-unfolding-final/response_matrices/response_matrix_AA_cent_0_r03_PRIMER1_PRIOR.root", "READ");
    auto *f22 = TFile::Open("/home/tmengel/PPG14/version1/dijet-unfolding-final/unfolding_hists/unfolding_hists_AA_cent_0_r03_PRIMER1_PRIOR.root", "READ");

    auto * h11 = (TH1D*) f11->Get("h_truth_flat_pt1pt2");
    auto * h12 = (TH1D*) f12->Get("h_flat_unfold_pt1pt2_1");
    h11 -> SetName("h11");
    h12 -> SetName("h12");
    h11 -> Scale( 1.0 / h11->Integral() );
    h12 -> Scale( 1.0 / h12->Integral() );
    h11 -> SetLineColor(kRed);
    h12 -> SetLineColor(kRed);
    auto * h21 = (TH1D*) f21->Get("h_truth_flat_pt1pt2");
    auto * h22 = (TH1D*) f22->Get("h_flat_unfold_pt1pt2_1");
    h21 -> SetName("h21");
    h22 -> SetName("h22");
    h21 -> Scale( 1.0 / h21->Integral() );
    h22 -> Scale( 1.0 / h22->Integral() );
    h21 -> SetLineColor(kBlue);
    h22 -> SetLineColor(kBlue);

    auto ratio = (TH1D*) h12->Clone("ratio");
    ratio -> Divide(h11);
    ratio -> SetLineColor(kRed);
    auto ratio2 = (TH1D*) h22->Clone("ratio2");
    ratio2 -> Divide(h21);
    ratio2 -> SetLineColor(kBlue);

    PlotUtils::set_sphenix_style();
    auto * c1 = new TCanvas();
    c1 -> cd();
    // h11 -> Draw("hist");
    // h21 -> Draw("hist same");
    // h12 -> Draw("hist same");
    // h22 -> Draw("hist same");
    ratio -> Draw("hist");
    ratio2 -> Draw("hist same");
    PlotUtils::myText(0.2, 0.8,kRed, "v0");
    PlotUtils::myText(0.2, 0.7,kBlue, "v1");
    c1 -> SetLogy();
    c1 -> Print("quick_comp.pdf");

}
