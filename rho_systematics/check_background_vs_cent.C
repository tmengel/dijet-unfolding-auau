// Draws (jet_unsub_pT - jet_pT) vs centrality for two input files side by side,
// with the makeUnfoldingHists.C background-quality cut fcut(cent) overlaid.
//
// Usage:
//   root -l -b -q 'check_background_vs_cent.C("fileA.root", "fileB.root")'

#include "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/PlotUtils.h"

const std::string & tree_crit = "jet_E>0 && jet_accept_eta==1 && jet_pT > 5 && fabs(zvrtx)<60";

// void draw_v001(TPad *pad, const char *path, const char *label, Long64_t nentries)
// {
//   pad->cd();
//   pad->SetRightMargin(0.14);
//   pad->SetLogz();

//   TFile *f = TFile::Open(path, "read");
//   TTree *t = (TTree *)f->Get("T");
//   Long64_t n = t->GetEntries();
//   if (nentries > n) nentries = n;

//   TH2D *h = new TH2D(Form("h_%s", label), Form("%s;centrality;jet_unsub_pT - jet_pT  [GeV]", label),
//                       100, 0, 100, 250, -25, 100);
//   t->Draw(Form("(jet_unsub_pT-jet_pT):cent>>h_%s", label), "jet_pT>7 && fabs(zvrtx)<60", "goff", nentries);
//   h->SetStats(0);
//   h->Draw("COLZ");

//   // same fcut(cent) used as the background-quality cut in makeUnfoldingHists.C
//   TF1 *fcut = new TF1(Form("fcut_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
//   fcut->SetParameters(0.0, 40, 0.038);
//   fcut->SetLineColor(kRed);
//   fcut->SetLineWidth(2);
//   fcut->Draw("same");
// }
void draw_fit(TPad *pad, const char *path, const char *label, Long64_t nentries)
{
  pad->cd();
  pad->SetRightMargin(0.14);
  pad->SetLogz();

  TFile *f = TFile::Open(path, "read");
  TTree *t = (TTree *)f->Get("T");
  Long64_t n = t->GetEntries();
  if (nentries > n) nentries = n;

  TH2D *h = new TH2D(Form("h_%s", label), Form("%s;centrality;jet_unsub_pT - jet_pT  [GeV]", label), 50, 0, 100, 250, -25, 60);
  t->Draw(Form("(jet_comp_pT-jet_pT):cent>>h_%s", label), "jet_E>0 && jet_accept_eta==1 && jet_pT > 5 && fabs(zvrtx)<60", "goff");
 
  // make a profile 
  TProfile * p = h->ProfileX(Form("p_%s", label), 1, -1, "s");
  p->SetLineColor(kBlue);
  p->SetLineWidth(2);

  TGraph * g = new TGraph();
  // for every 2% central x slice, get the 1D proj, normalize, fit to gaus, extract mean and sigma, set point to mean + 3.5*sigma
  int xbin_60 = h->GetXaxis()->FindBin(60);
  for (int i = 1; i <= xbin_60; ++i) {
    TH1D *proj = h->ProjectionY(Form("proj_%s_%d", label, i), i, i);
    proj->Scale(1.0 / proj->Integral());
    TF1 *fit = new TF1(Form("fit_%s_%d", label, i), "gaus", -25, 60);
    proj->Fit(fit, "Q");
    proj->Fit(fit, "QR", "", fit->GetParameter(1) - 2*fit->GetParameter(2), fit->GetParameter(1) + 2*fit->GetParameter(2));
    double mu = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    g->SetPoint(i-1, h->GetXaxis()->GetBinCenter(i), mu + 3.5 * sigma);
  }
  g -> SetLineColor(kAzure-2);
  g -> SetLineWidth(2);

  h->SetStats(0);
  h->GetYaxis()->SetRangeUser(-5, 60);
  h->Draw("COLZ");
  p->Draw("same");
  g->Draw("L same");
  // create a second profile to show mu + 3.5 sigma
  TGraph *p2 = new TGraph();
  for (int i = 1; i <= p->GetNbinsX(); ++i) {
    double mu = p->GetBinContent(i);
    double sigma = p->GetBinError(i);
    p2->SetPoint(i-1, p->GetBinCenter(i), mu + 3.5 * sigma);
  }
  p2->SetLineColor(kRed);
  p2->SetLineWidth(2);
  // p2->Draw("L same");
  // fit p2 
  TF1 *f1 = new TF1(Form("f1_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
  g->Fit(f1, "R","",0, 60);
  f1->SetLineColor(kMagenta);
  f1->SetLineWidth(2);
  f1->Draw("same");
  
  // same fcut(cent) used as the background-quality cut in makeUnfoldingHists.C
  TF1 *fcut = new TF1(Form("fcut_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
  fcut->SetParameters(0.0, 40, 0.038);
  fcut->SetLineColor(kRed);
  fcut->SetLineWidth(2);
  fcut->SetLineStyle(2);
  fcut->Draw("same");

  PlotUtils::myText(0.2, 0.9, kBlack, label, 0.04);

  TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(p, "Profile", "l");
  leg->AddEntry(g, "Mean + 3.5 #sigma", "l");
  leg->AddEntry(f1, Form("p0 = %.2f, p1 = %.2f, p2 = %.3f", f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2)), "l");
  leg->AddEntry(fcut, Form("p0 = %.2f, p1 = %.2f, p2 = %.3f", fcut->GetParameter(0), fcut->GetParameter(1), fcut->GetParameter(2)), "l");
  leg->Draw("same");

}

void check_background_vs_cent(const char *fileA = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/input/data/data_v005_20260903_calibrated_merged.root",
                               const char *fileB = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/input/data/data_v004_20260821_calibrated_merged.root",
                               Long64_t nentries = 10000000,
                               const char *outpng = "background_vs_cent.png")
{

  PlotUtils::set_sphenix_style();
  TCanvas *c = new TCanvas("c", "background vs centrality", 1600, 700);
  c->Divide(2, 1);

  draw_fit((TPad *)c->cd(1), fileA, "v005", nentries);
  draw_fit((TPad *)c->cd(2), fileB, "v004", nentries);

  c->SaveAs(outpng);
}
