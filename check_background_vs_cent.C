// Draws (jet_unsub_pT - jet_pT) vs centrality for two input files side by side,
// with the makeUnfoldingHists.C background-quality cut fcut(cent) overlaid.
//
// Usage:
//   root -l -b -q 'check_background_vs_cent.C("fileA.root", "fileB.root")'

void draw_one(TPad *pad, const char *path, const char *label, Long64_t nentries)
{
  pad->cd();
  pad->SetRightMargin(0.14);
  pad->SetLogz();

  TFile *f = TFile::Open(path, "read");
  TTree *t = (TTree *)f->Get("T");
  Long64_t n = t->GetEntries();
  if (nentries > n) nentries = n;

  TH2D *h = new TH2D(Form("h_%s", label), Form("%s;centrality;jet_unsub_pT - jet_pT  [GeV]", label),
                      100, 0, 100, 250, -25, 100);
  t->Draw(Form("(jet_unsub_pT-jet_pT):cent>>h_%s", label), "jet_pT>7 && fabs(zvrtx)<60", "goff", nentries);
  h->SetStats(0);
  h->Draw("COLZ");

  // same fcut(cent) used as the background-quality cut in makeUnfoldingHists.C
  TF1 *fcut = new TF1(Form("fcut_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
  fcut->SetParameters(0.0, 40, 0.038);
  fcut->SetLineColor(kRed);
  fcut->SetLineWidth(2);
  fcut->Draw("same");
}
void draw_other(TPad *pad, const char *path, const char *label, Long64_t nentries)
{
  pad->cd();
  pad->SetRightMargin(0.14);
  pad->SetLogz();

  TFile *f = TFile::Open(path, "read");
  TTree *t = (TTree *)f->Get("T");
  Long64_t n = t->GetEntries();
  if (nentries > n) nentries = n;

  TH2D *h = new TH2D(Form("h_%s", label), Form("%s;centrality;jet_unsub_pT - jet_pT  [GeV]", label), 50, 0, 100, 250, -25, 100);
  t->Draw(Form("(jet_comp_pT-jet_pT):cent>>h_%s", label), "jet_E>0&&jet_accept_eta==1&&jet_pT>5 && fabs(zvrtx)<60", "goff", nentries);
 
  // make a profile 
  TProfile * p = h->ProfileX(Form("p_%s", label), 1, -1, "s");
  p->SetLineColor(kBlue);
  p->SetLineWidth(2);

  h->SetStats(0);
  h->Draw("COLZ");
  p->Draw("same");
  // create a second profile to show mu + 3.5 sigma
  TGraph *p2 = new TGraph();
  for (int i = 1; i <= p->GetNbinsX(); ++i) {
    double mu = p->GetBinContent(i);
    double sigma = p->GetBinError(i);
    p2->SetPoint(i-1, p->GetBinCenter(i), mu + 3.5 * sigma);
  }
  p2->SetLineColor(kRed);
  p2->SetLineWidth(2);
  p2->Draw("L same");
  // fit p2 
  TF1 *f1 = new TF1(Form("f1_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
  p2->Fit(f1, "R","",0, 60);
  f1->SetLineColor(kMagenta);
  f1->SetLineWidth(2);
  f1->Draw("same");
  
  // same fcut(cent) used as the background-quality cut in makeUnfoldingHists.C
  TF1 *fcut = new TF1(Form("fcut_%s", label), "[0]+[1]*TMath::Exp(-[2]*x)", 0, 100);
  fcut->SetParameters(0.0, 40, 0.038);
  fcut->SetLineColor(kRed);
  fcut->SetLineWidth(2);
  fcut->Draw("same");
}

void check_background_vs_cent(const char *fileA = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root",
                               const char *fileB = "/home/tmengel/PPG14/rootfiles/data_v004_20260821_calibrated_merged.root",
                               Long64_t nentries = 5000000,
                               const char *outpng = "background_vs_cent.png")
{
  TCanvas *c = new TCanvas("c", "background vs centrality", 1600, 700);
  c->Divide(2, 1);

  draw_one((TPad *)c->cd(1), fileA, "v001", nentries);
  draw_other((TPad *)c->cd(2), fileB, "v004", nentries);

  c->SaveAs(outpng);
}
