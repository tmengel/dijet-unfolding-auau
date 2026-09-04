// Compare the cent0 background summary maps (h_pt1_pt2_Signal/ZYAM/ratio/sub,
// each a 15x15 TH2D in pT1 vs pT2) between this version's
// background_hists_AA_cent_0_r03.root and the v0 tune's. For each map, draws
// the current version, the v0 version, and their ratio side by side.
//
// Usage:
//   root -l -b -q plot_background_comparison.C

void plot_background_comparison(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(kBird);

  const char* file_current = "/home/tmengel/PPG14/version1/rho_jet_starting/unfolding_hists/background_hists_AA_cent_0_r03.root";
  const char* file_v0      = "/home/tmengel/PPG14/version0/v001_20260715/unfolding_hists/background_hists_AA_cent_0_r03.root";
  const char* outdir       = "plots";
  const char* pdfname      = "background_comparison.pdf";

  gSystem->mkdir(outdir, true);

  TFile *fc = TFile::Open(file_current, "READ");
  TFile *fv = TFile::Open(file_v0, "READ");
  if (!fc || fc->IsZombie() || !fv || fv->IsZombie()) {
    printf("ERROR: could not open input files\n");
    return;
  }

  std::vector<TString> names = {
    "h_pt1_pt2_Signal", "h_pt1_pt2_ZYAM", "h_pt1_pt2_ratio", "h_pt1_pt2_sub"
  };

  TString pdfpath = TString(outdir) + "/" + pdfname;

  for (size_t idx = 0; idx < names.size(); ++idx) {

    TString name = names[idx];

    TH2D *hcur = (TH2D*) fc->Get(name);
    TH2D *hv0  = (TH2D*) fv->Get(name);
    if (!hcur || !hv0) {
      printf("WARNING: %s missing in one of the files, skipping\n", name.Data());
      continue;
    }

    hcur = (TH2D*) hcur->Clone(name + "_current");
    hv0  = (TH2D*) hv0->Clone(name + "_v0");
    hcur->SetDirectory(0);
    hv0->SetDirectory(0);

    TH2D *hratio = (TH2D*) hcur->Clone(name + "_ratio");
    hratio->SetDirectory(0);
    hratio->Divide(hv0);

    double zmax = std::max(hcur->GetMaximum(), hv0->GetMaximum());
    double zmin = std::min(hcur->GetMinimum(), hv0->GetMinimum());
    hcur->GetZaxis()->SetRangeUser(zmin, zmax);
    hv0->GetZaxis()->SetRangeUser(zmin, zmax);

    hcur->GetXaxis()->SetTitle("p_{T,1} [GeV]");
    hcur->GetYaxis()->SetTitle("p_{T,2} [GeV]");
    hv0->GetXaxis()->SetTitle("p_{T,1} [GeV]");
    hv0->GetYaxis()->SetTitle("p_{T,2} [GeV]");
    hratio->GetXaxis()->SetTitle("p_{T,1} [GeV]");
    hratio->GetYaxis()->SetTitle("p_{T,2} [GeV]");
    hratio->GetZaxis()->SetTitle("current / v0");
    hratio->GetZaxis()->SetRangeUser(0.5, 1.5);

    TCanvas *c = new TCanvas(Form("c_%s", name.Data()), name, 1500, 500);
    c->Divide(3, 1);

    c->cd(1);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    hcur->Draw("COLZ");
    TLatex lt1; lt1.SetNDC(); lt1.SetTextSize(0.05);
    lt1.DrawLatex(0.15, 0.92, Form("%s  current", name.Data()));

    c->cd(2);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    hv0->Draw("COLZ");
    TLatex lt2; lt2.SetNDC(); lt2.SetTextSize(0.05);
    lt2.DrawLatex(0.15, 0.92, Form("%s  v0", name.Data()));

    c->cd(3);
    gPad->SetRightMargin(0.16);
    gPad->SetLeftMargin(0.13);
    hratio->Draw("COLZ");
    TLatex lt3; lt3.SetNDC(); lt3.SetTextSize(0.05);
    lt3.DrawLatex(0.15, 0.92, Form("%s  current/v0", name.Data()));

    c->cd();
    c->SaveAs(Form("%s/%s.png", outdir, name.Data()));

    if (idx == 0)                    c->Print(pdfpath + "(");
    else if (idx == names.size() - 1) c->Print(pdfpath + ")");
    else                               c->Print(pdfpath);
  }

  printf("Wrote %s and per-histogram PNGs in %s/\n", pdfpath.Data(), outdir);
}
