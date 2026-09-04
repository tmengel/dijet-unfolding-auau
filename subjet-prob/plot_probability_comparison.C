// Compare the reconstruction-probability correction histograms between
// this version's probability_hists_AA_r03.root and the v0 tune's.
// Only h_pt2_bin_log_correction_0..3 exist by the same name in both files;
// each is drawn overlaid on top with its ratio (current / v0) underneath.
//
// Usage:
//   root -l -b -q plot_probability_comparison.C

void plot_probability_comparison(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const char* file_current = "/home/tmengel/PPG14/version1/rho_jet_starting/unfolding_hists/probability_hists_AA_r03.root";
  const char* file_v0      = "/home/tmengel/PPG14/version0/v001_20260715/unfolding_hists/probability_hists_AA_r03.root";
  const char* outdir       = "plots";
  const char* pdfname      = "probability_comparison.pdf";

  gSystem->mkdir(outdir, true);

  TFile *fc = TFile::Open(file_current, "READ");
  TFile *fv = TFile::Open(file_v0, "READ");
  if (!fc || fc->IsZombie() || !fv || fv->IsZombie()) {
    printf("ERROR: could not open input files\n");
    return;
  }

  const int nbins = 4; // h_pt2_bin_log_correction_0 .. 3
  TString pdfpath = TString(outdir) + "/" + pdfname;

  for (int ib = 0; ib < nbins; ++ib) {

    TString name = Form("h_pt2_bin_log_correction_%d", ib);

    TH1D *hcur = (TH1D*) fc->Get(name);
    TH1D *hv0  = (TH1D*) fv->Get(name);
    if (!hcur || !hv0) {
      printf("WARNING: %s missing in one of the files, skipping\n", name.Data());
      continue;
    }

    hcur = (TH1D*) hcur->Clone(name + "_current");
    hv0  = (TH1D*) hv0->Clone(name + "_v0");
    hcur->SetDirectory(0);
    hv0->SetDirectory(0);

    hcur->SetLineWidth(0);
    hcur->SetMarkerStyle(20);
    hcur->SetMarkerColor(kAzure+2);
    hcur->SetMarkerSize(1.1);

    hv0->SetLineWidth(0);
    hv0->SetMarkerStyle(24);
    hv0->SetMarkerColor(kRed+1);
    hv0->SetMarkerSize(1.1);

    double ymax = std::max(hcur->GetMaximum(), hv0->GetMaximum()) * 1.3;
    double xmin = hcur->GetXaxis()->GetXmin();
    double xmax = hcur->GetXaxis()->GetXmax();

    TCanvas *c = new TCanvas(Form("c_%s", name.Data()), name, 700, 700);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.32);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.32);
    pad2->SetLeftMargin(0.14);
    pad2->Draw();

    // --- top pad: overlaid histograms ---
    pad1->cd();
    TH1F *hframe = pad1->DrawFrame(xmin, 0, xmax, ymax);
    hframe->GetYaxis()->SetTitle("reconstruction probability");
    hframe->GetYaxis()->SetTitleSize(0.06);
    hframe->GetYaxis()->SetTitleOffset(1.05);
    hframe->GetYaxis()->SetLabelSize(0.045);
    hframe->GetXaxis()->SetLabelSize(0);
    hcur->Draw("PE SAME");
    hv0->Draw("PE SAME");

    TLatex lt;
    lt.SetNDC();
    lt.SetTextSize(0.06);
    lt.DrawLatex(0.18, 0.85, Form("bin %d", ib));

    TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(hcur, "current", "p");
    leg->AddEntry(hv0, "v0", "p");
    leg->Draw();

    // --- bottom pad: ratio ---
    pad2->cd();
    TH1D *hratio = (TH1D*) hcur->Clone(name + "_ratio");
    hratio->SetDirectory(0);
    hratio->Divide(hv0);
    hratio->SetLineWidth(0);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);

    double rmin = 1e9, rmax = -1e9;
    for (int i = 1; i <= hratio->GetNbinsX(); ++i) {
      double c_ = hratio->GetBinContent(i);
      if (c_ == 0) continue;
      rmin = std::min(rmin, c_);
      rmax = std::max(rmax, c_);
    }
    if (rmin > rmax) { rmin = 0.9; rmax = 1.1; }
    double rpad = 0.15 * (rmax - rmin + 1e-6);

    TH1F *hframe2 = pad2->DrawFrame(xmin, rmin - rpad, xmax, rmax + rpad);
    hframe2->GetYaxis()->SetTitle("current / v0");
    hframe2->GetXaxis()->SetTitle(hcur->GetXaxis()->GetTitle());
    hframe2->GetYaxis()->SetTitleSize(0.11);
    hframe2->GetXaxis()->SetTitleSize(0.11);
    hframe2->GetYaxis()->SetTitleOffset(0.55);
    hframe2->GetXaxis()->SetTitleOffset(1.15);
    hframe2->GetYaxis()->SetLabelSize(0.09);
    hframe2->GetXaxis()->SetLabelSize(0.09);
    hframe2->GetYaxis()->SetNdivisions(505);

    hratio->Draw("PE SAME");

    TLine *line1 = new TLine(xmin, 1.0, xmax, 1.0);
    line1->SetLineStyle(3);
    line1->SetLineColor(kGray+2);
    line1->Draw();

    c->cd();
    c->SaveAs(Form("%s/%s.png", outdir, name.Data()));

    if (ib == 0)              c->Print(pdfpath + "(");
    else if (ib == nbins - 1) c->Print(pdfpath + ")");
    else                       c->Print(pdfpath);
  }

  printf("Wrote %s and per-bin PNGs in %s/\n", pdfpath.Data(), outdir);
}
