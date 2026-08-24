// Compare JER smear functions between the current tune and the v0 tune.
// For every TF1 present in both files, draws the two curves overlaid on
// top with their ratio (current / v0) underneath, and saves one page per
// function into a multi-page PDF plus individual PNGs.
//
// Usage:
//   root -l -b -q plot_jer_comparison.C

void plot_jer_comparison(){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  const char* file_current = "jer_smear_functions.root";
  const char* file_v0      = "jer_smear_functions_v0.root";
  const char* outdir       = "plots";
  const char* pdfname      = "jer_smear_comparison.pdf";

  gSystem->mkdir(outdir, true);

  TFile *fc = TFile::Open(file_current, "READ");
  TFile *fv = TFile::Open(file_v0, "READ");
  if (!fc || fc->IsZombie() || !fv || fv->IsZombie()) {
    printf("ERROR: could not open input files\n");
    return;
  }

  // function base names common to both files
  const char* types[3] = {"nominal", "positive", "negative"};
  const int nbins = 4; // f_<type>_0 .. f_<type>_3

  double xmin = 5, xmax = 100; // plotting range (GeV)
  int npts = 200;

  // collect the (type, bin, name) triples that exist in both files first,
  // so we know which one is last (needed to close the multi-page pdf)
  struct FnInfo { TString type; int bin; TString name; };
  std::vector<FnInfo> fns;
  for (int it = 0; it < 3; ++it) {
    for (int ib = 0; ib < nbins; ++ib) {
      TString name = Form("f_%s_%d", types[it], ib);
      if (fc->Get(name) && fv->Get(name)) fns.push_back({types[it], ib, name});
      else printf("WARNING: %s missing in one of the files, skipping\n", name.Data());
    }
  }

  TString pdfpath = TString(outdir) + "/" + pdfname;

  for (size_t idx = 0; idx < fns.size(); ++idx) {
      TString name = fns[idx].name;
      TString type = fns[idx].type;
      int ib = fns[idx].bin;

      TF1 *fcur = (TF1*) fc->Get(name);
      TF1 *fv0  = (TF1*) fv->Get(name);

      fcur = (TF1*) fcur->Clone(name + "_current");
      fv0  = (TF1*) fv0->Clone(name + "_v0");
      fcur->SetRange(xmin, xmax);
      fv0->SetRange(xmin, xmax);

      fcur->SetLineColor(kAzure+2);
      fcur->SetLineWidth(3);
      fcur->SetLineStyle(1);

      fv0->SetLineColor(kRed+1);
      fv0->SetLineWidth(3);
      fv0->SetLineStyle(2);

      // find a sensible y-range for the top pad
      double ymax = 0;
      for (int i = 0; i <= npts; ++i) {
        double x = xmin + (xmax - xmin) * i / npts;
        ymax = std::max({ymax, fcur->Eval(x), fv0->Eval(x)});
      }
      ymax *= 1.3;

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

      // --- top pad: overlaid functions ---
      pad1->cd();
      TH1F *hframe = pad1->DrawFrame(xmin, 0, xmax, ymax);
      hframe->GetYaxis()->SetTitle("#sigma(p_{T}) / p_{T}");
      hframe->GetYaxis()->SetTitleSize(0.06);
      hframe->GetYaxis()->SetTitleOffset(1.05);
      hframe->GetYaxis()->SetLabelSize(0.045);
      hframe->GetXaxis()->SetLabelSize(0);
      fcur->Draw("L SAME");
      fv0->Draw("L SAME");

      TLatex lt;
      lt.SetNDC();
      lt.SetTextSize(0.06);
      lt.DrawLatex(0.18, 0.85, Form("%s, bin %d", type.Data(), ib));

      TLegend *leg = new TLegend(0.55, 0.72, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.05);
      leg->AddEntry(fcur, "current", "l");
      leg->AddEntry(fv0, "v0", "l");
      leg->Draw();

      // --- bottom pad: ratio ---
      pad2->cd();
      TGraph *gr = new TGraph(npts + 1);
      double rmin = 1e9, rmax = -1e9;
      for (int i = 0; i <= npts; ++i) {
        double x = xmin + (xmax - xmin) * i / npts;
        double yv0 = fv0->Eval(x);
        double r = (yv0 != 0) ? fcur->Eval(x) / yv0 : 0;
        gr->SetPoint(i, x, r);
        rmin = std::min(rmin, r);
        rmax = std::max(rmax, r);
      }
      double rpad = 0.1 * (rmax - rmin + 1e-6);

      TH1F *hframe2 = pad2->DrawFrame(xmin, rmin - rpad, xmax, rmax + rpad);
      hframe2->GetYaxis()->SetTitle("current / v0");
      hframe2->GetXaxis()->SetTitle("p_{T} [GeV]");
      hframe2->GetYaxis()->SetTitleSize(0.11);
      hframe2->GetXaxis()->SetTitleSize(0.11);
      hframe2->GetYaxis()->SetTitleOffset(0.55);
      hframe2->GetXaxis()->SetTitleOffset(1.15);
      hframe2->GetYaxis()->SetLabelSize(0.09);
      hframe2->GetXaxis()->SetLabelSize(0.09);
      hframe2->GetYaxis()->SetNdivisions(505);

      gr->SetLineColor(kBlack);
      gr->SetLineWidth(2);
      gr->Draw("L SAME");

      TLine *line1 = new TLine(xmin, 1.0, xmax, 1.0);
      line1->SetLineStyle(3);
      line1->SetLineColor(kGray+2);
      line1->Draw();

      c->cd();
      c->SaveAs(Form("%s/%s.png", outdir, name.Data()));

      if (idx == 0)                    c->Print(pdfpath + "(");
      else if (idx == fns.size() - 1)  c->Print(pdfpath + ")");
      else                              c->Print(pdfpath);
  }

  printf("Wrote %s and per-function PNGs in %s/\n", pdfpath.Data(), outdir);
}
