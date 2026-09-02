// Compares basic jet kinematics (pT, unsub pT, eta, phi, multiplicity) between
// two input files, overlaid. Produces both a normalized (unit-area, for shape
// comparison) and a raw/unnormalized (for comparing absolute yields) version.
//
// Usage:
//   root -l -b -q 'compare_jet_kinematics.C("fileA.root", "fileB.root")'

void draw_overlay(TPad *pad, TH1 *hA, TH1 *hB, const char *nameA, const char *nameB, bool logy, bool normalize)
{
  pad->cd();
  if (logy) pad->SetLogy();

  if (normalize) {
    if (hA->Integral() > 0) hA->Scale(1.0 / hA->Integral());
    if (hB->Integral() > 0) hB->Scale(1.0 / hB->Integral());
    hA->GetYaxis()->SetTitle("a.u.");
  } else {
    hA->GetYaxis()->SetTitle("counts");
  }

  hA->SetLineColor(kBlue + 1);
  hB->SetLineColor(kRed + 1);
  hA->SetLineWidth(2);
  hB->SetLineWidth(2);
  hA->SetStats(0);

  double ymax = std::max(hA->GetMaximum(), hB->GetMaximum());
  hA->SetMaximum(ymax * (logy ? 5 : 1.3));
  if (logy) hA->SetMinimum(normalize ? 1e-6 : 0.5);

  hA->Draw("hist");
  hB->Draw("hist same");

  TLegend *leg = new TLegend(0.62, 0.75, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->AddEntry(hA, nameA, "l");
  leg->AddEntry(hB, nameB, "l");
  leg->Draw();
}

void draw_all(const char *outpng, bool normalize,
              TH1D *hptA, TH1D *hptB, TH1D *hunsubA, TH1D *hunsubB,
              TH1D *hetaA, TH1D *hetaB, TH1D *hphiA, TH1D *hphiB,
              TH1D *hmultA, TH1D *hmultB,
              const char *nameA, const char *nameB)
{
  TCanvas *c = new TCanvas(Form("c_%d", normalize), "jet kinematics comparison", 1800, 1000);
  c->Divide(3, 2);
  draw_overlay((TPad *)c->cd(1), (TH1*)hptA->Clone(), (TH1*)hptB->Clone(), nameA, nameB, true, normalize);
  draw_overlay((TPad *)c->cd(2), (TH1*)hunsubA->Clone(), (TH1*)hunsubB->Clone(), nameA, nameB, true, normalize);
  draw_overlay((TPad *)c->cd(3), (TH1*)hmultA->Clone(), (TH1*)hmultB->Clone(), nameA, nameB, true, normalize);
  draw_overlay((TPad *)c->cd(4), (TH1*)hetaA->Clone(), (TH1*)hetaB->Clone(), nameA, nameB, false, normalize);
  draw_overlay((TPad *)c->cd(5), (TH1*)hphiA->Clone(), (TH1*)hphiB->Clone(), nameA, nameB, false, normalize);
  c->SaveAs(outpng);
}

void compare_jet_kinematics(const char *fileA = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root",
                             const char *fileB = "/home/tmengel/PPG14/rootfiles/v004_20260821/data_v004_20260821_calibrated_merged.root",
                             const char *nameA = "v001",
                             const char *nameB = "v004",
                             Long64_t nentries = 5000000,
                             const char *outbase = "jet_kinematics_compare")
{
  const char *sel = "jet_pT>7 && fabs(zvrtx)<60";

  TFile *fA = TFile::Open(fileA, "read");
  TTree *tA = (TTree *)fA->Get("T");
  Long64_t nA = std::min(nentries, tA->GetEntries());

  TFile *fB = TFile::Open(fileB, "read");
  TTree *tB = (TTree *)fB->Get("T");
  Long64_t nB = std::min(nentries, tB->GetEntries());

  TH1D *hptA = new TH1D("hptA", ";jet_pT [GeV]", 100, 0, 60);
  TH1D *hptB = new TH1D("hptB", ";jet_pT [GeV]", 100, 0, 60);
  tA->Draw("jet_pT>>hptA", sel, "goff", nA);
  tB->Draw("jet_pT>>hptB", sel, "goff", nB);

  TH1D *hunsubA = new TH1D("hunsubA", ";jet_unsub_pT [GeV]", 150, 0, 150);
  TH1D *hunsubB = new TH1D("hunsubB", ";jet_unsub_pT [GeV]", 150, 0, 150);
  tA->Draw("jet_unsub_pT>>hunsubA", sel, "goff", nA);
  tB->Draw("jet_unsub_pT>>hunsubB", sel, "goff", nB);

  TH1D *hetaA = new TH1D("hetaA", ";jet_eta", 80, -1.2, 1.2);
  TH1D *hetaB = new TH1D("hetaB", ";jet_eta", 80, -1.2, 1.2);
  tA->Draw("jet_eta>>hetaA", sel, "goff", nA);
  tB->Draw("jet_eta>>hetaB", sel, "goff", nB);

  TH1D *hphiA = new TH1D("hphiA", ";jet_phi", 64, -3.2, 3.2);
  TH1D *hphiB = new TH1D("hphiB", ";jet_phi", 64, -3.2, 3.2);
  tA->Draw("jet_phi>>hphiA", sel, "goff", nA);
  tB->Draw("jet_phi>>hphiB", sel, "goff", nB);

  TH1D *hmultA = new TH1D("hmultA", ";jets per event (pT>7)", 15, 0, 15);
  TH1D *hmultB = new TH1D("hmultB", ";jets per event (pT>7)", 15, 0, 15);
  tA->Draw("Sum$(jet_pT>7)>>hmultA", "fabs(zvrtx)<60", "goff", nA);
  tB->Draw("Sum$(jet_pT>7)>>hmultB", "fabs(zvrtx)<60", "goff", nB);

  draw_all(Form("%s_norm.png", outbase), true,
           hptA, hptB, hunsubA, hunsubB, hetaA, hetaB, hphiA, hphiB, hmultA, hmultB, nameA, nameB);
  draw_all(Form("%s_raw.png", outbase), false,
           hptA, hptB, hunsubA, hunsubB, hetaA, hetaB, hphiA, hphiB, hmultA, hmultB, nameA, nameB);
}
