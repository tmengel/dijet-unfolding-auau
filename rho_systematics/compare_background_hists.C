// Compare the histograms in two background-histogram ROOT files.
// Run with:
//   root -l -b -q 'compare_background_hists.C()'
// The result is a multipage PDF in comparison_plots/.

#include <map>
#include <set>
#include <string>

#include "TCanvas.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"

namespace {
void collectHistograms(TDirectory *directory, const std::string &prefix,
                       std::map<std::string, TH1 *> &histograms) {
  TIter next(directory->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next())) {
    const std::string name = key->GetName();
    const std::string path = prefix.empty() ? name : prefix + "/" + name;
    TClass *klass = TClass::GetClass(key->GetClassName());
    if (klass && klass->InheritsFrom(TDirectory::Class())) {
      collectHistograms(static_cast<TDirectory *>(key->ReadObj()), path, histograms);
    } else if (klass && klass->InheritsFrom(TH1::Class())) {
      // Ignore older cycles of the same object.
      TH1 *hist = static_cast<TH1 *>(key->ReadObj());
      if (!histograms.count(path)) histograms[path] = hist;
    }
  }
}

void labelPad(const char *text) {
  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.035);
  label.DrawLatex(0.13, 0.93, text);
}

bool isEmpty(const TH1 *hist) {
  for (int x = 0; x <= hist->GetNbinsX() + 1; ++x)
    for (int y = 0; y <= hist->GetNbinsY() + 1; ++y)
      for (int z = 0; z <= hist->GetNbinsZ() + 1; ++z)
        if (hist->GetBinContent(hist->GetBin(x, y, z)) != 0.) return false;
  return true;
}
} // namespace

void compare_background_hists(
    const char *newFile = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/unfolding_hists/background_hists_AA_cent_0_r03.root",
    const char *oldFile = "/home/tmengel/PPG14/version0/v001_20260715/unfolding_hists/background_hists_AA_cent_0_r03.root",
    const char *output = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/comparison_plots/background_hists_AA_cent_0_r03_comparison.pdf") {
  TFile current(newFile, "READ");
  TFile previous(oldFile, "READ");
  if (current.IsZombie() || previous.IsZombie()) {
    Error("compare_background_hists", "Could not open one or both input files.");
    return;
  }

  std::map<std::string, TH1 *> newHists, oldHists;
  collectHistograms(&current, "", newHists);
  collectHistograms(&previous, "", oldHists);

  std::set<std::string> names;
  for (const auto &entry : newHists) names.insert(entry.first);
  for (const auto &entry : oldHists) names.insert(entry.first);
  if (names.empty()) {
    Warning("compare_background_hists", "No histograms found.");
    return;
  }

  gSystem->mkdir(gSystem->DirName(output), true);
  gStyle->SetOptStat(0);
  TCanvas canvas("canvas", "background histogram comparison", 1200, 900);
  bool first = true;

  for (const auto &name : names) {
    TH1 *newHist = newHists.count(name) ? newHists[name] : nullptr;
    TH1 *oldHist = oldHists.count(name) ? oldHists[name] : nullptr;
    if (newHist && oldHist && isEmpty(newHist) && isEmpty(oldHist)) continue;
    TH1 *ratio1D = nullptr;
    TH2 *ratio2D = nullptr;
    TLegend *legend = nullptr;
    TLine *unity = nullptr;
    canvas.Clear();

    if (!newHist || !oldHist) {
      canvas.cd();
      labelPad((name + (newHist ? " (missing in previous file)" : " (missing in current file)")).c_str());
    } else if (newHist->GetDimension() == 1 && oldHist->GetDimension() == 1) {
      canvas.Divide(1, 2);
      canvas.cd(1);
      newHist->SetLineColor(kBlue + 1); newHist->SetMarkerColor(kBlue + 1);
      oldHist->SetLineColor(kRed + 1);  oldHist->SetMarkerColor(kRed + 1);
      newHist->SetLineWidth(2); oldHist->SetLineWidth(2);
      newHist->SetTitle(name.c_str());
      newHist->Draw("E"); oldHist->Draw("E SAME");
      legend = new TLegend(0.68, 0.74, 0.88, 0.88);
      legend->AddEntry(newHist, "current (version1)", "lep");
      legend->AddEntry(oldHist, "previous (version0)", "lep");
      legend->Draw();
      canvas.cd(2);
      ratio1D = static_cast<TH1 *>(newHist->Clone(("ratio_" + name).c_str()));
      ratio1D->SetDirectory(nullptr); ratio1D->Divide(oldHist);
      ratio1D->SetTitle("current / previous"); ratio1D->GetYaxis()->SetTitle("ratio");
      ratio1D->SetMinimum(0.5); ratio1D->SetMaximum(1.5); ratio1D->Draw("E");
      unity = new TLine(ratio1D->GetXaxis()->GetXmin(), 1., ratio1D->GetXaxis()->GetXmax(), 1.);
      unity->SetLineStyle(2); unity->Draw();
    } else if (newHist->GetDimension() == 2 && oldHist->GetDimension() == 2) {
      canvas.Divide(3, 1);
      canvas.cd(1); newHist->SetTitle((name + " (current)").c_str()); newHist->Draw("COLZ");
      canvas.cd(2); oldHist->SetTitle((name + " (previous)").c_str()); oldHist->Draw("COLZ");
      canvas.cd(3);
      ratio2D = static_cast<TH2 *>(newHist->Clone(("ratio_" + name).c_str()));
      ratio2D->SetDirectory(nullptr); ratio2D->Divide(oldHist);
      ratio2D->SetTitle((name + " (current / previous)").c_str()); ratio2D->Draw("COLZ");
    } else {
      canvas.cd();
      labelPad((name + " (different dimensions)").c_str());
    }

    canvas.Print((std::string(output) + (first ? "(" : "")).c_str());
    first = false;
    delete unity;
    delete legend;
    delete ratio1D;
    delete ratio2D;
  }
  canvas.Print((std::string(output) + ")").c_str());
  Info("compare_background_hists", "Wrote %s", output);
}
