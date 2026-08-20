#include <algorithm>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"

#include "dlUtility.h"
#include "read_binning.h"

namespace
{
std::unique_ptr<TH1D> normalizedClone(TFile *file, const char *name,
                                      const char *newName)
{
  TH1D *source = file ? dynamic_cast<TH1D*>(file->Get(name)) : nullptr;
  if (!source || source->Integral() <= 0) return nullptr;
  auto histogram = std::unique_ptr<TH1D>(static_cast<TH1D*>(source->Clone(newName)));
  histogram->SetDirectory(nullptr);
  histogram->Scale(1.0/histogram->Integral(), "width");
  return histogram;
}

void drawValidation(const std::string &quantity, const std::string &axisTitle,
                    TH1D *data, TH1D *before, TH1D *after,
                    const std::string &system, const int coneSize,
                    const std::string &basePath)
{
  TCanvas *canvas = new TCanvas(Form("c_reweight_%s", quantity.c_str()), "", 650, 750);
  dlutility::ratioPanelCanvas(canvas, 0.4);
  canvas->cd(1);
  data->SetTitle(Form(";%s;normalized counts", axisTitle.c_str()));
  data->SetMarkerStyle(20);
  data->SetMarkerColor(kBlack);
  data->SetLineColor(kBlack);
  before->SetLineColor(kRed + 1);
  before->SetMarkerColor(kRed + 1);
  before->SetMarkerStyle(24);
  after->SetLineColor(kBlue + 1);
  after->SetMarkerColor(kBlue + 1);
  after->SetMarkerStyle(25);
  data->SetMaximum(1.25*std::max({data->GetMaximum(), before->GetMaximum(),
                                  after->GetMaximum()}));
  data->Draw("E1");
  before->Draw("same E1");
  after->Draw("same E1");
  TLegend legend(0.52, 0.66, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.AddEntry(data, "Au+Au data", "pl");
  legend.AddEntry(before, "simulation before weights", "pl");
  legend.AddEntry(after, "simulation after weights", "pl");
  legend.Draw();

  canvas->cd(2);
  std::unique_ptr<TH1D> beforeRatio(static_cast<TH1D*>(before->Clone(
    Form("h_%s_before_over_data", quantity.c_str()))));
  std::unique_ptr<TH1D> afterRatio(static_cast<TH1D*>(after->Clone(
    Form("h_%s_after_over_data", quantity.c_str()))));
  beforeRatio->SetDirectory(nullptr);
  afterRatio->SetDirectory(nullptr);
  beforeRatio->Divide(data);
  afterRatio->Divide(data);
  beforeRatio->SetTitle(Form(";%s;simulation / data", axisTitle.c_str()));
  beforeRatio->SetMinimum(0.0);
  beforeRatio->SetMaximum(2.0);
  beforeRatio->Draw("E1");
  afterRatio->Draw("same E1");
  TLine unity(data->GetXaxis()->GetXmin(), 1.0,
              data->GetXaxis()->GetXmax(), 1.0);
  unity.SetLineStyle(2);
  unity.Draw("same");
  canvas->Print(Form("%s/unfolding_plots/reweight_validation_%s_%s_r%02d.png",
                     basePath.c_str(), quantity.c_str(), system.c_str(), coneSize));
  canvas->Print(Form("%s/unfolding_plots/reweight_validation_%s_%s_r%02d.pdf",
                     basePath.c_str(), quantity.c_str(), system.c_str(), coneSize));
  delete canvas;
}
}

void validateReweighting_AA(const int cone_size = 3,
                            const int centrality_bin = 0,
                            const std::string configfile = "binning_AA.config")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  read_binning rb(configfile);
  const std::string system = "AA_cent_" + std::to_string(centrality_bin);
  const TString beforePath = Form(
    "%s/response_matrices/response_matrix_%s_r%02d_PRIMER1_nominal.root",
    rb.get_code_location().c_str(), system.c_str(), cone_size);
  const TString afterPath = Form(
    "%s/response_matrices/response_matrix_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.c_str(), cone_size);
  const TString dataPath = Form(
    "%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root",
    rb.get_code_location().c_str(), system.c_str(), cone_size);
  std::unique_ptr<TFile> beforeFile(TFile::Open(beforePath, "READ"));
  std::unique_ptr<TFile> afterFile(TFile::Open(afterPath, "READ"));
  std::unique_ptr<TFile> dataFile(TFile::Open(dataPath, "READ"));
  if (!beforeFile || beforeFile->IsZombie() || !afterFile || afterFile->IsZombie() ||
      !dataFile || dataFile->IsZombie())
    {
      std::cerr << "Cannot open before/after/data reweighting validation inputs for "
                << system << std::endl;
      return;
    }

  struct Quantity { const char *name; const char *axis; };
  const Quantity quantities[] = {
    {"h_mbd_vertex", "z_{vtx} [cm]"},
    {"h_centrality", "centrality [%]"},
    {"h_sumeT", "#Sigma E_{T}"}
  };
  std::unique_ptr<TFile> output(TFile::Open(Form(
    "%s/unfolding_hists/reweight_validation_%s_r%02d.root",
    rb.get_code_location().c_str(), system.c_str(), cone_size), "RECREATE"));
  for (const auto &quantity : quantities)
    {
      const std::string shortName = quantity.name == std::string("h_mbd_vertex")
        ? "vertex" : quantity.name == std::string("h_centrality")
        ? "centrality" : "sumeT";
      auto data = normalizedClone(dataFile.get(), quantity.name,
                                  Form("h_%s_data", shortName.c_str()));
      auto before = normalizedClone(beforeFile.get(), quantity.name,
                                    Form("h_%s_sim_before", shortName.c_str()));
      auto after = normalizedClone(afterFile.get(), quantity.name,
                                   Form("h_%s_sim_after", shortName.c_str()));
      if (!data || !before || !after)
        {
          std::cerr << "Missing or empty " << quantity.name
                    << " in reweighting validation inputs for " << system << std::endl;
          continue;
        }
      drawValidation(shortName, quantity.axis, data.get(), before.get(), after.get(),
                     system, cone_size, rb.get_code_location());
      output->cd();
      data->Write();
      before->Write();
      after->Write();
    }
  output->Write();
}
