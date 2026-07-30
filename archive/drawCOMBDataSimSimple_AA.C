#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"

#include "dlUtility.h"
#include "read_binning.h"

namespace
{
  std::set<std::string> parseCurves(const std::string &curveList)
  {
    std::set<std::string> result;
    std::stringstream stream(curveList);
    std::string curve;
    while (std::getline(stream, curve, ','))
    {
      const std::size_t first = curve.find_first_not_of(" \t");
      const std::size_t last = curve.find_last_not_of(" \t");
      if (first != std::string::npos){ result.insert(curve.substr(first, last - first + 1));}
    }
    return result;
  }
  constexpr double kFlowFitLow = 0.0;
  constexpr double kFlowFitHigh = 2.5;
  constexpr double kNormalizationLow = 0.8;
  constexpr double kNormalizationHigh = 2.5;
  struct Curve
  {
    std::unique_ptr<TH1D> histogram;
    std::string label;
    int color = kBlack;
    int style = 1;
    std::string selection;
  };
  struct Function
  {
    std::unique_ptr<TF1> function;
    std::string label;
    int color = kBlack;
    int style = 1;
    std::string selection;
  };
}

// Available curve names:
//   data_raw, data_background, data_nominal, data_down, data_up,
//   sim_raw,  sim_background,  sim_nominal,  sim_down,  sim_up
//
// scale_sim_to_data=true scales every displayed simulation component by the
// raw data/simulation integral ratio in 0.8 <= dphi < 2.5.  The cached weighted
// counts are never modified.  Set it false to compare absolute weighted counts.
void drawCOMBDataSimSimple_AA(
  const int cone_size = 3, const int centrality_bin = 1,
  const std::string curves = "data_raw,data_nominal,sim_raw,sim_nominal",
  const std::string funcs = "data_flow_fit_nominal,data_flow_fit_COMBDown,data_flow_fit_COMBUp",
  const bool scale_sim_to_data = true,
  const std::string config = "configs/binning_AA.config"
)
{
  read_binning rb(config);
  std::unique_ptr<float[]> ptBins(new float[rb.get_nbins() + 1]);
  rb.get_pt_bins(ptBins.get());

  const TString directory = Form("%s/dphi_plots", rb.get_code_location().c_str());
  const TString dataPath = Form("%s/dphi_COMB_modulation_AA_cent_%d_r%02d.root", directory.Data(), centrality_bin, cone_size);
  const TString simPath = Form( "%s/dphi_COMB_modulation_sim_combined_AA_cent_%d_r%02d.root", directory.Data(), centrality_bin, cone_size);
  std::unique_ptr<TFile> dataFile(TFile::Open(dataPath, "READ"));
  std::unique_ptr<TFile> simFile(TFile::Open(simPath, "READ"));

  if (!dataFile || dataFile->IsZombie() || !simFile || simFile->IsZombie())
  {
    std::cerr << "Cannot open data and simulation component caches" << std::endl;
    return;
  }

  const auto selected = parseCurves(curves);
  const std::set<std::string> valid = {
    "data_pairs", "data_eta", "data_fit", "data_background", "data_background_COMBDown", "data_background_COMBUp",
    "data_nominal", "data_down", "data_up",
    "data_signal_region", "data_normalization_region",
    "data_background_down_alt", "data_background_up_alt", "data_background_down_alt2", "data_background_up_alt2",
    "sim_pairs", "sim_eta", "sim_fit", "sim_background", "sim_background_COMBDown", "sim_background_COMBUp",
    "sim_nominal", "sim_down", "sim_up",
    "sim_signal_region", "sim_normalization_region",
    "sim_truth_pairs", "sim_truth_eta"
  };
  const auto selectedFuncs = parseCurves(funcs);
  const std::set<std::string> validFuncs = {
    "data_flow_fit_nominal", "data_flow_fit_COMBDown", "data_flow_fit_COMBUp",
    "sim_flow_fit_nominal", "sim_flow_fit_COMBDown", "sim_flow_fit_COMBUp",
    "data_flow_fit_down_alt", "data_flow_fit_up_alt", "data_flow_fit_down_alt2", "data_flow_fit_up_alt2"
  };
  
  for (const auto &name : selected)
  { 
    if (!valid.count(name))
    {
      std::cerr << "Unknown curve '" << name << "'. Available curves:" << std::endl;
      for (const auto &choice : valid) {std::cerr << "  " << choice << std::endl;}
      return;
    }
    // std::cout << "Selected curve: " << name << std::endl;
  }
  if (selected.empty())
  {
    std::cerr << "No curves selected" << std::endl;
    return;
  }
  // std::cout << "Selected curves size: " << selected.size() << std::endl;

  TH1D *dataRaw = dynamic_cast<TH1D*>(dataFile->Get("h_dphi_pairs"));
  TH1D *simRaw = dynamic_cast<TH1D*>(simFile->Get("h_dphi_pairs"));
  if (!dataRaw || !simRaw) return;
  const int firstNormBin = dataRaw->FindBin(0.8);
  const int lastNormBin = dataRaw->FindBin(std::nextafter(2.5, 0.8));
  const double simNorm = simRaw->Integral(firstNormBin, lastNormBin);
  const double simulationScale = scale_sim_to_data && simNorm > 0 ? dataRaw->Integral(firstNormBin, lastNormBin)/simNorm : 1.0;

  struct Definition
  {
    const char *selection;
    TFile *file;
    const char *object;
    const char *label;
    int color;
    int style;
    bool simulation;
  };

  // data
  // h_dphi_pairs = exclusive
  // h_dphi_eta_separated = exclusive w/ |deta| > 0.8
  // h_dphi_eta_separated_fit = exclusive w/ |deta| > 0.8, fit from 0-2.5
  // h_flow_background_nominal = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5
  // h_flow_background_COMBDown = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5, with COMBDown variation
  // h_flow_background_COMBUp = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5, with COMBUp variation
  // h_dphi_nominal_subtracted = exclusive pairs - h
  // h_dphi_COMBDown_subtracted = exclusive pairs - h_flow_background_COMBDown
  // h_dphi_COMBUp_subtracted = exclusive pairs - h_flow_background_COMBUp
  // h_signal_region = exclusive pairs in 0 <= dphi < 0.8
  // h_flow_normalization_region = exclusive pairs in 0.8 <= dphi < 2.5

  // sim 
  // h_dphi_pairs = exclusive
  // h_dphi_eta_separated = exclusive w/ |deta| > 0.8
  // h_dphi_eta_separated_fit = exclusive w/ |deta| > 0.8, fit from 0-2.5
  // h_flow_background_nominal = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5
  // h_flow_background_COMBDown = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5, with COMBDown variation
  // h_flow_background_COMBUp = fitted modulation from h_dphi_eta_separated_fit scaled to yield of exclusive pairs in 0.8 <= dphi < 2.5, with COMBUp variation
  // h_dphi_nominal_subtracted_sim = exclusive pairs - h
  // h_dphi_COMBDown_subtracted_sim = exclusive pairs - h_flow_background_COMBDown
  // h_dphi_COMBUp_subtracted_sim = exclusive pairs - h_flow_background_COMBUp
  // h_signal_region = exclusive pairs in 0 <= dphi < 0.8
  // h_flow_normalization_region = exclusive pairs in 0.8 <= dphi < 2.5
  // h_flow_pairs_truth = exclusive pairs from truth-level simulation
  // h_flow_pairs_eta_separated_truth = exclusive pairs from truth-level simulation w/ |deta| > 0.8



  const std::vector<Definition> definitions = { 
    {"data_pairs", dataFile.get(), "h_dphi_pairs", "Exclusive", kBlack, 1, false},
    {"data_eta", dataFile.get(), "h_dphi_eta_separated", "Inclusive |#Delta#eta| > 0.8", kAzure-2, 1, false},
    {"data_fit", dataFile.get(), "h_dphi_eta_separated_fit", "Inclusive |#Delta#eta| > 0.8", kAzure-2, 1, false},
    {"data_background", dataFile.get(), "h_flow_background_nominal", "Scaled fit", kRed, 1, false},
    {"data_background_COMBDown", dataFile.get(), "h_flow_background_COMBDown", "F(1.3v_{2,2},1.3v_{3,3})", kRed, 3, false},
    {"data_background_COMBUp", dataFile.get(), "h_flow_background_COMBUp", "F(0.7v_{2,2},0.7v_{3,3})", kRed, 2, false},
    {"data_nominal", dataFile.get(), "h_dphi_nominal_subtracted", "Data nominal subtracted", kBlue, 1, false},
    {"data_down", dataFile.get(), "h_dphi_COMBDown_subtracted", "Data COMBDown subtracted", kBlue, 2, false},
    {"data_up", dataFile.get(), "h_dphi_COMBUp_subtracted", "Data COMBUp subtracted", kBlue, 3, false},
    {"data_signal_region", dataFile.get(), "h_signal_region", "Data signal region (0 <= #Delta#phi < 0.8)", kGreen + 2, 1, false},
    {"data_normalization_region", dataFile.get(), "h_flow_normalization_region", "Data normalization region (0.8 <= #Delta#phi < 2.5)", kYellow + 2, 2, false},
    {"data_background_down_alt", dataFile.get(), "h_flow_background_COMBDown_alt", "F(1.5v_{2,2},0.5v_{3,3})", kRed+1, 3, false},
    {"data_background_up_alt", dataFile.get(), "h_flow_background_COMBUp_alt", "F(1.5v_{2,2},1.5v_{3,3})", kRed+2, 2, false},
    {"data_background_down_alt2", dataFile.get(), "h_flow_background_COMBDown_alt2", "F(0.5v_{2,2},1.5v_{3,3})", kRed+3, 3, false},
    {"data_background_up_alt2", dataFile.get(), "h_flow_background_COMBUp_alt2", "F(0.5v_{2,2},0.5v_{3,3})", kRed+4, 2, false},
    {"sim_pairs", simFile.get(), "h_dphi_pairs", "Exclusive", kBlack, 1, true},
    {"sim_eta", simFile.get(), "h_dphi_eta_separated", "Inclusive |#Delta#eta| > 0.8", kAzure-2, 2, true},
    {"sim_fit", simFile.get(), "h_dphi_eta_separated_fit", "Inclusive |#Delta#eta| > 0.8", kAzure-2, 3, true},
    {"sim_background", simFile.get(), "h_flow_background_nominal", "Simulation background (nominal)", kRed + 1, 1, true},
    {"sim_background_COMBDown", simFile.get(), "h_flow_background_COMBDown", "Simulation background (COMBDown)", kRed + 1, 2, true},
    {"sim_background_COMBUp", simFile.get(), "h_flow_background_COMBUp", "Simulation background (COMBUp)", kRed + 1, 3, true},
    {"sim_nominal", simFile.get(), "h_dphi_nominal_subtracted_sim", "Simulation nominal subtracted", kBlue + 1, 1, true},
    {"sim_down", simFile.get(), "h_dphi_down_subtracted_sim", "Simulation COMBDown subtracted", kBlue + 1, 2, true},
    {"sim_up", simFile.get(), "h_dphi_up_subtracted_sim", "Simulation COMBUp subtracted", kBlue + 1, 3, true},
    {"sim_signal_region", simFile.get(), "h_signal_region", "Simulation signal region (0 <= #Delta#phi < 0.8)", kGreen + 3, 1, true},
    {"sim_normalization_region", simFile.get(), "h_flow_normalization_region_sim", "Simulation normalization region (0.8 <= #Delta#phi < 2.5)", kGreen + 3, 2, true},
    {"sim_truth_pairs", simFile.get(), "h_dphi_pairs_truth", "Simulation truth raw", kBlack, 1, true},
    {"sim_truth_eta", simFile.get(), "h_dphi_eta_separated_truth", "Simulation truth |#Delta#eta| > 0.8", kBlack, 2, true}
  };


  const std::vector<Definition> f1_definitions = { 
    {"data_flow_fit", dataFile.get(), "global_flow_fit", "Data flow fit", kOrange + 7, 1, false},
    {"data_flow_fit_nominal", dataFile.get(), "f_flow_fit_0", "Data flow fit (nominal)", kBlack, 1, false},
    {"data_flow_fit_COMBDown", dataFile.get(), "f_flow_fit_1", "Data flow fit (COMBDown)", kBlack, 2, false},
    {"data_flow_fit_COMBUp", dataFile.get(), "f_flow_fit_2", "Data flow fit (COMBUp)", kBlack, 3, false},
    {"data_flow_fit_down_alt", dataFile.get(), "f_flow_fit_alt_1", "Data flow fit (COMBDown alt)", kBlack, 2, false},
    {"data_flow_fit_up_alt", dataFile.get(), "f_flow_fit_alt_0", "Data flow fit (COMBUp alt)", kBlack, 3, false},
    {"data_flow_fit_down_alt2", dataFile.get(), "f_flow_fit_alt_3", "Data flow fit (COMBDown alt2)", kBlack, 2, false},
    {"data_flow_fit_up_alt2", dataFile.get(), "f_flow_fit_alt_2", "Data flow fit (COMBUp alt2)", kBlack, 3, false},
    {"sim_flow_fit", simFile.get(), "f_flow_fit", "Simulation flow fit", kOrange + 7, 3, true},
    {"sim_flow_fit_nominal", simFile.get(), "f_flow_fit_0_sim", "Simulation flow fit (nominal)", kBlack, 1, true},
    {"sim_flow_fit_COMBDown", simFile.get(), "f_flow_fit_1_sim", "Simulation flow fit (COMBDown)", kBlack, 2, true},
    {"sim_flow_fit_COMBUp", simFile.get(), "f_flow_fit_2_sim", "Simulation flow fit (COMBUp)", kBlack, 3, true}
  };
  
  const std::vector<std::string> draw_as_lines = {
    "data_background", "data_background_COMBDown",
    "data_background_COMBUp",
    "data_fit",
    "data_background_down_alt", 
    "data_background_up_alt", 
    "data_background_down_alt2", 
    "data_background_up_alt2",
  };
  
  std::vector<Curve> displayed;
  std::vector<bool> isSimulation;
  std::vector<bool> isregion;
  double maximum = 0;
  double minimum = 0;
  for (const auto &definition : definitions)
  {
    if (!selected.count(definition.selection)) continue;
    TH1D *source = dynamic_cast<TH1D*>(definition.file->Get(definition.object));
    if (!source)
    {
      std::cerr << "Missing " << definition.object << std::endl;
      return;
    }
    // std::cout << "Selected curve: " << definition.selection << std::endl;
    Curve curve;
    curve.histogram.reset(static_cast<TH1D*>(source->Clone(Form(
      "h_overlay_%s", definition.selection))));
    curve.histogram->SetDirectory(nullptr);
    if (definition.simulation) curve.histogram->Scale(simulationScale);
    else curve.histogram->SetMarkerStyle(20);
    curve.selection = definition.selection;
    const std::string selectionName(definition.selection);
    if (selectionName == "data_signal_region" ||
        selectionName == "data_normalization_region")
      isregion.push_back(true);
    else
      isregion.push_back(false);

    if (definition.simulation ) isSimulation.push_back(true);
    else isSimulation.push_back(false);
    curve.label = definition.label;
    if (definition.simulation && scale_sim_to_data)
      curve.label += " (Sim)";
    curve.color = definition.color;
    curve.style = definition.style;
    curve.histogram->SetLineColor(curve.color);
    curve.histogram->SetLineWidth(2);
    curve.histogram->SetLineStyle(curve.style);
    curve.histogram->SetMarkerColor(curve.color);
    // maximum = std::max(maximum, curve.histogram->GetMaximum());
    if (maximum == 0) maximum = std::max(maximum, curve.histogram->GetMaximum());;
    minimum = std::min(minimum, curve.histogram->GetMinimum());
    displayed.push_back(std::move(curve));
  }
  // std::cout << "Displayed curves size: " << displayed.size() << std::endl;
  // for (const auto &curve : displayed)
  // {
  //   std::cout << "Curve: " << curve.selection << std::endl;
  // }

  // get the sim eta separted and data eta separated histograms to get the v22 and v33 values
  // auto data_eta =(TH1D*) dynamic_cast<TH1D*>(dataFile->Get("h_dphi_eta_separated")) -> Clone("data_eta_tmp");
  // auto sim_eta = (TH1D*)dynamic_cast<TH1D*>(simFile->Get("h_dphi_eta_separated")) -> Clone("sim_eta_tmp");
  // auto * fit = new TF1("fit", "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x))", 0, 2.5);
  // data_eta->Fit(fit, "0RlQ","", 0, 2.5);
  // double vv22_data = fit->GetParameter(1);
  // double vv33_data = fit->GetParameter(2);
  // double vv22_data_err = fit->GetParError(1);
  // double vv33_data_err = fit->GetParError(2);
  // sim_eta->Fit(fit, "0RlQ","", 0, 2.5);
  // double vv22_sim = fit->GetParameter(1);
  // double vv33_sim = fit->GetParameter(2);
  // double vv22_sim_err = fit->GetParError(1);
  // double vv33_sim_err = fit->GetParError(2);
  // std::cout << "Data v22 = " << vv22_data << " +/- " << vv22_data_err << std::endl;
  // std::cout << "Data v33 = " << vv33_data << " +/- " << vv33_data_err << std::endl;
  // std::cout << "Sim v22 = " << vv22_sim << " +/- " << vv22_sim_err << std::endl;
  // std::cout << "Sim v33 = " << vv33_sim << " +/- " << vv33_sim_err << std::endl;  
  // std::cout << "Data v22/v33 = " << vv22_data/vv22_sim << std::endl;
  // std::cout << "Sim v22/v33 = " << vv22_sim/vv33_sim << std::endl;
  // delete fit;

  std::vector<Function> displayedFuncs;
  for (const auto &definition : f1_definitions)
  {
    if (!selectedFuncs.count(definition.selection)) continue;
    TF1 *source = dynamic_cast<TF1*>(definition.file->Get(definition.object));
    if (!source)
      {
        std::cerr << "Missing " << definition.object << std::endl;
        return;
      }
    Function func;
    func.function.reset(static_cast<TF1*>(source->Clone(Form(
      "f_overlay_%s", definition.selection))));
    // func.function->SetDirectory(nullptr);
    func.selection = definition.selection;
    if (definition.simulation) func.function->SetParameter(0, func.function->GetParameter(0)*simulationScale);
    func.label = definition.label;
    func.color = definition.color;
    func.style = definition.style;
    func.function->SetLineColor(func.color);
    func.function->SetLineWidth(2);
    func.function->SetLineStyle(func.style);
    displayedFuncs.push_back(std::move(func));
  }


  std::vector<double> v22_data;
  std::vector<double> v33_data;
  std::vector<double> v22_data_err;
  std::vector<double> v33_data_err;
  std::vector<double> v22_sim;
  std::vector<double> v33_sim;
  std::vector<double> v22_sim_err;  
  std::vector<double> v33_sim_err;
  for (const auto &name : {"global_flow_fit"})
  {
    TF1 *fit = dynamic_cast<TF1*>(dataFile->Get(name));
    if (!fit)
      {
        std::cerr << "Missing " << name << std::endl;
        return;
      }
      v22_data.push_back(fit->GetParameter(1));
      v33_data.push_back(fit->GetParameter(2));
      v22_data_err.push_back(fit->GetParError(1));
      v33_data_err.push_back(fit->GetParError(2));
  }
  for (const auto &name : {"f_flow_fit" })
  {
    TF1 *fit = dynamic_cast<TF1*>(simFile->Get(name));
    if (!fit)
      {
        std::cerr << "Missing " << name << std::endl;
        return;
      }
      v22_sim.push_back(fit->GetParameter(1));
      v33_sim.push_back(fit->GetParameter(2));
      v22_sim_err.push_back(fit->GetParError(1));
      v33_sim_err.push_back(fit->GetParError(2));
  }

  std::cout << "Data v22 = " << v22_data[0] << " +/- " << v22_data_err[0] << std::endl;
  std::cout << "Data v33 = " << v33_data[0] << " +/- " << v33_data_err[0] << std::endl;
  std::cout << "Sim v22 = " << v22_sim[0] << " +/- " << v22_sim_err[0] << std::endl;
  std::cout << "Sim v33 = " << v33_sim[0] << " +/- " << v33_sim_err[0] << std::endl;  
  std::cout << "Data v22/v33 = " << v22_data[0]/v22_sim[0] << std::endl;
  std::cout << "Sim v22/v33 = " << v33_data[0]/v33_sim[0] << std::endl;  


  

  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  displayed.front().histogram->SetTitle(";#Delta#phi;Counts");
  displayed.front().histogram->SetMaximum(1.05*maximum);
  displayed.front().histogram->SetMinimum(std::min(-0.05*maximum, 1.15*minimum));
  TCanvas canvas("c_comb_data_sim_overlay", "c_comb_data_sim_overlay", 700, 600);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.04);
  canvas.SetBottomMargin(0.12);
  // for (std::size_t index = 0; index < displayed.size(); ++index)
  //   displayed[index].histogram->Draw(index == 0 ? "hist" : "hist same");
  

  for (std::size_t index = 0; index < displayed.size(); ++index)
  {
    // std::cout << "Drawing " << displayed[index].label << " as " << (isSimulation[index] ? "simulation" : "data") << std::endl;
    if (isregion[index])
    {
      // shadeen fill 
      displayed[index].histogram->SetFillColorAlpha(displayed[index].color, 0.12);
      displayed[index].histogram->SetLineColor(displayed[index].color);
      displayed[index].histogram->SetLineStyle(3);
      displayed[index].histogram->SetLineWidth(2);
      displayed[index].histogram->Draw("hist same");
      continue;
    }

    if (isSimulation[index])
    {
      // std::cout << "Drawing " << definitions[index].selection << " as simulation" << std::endl;
      displayed[index].histogram->Draw(index == 0 ? "hist" : "hist same");
    }
    else if (draw_as_lines.end() != std::find(draw_as_lines.begin(), draw_as_lines.end(), displayed[index].selection))
    { 
      // std::cout << "Drawing " << definitions[index].selection << " as line" << std::endl;
      displayed[index].histogram->Draw(index == 0 ? "hist" : "hist same");
    }
    else 
    {
      // std::cout << "Drawing " << definitions[index].selection << " as points" << std::endl;
      displayed[index].histogram->Draw(index == 0 ? "p" : "p same");
    }
  }
    // displayed[index].histogram->Draw(index == 0 ? "hist" : "hist same");
  for (std::size_t index = 0; index < displayedFuncs.size(); ++index)
    displayedFuncs[index].function->Draw("same");

  float centralityBins[5] = {0};
  rb.get_centrality_bins(centralityBins);
  dlutility::DrawSPHENIX(0.18, 0.91, 0.040);
  dlutility::drawText(Form("p_{T,1} > %.1f, p_{T,2} > %.1f GeV",
                           rb.get_reco_leading_cut(), rb.get_reco_subleading_cut()),
                      0.18, 0.79, 0, kBlack, 0.033);
  dlutility::drawText(Form("%.0f - %.0f %%", centralityBins[centrality_bin],
                           centralityBins[centrality_bin + 1]),

                      0.18, 0.73, 0, kBlack, 0.034);
  TLegend legend(0.55, 0.66, 0.9, 0.9);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.026);
    // get vn

  for (auto &curve : displayed)
    // legend.AddEntry(curve.histogram.get(), curve.label.c_str(), "l");
  {
    // std::cout << "Adding legend entry for " << curve.label << std::endl;
    if (isregion[&curve - &displayed[0]])
      continue;
    else if (isSimulation[&curve - &displayed[0]])
      legend.AddEntry(curve.histogram.get(), curve.label.c_str(), "l");
    else if (draw_as_lines.end() != std::find(draw_as_lines.begin(), draw_as_lines.end(), displayed[&curve - &displayed[0]].selection))
      legend.AddEntry(curve.histogram.get(), curve.label.c_str(), "l");
    else
      legend.AddEntry(curve.histogram.get(), curve.label.c_str(), "pe");
  }
  for (auto &func : displayedFuncs)
    legend.AddEntry(func.function.get(), func.label.c_str(), "l");
  legend.Draw();

  dlutility::drawText(Form("Data v_{2,2} = %.3f , v_{3,3} = %.3f", v22_data[0], v33_data[0]), 0.18, 0.66, 0, kBlack, 0.028);
  dlutility::drawText(Form("Sim v_{2,2} = %.3f , v_{3,3} = %.3f", v22_sim[0], v33_sim[0]), 0.18, 0.62, 0, kBlack, 0.028);

  const TString output = Form("%s/dphi_COMB_data_sim_overlay_AA_cent_%d_r%02d.pdf",
                              directory.Data(), centrality_bin, cone_size);
  canvas.SaveAs(output);
  std::cout << "Wrote " << output << " using simulation scale "
            << simulationScale << std::endl;
}
