#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TTree.h"

#include "read_binning.h"

namespace
{
struct Jet
{
  int id = -1;
  float pt = 0.0F;
  float eta = 0.0F;
  float phi = 0.0F;
};

using Match = std::pair<Jet, Jet>;

float deltaPhi(const float first, const float second)
{
  float difference = first - second;
  while (difference > TMath::Pi()) difference -= 2.0F*TMath::Pi();
  while (difference < -TMath::Pi()) difference += 2.0F*TMath::Pi();
  return std::fabs(difference);
}

float deltaR(const Jet& first, const Jet& second)
{
  return std::hypot(first.eta - second.eta,
                    deltaPhi(first.phi, second.phi));
}

// This intentionally preserves the previous analysis matching rule: truth
// jets and reco jets are pT ordered, and the first unused reco jet inside the
// legacy dR threshold is selected.  For cone_size=3 the historical threshold
// is dR < 3, not dR < 0.3.
std::vector<Match> matchJets(const std::vector<Jet>& truth,
                             const std::vector<Jet>& reco,
                             const float maximumDeltaR)
{
  std::vector<Match> result;
  std::vector<bool> used(reco.size(), false);
  for (const auto& truthJet : truth)
    {
      for (std::size_t recoIndex = 0; recoIndex < reco.size(); ++recoIndex)
        {
          if (used[recoIndex]) continue;
          if (deltaR(truthJet, reco[recoIndex]) >= maximumDeltaR) continue;
          result.emplace_back(truthJet, reco[recoIndex]);
          used[recoIndex] = true;
          break;
        }
    }
  return result;
}

const Match* findTruthMatch(const std::vector<Match>& matches, const int truthId)
{
  const auto found = std::find_if(matches.begin(), matches.end(),
    [truthId](const Match& match) { return match.first.id == truthId; });
  return found == matches.end() ? nullptr : &*found;
}

struct Output
{
  std::unique_ptr<TFile> file;
  std::unique_ptr<TNtuple> matches;
  std::unique_ptr<TNtuple> stats;
  std::array<std::unique_ptr<TEfficiency>, 10> truthEfficiency;
  std::array<std::unique_ptr<TEfficiency>, 10> dijetEfficiency;
  std::array<std::unique_ptr<TH2D>, 10> response;
  std::array<std::unique_ptr<TProfile>, 10> responseProfile;
};
}

void makeMatchedTreesInclusiveAuAu(
  const int cone_size = 3,
  const std::string input_path =
    "/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root",
  const std::string output_directory = "",
  const std::string configfile = "binning_AA.config",
  const std::string suffix = "ProdA_2024-00000030_sumeT")
{
  if (cone_size != 3)
    {
      std::cerr << "The supplied simulation contains R=0.3 jets only." << std::endl;
      return;
    }

  read_binning rb(configfile.c_str());
  const std::string outputDirectory = output_directory.empty()
    ? rb.get_tntuple_location() : output_directory;
  std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
  TTree* tree = input && !input->IsZombie()
    ? dynamic_cast<TTree*>(input->Get("T")) : nullptr;
  if (!tree)
    {
      std::cerr << "Cannot read tree T from " << input_path << std::endl;
      return;
    }

  const std::vector<std::string> required = {
    "jetid", "zvrtx", "cent", "sumeT", "jet_pT", "jet_E",
    "jet_unsub_E", "jet_eta", "jet_phi", "truth_jet_pT",
    "truth_jet_eta", "truth_jet_phi"};
  for (const auto& branch : required)
    if (!tree->GetBranch(branch.c_str()))
      {
        std::cerr << "Missing branch " << branch << " in " << input_path
                  << std::endl;
        return;
      }

  tree->SetBranchStatus("*", 0);
  for (const auto& branch : required) tree->SetBranchStatus(branch.c_str(), 1);
  tree->SetCacheSize(256*1024*1024);
  tree->AddBranchToCache("*", true);

  int jetId = 0;
  int centrality = 0;
  float vertex = 0.0F;
  float sumEt = 0.0F;
  std::vector<float>* truthPt = nullptr;
  std::vector<float>* truthEta = nullptr;
  std::vector<float>* truthPhi = nullptr;
  std::vector<float>* recoPt = nullptr;
  std::vector<float>* recoE = nullptr;
  std::vector<float>* recoUnsubE = nullptr;
  std::vector<float>* recoEta = nullptr;
  std::vector<float>* recoPhi = nullptr;
  tree->SetBranchAddress("jetid", &jetId);
  tree->SetBranchAddress("zvrtx", &vertex);
  tree->SetBranchAddress("cent", &centrality);
  tree->SetBranchAddress("sumeT", &sumEt);
  tree->SetBranchAddress("truth_jet_pT", &truthPt);
  tree->SetBranchAddress("truth_jet_eta", &truthEta);
  tree->SetBranchAddress("truth_jet_phi", &truthPhi);
  tree->SetBranchAddress("jet_pT", &recoPt);
  tree->SetBranchAddress("jet_E", &recoE);
  tree->SetBranchAddress("jet_unsub_E", &recoUnsubE);
  tree->SetBranchAddress("jet_eta", &recoEta);
  tree->SetBranchAddress("jet_phi", &recoPhi);

  const std::array<int, 3> sampleIds = {10, 20, 30};
  std::array<Output, 3> outputs;
  for (std::size_t sample = 0; sample < sampleIds.size(); ++sample)
    {
      // const std::string path = Form(
        // "%s/TREE_MATCH_r0%d_v15_%d_new_ProdA_2024-00000030_sumeT.root",
        // outputDirectory.c_str(), cone_size, sampleIds[sample]);
      const std::string path = Form(
        "%s/TREE_MATCH_r0%d_%s_%d.root",
        outputDirectory.c_str(), cone_size, suffix.c_str(), sampleIds[sample]);
      outputs[sample].file.reset(TFile::Open(path.c_str(), "RECREATE"));
      if (!outputs[sample].file || outputs[sample].file->IsZombie())
        {
          std::cerr << "Cannot create " << path << std::endl;
          return;
        }
      outputs[sample].matches = std::make_unique<TNtuple>(
        "tn_match", "matched truth and reco",
        "maxpttruth:pt1_truth:pt2_truth:dphi_truth:pt1_reco:pt2_reco:"
        "dphi_reco:matched:nrecojets:centrality:mbd_vertex:sumeT");
      outputs[sample].stats = std::make_unique<TNtuple>(
        "tn_stats", "run stats", "nevents");
      outputs[sample].matches->SetDirectory(nullptr);
      outputs[sample].stats->SetDirectory(nullptr);
      for (int centralityBin = 0; centralityBin < 10; ++centralityBin)
        {
          outputs[sample].truthEfficiency[centralityBin] =
            std::make_unique<TEfficiency>(
              Form("he_pt_truth_%d", centralityBin), "", 25, 0, 50);
          outputs[sample].dijetEfficiency[centralityBin] =
            std::make_unique<TEfficiency>(
              Form("he_pt_dijet_%d", centralityBin), "", 25, 0, 50);
          outputs[sample].response[centralityBin] = std::make_unique<TH2D>(
            Form("h2_ptt_ptrptt_%d", centralityBin),
            ";p_{T}^{truth};p_{T}^{reco}/p_{T}^{truth}",
            25, 0, 50, 120, 0, 1.2);
          outputs[sample].responseProfile[centralityBin] =
            std::make_unique<TProfile>(
              Form("hp_ptt_ptrptt_%d", centralityBin), "", 25, 0, 50, "s");
          outputs[sample].truthEfficiency[centralityBin]->SetDirectory(nullptr);
          outputs[sample].dijetEfficiency[centralityBin]->SetDirectory(nullptr);
          outputs[sample].response[centralityBin]->SetDirectory(nullptr);
          outputs[sample].responseProfile[centralityBin]->SetDirectory(nullptr);
        }
    }

  Long64_t perSampleLimit = 0;
  if (const char* value = std::getenv("DIJET_MAX_SIM_EVENTS"))
    perSampleLimit = std::max(0LL, std::atoll(value));
  std::array<Long64_t, 3> processed = {0, 0, 0};
  std::array<Long64_t, 3> written = {0, 0, 0};
  std::vector<Jet> truthJets;
  std::vector<Jet> recoJets;
  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
    {
      tree->GetEntry(entry);
      const int sample = jetId == 10 ? 0 : jetId == 20 ? 1 : jetId == 30 ? 2 : -1;
      if (sample < 0) continue;
      if (perSampleLimit > 0 && processed[sample] >= perSampleLimit) continue;
      ++processed[sample];
      if (entry % 100000 == 0)
        std::cout << "Simulation event " << entry << " / " << entries << "\r"
                  << std::flush;
      if (std::fabs(vertex) > 60.0F) continue;
      if (!truthPt || !truthEta || !truthPhi || !recoPt || !recoE ||
          !recoUnsubE || !recoEta || !recoPhi) continue;
      if (truthPt->size() != truthEta->size() ||
          truthPt->size() != truthPhi->size() ||
          recoPt->size() != recoE->size() ||
          recoPt->size() != recoUnsubE->size() ||
          recoPt->size() != recoEta->size() ||
          recoPt->size() != recoPhi->size() || truthPt->empty()) continue;

      const float maximumTruth = *std::max_element(truthPt->begin(), truthPt->end());
      truthJets.clear();
      recoJets.clear();
      for (std::size_t jet = 0; jet < truthPt->size(); ++jet)
        if (truthPt->at(jet) >= 3.0F && std::fabs(truthEta->at(jet)) <= 0.8F)
          truthJets.push_back({static_cast<int>(jet), truthPt->at(jet),
                               truthEta->at(jet), truthPhi->at(jet)});
      for (std::size_t jet = 0; jet < recoPt->size(); ++jet)
        if (recoPt->at(jet) >= 3.0F && recoE->at(jet) >= 0.0F &&
            recoUnsubE->at(jet) >= 0.0F && std::fabs(recoEta->at(jet)) <= 0.8F)
          recoJets.push_back({static_cast<int>(jet), recoPt->at(jet),
                              recoEta->at(jet), recoPhi->at(jet)});
      const auto descendingPt = [](const Jet& first, const Jet& second)
        { return first.pt > second.pt; };
      std::sort(truthJets.begin(), truthJets.end(), descendingPt);
      std::sort(recoJets.begin(), recoJets.end(), descendingPt);
      if (truthJets.size() < 2) continue;

      const auto matches = matchJets(
        truthJets, recoJets, static_cast<float>(cone_size));
      const int centralityBin = std::clamp(centrality/10, 0, 9);
      for (const auto& truthJet : truthJets)
        outputs[sample].truthEfficiency[centralityBin]->Fill(
          findTruthMatch(matches, truthJet.id) != nullptr, truthJet.pt);
      for (const auto& match : matches)
        if (match.first.pt > 0.0F)
          {
            outputs[sample].response[centralityBin]->Fill(
              match.first.pt, match.second.pt/match.first.pt);
            outputs[sample].responseProfile[centralityBin]->Fill(
              match.first.pt, match.second.pt/match.first.pt);
          }

      const Match* leading = findTruthMatch(matches, truthJets[0].id);
      const Match* subleading = findTruthMatch(matches, truthJets[1].id);
      const bool matched = leading && subleading;
      outputs[sample].dijetEfficiency[centralityBin]->Fill(
        matched, truthJets[0].pt);
      outputs[sample].matches->Fill(
        maximumTruth, truthJets[0].pt, truthJets[1].pt,
        deltaPhi(truthJets[0].phi, truthJets[1].phi),
        matched ? leading->second.pt : 0.0F,
        matched ? subleading->second.pt : 0.0F,
        matched ? deltaPhi(leading->second.phi, subleading->second.phi) : 0.0F,
        matched ? 1.0F : 0.0F,
        matched ? static_cast<float>(recoJets.size()) : 0.0F,
        static_cast<float>(centrality), vertex, sumEt);
      ++written[sample];
    }
  std::cout << std::endl;

  for (std::size_t sample = 0; sample < outputs.size(); ++sample)
    {
      outputs[sample].file->cd();
      outputs[sample].stats->Fill(static_cast<float>(processed[sample]));
      outputs[sample].matches->Write();
      outputs[sample].stats->Write();
      for (int centralityBin = 0; centralityBin < 10; ++centralityBin)
        {
          outputs[sample].truthEfficiency[centralityBin]->Write();
          outputs[sample].dijetEfficiency[centralityBin]->Write();
          outputs[sample].response[centralityBin]->Write();
          outputs[sample].responseProfile[centralityBin]->Write();
        }
      std::cout << "jetid " << sampleIds[sample] << ": processed "
                << processed[sample] << ", wrote " << written[sample]
                << " tn_match rows to " << outputs[sample].file->GetName()
                << std::endl;
      outputs[sample].file->Close();
    }
}
