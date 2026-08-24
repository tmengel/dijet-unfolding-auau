
#ifndef _MAKEMATCHEDTREESTAGGEDAUAU_C_
#define _MAKEMATCHEDTREESTAGGEDAUAU_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

namespace
{

struct DijetCand
{
  int id = -1;
  float pt = 0.0F;
  float eta = 0.0F;
  float phi = 0.0F;
  int flavor = -999; // truth jets only
};

using Match = std::pair<DijetCand, DijetCand>;

// greedy nearest-neighbor matching: truth jets in descending-pT order each
// claim the closest not-yet-used reco jet inside maxDr, if any.
std::vector<Match> matchJets(const std::vector<DijetCand>& truth, const std::vector<DijetCand>& reco, const float maxDr)
{
  std::vector<Match> result;
  std::vector<bool> used(reco.size(), false);
  for (const auto& truthJet : truth)
  {
    int best = -1;
    float bestDr = maxDr;
    for (std::size_t i = 0; i < reco.size(); ++i)
    {
      if (used[i]) continue;
      const float dr = AnaUtils::calc_dr(truthJet.eta, truthJet.phi, reco[i].eta, reco[i].phi);
      if (dr < bestDr) { bestDr = dr; best = static_cast<int>(i); }
    }
    if (best >= 0) { result.emplace_back(truthJet, reco[best]); used[best] = true; }
  }
  return result;
}

const Match* findTruthMatch(const std::vector<Match>& matches, const int truthId)
{
  const auto found = std::find_if(matches.begin(), matches.end(), [truthId](const Match& match) { return match.first.id == truthId; });
  return found == matches.end() ? nullptr : &*found;
}

enum FlavorCategory { 
  kUnmatchedFlavor = 0, 
  kQQ = 1, 
  kQG = 2, 
  kGG = 3, 
  kOtherFlavor = 4 
};

FlavorCategory categorizeFlavor(const int flavor1, const int flavor2)
{
  const bool isQuark1 = (flavor1 >= 1 && flavor1 <= 6);
  const bool isQuark2 = (flavor2 >= 1 && flavor2 <= 6);
  const bool isGluon1 = (flavor1 == 21);
  const bool isGluon2 = (flavor2 == 21);
  if (isQuark1 && isQuark2) return kQQ;
  if ((isQuark1 && isGluon2) || (isGluon1 && isQuark2)) return kQG;
  if (isGluon1 && isGluon2) return kGG;
  if (flavor1 <= 0 || flavor2 <= 0) return kUnmatchedFlavor;
  return kOtherFlavor;
}

// meaning of sublead_flag in the output ntuple, filled only when the
// fiducial-selected truth subleading jet has no matched reco jet:
//   kNoRecoCandidate      - the event has no second reco jet inside the
//                           fiducial cut at all.
//   kRealElsewhere        - the true physical subleading jet (the 2nd
//                           highest-pT truth jet in the *full*, unfiltered
//                           truth jet list -- acceptance ignored) is not
//                           what the observed 2nd-highest-pT reco jet
//                           reconstructs: either that reco jet matches a
//                           different truth jet (its real reco partner got
//                           claimed elsewhere by the greedy matcher, e.g. by
//                           the leading jet) or matches nothing at all.
//   kTruthOutsideFiducial - the observed 2nd-highest-pT reco jet DOES
//                           dR-match the true physical subleading jet, but
//                           that truth jet itself fails the fiducial
//                           pt/|eta| cut, so it was never a candidate for
//                           "subleading" in the fiducial-selected truth list.
enum SubleadFlag { 
  kSubleadMatched = -1, 
  kNoRecoCandidate = 0, 
  kRealElsewhere = 1, 
  kTruthOutsideFiducial = 2 
};

struct Category
{
  std::unique_ptr<TFile> file;
  std::unique_ptr<TNtuple> matches;
  std::unique_ptr<TNtuple> stats;
};

} // namespace

// Reads the raw per-event simulation trees (AnaTreev1 "CALO_tree" output,
// listed by a catalog .list file, e.g.
// /sphenix/user/tmengel/dijet-ana-auau/catalog/jet30_hijing_scaled/jet30_hijing_scaled-000.list)
// which carry the *full*, unfiltered truth jet list per event (unlike the
// flat trees written by match.C, which only retain truth jets that already
// passed a selection cut). For each event, ranks the fiducial-selected truth
// jets (pt >= min_pt, |eta| <= fiducial_eta) by pT to define the leading and
// subleading truth jet, matches them to the fiducial-selected jets of one
// reco collection (reco_type_name: "rho" or "sub1") via greedy nearest-
// neighbor dR matching, and writes the result to four TNtuple output files:
// an inclusive one with every leading+subleading truth dijet, and qq/qg/gg
// ones restricted to dijets whose leading+subleading truth jets both carry a
// hard-parton flavor code matching that tag. Dijets with an unmatched-to-
// parton or otherwise uncategorizable flavor combination go to the
// inclusive file only. See the SubleadFlag comment above for the meaning of
// the extra sublead_flag column.
void makeMatchedTreesTaggedAuAu(
  const std::string infile = "/sphenix/user/tmengel/dijet-ana-auau/catalog/jet30_hijing_scaled/jet30_hijing_scaled-000.list",
  const std::string output_directory = "",
  const std::string suffix = "08_23_2026_v001",
  const std::string reco_type_name = "rho",
  const float min_pt = 3.0F,
  const float fiducial_eta = 0.8F,
  const float max_vertex = 60.0F,
  const float max_dr_factor = 0.75F
)
{
  
  if (reco_type_name != "rho" && reco_type_name != "sub1")
  {
    std::cerr << "reco_type_name must be \"rho\" or \"sub1\", got \"" << reco_type_name << "\"." << std::endl;
    return;
  }
  const bool useRho = reco_type_name == "rho";

  auto* chain = new TChain("T");
  const auto files = AnaUtils::getFilelist(infile, ".root");
  for (const auto& file : files) chain->Add(file.c_str());

  const Long64_t nentries = chain->GetEntries();
  if (nentries < 1)
  {
    std::cerr << "Error: no entries in TChain built from " << infile << std::endl;
    return;
  }
  std::cout << "Total entries: " << nentries << " (" << files.size() << " files)" << std::endl;

  chain->SetBranchStatus("*", false);

  int eventId = 0;
  float zvrtx = 0.0F;
  int centrality = -1;
  float psi2 = 0.0F;
  float psi3 = 0.0F;

  static const int NETA = 24;
  static const int NPHI = 64;
  static float cemcTowerE[NETA][NPHI];
  static int cemcTowerIsGood[NETA][NPHI];
  static float ihcalTowerE[NETA][NPHI];
  static int ihcalTowerIsGood[NETA][NPHI];
  static float ohcalTowerE[NETA][NPHI];
  static int ohcalTowerIsGood[NETA][NPHI];

  float truthJetR = 0.0F;
  float truthJetMaxR04pT = 0.0F;
  std::vector<float>* truthPt = nullptr;
  std::vector<float>* truthEta = nullptr;
  std::vector<float>* truthPhi = nullptr;
  std::vector<int>* truthFlavor = nullptr;

  float recoJetR = 0.0F;
  std::vector<float>* recoPt = nullptr;
  std::vector<float>* recoEta = nullptr;
  std::vector<float>* recoPhi = nullptr;
  std::vector<float>* recoE = nullptr;

  const auto enable = [&](const char* name, void* addr)
  {
    chain->SetBranchStatus(name, true);
    chain->SetBranchAddress(name, addr);
  };

  enable("event_id", &eventId);
  enable("zvrtx", &zvrtx);
  enable("cent", &centrality);
  enable("psi2", &psi2);
  enable("psi3", &psi3);

  enable("cemc_tower_E", cemcTowerE);
  enable("cemc_tower_isgood", cemcTowerIsGood);
  enable("ihcal_tower_E", ihcalTowerE);
  enable("ihcal_tower_isgood", ihcalTowerIsGood);
  enable("ohcal_tower_E", ohcalTowerE);
  enable("ohcal_tower_isgood", ohcalTowerIsGood);

  enable("truth_jet_maxpT_r04", &truthJetMaxR04pT);
  enable("truth_jet_R", &truthJetR);
  enable("truth_jet_pT", &truthPt);
  enable("truth_jet_eta", &truthEta);
  enable("truth_jet_phi", &truthPhi);
  enable("truth_jet_flavor", &truthFlavor);

  if (useRho)
  {
    enable("rho_jet_R", &recoJetR);
    enable("rho_jet_pT", &recoPt);
    enable("rho_jet_eta", &recoEta);
    enable("rho_jet_phi", &recoPhi);
    enable("rho_jet_E", &recoE);
  }
  else
  {
    enable("sub1_jet_R", &recoJetR);
    enable("sub1_jet_pT", &recoPt);
    enable("sub1_jet_eta", &recoEta);
    enable("sub1_jet_phi", &recoPhi);
    enable("sub1_jet_E", &recoE);
  }

  std::string base = infile;
  const auto slash = base.find_last_of('/');
  if (slash != std::string::npos) base = base.substr(slash + 1);
  const auto dot = base.rfind(".list");
  if (dot != std::string::npos) base = base.substr(0, dot);

  const std::array<std::string, 4> tagNames = {"inclusive", "qq", "qg", "gg"};
  std::array<Category, 4> categories;
  for (std::size_t tag = 0; tag < tagNames.size(); ++tag)
  {
    const std::string path = Form("%s/TNTUPLE_DIJET_%s_%s_%s.root", output_directory.c_str(), tagNames[tag].c_str(), base.c_str(), suffix.c_str());
    categories[tag].file.reset(TFile::Open(path.c_str(), "RECREATE"));
    if (!categories[tag].file || categories[tag].file->IsZombie())
    {
      std::cerr << "Cannot create " << path << std::endl;
      return;
    }
    
    categories[tag].matches = std::make_unique<TNtuple>(
        "tn_match", "matched truth and reco",
        "maxptruthr04:maxpttruth:pt1_truth:pt2_truth:dphi_truth:pt1_reco:pt2_reco:"
        "dphi_reco:matched:nrecojets:centrality:mbd_vertex:sumeT:psi2:psi3:sublead_flag");
      categories[tag].stats = std::make_unique<TNtuple>("tn_stats", "run stats", "nevents:ndijets");
      categories[tag].matches->SetDirectory(nullptr);
      categories[tag].stats->SetDirectory(nullptr);
    }

    Long64_t nEventsSeen = 0;
    Long64_t nEventsWithDijet = 0;
    std::array<Long64_t, 4> nDijetsWritten = {0, 0, 0, 0};

    std::vector<DijetCand> fullTruth;
    std::vector<DijetCand> truthSel;
    std::vector<DijetCand> recoSel;

    for (Long64_t entry = 0; entry < nentries; ++entry)
    {
      chain->GetEntry(entry);
      if (nentries >= 10 && entry % (nentries/10) == 0 && entry > 0) std::cout << "Entry " << entry << " / " << nentries << "\r" << std::flush;

      ++nEventsSeen;
      if (std::fabs(zvrtx) > max_vertex) continue;
      if (!truthPt || !truthEta || !truthPhi || !truthFlavor || !recoPt || !recoEta || !recoPhi || !recoE) continue;
      if (truthPt->size() != truthEta->size() 
        || truthPt->size() != truthPhi->size() 
        || truthPt->size() != truthFlavor->size() 
        || recoPt->size() != recoEta->size() 
        || recoPt->size() != recoPhi->size() 
        || recoPt->size() != recoE->size()
      ){ continue; }

      fullTruth.clear();
      for (std::size_t i = 0; i < truthPt->size(); ++i)
      {
        fullTruth.push_back( {static_cast<int>(i), truthPt->at(i), truthEta->at(i), truthPhi->at(i), truthFlavor->at(i)} );
      }
      if (fullTruth.size() < 2) continue;
      
      const auto descendingPt = [](const DijetCand& a, const DijetCand& b) { return a.pt > b.pt; };
      std::sort(fullTruth.begin(), fullTruth.end(), descendingPt);
      const DijetCand& realLeading = fullTruth[0];
      const DijetCand& realSubleading = fullTruth[1];

      truthSel.clear();
      for (const auto& jet : fullTruth)
      {  
        if (jet.pt >= min_pt && std::fabs(jet.eta) <= fiducial_eta) truthSel.push_back(jet);
      }
      if (truthSel.size() < 2) continue;

      recoSel.clear();
      for (std::size_t i = 0; i < recoPt->size(); ++i)
      {
        if ( 
          recoPt->at(i) >= min_pt 
          && recoE->at(i) >= 0.0F
          && std::fabs(recoEta->at(i)) <= fiducial_eta
        )
        { 
          recoSel.push_back({static_cast<int>(i), recoPt->at(i), recoEta->at(i), recoPhi->at(i)});
        }
      }  
      std::sort(recoSel.begin(), recoSel.end(), descendingPt);

      ++nEventsWithDijet;

      const float sumEt = AnaUtils::calc_sumeT(AnaUtils::CEMC, zvrtx, cemcTowerE, cemcTowerIsGood)
        + AnaUtils::calc_sumeT(AnaUtils::HCALIN, zvrtx, ihcalTowerE, ihcalTowerIsGood)
        + AnaUtils::calc_sumeT(AnaUtils::HCALOUT, zvrtx, ohcalTowerE, ohcalTowerIsGood);

      const float maxDr = max_dr_factor * recoJetR;
      const auto matches = matchJets(truthSel, recoSel, maxDr);

      const DijetCand& lead = truthSel[0];
      const DijetCand& sub = truthSel[1];
      const Match* leadMatch = findTruthMatch(matches, lead.id);
      const Match* subMatch = findTruthMatch(matches, sub.id);
      const bool dijetMatched = leadMatch && subMatch;

      const float pt1Reco = leadMatch ? leadMatch->second.pt : 0.0F;
      const float pt2Reco = subMatch ? subMatch->second.pt : 0.0F;
      const float dphiTruth = AnaUtils::dphi_wrap(lead.phi, sub.phi);
      const float dphiReco = dijetMatched ? AnaUtils::dphi_wrap(leadMatch->second.phi, subMatch->second.phi) : 0.0F;

      int subFlag = kSubleadMatched;
      if (!subMatch)
      {
        if (recoSel.size() < 2) subFlag = kNoRecoCandidate;
        else
        {
          const DijetCand& observedSub = recoSel[1];
          const float drToReal = AnaUtils::calc_dr(realSubleading.eta, realSubleading.phi, observedSub.eta, observedSub.phi);
          if (drToReal >= maxDr) subFlag = kRealElsewhere;
          else
          { 
            const bool realSubInFiducial =
              realSubleading.pt >= min_pt && std::fabs(realSubleading.eta) <= fiducial_eta;
              subFlag = realSubInFiducial ? kRealElsewhere : kTruthOutsideFiducial;
          }
        }
      }

      const float nRecoJets = dijetMatched ? static_cast<float>(recoSel.size()) : 0.0F;

      // maxpttruth is the highest-pT truth jet in the whole event (no
      // fiducial cut at all), same convention as the legacy ntuple -- it can
      // exceed pt1_truth (the highest-pT *fiducial-selected* truth jet) when
      // an even harder truth jet exists outside the fiducial acceptance.
      const std::array<float, 16> values = {
        truthJetMaxR04pT,
        realLeading.pt, lead.pt, sub.pt, dphiTruth,
        pt1Reco, pt2Reco, dphiReco,
        dijetMatched ? 1.0F : 0.0F, nRecoJets,
        static_cast<float>(centrality), zvrtx, sumEt,
        psi2, psi3,
        static_cast<float>(subFlag)
      };

      categories[0].matches->Fill(values.data());
      ++nDijetsWritten[0];

      const FlavorCategory category = categorizeFlavor(lead.flavor, sub.flavor);
      int tagIndex = -1;
      if (category == kQQ) tagIndex = 1;
      else if (category == kQG) tagIndex = 2;
      else if (category == kGG) tagIndex = 3;
      if (tagIndex >= 0)
        {
          categories[tagIndex].matches->Fill(values.data());
          ++nDijetsWritten[tagIndex];
        }
    }
  std::cout << std::endl;

  for (std::size_t tag = 0; tag < categories.size(); ++tag)
    {
      categories[tag].file->cd();
      categories[tag].stats->Fill(static_cast<float>(nEventsSeen), static_cast<float>(nDijetsWritten[tag]));
      categories[tag].matches->Write();
      categories[tag].stats->Write();
      std::cout << tagNames[tag] << ": wrote " << nDijetsWritten[tag] << " / " << nEventsWithDijet
                << " dijets to " << categories[tag].file->GetName() << std::endl;
      categories[tag].file->Close();
    }
}

#endif
