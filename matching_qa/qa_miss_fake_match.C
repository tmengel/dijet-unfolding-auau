// QA on the miss/fake/matched scheme used by createResponse_noempty_AA.cxx
// (lines ~792-968) and the dR-matching algorithm that produces the
// "matched" flag it reads (makeMatchedTreesInclusiveAuAu.C), run directly
// against the raw production tree "T" in
//   /home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root
//
// Reference definitions (createResponse_noempty_AA.cxx:792-799):
//   truth_good = maxit >= sample_boundary[1] && minit >= truth_subleading_cut
//                && dphi_truth >= dphicuttruth
//     NOTE: the code comment there says "truth_leading_cut" but the actual
//     leading-truth threshold used is sample_boundary[1] = 20.85 GeV, not
//     truth_leading_cut = 14.49 GeV. That is a real discrepancy between the
//     comment and the code, not a QA artifact -- flagged again below.
//   reco_good  = maxi  >= reco_leading_cut && mini  >= reco_subleading_cut
//                && dphi_reco  >= dphicut
//   skip_pair         = !truth_good && !reco_good
//   miss_pair         =  truth_good && !reco_good
//   fake_pair         = !truth_good &&  reco_good && false   // disabled
//   real_pair         =  truth_good &&  reco_good && matched
//   missed_real_pair  =  truth_good &&  reco_good && !matched && false // disabled
//
// i.e. any event with truth_good && reco_good && !matched is currently
// dropped on the floor -- filed under none of the five categories.
//
// matchJets() (makeMatchedTreesInclusiveAuAu.C:52-70): truth jets processed
// leading-to-trailing in pT order; each claims the first not-yet-used reco
// jet within dR < cone_size = 3.0 (a very loose "same hemisphere-ish" cut,
// NOT dR < 0.3). This greedy, order-dependent assignment is exactly the
// mechanism that can let jet: truthJets[0] (leading) grab a reco jet that
// was geometrically the better match for truthJets[2] (a third jet), pushing
// truthJets[1] (subleading) onto a worse partner or truthJets[2] onto
// nothing at all.
//
// Truth/reco fiducial cuts reproduced verbatim from
// makeMatchedTreesInclusiveAuAu.C:236-244:
//   truth jet:  truth_jet_pT >= 3.0  && |truth_jet_eta| <= 0.8
//   reco  jet:  jet_pT >= 3.0 && jet_E >= 0 && jet_unsub_E >= 0 && |jet_eta| <= 0.8
//   event:      |zvrtx| <= 60, and >=2 fiducial truth jets required to even
//               write a tn_match row (events with <2 pass truth jets are
//               silently absent from tn_match entirely -- they can never
//               become "fake" there even if fake_pair were enabled).
//
// This macro answers two specific questions against that scheme:
//
// Q1: of the reco-good-but-truth-bad ("would-be fake") events, how many are
//     that way only because the truth fiducial cut threw out a real second
//     truth jet -- i.e. how many would flip to a real match if the truth
//     pT/eta fiducial window were widened?
//
// Q2: for real (matched) pairs, how often is the reco jet assigned to the
//     truth-SUBLEADING slot actually a better geometric match to a THIRD
//     (or higher) truth jet than to its own official partner -- i.e. how
//     often does a third truth jet's rightful reco partner get stolen by/
//     assigned into the subleading slot by the greedy matcher?

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLegend.h"

#include "../read_binning.h"

namespace {

struct Jet { int id; float pt, eta, phi; };
using Match = std::pair<Jet, Jet>;

float deltaPhi(float a, float b)
{
    float d = a - b;
    while (d > TMath::Pi()) d -= 2.0F * TMath::Pi();
    while (d < -TMath::Pi()) d += 2.0F * TMath::Pi();
    return std::fabs(d);
}
float deltaR(const Jet &a, const Jet &b) { return std::hypot(a.eta - b.eta, deltaPhi(a.phi, b.phi)); }

// Verbatim copy of makeMatchedTreesInclusiveAuAu.C's matchJets().
std::vector<Match> matchJets(const std::vector<Jet> &truth, const std::vector<Jet> &reco, float maxDR)
{
    std::vector<Match> result;
    std::vector<bool> used(reco.size(), false);
    for (const auto &tj : truth)
    {
        for (std::size_t r = 0; r < reco.size(); ++r)
        {
            if (used[r]) continue;
            if (deltaR(tj, reco[r]) >= maxDR) continue;
            result.emplace_back(tj, reco[r]);
            used[r] = true;
            break;
        }
    }
    return result;
}

const Match *findTruthMatch(const std::vector<Match> &matches, int truthId)
{
    auto it = std::find_if(matches.begin(), matches.end(), [truthId](const Match &m) { return m.first.id == truthId; });
    return it == matches.end() ? nullptr : &*it;
}

} // namespace

void qa_miss_fake_match(const std::string &infile = "/home/tmengel/PPG14/rootfiles/v001_20260720/hijing_rho_jet.root",
                         const std::string &configfile = "configs/binning_AA.config")
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("matching_qa/plots", true);

    read_binning rb(configfile.c_str());
    const int nbins = rb.get_nbins();
    std::vector<float> ipt_bins(nbins + 1), ixj_bins(nbins + 1);
    rb.get_pt_bins(ipt_bins.data());
    rb.get_xj_bins(ixj_bins.data());

    const double sample_boundary1 = rb.get_sample_boundary(1); // "truth leading" cut actually used
    const double truth_leading_cut_commented = rb.get_truth_leading_cut(); // what the comment claims
    const double truth_subleading_cut = rb.get_truth_subleading_cut();
    const double reco_leading_cut = rb.get_reco_leading_cut();
    const double reco_subleading_cut = rb.get_reco_subleading_cut();
    const double dphicut = rb.get_dphicut();

    std::cout << "sample_boundary[1] (actual truth-leading threshold used): " << sample_boundary1 << std::endl;
    std::cout << "truth_leading_cut (what the code COMMENT claims is used): " << truth_leading_cut_commented << std::endl;
    std::cout << "truth_subleading_cut: " << truth_subleading_cut << std::endl;
    std::cout << "reco_leading_cut: " << reco_leading_cut << " reco_subleading_cut: " << reco_subleading_cut << std::endl;
    std::cout << "dphicut: " << dphicut << std::endl;

    TFile *f = TFile::Open(infile.c_str(), "READ");
    TTree *t = (TTree *) f->Get("T");

    int jetid = 0, centrality = 0;
    float vertex = 0;
    std::vector<float> *truthPt = nullptr, *truthEta = nullptr, *truthPhi = nullptr;
    std::vector<float> *recoPt = nullptr, *recoE = nullptr, *recoUnsubE = nullptr, *recoEta = nullptr, *recoPhi = nullptr;
    t->SetBranchStatus("*", 0);
    for (const char *b : {"jetid", "zvrtx", "cent", "truth_jet_pT", "truth_jet_eta", "truth_jet_phi",
                           "jet_pT", "jet_E", "jet_unsub_E", "jet_eta", "jet_phi"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("jetid", &jetid);
    t->SetBranchAddress("zvrtx", &vertex);
    t->SetBranchAddress("cent", &centrality);
    t->SetBranchAddress("truth_jet_pT", &truthPt);
    t->SetBranchAddress("truth_jet_eta", &truthEta);
    t->SetBranchAddress("truth_jet_phi", &truthPhi);
    t->SetBranchAddress("jet_pT", &recoPt);
    t->SetBranchAddress("jet_E", &recoE);
    t->SetBranchAddress("jet_unsub_E", &recoUnsubE);
    t->SetBranchAddress("jet_eta", &recoEta);
    t->SetBranchAddress("jet_phi", &recoPhi);
    t->SetCacheSize(256 * 1024 * 1024);
    t->AddBranchToCache("*", kTRUE);

    // ---- Set A: scheme population counts ----
    Long64_t n_skip = 0, n_miss = 0, n_wouldbefake = 0, n_real = 0, n_wouldbemissedreal = 0;
    Long64_t n_events_used = 0;

    TH1D *h_truth_lead_miss = new TH1D("h_truth_lead_miss", ";truth leading jet p_{T} [GeV];counts", 40, 0, 100);
    TH1D *h_truth_lead_real = new TH1D("h_truth_lead_real", ";truth leading jet p_{T} [GeV];counts", 40, 0, 100);
    TH1D *h_ntruthfid = new TH1D("h_ntruthfid", ";# fiducial truth jets (p_{T}#geq3, |#eta|#leq0.8);counts", 8, 0, 8);

    // ---- Set B: fiducial-cut recovery (Q1) ----
    // Would-be-fake events specifically caused by <2 fiducial truth jets.
    Long64_t n_wouldbefake_truthlt2 = 0;
    Long64_t n_wouldbefake_truthge2_kinematic = 0; // >=2 fiducial truth jets exist, but truth_good kinematics fail
    TH1D *h_missing_truth_pt = new TH1D("h_missing_truth_pt", ";next-best truth jet p_{T} [GeV] (fails nominal fiducial cut);counts", 30, 0, 3);
    TH1D *h_missing_truth_eta = new TH1D("h_missing_truth_eta", ";next-best truth jet |#eta| (fails nominal fiducial cut);counts", 30, 0.8, 1.4);
    TH1D *h_missing_truth_reason = new TH1D("h_missing_truth_reason", ";;count", 3, 0, 3);
    // recovery counts under a few candidate widened fiducial definitions
    struct WideCut { std::string label; float ptmin; float etamax; };
    const std::vector<WideCut> wideCuts = {
        {"pT>=3, |#eta|#leq1.1", 3.0F, 1.1F},
        {"pT>=1, |#eta|#leq0.8", 1.0F, 0.8F},
        {"pT>=1, |#eta|#leq1.1", 1.0F, 1.1F},
    };
    std::vector<Long64_t> n_recovered_2truth(wideCuts.size(), 0);      // now has >=2 fiducial truth jets
    std::vector<Long64_t> n_recovered_truthgood(wideCuts.size(), 0);   // AND that pair also passes truth_good kinematics
    std::vector<Long64_t> n_recovered_realmatch(wideCuts.size(), 0);   // AND dR-matches to the reco pair (truly becomes real_pair)

    // ---- Set C: third-jet contamination of the subleading match (Q2) ----
    Long64_t n_matched_with_3rd_truth = 0;
    Long64_t n_thirdjet_closer = 0;       // reco "subleading" jet is dR-closer to truth[2] than to truth[1]
    Long64_t n_third_truth_unmatched = 0; // truth[2] itself found no reco partner at all
    Long64_t n_subleading_rank_not2 = 0;  // reco jet assigned to subleading isn't recoNominal's actual #2 by pT
    TH1D *h_dR_diff = new TH1D("h_dR_diff", ";dR(reco_{sub}, truth_{3rd}) - dR(reco_{sub}, truth_{sub});counts", 60, -3, 3);
    TH1D *h_subleading_rank = new TH1D("h_subleading_rank", ";reco rank (by p_{T}) of jet assigned to truth-subleading;counts", 8, 1, 9);
    TH2D *h_thirdjet_closer_vs_pt = new TH2D("h_thirdjet_closer_vs_pt", ";truth subleading p_{T} [GeV];third-jet-closer?", 20, 0, 40, 2, 0, 2);

    const Long64_t entries = t->GetEntries();
    std::cout << "raw entries: " << entries << std::endl;

    for (Long64_t i = 0; i < entries; i++)
    {
        t->GetEntry(i);
        if (i % (entries / 10 == 0 ? 1 : entries / 10) == 0)
            std::cout << i << " / " << entries << "\r" << std::flush;

        if (std::fabs(vertex) > 60.0F) continue;
        if (!truthPt || truthPt->empty()) continue;
        if (truthPt->size() != truthEta->size() || truthPt->size() != truthPhi->size()) continue;
        if (recoPt->size() != recoE->size() || recoPt->size() != recoUnsubE->size() ||
            recoPt->size() != recoEta->size() || recoPt->size() != recoPhi->size()) continue;

        // full unfiltered truth list, kept for the "widen the cut" study
        std::vector<Jet> truthAll;
        for (std::size_t j = 0; j < truthPt->size(); j++)
            truthAll.push_back({(int) j, truthPt->at(j), truthEta->at(j), truthPhi->at(j)});

        std::vector<Jet> truthNom, recoNom;
        for (const auto &tj : truthAll)
            if (tj.pt >= 3.0F && std::fabs(tj.eta) <= 0.8F) truthNom.push_back(tj);
        for (std::size_t j = 0; j < recoPt->size(); j++)
            if (recoPt->at(j) >= 3.0F && recoE->at(j) >= 0.0F && recoUnsubE->at(j) >= 0.0F && std::fabs(recoEta->at(j)) <= 0.8F)
                recoNom.push_back({(int) j, recoPt->at(j), recoEta->at(j), recoPhi->at(j)});

        const auto byPt = [](const Jet &a, const Jet &b) { return a.pt > b.pt; };
        std::sort(truthAll.begin(), truthAll.end(), byPt);
        std::sort(truthNom.begin(), truthNom.end(), byPt);
        std::sort(recoNom.begin(), recoNom.end(), byPt);

        h_ntruthfid->Fill(std::min<int>(truthNom.size(), 7));

        // reco_good is purely reco-side: needs >=2 fiducial reco jets passing the thresholds
        bool reco_good = false;
        float recoLeadPt = 0, recoSubPt = 0, dphi_reco = 0;
        if (recoNom.size() >= 2)
        {
            recoLeadPt = recoNom[0].pt;
            recoSubPt = recoNom[1].pt;
            dphi_reco = deltaPhi(recoNom[0].phi, recoNom[1].phi);
            reco_good = (recoLeadPt >= reco_leading_cut && recoSubPt >= reco_subleading_cut && dphi_reco >= dphicut);
        }

        bool truth_good = false;
        float truthLeadPt = 0, truthSubPt = 0, dphi_truth = 0;
        if (truthNom.size() >= 2)
        {
            truthLeadPt = truthNom[0].pt;
            truthSubPt = truthNom[1].pt;
            dphi_truth = deltaPhi(truthNom[0].phi, truthNom[1].phi);
            truth_good = (truthLeadPt >= sample_boundary1 && truthSubPt >= truth_subleading_cut && dphi_truth >= dphicut);
        }

        const bool skip_pair = !truth_good && !reco_good;
        if (skip_pair) { n_skip++; continue; }

        n_events_used++;

        // ---- matching (only meaningful/computed by production when truthNom.size()>=2) ----
        bool matched = false;
        std::vector<Match> matches;
        if (truthNom.size() >= 2)
        {
            matches = matchJets(truthNom, recoNom, 3.0F);
            matched = findTruthMatch(matches, truthNom[0].id) != nullptr &&
                      findTruthMatch(matches, truthNom[1].id) != nullptr;
        }

        const bool miss_pair = truth_good && !reco_good;
        const bool wouldbefake_pair = !truth_good && reco_good;
        const bool real_pair = truth_good && reco_good && matched;
        const bool wouldbemissedreal_pair = truth_good && reco_good && !matched;

        if (miss_pair) { n_miss++; h_truth_lead_miss->Fill(truthLeadPt); }
        if (real_pair) { n_real++; h_truth_lead_real->Fill(truthLeadPt); }
        if (wouldbemissedreal_pair) n_wouldbemissedreal++;

        // ================= Q1: fiducial-cut-driven fakes =================
        if (wouldbefake_pair)
        {
            n_wouldbefake++;
            if (truthNom.size() < 2)
            {
                n_wouldbefake_truthlt2++;
                h_missing_truth_reason->Fill(0.5); // "fewer than 2 fiducial truth jets"

                // find, among ALL truth jets, the best dR match to each reco jet
                // used in the reco-good pair, restricted to jets that FAIL the
                // nominal fiducial cut (so they're the "recoverable" candidates)
                for (const Jet &rj : {recoNom[0], recoNom[1]})
                {
                    float bestDR = 1e9; const Jet *best = nullptr;
                    for (const auto &tj : truthAll)
                    {
                        if (tj.pt >= 3.0F && std::fabs(tj.eta) <= 0.8F) continue; // already fiducial, not "missing"
                        const float dr = deltaR(rj, tj);
                        if (dr < bestDR) { bestDR = dr; best = &tj; }
                    }
                    if (best && bestDR < 3.0F)
                    {
                        h_missing_truth_pt->Fill(best->pt);
                        h_missing_truth_eta->Fill(std::fabs(best->eta));
                    }
                }

                // recovery check under each candidate widened fiducial cut
                for (std::size_t w = 0; w < wideCuts.size(); w++)
                {
                    std::vector<Jet> truthWide;
                    for (const auto &tj : truthAll)
                        if (tj.pt >= wideCuts[w].ptmin && std::fabs(tj.eta) <= wideCuts[w].etamax)
                            truthWide.push_back(tj);
                    std::sort(truthWide.begin(), truthWide.end(), byPt);
                    if (truthWide.size() < 2) continue;
                    n_recovered_2truth[w]++;

                    const float wLead = truthWide[0].pt, wSub = truthWide[1].pt;
                    const float wDphi = deltaPhi(truthWide[0].phi, truthWide[1].phi);
                    const bool wTruthGood = (wLead >= sample_boundary1 && wSub >= truth_subleading_cut && wDphi >= dphicut);
                    if (!wTruthGood) continue;
                    n_recovered_truthgood[w]++;

                    const auto wMatches = matchJets(truthWide, recoNom, 3.0F);
                    const bool wMatched = findTruthMatch(wMatches, truthWide[0].id) != nullptr &&
                                          findTruthMatch(wMatches, truthWide[1].id) != nullptr;
                    if (wMatched) n_recovered_realmatch[w]++;
                }
            }
            else
            {
                n_wouldbefake_truthge2_kinematic++;
                h_missing_truth_reason->Fill(1.5); // ">=2 fiducial truth jets, but kinematic cuts fail"
            }
        }

        // ================= Q2: third-jet contamination =================
        if (real_pair && truthNom.size() >= 3)
        {
            n_matched_with_3rd_truth++;
            const Match *subMatch = findTruthMatch(matches, truthNom[1].id);
            const Jet &recoSub = subMatch->second;

            const float dR_to_sub = deltaR(recoSub, truthNom[1]);
            const float dR_to_third = deltaR(recoSub, truthNom[2]);
            h_dR_diff->Fill(dR_to_third - dR_to_sub);
            const bool thirdCloser = dR_to_third < dR_to_sub;
            if (thirdCloser) n_thirdjet_closer++;
            h_thirdjet_closer_vs_pt->Fill(truthNom[1].pt, thirdCloser ? 1 : 0);

            if (findTruthMatch(matches, truthNom[2].id) == nullptr) n_third_truth_unmatched++;

            // reco-pT rank of the jet assigned to the subleading slot
            int rank = 1;
            for (const auto &rj : recoNom) { if (rj.pt > recoSub.pt) rank++; }
            h_subleading_rank->Fill(std::min(rank, 8));
            if (rank != 2) n_subleading_rank_not2++;
        }
    }
    std::cout << std::endl;

    // ---------------------------------------------------------------------
    std::cout << "\n==================== Set A: scheme population ====================\n";
    std::cout << "events considered (skip excluded): " << n_events_used << std::endl;
    std::cout << "  skip_pair                 : " << n_skip << std::endl;
    std::cout << "  miss_pair                 : " << n_miss << std::endl;
    std::cout << "  would-be fake_pair        : " << n_wouldbefake << "  (currently disabled -> dropped)" << std::endl;
    std::cout << "  real_pair                 : " << n_real << std::endl;
    std::cout << "  would-be missed_real_pair : " << n_wouldbemissedreal << "  (currently disabled -> dropped)" << std::endl;
    std::cout << "  [note] truth_good uses sample_boundary[1]=" << sample_boundary1
              << " as the leading-truth threshold, NOT truth_leading_cut=" << truth_leading_cut_commented
              << " that the code comment claims." << std::endl;

    std::cout << "\n==================== Q1: fiducial-cut-driven fakes ====================\n";
    std::cout << "would-be fake events                          : " << n_wouldbefake << std::endl;
    std::cout << "  ...caused by <2 fiducial truth jets          : " << n_wouldbefake_truthlt2
              << " (" << (n_wouldbefake > 0 ? 100.0 * n_wouldbefake_truthlt2 / n_wouldbefake : 0) << "%)" << std::endl;
    std::cout << "  ...caused by kinematic cuts (>=2 truth jets) : " << n_wouldbefake_truthge2_kinematic
              << " (" << (n_wouldbefake > 0 ? 100.0 * n_wouldbefake_truthge2_kinematic / n_wouldbefake : 0) << "%)" << std::endl;
    for (std::size_t w = 0; w < wideCuts.size(); w++)
    {
        std::cout << "  widened cut [" << wideCuts[w].label << "]: "
                  << "now has >=2 truth jets: " << n_recovered_2truth[w]
                  << ", also passes truth_good: " << n_recovered_truthgood[w]
                  << ", also dR-matches (-> real_pair): " << n_recovered_realmatch[w]
                  << "  (recovery rate " << (n_wouldbefake_truthlt2 > 0 ? 100.0 * n_recovered_realmatch[w] / n_wouldbefake_truthlt2 : 0)
                  << "% of the <2-truth-jet fakes)" << std::endl;
    }

    std::cout << "\n==================== Q2: third-jet contamination of subleading match ====================\n";
    std::cout << "real_pair events with a 3rd+ fiducial truth jet present : " << n_matched_with_3rd_truth << std::endl;
    std::cout << "  ...where reco_sub is dR-CLOSER to truth[2] than truth[1] (\"stolen\") : " << n_thirdjet_closer
              << " (" << (n_matched_with_3rd_truth > 0 ? 100.0 * n_thirdjet_closer / n_matched_with_3rd_truth : 0) << "%)" << std::endl;
    std::cout << "  ...where truth[2] itself ended up completely unmatched              : " << n_third_truth_unmatched
              << " (" << (n_matched_with_3rd_truth > 0 ? 100.0 * n_third_truth_unmatched / n_matched_with_3rd_truth : 0) << "%)" << std::endl;
    std::cout << "reco-pT rank of the jet assigned to the subleading slot != 2 : " << n_subleading_rank_not2
              << " / " << n_real << " real_pair events" << std::endl;

    // ---------------------------------------------------------------------
    // plots
    {
        TCanvas c("c1", "c1", 700, 500);
        h_truth_lead_real->SetLineColor(kBlack); h_truth_lead_real->SetLineWidth(2);
        h_truth_lead_miss->SetLineColor(kRed + 1); h_truth_lead_miss->SetLineWidth(2);
        h_truth_lead_real->SetTitle("");
        h_truth_lead_real->Draw("hist");
        h_truth_lead_miss->Draw("hist same");
        TLegend leg(0.55, 0.7, 0.88, 0.85);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(h_truth_lead_real, "real_pair", "l");
        leg.AddEntry(h_truth_lead_miss, "miss_pair", "l");
        leg.Draw();
        c.Print("matching_qa/plots/qa_truth_lead_miss_vs_real.pdf");
    }
    {
        TCanvas c("c2", "c2", 700, 500);
        h_ntruthfid->Draw("hist");
        c.Print("matching_qa/plots/qa_n_fiducial_truth_jets.pdf");
    }
    {
        TCanvas c("c3", "c3", 700, 500);
        h_missing_truth_pt->SetLineColor(kBlue + 1); h_missing_truth_pt->SetLineWidth(2);
        h_missing_truth_pt->Draw("hist");
        c.Print("matching_qa/plots/qa_missing_truth_pt.pdf");
    }
    {
        TCanvas c("c4", "c4", 700, 500);
        h_missing_truth_eta->SetLineColor(kBlue + 1); h_missing_truth_eta->SetLineWidth(2);
        h_missing_truth_eta->Draw("hist");
        c.Print("matching_qa/plots/qa_missing_truth_eta.pdf");
    }
    {
        TCanvas c("c5", "c5", 700, 500);
        h_dR_diff->SetLineColor(kMagenta + 1); h_dR_diff->SetLineWidth(2);
        h_dR_diff->Draw("hist");
        c.Print("matching_qa/plots/qa_dR_diff_thirdjet.pdf");
    }
    {
        TCanvas c("c6", "c6", 700, 500);
        h_subleading_rank->SetLineColor(kGreen + 2); h_subleading_rank->SetLineWidth(2);
        h_subleading_rank->Draw("hist");
        c.Print("matching_qa/plots/qa_subleading_reco_rank.pdf");
    }

    TFile fout("matching_qa/matching_qa.root", "RECREATE");
    h_truth_lead_miss->Write(); h_truth_lead_real->Write(); h_ntruthfid->Write();
    h_missing_truth_pt->Write(); h_missing_truth_eta->Write(); h_missing_truth_reason->Write();
    h_dR_diff->Write(); h_subleading_rank->Write(); h_thirdjet_closer_vs_pt->Write();
    fout.Close();

    std::cout << "\nWrote matching_qa/plots/*.pdf and matching_qa/matching_qa.root" << std::endl;
}
