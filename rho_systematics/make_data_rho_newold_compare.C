// Comparison of the NEW vs OLD rho-subtracted AuAu jet production, r=0.3.
//
//   new: input/data/data_v004_20260821_calibrated_merged.root  (296.56M events)
//   old: input/data/run2auau_rho_jet.root                      (295.20M events)
//
// Both are the raw per-event production tree "T" of anti-kt R=0.3 jets
// reconstructed on rho/UE-subtracted towers (AntiKt_TowerInfo_r03_Rho2 in
// macros/data/Fun4All_Dijets_AuAu.C). The new production adds the eta-shape
// rho calibration (SubtractTowersRhov1::set_etaCalib_directPath) and a
// MinBias/trigger event selection, so this macro is the data-level analogue
// of make_data_rho_compare.C but with "rho vs no-rho" replaced by
// "new rho vs old rho".
//
// Branch names differ between the two productions:
//   mbd charge   : mbd_q                       <-> mbdQ
//   rho (CEMC)   : rho_val_TowerRho_MULT_CEMC  <-> rho_cemc
//   rho (HCALIN) : rho_val_TowerRho_MULT_HCALIN<-> rho_hcalin
//   rho (HCALOUT): rho_val_TowerRho_MULT_HCALOUT<->rho_hcalout
// The jet vectors (jet_pT/eta/phi/E/unsub_pT/unsub_E) and cent/zvrtx/sumeT
// have the same names in both.
//
// Produces, for both productions:
//   - event-level: centrality, z-vertex, sum-eT, mbd charge, rho, jets/event
//   - per-jet: pT (all and with the jet_E>=0 && jet_unsub_E>=0 cut), eta, phi,
//     unsubtracted pT, subtracted energy fraction, pT_unsub - pT_sub
//   - the two selections that actually gate the unfolding input: the E>=0
//     positivity cut and the UE fake-jet cut
//     pt_unsub - pt_sub > fcut(centrality)  (makeDataTreeAuAu.C:152-193 /
//     getBackground.C:197-233 / anaConf.h:309-313)
//   - <rho> and <jets/event> vs centrality profiles
//   - a printed + written text table of event/jet statistics
//
// Spectra are normalised per event before overlaying, since the two samples
// do not have identical event counts.
//
// Usage (from macros/unfolding/dijet-unfolding-auau):
//   root -l -b -q 'rho_systematics/make_data_rho_newold_compare.C(1)'
// The argument is an event stride (1 = full sample).

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TStopwatch.h"

namespace {

const double reco_pt_min = 8.0;
const double abs_eta_max = 0.8;

struct Sample
{
    std::string tag;         // "new" / "old"
    std::string label;       // legend label
    std::string file;
    std::string br_mbd;
    std::string br_rho_cemc;
    std::string br_rho_hcalin;
    std::string br_rho_hcalout;
    bool has_minbias = false;
};

struct Hists
{
    // ---- event level ----
    TH1D *h_cent = nullptr;
    TH1D *h_zvrtx = nullptr;
    TH1D *h_sumeT = nullptr;
    TH1D *h_mbdq = nullptr;
    TH1D *h_rho_cemc = nullptr;
    TH1D *h_rho_hcalin = nullptr;
    TH1D *h_rho_hcalout = nullptr;
    TH1D *h_njet_all = nullptr;   // jets/event, no selection
    TH1D *h_njet_sel = nullptr;   // jets/event, pT>8, |eta|<0.8
    TH1D *h_lead_pt = nullptr;    // leading selected jet pT
    TProfile *p_rho_vs_cent = nullptr;
    TProfile *p_njet_vs_cent = nullptr;

    // ---- per jet ----
    TH1D *h_pt_all = nullptr;
    TH1D *h_pt_epos = nullptr;
    TH1D *h_pt_unsub = nullptr;
    TH1D *h_eta = nullptr;
    TH1D *h_phi = nullptr;
    TH1D *h_efrac = nullptr;           // (E_unsub - E)/E_unsub
    TH1D *h_ptexcess = nullptr;        // pT_unsub - pT_sub
    TH1D *h_uecut_pt_before = nullptr;
    TH1D *h_uecut_pt_after = nullptr;
    TProfile *p_efrac_vs_pt = nullptr;
    TProfile *p_ptexcess_vs_cent = nullptr;

    // ---- counters ----
    Long64_t n_events = 0;
    Long64_t n_events_minbias = 0;
    Long64_t n_jets_raw = 0;       // all jets in tree
    Long64_t n_jets_sel = 0;       // pT>8, |eta|<0.8
    Long64_t n_jets_epos = 0;      // + E>=0 && E_unsub>=0
    Long64_t n_jets_uepass = 0;    // + UE fake-jet cut
    Long64_t n_ev_ge1jet = 0;
    Long64_t n_ev_ge2jet = 0;
    double   sum_pt_sel = 0;
};

Hists make_hists(const std::string &tag)
{
    Hists h;
    const std::string s = "_" + tag;

    h.h_cent    = new TH1D(("h_cent" + s).c_str(),    ";centrality [%];events", 100, 0, 100);
    h.h_zvrtx   = new TH1D(("h_zvrtx" + s).c_str(),   ";z_{vtx} [cm];events", 120, -60, 60);
    h.h_sumeT   = new TH1D(("h_sumeT" + s).c_str(),   ";#sum E_{T} [GeV];events", 250, -500, 2000);
    h.h_mbdq    = new TH1D(("h_mbdq" + s).c_str(),    ";MBD charge;events", 300, 0, 3000);
    h.h_rho_cemc    = new TH1D(("h_rho_cemc" + s).c_str(),    ";#rho_{CEMC} [GeV];events", 200, 0, 1.0);
    h.h_rho_hcalin  = new TH1D(("h_rho_hcalin" + s).c_str(),  ";#rho_{HCALIN} [GeV];events", 200, 0, 1.0);
    h.h_rho_hcalout = new TH1D(("h_rho_hcalout" + s).c_str(), ";#rho_{HCALOUT} [GeV];events", 200, 0, 1.0);
    h.h_njet_all = new TH1D(("h_njet_all" + s).c_str(), ";jets / event (all);events", 30, -0.5, 29.5);
    h.h_njet_sel = new TH1D(("h_njet_sel" + s).c_str(), ";jets / event (p_{T}>8, |#eta|<0.8);events", 15, -0.5, 14.5);
    h.h_lead_pt  = new TH1D(("h_lead_pt" + s).c_str(),  ";leading jet p_{T} [GeV];events", 100, 0, 100);
    h.p_rho_vs_cent  = new TProfile(("p_rho_vs_cent" + s).c_str(),  ";centrality [%];#LT#rho_{CEMC}#GT [GeV]", 100, 0, 100);
    h.p_njet_vs_cent = new TProfile(("p_njet_vs_cent" + s).c_str(), ";centrality [%];#LTjets/event#GT (p_{T}>8)", 100, 0, 100);

    h.h_pt_all   = new TH1D(("h_pt_all" + s).c_str(),   ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_pt_epos  = new TH1D(("h_pt_epos" + s).c_str(),  ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_pt_unsub = new TH1D(("h_pt_unsub" + s).c_str(), ";unsubtracted jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_eta      = new TH1D(("h_eta" + s).c_str(),      ";reco jet #eta;counts", 40, -1.0, 1.0);
    h.h_phi      = new TH1D(("h_phi" + s).c_str(),      ";reco jet #phi;counts", 40, -3.2, 3.2);
    h.h_efrac    = new TH1D(("h_efrac" + s).c_str(),    ";(E_{unsub} - E) / E_{unsub};counts", 100, -0.2, 1.0);
    h.h_ptexcess = new TH1D(("h_ptexcess" + s).c_str(), ";p_{T}^{unsub} - p_{T}^{sub} [GeV];counts", 120, -10, 50);
    h.h_uecut_pt_before = new TH1D(("h_uecut_pt_before" + s).c_str(), ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_uecut_pt_after  = new TH1D(("h_uecut_pt_after" + s).c_str(),  ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.p_efrac_vs_pt      = new TProfile(("p_efrac_vs_pt" + s).c_str(),      ";reco jet p_{T} [GeV];#LT(E_{unsub}-E)/E_{unsub}#GT", 100, 0, 100);
    h.p_ptexcess_vs_cent = new TProfile(("p_ptexcess_vs_cent" + s).c_str(), ";centrality [%];#LTp_{T}^{unsub} - p_{T}^{sub}#GT [GeV]", 100, 0, 100);

    for (TH1 *hh : {(TH1*)h.h_cent, (TH1*)h.h_zvrtx, (TH1*)h.h_sumeT, (TH1*)h.h_mbdq,
                    (TH1*)h.h_rho_cemc, (TH1*)h.h_rho_hcalin, (TH1*)h.h_rho_hcalout,
                    (TH1*)h.h_njet_all, (TH1*)h.h_njet_sel, (TH1*)h.h_lead_pt,
                    (TH1*)h.h_pt_all, (TH1*)h.h_pt_epos, (TH1*)h.h_pt_unsub, (TH1*)h.h_eta,
                    (TH1*)h.h_phi, (TH1*)h.h_efrac, (TH1*)h.h_ptexcess,
                    (TH1*)h.h_uecut_pt_before, (TH1*)h.h_uecut_pt_after})
    { hh->SetDirectory(nullptr); hh->Sumw2(); }
    for (TProfile *pp : {h.p_rho_vs_cent, h.p_njet_vs_cent, h.p_efrac_vs_pt, h.p_ptexcess_vs_cent})
        pp->SetDirectory(nullptr);

    return h;
}

Hists fill_sample(const Sample &s, Long64_t stride)
{
    Hists h = make_hists(s.tag);

    TFile *f = TFile::Open(s.file.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << s.file << std::endl; return h; }
    TTree *t = (TTree*) f->Get("T");
    if (!t) { std::cerr << "No tree T in " << s.file << std::endl; return h; }

    std::vector<float> *jet_pT = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *jet_E = nullptr, *jet_unsub_E = nullptr, *jet_unsub_pT = nullptr;
    int centrality = 0, is_minbias = 1;
    float zvrtx = 0, sumeT = 0, mbdq = 0;
    float rho_cemc = 0, rho_hcalin = 0, rho_hcalout = 0;

    t->SetBranchStatus("*", 0);
    std::vector<std::string> on = {"jet_pT", "jet_eta", "jet_phi", "jet_E", "jet_unsub_E", "jet_unsub_pT",
                                   "cent", "zvrtx", "sumeT",
                                   s.br_mbd, s.br_rho_cemc, s.br_rho_hcalin, s.br_rho_hcalout};
    if (s.has_minbias) on.push_back("is_minbias");
    for (const std::string &b : on) t->SetBranchStatus(b.c_str(), 1);

    t->SetBranchAddress("jet_pT", &jet_pT);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);
    t->SetBranchAddress("jet_E", &jet_E);
    t->SetBranchAddress("jet_unsub_E", &jet_unsub_E);
    t->SetBranchAddress("jet_unsub_pT", &jet_unsub_pT);
    t->SetBranchAddress("cent", &centrality);
    t->SetBranchAddress("zvrtx", &zvrtx);
    t->SetBranchAddress("sumeT", &sumeT);
    t->SetBranchAddress(s.br_mbd.c_str(), &mbdq);
    t->SetBranchAddress(s.br_rho_cemc.c_str(), &rho_cemc);
    t->SetBranchAddress(s.br_rho_hcalin.c_str(), &rho_hcalin);
    t->SetBranchAddress(s.br_rho_hcalout.c_str(), &rho_hcalout);
    if (s.has_minbias) t->SetBranchAddress("is_minbias", &is_minbias);

    t->SetCacheSize(256 * 1024 * 1024);
    t->AddBranchToCache("*", kTRUE);

    // Same UE/fake-jet rejection cut used in makeDataTreeAuAu.C / getBackground.C / anaConf.h.
    TF1 fcut("fcut", "[0]+[1]*TMath::Exp(-[2]*x)", 0.0, 100.0);
    fcut.SetParameters(0.0, 40, 0.038);

    const Long64_t total_entries = t->GetEntries();
    const Long64_t entries = total_entries / stride;
    std::cout << "[" << s.tag << "] " << s.file << "\n[" << s.tag << "] raw entries: " << total_entries
              << " (stride " << stride << " -> " << entries << " events)" << std::endl;

    TStopwatch sw; sw.Start();
    for (Long64_t k = 0; k < entries; k++)
    {
        t->GetEntry(k * stride);
        if (k % (entries/20 == 0 ? 1 : entries/20) == 0)
            std::cout << "[" << s.tag << "] " << k << " / " << entries
                      << "  (" << std::fixed << std::setprecision(0) << sw.RealTime() << " s)"
                      << "          \r" << std::flush, sw.Continue();

        h.n_events++;
        if (s.has_minbias && is_minbias) h.n_events_minbias++;

        h.h_cent->Fill(centrality);
        h.h_zvrtx->Fill(zvrtx);
        h.h_sumeT->Fill(sumeT);
        h.h_mbdq->Fill(mbdq);
        h.h_rho_cemc->Fill(rho_cemc);
        h.h_rho_hcalin->Fill(rho_hcalin);
        h.h_rho_hcalout->Fill(rho_hcalout);
        h.p_rho_vs_cent->Fill(centrality, rho_cemc);

        h.n_jets_raw += (Long64_t) jet_pT->size();
        h.h_njet_all->Fill((double) jet_pT->size());

        const double cut_value = fcut.Eval(centrality);

        int njet_sel = 0;
        double lead_pt = -1;

        for (size_t j = 0; j < jet_pT->size(); j++)
        {
            const float pt = jet_pT->at(j);
            const float eta = jet_eta->at(j);
            if (pt < reco_pt_min || std::fabs(eta) > abs_eta_max) continue;

            njet_sel++;
            h.n_jets_sel++;
            h.sum_pt_sel += pt;
            if (pt > lead_pt) lead_pt = pt;

            h.h_pt_all->Fill(pt);

            const bool e_positive = (jet_E->at(j) >= 0) && (jet_unsub_E->at(j) >= 0);
            if (!e_positive) continue;
            h.n_jets_epos++;

            h.h_pt_epos->Fill(pt);
            h.h_pt_unsub->Fill(jet_unsub_pT->at(j));
            h.h_eta->Fill(eta);
            h.h_phi->Fill(jet_phi->at(j));

            const float unsub_E = jet_unsub_E->at(j);
            if (unsub_E > 0)
            {
                const double ef = (unsub_E - jet_E->at(j)) / unsub_E;
                h.h_efrac->Fill(ef);
                h.p_efrac_vs_pt->Fill(pt, ef);
            }

            const float pt_unsub_excess = jet_unsub_pT->at(j) - pt;
            h.h_ptexcess->Fill(pt_unsub_excess);
            h.p_ptexcess_vs_cent->Fill(centrality, pt_unsub_excess);

            h.h_uecut_pt_before->Fill(pt);
            if (pt_unsub_excess <= cut_value)
            {
                h.h_uecut_pt_after->Fill(pt);
                h.n_jets_uepass++;
            }
        }

        h.h_njet_sel->Fill(njet_sel);
        h.p_njet_vs_cent->Fill(centrality, njet_sel);
        if (njet_sel >= 1) { h.n_ev_ge1jet++; h.h_lead_pt->Fill(lead_pt); }
        if (njet_sel >= 2) h.n_ev_ge2jet++;
    }
    sw.Stop();
    std::cout << "\n[" << s.tag << "] done in " << sw.RealTime() << " s" << std::endl;
    f->Close();
    return h;
}

void draw_overlay(TH1 *h_new, TH1 *h_old, const std::string &xtitle, const std::string &outpdf,
                  bool logy = true, double ratio_lo = 0.5, double ratio_hi = 1.5,
                  const std::string &ytitle = "")
{
    h_new->SetLineColor(kRed + 1);  h_new->SetLineWidth(2); h_new->SetMarkerColor(kRed + 1);
    h_old->SetLineColor(kBlack);    h_old->SetLineWidth(2); h_old->SetMarkerColor(kBlack);

    TCanvas c("c", "c", 600, 700);
    TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1.0);
    TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.3);
    p1->SetBottomMargin(0.02); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.3);
    p1->Draw(); p2->Draw();

    p1->cd();
    if (logy) p1->SetLogy();
    h_old->GetXaxis()->SetLabelSize(0);
    h_old->SetTitle("");
    if (!ytitle.empty()) h_old->GetYaxis()->SetTitle(ytitle.c_str());
    h_old->Draw("hist e");
    h_new->Draw("hist e same");
    TLegend leg(0.52, 0.68, 0.88, 0.86);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(h_old, "old #rho (run2auau_rho_jet)", "l");
    leg.AddEntry(h_new, "new #rho (v004 calibrated)", "l");
    leg.Draw();

    p2->cd();
    TH1 *ratio = (TH1*) h_new->Clone("ratio_tmp");
    ratio->SetDirectory(nullptr);
    ratio->Divide(h_old);
    ratio->SetTitle((";" + xtitle + ";new / old").c_str());
    ratio->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetTitleSize(0.11); ratio->GetXaxis()->SetLabelSize(0.09);
    ratio->GetYaxis()->SetTitleSize(0.10); ratio->GetYaxis()->SetLabelSize(0.09);
    ratio->GetYaxis()->SetTitleOffset(0.45);
    ratio->SetMarkerStyle(20); ratio->SetMarkerSize(0.7);
    ratio->Draw("ep");
    TLine line(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
    line.SetLineStyle(2);
    line.Draw();

    c.Print(outpdf.c_str());
    delete ratio;
}

// per-event-normalised clone, so the two samples' different event counts cancel
TH1D * norm_clone(TH1D *h, Long64_t nev, const char *suffix)
{
    TH1D *c = (TH1D*) h->Clone((std::string(h->GetName()) + suffix).c_str());
    c->SetDirectory(nullptr);
    if (nev > 0) c->Scale(1.0 / (double) nev);
    c->GetYaxis()->SetTitle("counts / event");
    return c;
}

std::string stats_table(const Hists &n, const Hists &o)
{
    std::ostringstream ss;
    ss << std::fixed;
    auto row = [&](const std::string &name, double vn, double vo, int prec) {
        ss << std::setw(38) << std::left << name << std::right
           << std::setw(18) << std::setprecision(prec) << vn
           << std::setw(18) << std::setprecision(prec) << vo
           << std::setw(12) << std::setprecision(4)
           << (vo != 0 ? vn / vo : 0.0) << "\n";
    };

    ss << "\n" << std::string(86, '=') << "\n";
    ss << std::setw(38) << std::left << "QUANTITY" << std::right
       << std::setw(18) << "NEW (v004)" << std::setw(18) << "OLD (run2auau)" << std::setw(12) << "new/old" << "\n";
    ss << std::string(86, '-') << "\n";

    ss << "-- events --\n";
    row("events processed", n.n_events, o.n_events, 0);
    if (n.n_events_minbias > 0)
        row("  of which is_minbias==1 (new only)", n.n_events_minbias, 0, 0);
    row("<centrality> [%]", n.h_cent->GetMean(), o.h_cent->GetMean(), 3);
    row("<z_vtx> [cm]", n.h_zvrtx->GetMean(), o.h_zvrtx->GetMean(), 3);
    row("RMS z_vtx [cm]", n.h_zvrtx->GetRMS(), o.h_zvrtx->GetRMS(), 3);
    row("<sum E_T> [GeV]  (*)", n.h_sumeT->GetMean(), o.h_sumeT->GetMean(), 3);
    row("<MBD charge>  (*)", n.h_mbdq->GetMean(), o.h_mbdq->GetMean(), 3);
    row("<rho_CEMC> [GeV]", n.h_rho_cemc->GetMean(), o.h_rho_cemc->GetMean(), 5);
    row("<rho_HCALIN> [GeV]", n.h_rho_hcalin->GetMean(), o.h_rho_hcalin->GetMean(), 5);
    row("<rho_HCALOUT> [GeV]", n.h_rho_hcalout->GetMean(), o.h_rho_hcalout->GetMean(), 5);

    ss << "-- jets --\n";
    row("jets in tree (all)", n.n_jets_raw, o.n_jets_raw, 0);
    row("jets/event (all)", n.n_events ? (double)n.n_jets_raw/n.n_events : 0,
                            o.n_events ? (double)o.n_jets_raw/o.n_events : 0, 5);
    row("jets pT>8, |eta|<0.8", n.n_jets_sel, o.n_jets_sel, 0);
    row("  per event", n.n_events ? (double)n.n_jets_sel/n.n_events : 0,
                       o.n_events ? (double)o.n_jets_sel/o.n_events : 0, 6);
    row("  <pT> [GeV]", n.n_jets_sel ? n.sum_pt_sel/n.n_jets_sel : 0,
                        o.n_jets_sel ? o.sum_pt_sel/o.n_jets_sel : 0, 4);
    row("jets + E>=0 & E_unsub>=0", n.n_jets_epos, o.n_jets_epos, 0);
    row("  E>=0 pass fraction", n.n_jets_sel ? (double)n.n_jets_epos/n.n_jets_sel : 0,
                                o.n_jets_sel ? (double)o.n_jets_epos/o.n_jets_sel : 0, 5);
    row("jets + UE fake-jet cut", n.n_jets_uepass, o.n_jets_uepass, 0);
    row("  UE cut pass fraction", n.n_jets_epos ? (double)n.n_jets_uepass/n.n_jets_epos : 0,
                                  o.n_jets_epos ? (double)o.n_jets_uepass/o.n_jets_epos : 0, 5);
    row("  overall pass fraction", n.n_jets_sel ? (double)n.n_jets_uepass/n.n_jets_sel : 0,
                                   o.n_jets_sel ? (double)o.n_jets_uepass/o.n_jets_sel : 0, 5);

    ss << "-- event yields --\n";
    row("events with >=1 selected jet", n.n_ev_ge1jet, o.n_ev_ge1jet, 0);
    row("  fraction", n.n_events ? (double)n.n_ev_ge1jet/n.n_events : 0,
                      o.n_events ? (double)o.n_ev_ge1jet/o.n_events : 0, 6);
    row("events with >=2 selected jets", n.n_ev_ge2jet, o.n_ev_ge2jet, 0);
    row("  fraction", n.n_events ? (double)n.n_ev_ge2jet/n.n_events : 0,
                      o.n_events ? (double)o.n_ev_ge2jet/o.n_events : 0, 6);
    row("<leading jet pT> [GeV]", n.h_lead_pt->GetMean(), o.h_lead_pt->GetMean(), 4);

    ss << "-- subtraction --\n";
    row("<(E_unsub-E)/E_unsub>", n.h_efrac->GetMean(), o.h_efrac->GetMean(), 5);
    row("<pT_unsub - pT_sub> [GeV]", n.h_ptexcess->GetMean(), o.h_ptexcess->GetMean(), 5);
    row("RMS pT_unsub - pT_sub [GeV]", n.h_ptexcess->GetRMS(), o.h_ptexcess->GetRMS(), 5);

    ss << std::string(86, '=') << "\n";
    ss << "\n(*) NOT directly comparable between the two productions:\n"
       << "    sumeT   : new = sum E_T over RHO-SUBTRACTED towers (fluctuates about 0, can be\n"
       << "              negative; components in sumeT_cemc/ihcal/ohcal). old = raw sum E_T.\n"
       << "    mbd_q   : new mbd_q and old mbdQ use different normalisations.\n"
       << "    cent, zvrtx and the rho values ARE on the same footing.\n"
       << "\n(!) jet_unsub_pT / jet_unsub_E changed meaning between the productions:\n"
       << "    new = true unsubtracted (raw tower) jet, so pT_unsub - pT_sub is the full UE\n"
       << "          (tens of GeV in central events).\n"
       << "    old = recomputed from already-subtracted constituents, so pT_unsub is BELOW\n"
       << "          pT_sub and pT_unsub - pT_sub is small and mostly negative.\n"
       << "    => the UE fake-jet cut pT_unsub - pT_sub > 40*exp(-0.038*cent), tuned on the\n"
       << "       old/Sub1 convention, is meaningless on the new sample as written and must\n"
       << "       be retuned before it is applied to the new production. The rows above are\n"
       << "       reported as-is to quantify that, not as a physics comparison.\n";
    return ss.str();
}

} // namespace

void make_data_rho_newold_compare(Long64_t stride = 1)
{
    gStyle->SetOptStat(0);
    const std::string outdir = "rho_systematics/data_newold";
    gSystem->mkdir((outdir + "/plots").c_str(), true);

    const std::string in_dir = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/input/data/";

    Sample s_new;
    s_new.tag = "new"; s_new.label = "new #rho (v004)";
    s_new.file = in_dir + "data_v004_20260821_calibrated_merged.root";
    s_new.br_mbd = "mbd_q";
    s_new.br_rho_cemc    = "rho_val_TowerRho_MULT_CEMC";
    s_new.br_rho_hcalin  = "rho_val_TowerRho_MULT_HCALIN";
    s_new.br_rho_hcalout = "rho_val_TowerRho_MULT_HCALOUT";
    s_new.has_minbias = true;

    Sample s_old;
    s_old.tag = "old"; s_old.label = "new #rho (v005)";
    s_old.file = in_dir + "run2auau_rho_jet.root";
    s_old.br_mbd = "mbdQ";
    s_old.br_rho_cemc    = "rho_cemc";
    s_old.br_rho_hcalin  = "rho_hcalin";
    s_old.br_rho_hcalout = "rho_hcalout";
    s_old.has_minbias = false;

    Hists hn = fill_sample(s_new, stride);
    Hists ho = fill_sample(s_old, stride);

    const std::string p = outdir + "/plots/";

    // ---- event level (per-event normalised) ----
    draw_overlay(norm_clone(hn.h_cent, hn.n_events, "_n"), norm_clone(ho.h_cent, ho.n_events, "_n"),
                 "centrality [%]", p + "ev_centrality.pdf", false, 0.8, 1.2, "events / event");
    draw_overlay(norm_clone(hn.h_zvrtx, hn.n_events, "_n"), norm_clone(ho.h_zvrtx, ho.n_events, "_n"),
                 "z_{vtx} [cm]", p + "ev_zvertex.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_sumeT, hn.n_events, "_n"), norm_clone(ho.h_sumeT, ho.n_events, "_n"),
                 "#sum E_{T} [GeV]", p + "ev_sumeT.pdf", true, 0.5, 1.5, "events / event");
    draw_overlay(norm_clone(hn.h_mbdq, hn.n_events, "_n"), norm_clone(ho.h_mbdq, ho.n_events, "_n"),
                 "MBD charge", p + "ev_mbd_charge.pdf", true, 0.5, 1.5, "events / event");
    draw_overlay(norm_clone(hn.h_rho_cemc, hn.n_events, "_n"), norm_clone(ho.h_rho_cemc, ho.n_events, "_n"),
                 "#rho_{CEMC} [GeV]", p + "ev_rho_cemc.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_rho_hcalin, hn.n_events, "_n"), norm_clone(ho.h_rho_hcalin, ho.n_events, "_n"),
                 "#rho_{HCALIN} [GeV]", p + "ev_rho_hcalin.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_rho_hcalout, hn.n_events, "_n"), norm_clone(ho.h_rho_hcalout, ho.n_events, "_n"),
                 "#rho_{HCALOUT} [GeV]", p + "ev_rho_hcalout.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_njet_all, hn.n_events, "_n"), norm_clone(ho.h_njet_all, ho.n_events, "_n"),
                 "jets / event (all)", p + "ev_njet_all.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_njet_sel, hn.n_events, "_n"), norm_clone(ho.h_njet_sel, ho.n_events, "_n"),
                 "jets / event (p_{T}>8, |#eta|<0.8)", p + "ev_njet_selected.pdf", true, 0.0, 2.0, "events / event");
    draw_overlay(norm_clone(hn.h_lead_pt, hn.n_events, "_n"), norm_clone(ho.h_lead_pt, ho.n_events, "_n"),
                 "leading jet p_{T} [GeV]", p + "ev_leading_jet_pt.pdf", true, 0.0, 2.0, "events / event");

    // ---- per jet (per-event normalised) ----
    draw_overlay(norm_clone(hn.h_pt_all, hn.n_events, "_n"), norm_clone(ho.h_pt_all, ho.n_events, "_n"),
                 "reco jet p_{T} [GeV]", p + "jet_reco_pt_all.pdf", true, 0.0, 2.0, "jets / event");
    draw_overlay(norm_clone(hn.h_pt_epos, hn.n_events, "_n"), norm_clone(ho.h_pt_epos, ho.n_events, "_n"),
                 "reco jet p_{T} [GeV] (E, E_{unsub} #geq 0)", p + "jet_reco_pt_Epositive.pdf", true, 0.0, 2.0, "jets / event");
    draw_overlay(norm_clone(hn.h_pt_unsub, hn.n_events, "_n"), norm_clone(ho.h_pt_unsub, ho.n_events, "_n"),
                 "unsubtracted jet p_{T} [GeV]", p + "jet_unsub_pt.pdf", true, 0.0, 2.0, "jets / event");
    draw_overlay(norm_clone(hn.h_eta, hn.n_events, "_n"), norm_clone(ho.h_eta, ho.n_events, "_n"),
                 "reco jet #eta", p + "jet_reco_eta.pdf", false, 0.0, 2.0, "jets / event");
    draw_overlay(norm_clone(hn.h_phi, hn.n_events, "_n"), norm_clone(ho.h_phi, ho.n_events, "_n"),
                 "reco jet #phi", p + "jet_reco_phi.pdf", false, 0.0, 2.0, "jets / event");
    draw_overlay(norm_clone(hn.h_efrac, hn.n_events, "_n"), norm_clone(ho.h_efrac, ho.n_events, "_n"),
                 "(E_{unsub} - E)/E_{unsub}", p + "jet_subtracted_E_fraction.pdf", true, 0.0, 3.0, "jets / event");
    draw_overlay(norm_clone(hn.h_ptexcess, hn.n_events, "_n"), norm_clone(ho.h_ptexcess, ho.n_events, "_n"),
                 "p_{T}^{unsub} - p_{T}^{sub} [GeV]", p + "jet_pt_excess.pdf", true, 0.0, 3.0, "jets / event");

    // ---- selection pass fractions (self-normalising, no event scaling needed) ----
    TH1D *fr_epos_new = (TH1D*) hn.h_pt_epos->Clone("h_Epos_passfrac_new");
    TH1D *fr_epos_old = (TH1D*) ho.h_pt_epos->Clone("h_Epos_passfrac_old");
    fr_epos_new->SetDirectory(nullptr); fr_epos_old->SetDirectory(nullptr);
    fr_epos_new->Divide(hn.h_pt_all); fr_epos_old->Divide(ho.h_pt_all);
    draw_overlay(fr_epos_new, fr_epos_old, "reco jet p_{T} [GeV]",
                 p + "sel_Epositive_pass_fraction.pdf", false, 0.8, 1.2, "E #geq 0 pass fraction");

    TH1D *fr_ue_new = (TH1D*) hn.h_uecut_pt_after->Clone("h_UEcut_passfrac_new");
    TH1D *fr_ue_old = (TH1D*) ho.h_uecut_pt_after->Clone("h_UEcut_passfrac_old");
    fr_ue_new->SetDirectory(nullptr); fr_ue_old->SetDirectory(nullptr);
    fr_ue_new->Divide(hn.h_uecut_pt_before); fr_ue_old->Divide(ho.h_uecut_pt_before);
    draw_overlay(fr_ue_new, fr_ue_old, "reco jet p_{T} [GeV]",
                 p + "sel_UEcut_pass_fraction.pdf", false, 0.8, 1.2, "UE fake-jet cut pass fraction");

    // spectrum actually entering the unfolding input (E>=0 and UE cut applied)
    draw_overlay(norm_clone(hn.h_uecut_pt_after, hn.n_events, "_n"), norm_clone(ho.h_uecut_pt_after, ho.n_events, "_n"),
                 "reco jet p_{T} [GeV] (all cuts)", p + "sel_final_pt_spectrum.pdf", true, 0.0, 2.0, "jets / event");

    // ---- profiles ----
    draw_overlay(hn.p_rho_vs_cent, ho.p_rho_vs_cent, "centrality [%]",
                 p + "prof_rho_cemc_vs_cent.pdf", false, 0.0, 2.0);
    draw_overlay(hn.p_njet_vs_cent, ho.p_njet_vs_cent, "centrality [%]",
                 p + "prof_njet_vs_cent.pdf", false, 0.0, 2.0);
    draw_overlay(hn.p_efrac_vs_pt, ho.p_efrac_vs_pt, "reco jet p_{T} [GeV]",
                 p + "prof_efrac_vs_pt.pdf", false, 0.0, 2.0);
    draw_overlay(hn.p_ptexcess_vs_cent, ho.p_ptexcess_vs_cent, "centrality [%]",
                 p + "prof_ptexcess_vs_cent.pdf", false, 0.0, 2.0);

    // ---- stats table ----
    const std::string table = stats_table(hn, ho);
    std::cout << table << std::endl;
    std::ofstream ftxt((outdir + "/data_rho_newold_stats.txt").c_str());
    ftxt << "new: " << s_new.file << "\nold: " << s_old.file << "\nstride: " << stride << "\n";
    ftxt << table;
    ftxt.close();

    // ---- output root ----
    TFile fout((outdir + "/data_rho_newold.root").c_str(), "RECREATE");
    for (const Hists *h : {&hn, &ho})
    {
        h->h_cent->Write(); h->h_zvrtx->Write(); h->h_sumeT->Write(); h->h_mbdq->Write();
        h->h_rho_cemc->Write(); h->h_rho_hcalin->Write(); h->h_rho_hcalout->Write();
        h->h_njet_all->Write(); h->h_njet_sel->Write(); h->h_lead_pt->Write();
        h->h_pt_all->Write(); h->h_pt_epos->Write(); h->h_pt_unsub->Write();
        h->h_eta->Write(); h->h_phi->Write(); h->h_efrac->Write(); h->h_ptexcess->Write();
        h->h_uecut_pt_before->Write(); h->h_uecut_pt_after->Write();
        h->p_rho_vs_cent->Write(); h->p_njet_vs_cent->Write();
        h->p_efrac_vs_pt->Write(); h->p_ptexcess_vs_cent->Write();
    }
    fr_epos_new->Write(); fr_epos_old->Write();
    fr_ue_new->Write(); fr_ue_old->Write();
    fout.Close();

    std::cout << "Wrote " << outdir << "/plots/*.pdf, "
              << outdir << "/data_rho_newold.root and "
              << outdir << "/data_rho_newold_stats.txt" << std::endl;
}
