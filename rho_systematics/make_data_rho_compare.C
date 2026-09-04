// Systematic comparison: rho-subtracted vs non-rho-subtracted jet
// reconstruction in AuAu DATA, r=0.3, mirroring make_sim_rho_compare.C.
//
// Reads the raw per-jet production tree "T" directly from
// /home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_[rho_]jet.root
// (295.2M events each -- same underlying events, reconstructed with and
// without UE/rho subtraction). No pthat-slice stitching is needed here:
// unlike the simulation, data is a single unweighted sample.
//
// Produces, for both rho and non-rho:
//   - reco jet pT, with and without the jet_E>=0 && jet_unsub_E>=0 cut
//   - eta, phi (E>=0-selected jets)
//   - subtracted energy fraction (jet_unsub_E - jet_E) / jet_unsub_E
//   - the "UE fake-jet" rejection cut used everywhere downstream
//     (pt_unsub - pt_sub > fcut(centrality), see makeDataTreeAuAu.C:152-193 /
//     getBackground.C:197-233 / anaConf.h:309-313), since that cut is a
//     background-affecting, rho-dependent selection that directly changes
//     which jets make it into the unfolding input.

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

const double reco_pt_min = 8.0;
const double abs_eta_max = 0.8;

struct PerJetHists
{
    TH1D *h_pt_all = nullptr;
    TH1D *h_pt_epos = nullptr;
    TH1D *h_eta = nullptr;
    TH1D *h_phi = nullptr;
    TH1D *h_efrac = nullptr;      // (unsub-sub)/unsub, E>=0 jets, unsub_E>0
    TH1D *h_uecut_pt_before = nullptr; // pT spectrum before the UE fake-jet cut
    TH1D *h_uecut_pt_after = nullptr;  // pT spectrum after it (same E>=0, eta selection)
};

PerJetHists make_hists(const std::string &tag)
{
    PerJetHists h;
    h.h_pt_all  = new TH1D(("h_pt_all_" + tag).c_str(),  ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_pt_epos = new TH1D(("h_pt_epos_" + tag).c_str(), ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_eta     = new TH1D(("h_eta_" + tag).c_str(),     ";reco jet #eta;counts", 40, -1.0, 1.0);
    h.h_phi     = new TH1D(("h_phi_" + tag).c_str(),     ";reco jet #phi;counts", 40, -3.2, 3.2);
    h.h_efrac   = new TH1D(("h_efrac_" + tag).c_str(),   ";(E_{unsub} - E) / E_{unsub};counts", 100, -0.2, 1.0);
    h.h_uecut_pt_before = new TH1D(("h_uecut_pt_before_" + tag).c_str(), ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    h.h_uecut_pt_after  = new TH1D(("h_uecut_pt_after_" + tag).c_str(),  ";reco jet p_{T} [GeV];counts", 100, 0, 100);
    for (TH1D *hh : {h.h_pt_all, h.h_pt_epos, h.h_eta, h.h_phi, h.h_efrac, h.h_uecut_pt_before, h.h_uecut_pt_after})
    { hh->SetDirectory(nullptr); hh->Sumw2(); }
    return h;
}

PerJetHists fill_per_jet(const std::string &tag, const std::string &raw_file, Long64_t stride)
{
    PerJetHists h = make_hists(tag);

    TFile *f = TFile::Open(raw_file.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << raw_file << std::endl; return h; }
    TTree *t = (TTree*) f->Get("T");

    std::vector<float> *jet_pT = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *jet_E = nullptr, *jet_unsub_E = nullptr, *jet_unsub_pT = nullptr;
    int centrality = 0;

    t->SetBranchStatus("*", 0);
    for (const char *b : {"jet_pT", "jet_eta", "jet_phi", "jet_E", "jet_unsub_E", "jet_unsub_pT", "cent"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("jet_pT", &jet_pT);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);
    t->SetBranchAddress("jet_E", &jet_E);
    t->SetBranchAddress("jet_unsub_E", &jet_unsub_E);
    t->SetBranchAddress("jet_unsub_pT", &jet_unsub_pT);
    t->SetBranchAddress("cent", &centrality);
    t->SetCacheSize(256 * 1024 * 1024);
    t->AddBranchToCache("*", kTRUE);

    // Same UE/fake-jet rejection cut used in makeDataTreeAuAu.C / getBackground.C / anaConf.h.
    TF1 fcut("fcut", "[0]+[1]*TMath::Exp(-[2]*x)", 0.0, 100.0);
    fcut.SetParameters(0.0, 40, 0.038);

    const Long64_t total_entries = t->GetEntries();
    const Long64_t entries = total_entries / stride;
    std::cout << "[" << tag << "] raw entries: " << total_entries << " (processing every " << stride
               << "th event, " << entries << " events)" << std::endl;

    for (Long64_t k = 0; k < entries; k++)
    {
        const Long64_t i = k * stride;
        t->GetEntry(i);
        if (k % (entries/10 == 0 ? 1 : entries/10) == 0)
            std::cout << "[" << tag << "] " << k << " / " << entries << "\r" << std::flush;

        const double cut_value = fcut.Eval(centrality);

        for (size_t j = 0; j < jet_pT->size(); j++)
        {
            const float pt = jet_pT->at(j);
            const float eta = jet_eta->at(j);
            if (pt < reco_pt_min || std::fabs(eta) > abs_eta_max) continue;

            h.h_pt_all->Fill(pt);

            const bool e_positive = (jet_E->at(j) >= 0) && (jet_unsub_E->at(j) >= 0);
            if (!e_positive) continue;

            h.h_pt_epos->Fill(pt);
            h.h_eta->Fill(eta);
            h.h_phi->Fill(jet_phi->at(j));

            const float unsub_E = jet_unsub_E->at(j);
            if (unsub_E > 0)
            {
                h.h_efrac->Fill((unsub_E - jet_E->at(j)) / unsub_E);
            }

            h.h_uecut_pt_before->Fill(pt);
            const float pt_unsub_excess = jet_unsub_pT->at(j) - pt;
            if (pt_unsub_excess <= cut_value)
            {
                h.h_uecut_pt_after->Fill(pt);
            }
        }
    }
    std::cout << std::endl;
    f->Close();
    return h;
}

void draw_overlay(TH1D *h_rho, TH1D *h_norho, const std::string &xtitle, const std::string &outpdf,
                   bool logy = true, double ratio_lo = 0.5, double ratio_hi = 1.5)
{
    h_rho->SetLineColor(kRed + 1);   h_rho->SetLineWidth(2);
    h_norho->SetLineColor(kBlack);   h_norho->SetLineWidth(2);

    TCanvas c("c", "c", 600, 700);
    TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1.0);
    TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.3);
    p1->SetBottomMargin(0.02); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.3);
    p1->Draw(); p2->Draw();

    p1->cd();
    if (logy) p1->SetLogy();
    h_norho->GetXaxis()->SetLabelSize(0);
    h_norho->SetTitle("");
    h_norho->Draw("hist e");
    h_rho->Draw("hist e same");
    TLegend leg(0.55, 0.65, 0.88, 0.85);
    leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(h_norho, "no #rho subtraction", "l");
    leg.AddEntry(h_rho, "#rho subtracted", "l");
    leg.Draw();

    p2->cd();
    TH1D *ratio = (TH1D*) h_rho->Clone("ratio_tmp");
    ratio->SetDirectory(nullptr);
    ratio->Divide(h_norho);
    ratio->SetTitle((";" + xtitle + ";#rho / no-#rho").c_str());
    ratio->GetYaxis()->SetRangeUser(ratio_lo, ratio_hi);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->SetMarkerStyle(20);
    ratio->Draw("ep");
    TLine line(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
    line.SetLineStyle(2);
    line.Draw();

    c.Print(outpdf.c_str());
    delete ratio;
}

} // namespace

void make_data_rho_compare(Long64_t stride = 1)
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("rho_systematics/data/plots", true);

    const std::string raw_dir = "/sphenix/user/tmengel/JetUESub-JSTG-TF03/macros/data/rootfiles";

    PerJetHists h_rho   = fill_per_jet("rho",   raw_dir + "run2auau_rho_jet.root", stride);
    PerJetHists h_norho = fill_per_jet("norho", raw_dir + "run2auau_jet.root", stride);

    draw_overlay(h_rho.h_pt_all,  h_norho.h_pt_all,  "reco jet p_{T} [GeV]",
                 "rho_systematics/data/plots/data_reco_pt_all.pdf");
    draw_overlay(h_rho.h_pt_epos, h_norho.h_pt_epos, "reco jet p_{T} [GeV] (E, E_{unsub} #geq 0)",
                 "rho_systematics/data/plots/data_reco_pt_Epositive.pdf");
    draw_overlay(h_rho.h_eta, h_norho.h_eta, "reco jet #eta",
                 "rho_systematics/data/plots/data_reco_eta.pdf", false, 0.7, 1.3);
    draw_overlay(h_rho.h_phi, h_norho.h_phi, "reco jet #phi",
                 "rho_systematics/data/plots/data_reco_phi.pdf", false, 0.7, 1.3);
    draw_overlay(h_rho.h_efrac, h_norho.h_efrac, "subtracted energy fraction",
                 "rho_systematics/data/plots/data_subtracted_E_fraction.pdf", true, 0.0, 5.0);

    TH1D *h_epos_frac_rho   = (TH1D*) h_rho.h_pt_epos->Clone("h_epos_frac_rho");
    TH1D *h_epos_frac_norho = (TH1D*) h_norho.h_pt_epos->Clone("h_epos_frac_norho");
    h_epos_frac_rho->Divide(h_rho.h_pt_all);
    h_epos_frac_norho->Divide(h_norho.h_pt_all);
    draw_overlay(h_epos_frac_rho, h_epos_frac_norho, "reco jet p_{T} [GeV]",
                 "rho_systematics/data/plots/data_Epositive_pass_fraction.pdf", false, 0.0, 2.0);

    TH1D *h_uecut_frac_rho   = (TH1D*) h_rho.h_uecut_pt_after->Clone("h_uecut_frac_rho");
    TH1D *h_uecut_frac_norho = (TH1D*) h_norho.h_uecut_pt_after->Clone("h_uecut_frac_norho");
    h_uecut_frac_rho->Divide(h_rho.h_uecut_pt_before);
    h_uecut_frac_norho->Divide(h_norho.h_uecut_pt_before);
    draw_overlay(h_uecut_frac_rho, h_uecut_frac_norho, "reco jet p_{T} [GeV]",
                 "rho_systematics/data/plots/data_UEcut_pass_fraction.pdf", false, 0.0, 1.5);

    TFile fout("rho_systematics/data/data_rho_vs_norho.root", "RECREATE");
    h_rho.h_pt_all->Write(); h_rho.h_pt_epos->Write(); h_rho.h_eta->Write(); h_rho.h_phi->Write();
    h_rho.h_efrac->Write(); h_rho.h_uecut_pt_before->Write(); h_rho.h_uecut_pt_after->Write();
    h_norho.h_pt_all->Write(); h_norho.h_pt_epos->Write(); h_norho.h_eta->Write(); h_norho.h_phi->Write();
    h_norho.h_efrac->Write(); h_norho.h_uecut_pt_before->Write(); h_norho.h_uecut_pt_after->Write();
    h_epos_frac_rho->Write(); h_epos_frac_norho->Write();
    h_uecut_frac_rho->Write(); h_uecut_frac_norho->Write();
    fout.Close();

    std::cout << "Wrote rho_systematics/data/plots/*.pdf and rho_systematics/data/data_rho_vs_norho.root" << std::endl;
}
