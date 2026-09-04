// Compare reco pT spectra of matched jets between the rho-subtracted and
// non-rho-subtracted "alljet" TNTUPLE samples, stitching jet10/jet20/jet30
// together by cross-section/n_events weight exactly as
// createResponse_noempty_AA.cxx does (see lines ~104-168 there).
//
//   scale_factor[i] = (n_events[jet30]/n_events[i]) * (cs_i / cs_30),  i = 10,20
//   scale_factor[jet30] = 1.0
//
// For each sample set (rho, no-rho) this fills one histogram per matched jet
// (both legs of every real pair, leading+subleading together, all three
// pthat slices weighted and added) using reco pT, then overlays the two
// stitched spectra and their ratio.

#include <iostream>
#include <string>

#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TStyle.h"

namespace {

// cross sections (pb), same values as createResponse_noempty_AA.cxx
const double cs_10 = 0.000003997;
const double cs_20 = 6.218e-8;
const double cs_30 = 2.505e-9;

// Fill a stitched (weight-combined) matched-jet reco pT spectrum from the
// three jetN files belonging to one sample set (e.g. all "rho" or all
// plain). Returns a new TH1D owned by the caller.
TH1D * stitched_reco_pt(const std::string &tag,
                         const std::string &f10, const std::string &f20, const std::string &f30,
                         int nbins, double lo, double hi)
{
    const char *fname[3] = { f10.c_str(), f20.c_str(), f30.c_str() };
    const double cs[3] = { cs_10, cs_20, cs_30 };

    TFile *fin[3];
    TNtuple *tn[3];
    float pt1_reco[3], pt2_reco[3], match[3];
    float n_events[3];

    for (int i = 0; i < 3; i++)
    {
        fin[i] = TFile::Open(fname[i], "READ");
        if (!fin[i] || fin[i]->IsZombie())
        {
            std::cerr << "Cannot open " << fname[i] << std::endl;
            return nullptr;
        }
        tn[i] = (TNtuple*) fin[i]->Get("tn_match");
        tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
        tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
        tn[i]->SetBranchAddress("matched", &match[i]);

        TNtuple *ts = (TNtuple*) fin[i]->Get("tn_stats");
        float b_n_events = 0;
        ts->SetBranchAddress("nevents", &b_n_events);
        ts->GetEntry(0);
        n_events[i] = b_n_events;
    }

    double scale_factor[3];
    scale_factor[0] = (n_events[2] / n_events[0]) * cs[0] / cs[2];
    scale_factor[1] = (n_events[2] / n_events[1]) * cs[1] / cs[2];
    scale_factor[2] = 1.0;

    std::cout << "[" << tag << "] n_events: " << n_events[0] << " " << n_events[1] << " " << n_events[2] << std::endl;
    std::cout << "[" << tag << "] scale_factor: " << scale_factor[0] << " " << scale_factor[1] << " " << scale_factor[2] << std::endl;

    TH1D *h = new TH1D(Form("h_matched_reco_pt_%s", tag.c_str()),
                        ";reco jet p_{T} [GeV];weighted counts", nbins, lo, hi);
    h->SetDirectory(nullptr);
    h->Sumw2();

    for (int isample = 0; isample < 3; isample++)
    {
        const Long64_t entries = tn[isample]->GetEntries();
        for (Long64_t i = 0; i < entries; i++)
        {
            tn[isample]->GetEntry(i);
            if (!match[isample]) continue;

            h->Fill(pt1_reco[isample], scale_factor[isample]);
            h->Fill(pt2_reco[isample], scale_factor[isample]);
        }
        fin[isample]->Close();
    }

    return h;
}

} // namespace

void compare_reco_pt_rho()
{
    gStyle->SetOptStat(0);

    const std::string dir = "rootfiles/";

    TH1D *h_rho = stitched_reco_pt("rho",
        dir + "TNTUPLE_DIJET_SIM_r03_jet10_hijing_alljet_rho_jet.root",
        dir + "TNTUPLE_DIJET_SIM_r03_jet20_hijing_alljet_rho_jet.root",
        dir + "TNTUPLE_DIJET_SIM_r03_jet30_hijing_alljet_rho_jet.root",
        100, 0, 100);

    TH1D *h_norho = stitched_reco_pt("norho",
        dir + "TNTUPLE_DIJET_SIM_r03_jet10_hijing_alljet_jet.root",
        dir + "TNTUPLE_DIJET_SIM_r03_jet20_hijing_alljet_jet.root",
        dir + "TNTUPLE_DIJET_SIM_r03_jet30_hijing_alljet_jet.root",
        100, 0, 100);

    if (!h_rho || !h_norho)
    {
        std::cerr << "Failed to build one or both stitched spectra." << std::endl;
        return;
    }

    h_rho->SetLineColor(kRed + 1);
    h_rho->SetLineWidth(2);
    h_norho->SetLineColor(kBlack);
    h_norho->SetLineWidth(2);

    TCanvas *c = new TCanvas("c_compare_reco_pt", "c_compare_reco_pt", 600, 700);
    TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1.0);
    TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.3);
    p1->SetBottomMargin(0.02);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetLogy();
    h_norho->GetXaxis()->SetLabelSize(0);
    h_norho->Draw("hist e");
    h_rho->Draw("hist e same");

    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_norho, "no #rho subtraction", "l");
    leg->AddEntry(h_rho, "#rho subtracted", "l");
    leg->Draw();

    p2->cd();
    TH1D *h_ratio = (TH1D*) h_rho->Clone("h_ratio_rho_over_norho");
    h_ratio->Divide(h_norho);
    h_ratio->SetTitle(";reco jet p_{T} [GeV];#rho / no-#rho");
    h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->SetMarkerStyle(20);
    h_ratio->Draw("ep");

    TLine *line1 = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line1->SetLineStyle(2);
    line1->Draw();

    c->Print("comparison_plots/compare_reco_pt_rho_vs_norho.pdf");

    TFile *fout = new TFile("comparison_plots/compare_reco_pt_rho_vs_norho.root", "RECREATE");
    h_rho->Write();
    h_norho->Write();
    h_ratio->Write();
    fout->Close();

    std::cout << "Wrote comparison_plots/compare_reco_pt_rho_vs_norho.pdf" << std::endl;
}
