// Systematic comparison: rho-subtracted vs non-rho-subtracted jet reconstruction
// in the HIJING+PYTHIA embedded simulation, r=0.3.
//
// Two input layers, both weight-stitched across the jet10/jet20/jet30 pthat
// slices using the same recipe as createResponse_noempty_AA.cxx
// (rootfiles/TNTUPLE_DIJET_SIM_r03_jetN_hijing_alljet_[rho_]jet.root, lines
// ~104-168 there):
//
//   scale_factor[i] = (n_events[jet30]/n_events[i]) * (cs_i / cs_30),  i = 10,20
//   scale_factor[jet30] = 1.0
//
// 1) Per-jet quantities (reco pT with/without the jet_E>=0 && jet_unsub_E>=0
//    positivity cut, eta, phi) come from the single MERGED raw production
//    tree "T" in rootfiles/../v001_20260720/hijing_alljet_[rho_]jet.root
//    (8.6M events = sum of the three TNTUPLE stats-tree nevents). Since this
//    raw tree has no per-slice tag, each event is routed to jet10/20/30 by
//    its leading truth-jet pT against the same sample_boundary edges
//    createResponse uses (read from configs/binning_AA.config), then scaled.
//
// 2) Matching efficiency vs truth pT comes from the already-computed
//    TEfficiency objects he_pt_truth_N (N = centrality bin 0-9, 10% wide,
//    see makeMatchedTreesInclusiveAuAu.C) stored in the per-jetN TNTUPLE
//    files. Passed/Total histograms are weight-summed across jetN and N
//    (matching the "no hard centrality cut on the MC" policy in
//    createResponse) and recombined into one TEfficiency per rho/non-rho set.

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

const double cs_10 = 0.000003997;
const double cs_20 = 6.218e-8;
const double cs_30 = 2.505e-9;

// From configs/binning_AA.config via read_binning::get_sample_boundary(),
// r=0.3 dijet binning (see _tmp_dump_config.C probe used to derive this).
const double sample_boundary[4] = {14.4914, 20.8505, 30.0, 100.0};

const double reco_pt_min = 8.0;   // read_binning::get_reco_pt_min_cut()
const double abs_eta_max = 0.8;   // read_binning::get_eta_cut_bkg() / get_abs_eta_acceptance(0.3)

struct SimSet
{
    std::string tag;
    std::string raw_file;
    std::string tntuple_jet10, tntuple_jet20, tntuple_jet30;
};

double n_events[3] = {0, 0, 0};

void load_n_events(const std::string &f10, const std::string &f20, const std::string &f30)
{
    const char *fname[3] = { f10.c_str(), f20.c_str(), f30.c_str() };
    for (int i = 0; i < 3; i++)
    {
        TFile *f = TFile::Open(fname[i], "READ");
        TNtuple *ts = (TNtuple*) f->Get("tn_stats");
        float b_n_events = 0;
        ts->SetBranchAddress("nevents", &b_n_events);
        ts->GetEntry(0);
        n_events[i] = b_n_events;
        f->Close();
    }
    std::cout << "n_events (jet10, jet20, jet30): " << n_events[0] << " " << n_events[1] << " " << n_events[2] << std::endl;
}

double scale_factor_for_maxpt(double maxpttruth, int &isample_out)
{
    static const double cs[3] = { cs_10, cs_20, cs_30 };
    for (int i = 0; i < 3; i++)
    {
        if (maxpttruth >= sample_boundary[i] && maxpttruth < sample_boundary[i + 1])
        {
            isample_out = i;
            if (i == 2) return 1.0;
            return (n_events[2] / n_events[i]) * (cs[i] / cs[2]);
        }
    }
    isample_out = -1;
    return -1.0;
}

struct PerJetHists
{
    TH1D *h_pt_all = nullptr;
    TH1D *h_pt_epos = nullptr;
    TH1D *h_eta = nullptr;
    TH1D *h_phi = nullptr;
};

PerJetHists make_hists(const std::string &tag)
{
    PerJetHists h;
    h.h_pt_all  = new TH1D(("h_pt_all_" + tag).c_str(),  ";reco jet p_{T} [GeV];weighted counts", 100, 0, 100);
    h.h_pt_epos = new TH1D(("h_pt_epos_" + tag).c_str(), ";reco jet p_{T} [GeV];weighted counts", 100, 0, 100);
    h.h_eta     = new TH1D(("h_eta_" + tag).c_str(),     ";reco jet #eta;weighted counts", 40, -1.0, 1.0);
    h.h_phi     = new TH1D(("h_phi_" + tag).c_str(),     ";reco jet #phi;weighted counts", 40, -3.2, 3.2);
    for (TH1D *hh : {h.h_pt_all, h.h_pt_epos, h.h_eta, h.h_phi}) { hh->SetDirectory(nullptr); hh->Sumw2(); }
    return h;
}

PerJetHists fill_per_jet(const std::string &tag, const std::string &raw_file)
{
    PerJetHists h = make_hists(tag);

    TFile *f = TFile::Open(raw_file.c_str(), "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << raw_file << std::endl; return h; }
    TTree *t = (TTree*) f->Get("T");

    std::vector<float> *jet_pT = nullptr, *jet_eta = nullptr, *jet_phi = nullptr;
    std::vector<float> *jet_E = nullptr, *jet_unsub_E = nullptr;
    std::vector<float> *truth_jet_pT = nullptr, *truth_jet_eta = nullptr;

    t->SetBranchStatus("*", 0);
    for (const char *b : {"jet_pT", "jet_eta", "jet_phi", "jet_E", "jet_unsub_E", "truth_jet_pT", "truth_jet_eta"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("jet_pT", &jet_pT);
    t->SetBranchAddress("jet_eta", &jet_eta);
    t->SetBranchAddress("jet_phi", &jet_phi);
    t->SetBranchAddress("jet_E", &jet_E);
    t->SetBranchAddress("jet_unsub_E", &jet_unsub_E);
    t->SetBranchAddress("truth_jet_pT", &truth_jet_pT);
    t->SetBranchAddress("truth_jet_eta", &truth_jet_eta);

    const Long64_t entries = t->GetEntries();
    std::cout << "[" << tag << "] raw entries: " << entries << std::endl;

    int n_out_of_range = 0;
    for (Long64_t i = 0; i < entries; i++)
    {
        t->GetEntry(i);
        if (i % (entries/10 == 0 ? 1 : entries/10) == 0)
            std::cout << "[" << tag << "] " << i << " / " << entries << "\r" << std::flush;

        double maxpttruth = -1;
        for (size_t k = 0; k < truth_jet_pT->size(); k++)
        {
            if (std::fabs(truth_jet_eta->at(k)) > abs_eta_max) continue;
            if (truth_jet_pT->at(k) > maxpttruth) maxpttruth = truth_jet_pT->at(k);
        }
        if (maxpttruth < 0) continue;

        int isample = -1;
        const double w = scale_factor_for_maxpt(maxpttruth, isample);
        if (isample < 0) { n_out_of_range++; continue; }

        for (size_t j = 0; j < jet_pT->size(); j++)
        {
            const float pt = jet_pT->at(j);
            const float eta = jet_eta->at(j);
            if (pt < reco_pt_min || std::fabs(eta) > abs_eta_max) continue;

            h.h_pt_all->Fill(pt, w);

            const bool e_positive = (jet_E->at(j) >= 0) && (jet_unsub_E->at(j) >= 0);
            if (e_positive)
            {
                h.h_pt_epos->Fill(pt, w);
                h.h_eta->Fill(eta, w);
                h.h_phi->Fill(jet_phi->at(j), w);
            }
        }
    }
    std::cout << std::endl;
    std::cout << "[" << tag << "] events with maxpttruth outside sample_boundary range: " << n_out_of_range << std::endl;
    f->Close();
    return h;
}

// Weight-sum he_pt_truth_N passed/total across jetN in {10,20,30} and
// centrality bin N in [0,9], using the same scale_factor recipe as above.
TEfficiency * combined_matching_efficiency(const std::string &tag,
                                            const std::string &f10, const std::string &f20, const std::string &f30)
{
    const char *fname[3] = { f10.c_str(), f20.c_str(), f30.c_str() };
    static const double cs[3] = { cs_10, cs_20, cs_30 };

    TH1D *passed_sum = nullptr;
    TH1D *total_sum = nullptr;

    for (int isample = 0; isample < 3; isample++)
    {
        TFile *f = TFile::Open(fname[isample], "READ");
        const double w = (isample == 2) ? 1.0 : (n_events[2] / n_events[isample]) * (cs[isample] / cs[2]);

        for (int N = 0; N <= 9; N++)
        {
            TEfficiency *e = (TEfficiency*) f->Get(Form("he_pt_truth_%d", N));
            if (!e) continue;
            TH1D *passed = (TH1D*) e->GetCopyPassedHisto();
            TH1D *total  = (TH1D*) e->GetCopyTotalHisto();
            passed->Scale(w);
            total->Scale(w);
            passed->SetDirectory(nullptr);
            total->SetDirectory(nullptr);

            if (!passed_sum)
            {
                passed_sum = (TH1D*) passed->Clone(("h_passed_sum_" + tag).c_str());
                total_sum  = (TH1D*) total->Clone(("h_total_sum_" + tag).c_str());
                passed_sum->SetDirectory(nullptr);
                total_sum->SetDirectory(nullptr);
            }
            else
            {
                passed_sum->Add(passed);
                total_sum->Add(total);
            }
            delete passed;
            delete total;
        }
        f->Close();
    }

    // TEfficiency requires the passed histogram bin content to be an integer
    // count under the default binomial statistic option; our weighted sums
    // are non-integer, so switch to the normal (Gaussian) approximation,
    // which TEfficiency supports for weighted entries.
    TEfficiency *eff = new TEfficiency(*passed_sum, *total_sum);
    eff->SetStatisticOption(TEfficiency::kFNormal);
    eff->SetName(("he_matching_eff_" + tag).c_str());
    return eff;
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

void make_sim_rho_compare()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("rho_systematics/sim/plots", true);

    const std::string tntuple_dir = "rootfiles/";
    const std::string raw_dir = "/home/tmengel/PPG14/rootfiles/v001_20260720/";

    load_n_events(tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet10_hijing_alljet_rho_jet.root",
                  tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet20_hijing_alljet_rho_jet.root",
                  tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet30_hijing_alljet_rho_jet.root");

    // ---- per-jet quantities: pT (all / E>=0), eta, phi ----
    PerJetHists h_rho   = fill_per_jet("rho",   raw_dir + "hijing_alljet_rho_jet.root");
    PerJetHists h_norho = fill_per_jet("norho", raw_dir + "hijing_alljet_jet.root");

    draw_overlay(h_rho.h_pt_all,  h_norho.h_pt_all,  "reco jet p_{T} [GeV]",
                 "rho_systematics/sim/plots/sim_reco_pt_all.pdf");
    draw_overlay(h_rho.h_pt_epos, h_norho.h_pt_epos, "reco jet p_{T} [GeV] (E, E_{unsub} #geq 0)",
                 "rho_systematics/sim/plots/sim_reco_pt_Epositive.pdf");
    draw_overlay(h_rho.h_eta, h_norho.h_eta, "reco jet #eta",
                 "rho_systematics/sim/plots/sim_reco_eta.pdf", false, 0.7, 1.3);
    draw_overlay(h_rho.h_phi, h_norho.h_phi, "reco jet #phi",
                 "rho_systematics/sim/plots/sim_reco_phi.pdf", false, 0.7, 1.3);

    // fraction of jets surviving the E>=0 cut, vs pT, as an explicit "with vs
    // without E>=0" efficiency-style comparison
    TH1D *h_epos_frac_rho   = (TH1D*) h_rho.h_pt_epos->Clone("h_epos_frac_rho");
    TH1D *h_epos_frac_norho = (TH1D*) h_norho.h_pt_epos->Clone("h_epos_frac_norho");
    h_epos_frac_rho->Divide(h_rho.h_pt_all);
    h_epos_frac_norho->Divide(h_norho.h_pt_all);
    draw_overlay(h_epos_frac_rho, h_epos_frac_norho, "reco jet p_{T} [GeV]",
                 "rho_systematics/sim/plots/sim_Epositive_pass_fraction.pdf", false, 0.0, 2.0);

    // ---- matching efficiency vs truth pT ----
    TEfficiency *eff_rho = combined_matching_efficiency("rho",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet10_hijing_alljet_rho_jet.root",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet20_hijing_alljet_rho_jet.root",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet30_hijing_alljet_rho_jet.root");
    TEfficiency *eff_norho = combined_matching_efficiency("norho",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet10_hijing_alljet_jet.root",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet20_hijing_alljet_jet.root",
        tntuple_dir + "TNTUPLE_DIJET_SIM_r03_jet30_hijing_alljet_jet.root");

    {
        TCanvas c("c_eff", "c_eff", 600, 700);
        TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1.0);
        TPad *p2 = new TPad("p2", "p2", 0, 0.0, 1, 0.3);
        p1->SetBottomMargin(0.02); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.3);
        p1->Draw(); p2->Draw();

        p1->cd();
        eff_norho->SetLineColor(kBlack); eff_norho->SetMarkerColor(kBlack); eff_norho->SetMarkerStyle(20);
        eff_rho->SetLineColor(kRed + 1); eff_rho->SetMarkerColor(kRed + 1); eff_rho->SetMarkerStyle(21);
        eff_norho->Draw("AP");
        gPad->Update();
        eff_norho->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 1.05);
        eff_norho->GetPaintedGraph()->GetXaxis()->SetLabelSize(0);
        eff_rho->Draw("P same");
        TLegend leg(0.55, 0.2, 0.88, 0.4);
        leg.SetBorderSize(0); leg.SetFillStyle(0);
        leg.AddEntry(eff_norho, "no #rho subtraction", "lep");
        leg.AddEntry(eff_rho, "#rho subtracted", "lep");
        leg.Draw();

        p2->cd();
        TH1D *h_ratio = (TH1D*) eff_rho->GetCopyTotalHisto()->Clone("h_eff_ratio");
        h_ratio->SetDirectory(nullptr);
        for (int b = 1; b <= h_ratio->GetNbinsX(); b++)
        {
            const double e1 = eff_rho->GetEfficiency(b);
            const double e2 = eff_norho->GetEfficiency(b);
            h_ratio->SetBinContent(b, e2 > 0 ? e1/e2 : 0);
            h_ratio->SetBinError(b, 0);
        }
        h_ratio->SetTitle(";truth jet p_{T} [GeV];#rho / no-#rho eff");
        h_ratio->GetYaxis()->SetRangeUser(0.9, 1.1);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->SetMarkerStyle(20);
        h_ratio->Draw("p");
        TLine line(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
        line.SetLineStyle(2);
        line.Draw();

        c.Print("rho_systematics/sim/plots/sim_matching_efficiency.pdf");
    }

    TFile fout("rho_systematics/sim/sim_rho_vs_norho.root", "RECREATE");
    h_rho.h_pt_all->Write(); h_rho.h_pt_epos->Write(); h_rho.h_eta->Write(); h_rho.h_phi->Write();
    h_norho.h_pt_all->Write(); h_norho.h_pt_epos->Write(); h_norho.h_eta->Write(); h_norho.h_phi->Write();
    h_epos_frac_rho->Write(); h_epos_frac_norho->Write();
    eff_rho->Write(); eff_norho->Write();
    fout.Close();

    std::cout << "Wrote rho_systematics/sim/plots/*.pdf and rho_systematics/sim/sim_rho_vs_norho.root" << std::endl;
}
