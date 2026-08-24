// Redraw the rho-vs-norho overlay PDFs from the already-saved histograms in
// sim_rho_vs_norho.root / data_rho_vs_norho.root, with the y-axis range set
// from the max of BOTH histograms (the original draw_overlay() in
// make_sim_rho_compare.C / make_data_rho_compare.C only auto-scaled off the
// first-drawn histogram, silently clipping off any bins where the second
// histogram ran higher -- eta/phi in particular were badly clipped).

#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

void draw_overlay(TH1D *h_rho, TH1D *h_norho, const std::string &xtitle, const std::string &outpdf,
                   bool logy, double ratio_lo, double ratio_hi)
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
    const double ymax = std::max(h_rho->GetMaximum(), h_norho->GetMaximum());
    h_norho->SetMaximum(logy ? ymax * 3 : ymax * 1.3);
    if (!logy) h_norho->SetMinimum(0);
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

void redraw_one(const std::string &rootfile, const std::string &tag_prefix, const std::string &outdir)
{
    TFile f(rootfile.c_str(), "READ");

    auto get = [&](const std::string &name) {
        TH1D *h = (TH1D*) f.Get(name.c_str());
        if (h) h->SetDirectory(nullptr);
        return h;
    };

    TH1D *pt_all_rho    = get("h_pt_all_rho"),    *pt_all_norho    = get("h_pt_all_norho");
    TH1D *pt_epos_rho   = get("h_pt_epos_rho"),   *pt_epos_norho   = get("h_pt_epos_norho");
    TH1D *eta_rho       = get("h_eta_rho"),       *eta_norho       = get("h_eta_norho");
    TH1D *phi_rho       = get("h_phi_rho"),       *phi_norho       = get("h_phi_norho");
    TH1D *eposfrac_rho  = get("h_epos_frac_rho"), *eposfrac_norho  = get("h_epos_frac_norho");

    draw_overlay(pt_all_rho, pt_all_norho, "reco jet p_{T} [GeV]",
                 outdir + "/" + tag_prefix + "_reco_pt_all.pdf", true, 0.5, 1.5);
    draw_overlay(pt_epos_rho, pt_epos_norho, "reco jet p_{T} [GeV] (E, E_{unsub} #geq 0)",
                 outdir + "/" + tag_prefix + "_reco_pt_Epositive.pdf", true, 0.5, 1.5);
    draw_overlay(eta_rho, eta_norho, "reco jet #eta",
                 outdir + "/" + tag_prefix + "_reco_eta.pdf", false, 0.5, 3.0);
    draw_overlay(phi_rho, phi_norho, "reco jet #phi",
                 outdir + "/" + tag_prefix + "_reco_phi.pdf", false, 0.5, 3.0);
    draw_overlay(eposfrac_rho, eposfrac_norho, "reco jet p_{T} [GeV]",
                 outdir + "/" + tag_prefix + "_Epositive_pass_fraction.pdf", false, 0.0, 2.0);

    TH1D *efrac_rho = get("h_efrac_rho");
    TH1D *efrac_norho = get("h_efrac_norho");
    if (efrac_rho && efrac_norho)
        draw_overlay(efrac_rho, efrac_norho, "subtracted energy fraction",
                     outdir + "/" + tag_prefix + "_subtracted_E_fraction.pdf", true, 0.0, 5.0);

    TH1D *uecut_rho = get("h_uecut_frac_rho");
    TH1D *uecut_norho = get("h_uecut_frac_norho");
    if (uecut_rho && uecut_norho)
        draw_overlay(uecut_rho, uecut_norho, "reco jet p_{T} [GeV]",
                     outdir + "/" + tag_prefix + "_UEcut_pass_fraction.pdf", false, 0.0, 1.5);
}

} // namespace

void redraw_overlays()
{
    gStyle->SetOptStat(0);
    redraw_one("rho_systematics/sim/sim_rho_vs_norho.root", "sim", "rho_systematics/sim/plots");
    redraw_one("rho_systematics/data/data_rho_vs_norho.root", "data", "rho_systematics/data/plots");
    std::cout << "Redrew overlay PDFs with corrected y-axis scaling." << std::endl;
}
