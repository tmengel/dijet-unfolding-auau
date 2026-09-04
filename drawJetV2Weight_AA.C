#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TParameter.h"
#include "TProfile.h"
#include "TStyle.h"

#include "dlUtility.h"
#include "PlotUtils.h"
#include "read_binning.h"

// Standalone QA plot: draws h_jetv2_weight, the per-truth-pair jet-v2
// cross-check weight 1 + jetv2_amp*(cos(2*dpsi2_lead) + cos(2*dpsi2_sub))
// that createResponse_exclusive_v2_AA.cxx applies when JETV2_SCALE != 1
// (see its "Debug: distribution of the jet-v2 cross-check weight" block and
// the JetV2 closure printout right after the main fill loop), plus, when
// available, h_jetv2_weight_vs_dpsi2_leg{1,2} -- <weight> profiled against
// each leg's own truth angle to the event plane. This second panel is the
// actual proof the injection is doing what it claims: a flat line there
// would mean the applied weight does not depend on dpsi2 at all, whereas
// the expected 1 + jetv2_amp*cos(2*dpsi2) modulation (overlaid using the
// jetv2_amp TParameter<double> the response file also carries) means the
// weight is doing exactly the v2 reweighting the cross-check is supposed to
// inject. Reads straight out of the response file
//   response_matrices/response_matrix_<system>_r0<cone>_<sys_name>.root,
// same file compareJetV2Xj_AA.C's response-side counterpart
// (drawResponse_AA.C) already reads for the same sys_name.
//
// The dpsi2 profiles and jetv2_amp are only present in response files built
// after they were added to createResponse_exclusive_v2_AA.cxx -- if missing
// (an older cached response file), this still draws the weight-distribution
// panel alone and says so on the console, rather than aborting.
void drawJetV2Weight_AA(
  const int cone_size = 3,
  const int centrality_bin = 0,
  const std::string sys_name = "JETV2",
  const std::string configfile = "/sphenix/user/tmengel/dijet-ana-auau/macros/unfolding/dijet-unfolding-auau/configs/binning_AA.config"
)
{
  PlotUtils::set_sphenix_style();
  gStyle->SetOptStat(0);

  read_binning rb(configfile.c_str());
  const std::string system_string = rb.get_system_string(centrality_bin);

  const TString path = Form("%s/response_matrices/response_matrix_%s_r%02d_%s.root",
                             rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str());
  TFile *f = TFile::Open(path, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "drawJetV2Weight_AA: cannot open " << path << std::endl;
    return;
  }

  TH1D *h = (TH1D *)f->Get("h_jetv2_weight");
  if (!h)
  {
    std::cerr << "drawJetV2Weight_AA: missing h_jetv2_weight in " << path << std::endl;
    f->Close();
    return;
  }
  h = (TH1D *)h->Clone(Form("h_jetv2_weight_%s", sys_name.c_str()));
  h->SetDirectory(nullptr);

  if (h->GetEntries() <= 0)
  {
    std::cerr << "drawJetV2Weight_AA: h_jetv2_weight in " << path
              << " is empty (JETV2_SCALE was 1, i.e. cross-check off for this file), nothing to plot." << std::endl;
    f->Close();
    return;
  }

  TProfile *p1 = (TProfile *)f->Get("h_jetv2_weight_vs_dpsi2_leg1");
  TProfile *p2 = (TProfile *)f->Get("h_jetv2_weight_vs_dpsi2_leg2");
  TParameter<double> *amp_param = (TParameter<double> *)f->Get("jetv2_amp");
  const bool haveDpsi2 = (p1 && p2 && amp_param);
  if (p1) { p1 = (TProfile *)p1->Clone(Form("h_jetv2_weight_vs_dpsi2_leg1_%s", sys_name.c_str())); p1->SetDirectory(nullptr); }
  if (p2) { p2 = (TProfile *)p2->Clone(Form("h_jetv2_weight_vs_dpsi2_leg2_%s", sys_name.c_str())); p2->SetDirectory(nullptr); }
  const double jetv2_amp = amp_param ? amp_param->GetVal() : 0.0;
  f->Close();

  if (!haveDpsi2)
  {
    std::cout << "drawJetV2Weight_AA: " << path
              << " has no h_jetv2_weight_vs_dpsi2_leg{1,2}/jetv2_amp (rebuild the response to get them) "
              << "-- drawing the weight distribution only." << std::endl;
  }

  dlutility::SetLineAtt(h, kMagenta + 1, 2, 1);
  h->SetFillColorAlpha(kMagenta + 1, 0.25);
  h->SetFillStyle(1001);
  dlutility::SetFont(h, 42, 0.05, 0.045, 0.045, 0.045);

  const TString label = Form("%s, R=0.%d, %s", system_string.c_str(), cone_size, sys_name.c_str());

  TCanvas *c = new TCanvas("cJetV2Weight", "cJetV2Weight", haveDpsi2 ? 1200 : 650, 600);
  if (haveDpsi2) { c->Divide(2, 1); }

  c->cd(1);
  gPad->SetLeftMargin(0.16);
  h->SetTitle(";jet-v_{2} cross-check weight;counts");
  h->Draw("hist");

  PlotUtils::myText(0.22, 0.88, kBlack, label.Data(), 0.035);
  PlotUtils::myText(0.22, 0.82, kBlack, Form("#LTw#GT = %.4f", h->GetMean()), 0.04);
  PlotUtils::myText(0.22, 0.77, kBlack, Form("RMS(w) = %.4f", h->GetRMS()), 0.04);
  PlotUtils::myText(0.22, 0.72, kBlack, Form("Entries = %.0f", h->GetEntries()), 0.04);

  if (haveDpsi2)
  {
    c->cd(2);
    gPad->SetLeftMargin(0.16);

    dlutility::SetLineAtt(p1, kAzure - 6, 2, 1);
    dlutility::SetMarkerAtt(p1, kAzure - 6, 1.1, 21);
    dlutility::SetLineAtt(p2, kRed + 1, 2, 1);
    dlutility::SetMarkerAtt(p2, kRed + 1, 1.1, 22);
    dlutility::SetFont(p1, 42, 0.05, 0.045, 0.045, 0.045);

    p1->SetTitle(";#Delta#psi_{2} (truth);#LTjet-v_{2} weight#GT");
    p1->SetMinimum(1.0 - 3.0 * jetv2_amp);
    p1->SetMaximum(1.0 + 3.0 * jetv2_amp);
    p1->Draw("pe");
    p2->Draw("pe same");

    PlotUtils::myText(0.22, 0.88, kBlack, label.Data(), 0.032);

    // Analytic proof-of-injection curve: <weight | dpsi2> = 1 + jetv2_amp*cos(2*dpsi2)
    // (the other leg's independent dpsi2 averages its cosine term to ~0 --
    // see createResponse_exclusive_v2_AA.cxx's pre-pass comment), evaluated
    // with the same jetv2_amp the response was actually built with. Data
    // points landing on this curve is the proof the reweighting is doing
    // what it claims; a flat set of points would mean it is not.
    TF1 *fcurve = new TF1("f_jetv2_weight_expected", "1.0 + [0]*cos(2.0*x)", 0.0, TMath::Pi() / 2.0);
    fcurve->SetParameter(0, jetv2_amp);
    fcurve->SetLineColor(kGray + 2);
    fcurve->SetLineStyle(2);
    fcurve->SetLineWidth(2);
    fcurve->Draw("same");

    TLegend *leg = new TLegend(0.18, 0.72, 0.55, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(p1, "Leg 1 (leading)", "pe");
    leg->AddEntry(p2, "Leg 2 (subleading)", "pe");
    leg->AddEntry(fcurve, "1 + amp #upoint cos(2#Delta#psi_{2})", "l");
    leg->Draw();
  }

  const TString outpath = Form("%s/unfolding_plots/jetv2_weight_%s_r%02d_%s.pdf",
                                rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str());
  c->Print(outpath);
  std::cout << "Wrote " << outpath << std::endl;
}
