#include <cmath>
#include "dlUtility.h"
#include "read_binning.h"

int color_sim = kRed - 2;
int color_data = kAzure - 6;

void getCentralityReweighting(const int cone_size = 4, const int centrality_bin = 9,  const std::string configfile = "binning_AA.config")
{

  bool ispp = ( centrality_bin < 0 );
  std::string system_string = (ispp?"pp":"AA_cent_" + std::to_string(centrality_bin));
  dlutility::SetyjPadStyle();
  read_binning rb(configfile.c_str());

  
  Int_t zyam_sys = rb.get_zyam_sys();
  Double_t flow_v22_scale = rb.get_flow_sys();
  Double_t flow_v33_scale = rb.get_flow_v33_sys();
  const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 ||
                        std::fabs(flow_v33_scale - 1.0) > 1e-6;
  Int_t inclusive_sys = rb.get_inclusive_sys();
  Int_t flavor_sys = rb.get_flavor_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  Int_t prior_sys = rb.get_prior_sys();
  
  std::string sys_name = "nominal";
  
  if (prior_sys)
    sys_name = "PRIOR";
  
  if (zyam_sys)
    sys_name = "ZYAM";

  if (flow_sys)
    sys_name = rb.get_flow_systematic_name();

  if (inclusive_sys)
    sys_name = "INCLUSIVE";

  if (flavor_sys == 1)
    sys_name = "QQ";

  if (flavor_sys == 2)
    sys_name = "QGGG";

  if (JER_sys < 0)
    sys_name = "negJER";

  if (JER_sys > 0)
    sys_name = "posJER";

  if (JES_sys < 0)
    sys_name = "negJES";

  if (JES_sys > 0)
    sys_name = "posJES";
    

  const TString sim_path = Form("%s/response_matrices/response_matrix_%s_r%02d_PRIMER1_%s.root",
                                rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str());
  TFile *f_sim = new TFile(sim_path,"r");
  if (!f_sim || f_sim->IsZombie())
    {
      std::cerr << "Cannot open primer response file " << sim_path << std::endl;
      return;
    }

  TH1D *h_mbd_sim = (TH1D*) f_sim->Get("h_mbd_vertex");
  TH1D *h_centrality_sim = (TH1D*) f_sim->Get("h_centrality");
  TH1D *h_sumeT_sim = (TH1D*) f_sim->Get("h_sumeT");
  if (!h_mbd_sim || !h_centrality_sim || (!ispp && !h_sumeT_sim))
    {
      std::cerr << "Primer response " << sim_path
                << " is missing one or more required reweighting histograms"
                << " (h_mbd_vertex, h_centrality, h_sumeT)." << std::endl;
      return;
    }
  h_mbd_sim->SetName("h_mbd_sim");
  h_mbd_sim->Rebin(5);
  h_centrality_sim->SetName("h_centrality_sim");
  if (h_sumeT_sim) h_sumeT_sim->SetName("h_sumeT_sim");

  const TString data_path = Form("%s/unfolding_hists/unfolding_hists_preload_%s_r%02d_nominal.root",
                                 rb.get_code_location().c_str(), system_string.c_str(), cone_size);
  TFile *f_data = new TFile(data_path,"r");
  if (!f_data || f_data->IsZombie())
    {
      std::cerr << "Cannot open data preload file " << data_path << std::endl;
      return;
    }

  TH1D *h_mbd_data = (TH1D*) f_data->Get("h_mbd_vertex");
  TH1D *h_centrality_data = (TH1D*) f_data->Get("h_centrality");
  TH1D *h_sumeT_data = (TH1D*) f_data->Get("h_sumeT");
  if (!h_mbd_data || !h_centrality_data || (!ispp && !h_sumeT_data))
    {
      std::cerr << "Data preload " << data_path
                << " is missing one or more required reweighting histograms"
                << " (h_mbd_vertex, h_centrality, h_sumeT)." << std::endl;
      return;
    }
  h_mbd_data->SetName("h_mbd_data");
  h_mbd_data->Rebin(5);
  h_centrality_data->SetName("h_centrality_data");
  if (h_sumeT_data) h_sumeT_data->SetName("h_sumeT_data");

  if (h_mbd_data->Integral() <= 0 || h_mbd_sim->Integral() <= 0 ||
      h_centrality_data->Integral() <= 0 || h_centrality_sim->Integral() <= 0)
    {
      std::cerr << "Cannot derive reweighting from empty data or simulation histograms in "
                << data_path << " and " << sim_path << std::endl;
      return;
    }

  h_mbd_data->Scale(1./h_mbd_data->Integral(),"width");
  h_mbd_sim->Scale(1./h_mbd_sim->Integral(),"width");

  h_centrality_data->Scale(1./h_centrality_data->Integral(),"width");
  h_centrality_sim->Scale(1./h_centrality_sim->Integral(),"width");

  TH1D *h_sumeT_compare = nullptr;
  if (h_sumeT_data && h_sumeT_sim && h_sumeT_data->Integral() > 0 && h_sumeT_sim->Integral() > 0)
    {
      h_sumeT_data->Scale(1./h_sumeT_data->Integral(),"width");
      h_sumeT_sim->Scale(1./h_sumeT_sim->Integral(),"width");
      h_sumeT_compare = (TH1D*) h_sumeT_data->Clone();
      h_sumeT_compare->SetName("h_sumeT_reweight");
      h_sumeT_compare->SetTitle(";#Sigma E_{T}; Data/MC");
      h_sumeT_compare->Divide(h_sumeT_sim);
    }
  else
    {
      std::cerr << "Warning: missing or empty h_sumeT in data/sim reweight inputs; no sumeT reweight file will be written." << std::endl;
    }

  TCanvas *c5 = new TCanvas("c5","c5", 500, 500);
  dlutility::ratioPanelCanvas(c5, 0.4);
  c5->cd(1);
  gPad->SetBottomMargin(0.1);
  dlutility::SetLineAtt(h_mbd_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_mbd_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_mbd_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_mbd_data, color_data, 0.5, 8);
  dlutility::SetFont(h_mbd_data, 42, 0.05);
  h_mbd_data->SetMinimum(0);
  h_mbd_data->SetMaximum(0.02);
  h_mbd_data->SetTitle("; z_{vtx}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dz_{vtx}} ");
  h_mbd_data->Draw("p");
  h_mbd_sim->Draw("same p");

  dlutility::DrawSPHENIX(0.7, 0.8);
  TLegend *leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_mbd_data,"Data");
  leg1->AddEntry(h_mbd_sim,"Sim Reco");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_compare = (TH1D*) h_mbd_data->Clone();
  h_compare->SetName("h_mbd_reweight");
  h_compare->SetTitle("; z_{vtx} [cm]; Data/MC");
  h_compare->Divide(h_mbd_sim);

  dlutility::SetLineAtt(h_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_compare, kBlack, 1, 8);

  dlutility::SetFont(h_compare, 42, 0.07);

  h_compare->SetMinimum(0.5);
  h_compare->SetMaximum(1.5);
  h_compare->Draw();
  TLine *linl = new TLine(-60, 1, 60 ,1);
  linl->SetLineColor(kBlack);
  linl->SetLineWidth(2);
  linl->SetLineStyle(4);
  linl->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print(Form("%s/unfolding_plots/datasim_mbd_%s_r%02d.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_mbd_%s_r%02d.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size));

  c5->cd(1);

  dlutility::SetLineAtt(h_centrality_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_centrality_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_centrality_sim, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_centrality_data, color_data, 0.5, 8);
  dlutility::SetFont(h_centrality_data, 42, 0.05);
  h_centrality_data->SetMinimum(0);
  h_centrality_data->SetMaximum(1.0);
  h_centrality_data->SetTitle("; N_{jet}; #frac{1}{N_{pair}}#frac{dN_{pair}}{d Centrality} ");
  h_centrality_data->Draw("p");
  h_centrality_sim->Draw("same p");

  dlutility::DrawSPHENIX(0.7, 0.8);
  leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_centrality_data,"Data");
  leg1->AddEntry(h_centrality_sim,"Sim Reco");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_centrality_compare = (TH1D*) h_centrality_data->Clone();
  h_centrality_compare->SetName("h_centrality_reweight");
  h_centrality_compare->SetTitle("; Centrality ; Data/MC");
  h_centrality_compare->Divide(h_centrality_sim);

  dlutility::SetLineAtt(h_centrality_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_centrality_compare, kBlack, 1, 8);

  dlutility::SetFont(h_centrality_compare, 42, 0.07);

  h_centrality_compare->SetMinimum(0.5);
  h_centrality_compare->SetMaximum(1.5);
  h_centrality_compare->Draw();
  linl = new TLine(-60, 1, 60 ,1);
  linl->SetLineColor(kBlack);
  linl->SetLineWidth(2);
  linl->SetLineStyle(4);
  linl->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print(Form("%s/unfolding_plots/datasim_centrality_%s_r%02d.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_centrality_%s_r%02d.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size));

  if (h_sumeT_compare)
    {
      c5->cd(1);
      dlutility::SetLineAtt(h_sumeT_sim, color_sim, 1, 1);
      dlutility::SetLineAtt(h_sumeT_data, color_data, 1, 1);
      dlutility::SetMarkerAtt(h_sumeT_sim, color_sim, 0.5, 8);
      dlutility::SetMarkerAtt(h_sumeT_data, color_data, 0.5, 8);
      dlutility::SetFont(h_sumeT_data, 42, 0.05);
      h_sumeT_data->SetMinimum(0);
      h_sumeT_data->SetTitle(";#Sigma E_{T}; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Sigma E_{T}} ");
      h_sumeT_data->Draw("p");
      h_sumeT_sim->Draw("same p");

      dlutility::DrawSPHENIX(0.7, 0.8);
      leg1 = new TLegend(0.22, 0.65, 0.65, 0.85);
      leg1->SetLineWidth(0);
      leg1->SetTextSize(0.04);
      leg1->SetTextFont(42);
      leg1->AddEntry(h_sumeT_data,"Data");
      leg1->AddEntry(h_sumeT_sim,"Sim Reco");
      leg1->Draw("same");

      c5->cd(2);
      gPad->SetTopMargin(0.05);
      dlutility::SetLineAtt(h_sumeT_compare, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_sumeT_compare, kBlack, 1, 8);
      dlutility::SetFont(h_sumeT_compare, 42, 0.07);
      h_sumeT_compare->SetMinimum(0.0);
      h_sumeT_compare->SetMaximum(2.5);
      h_sumeT_compare->Draw();
      linl = new TLine(h_sumeT_compare->GetXaxis()->GetXmin(), 1, h_sumeT_compare->GetXaxis()->GetXmax(), 1);
      linl->SetLineColor(kBlack);
      linl->SetLineWidth(2);
      linl->SetLineStyle(4);
      linl->Draw("same");
      c5->Print(Form("%s/unfolding_plots/datasim_sumeT_%s_r%02d.png", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
      c5->Print(Form("%s/unfolding_plots/datasim_sumeT_%s_r%02d.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size));
    }

  TFile *fout = new TFile(Form("%s/vertex/vertex_reweight_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"recreate");
  h_compare->Write();
  fout->Write();
  fout->Close();

  if (!ispp)
    {
      fout = new TFile(Form("%s/centrality/centrality_reweight_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"recreate");
      h_centrality_compare->Write();
      fout->Write();
      fout->Close();

      if (h_sumeT_compare)
        {
          fout = new TFile(Form("%s/sumeT/sumeT_reweight_%s_r%02d_%s.root", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()),"recreate");
          h_sumeT_compare->Write();
          fout->Write();
          fout->Close();
        }
    }
}
