
#include "dlUtility.h"

const bool NUCLEAR = false;
const bool useJES=true;

struct jet
{
  int id;
  float pt = 0;
  float ptraw = 0;
  float eta = 0;
  float eta_det = 0;
  float phi = 0;
  float t = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;

  void print()
  {
    std::cout << "Jet " << id << std::endl;
    std::cout << "    pt/eta/phi : " << pt << " / " << eta << " / " << phi << std::endl;
  };

};

const bool Debug = false;
const bool Nuclear = false;
const float dRcut = 0.75;
static float cone_size = .3;
const float truth_cut = 5;
const float reco_cut = 5;
const float etacut = 0.7;
const float isocut = 0.8;;
const bool use_jes = true;

float etacut_cone = 1.1;

const float dphicut = 3*TMath::Pi()/4.;

const float vertex_cut = 60;

float getDPHI(float phi1, float phi2);

std::vector<std::pair<struct jet, struct jet>>  match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets);

double getDR(struct jet j1, struct jet j2)
{ 
  double dphi = fabs(j1.phi - j2.phi);
  if (dphi > TMath::Pi())
    {
      dphi = 2*TMath::Pi() - dphi;
    }

  double dR = sqrt(TMath::Power(j1.eta - j2.eta, 2) + TMath::Power(dphi, 2));
  return dR;
}
bool check_dijet_reco(std::vector<struct jet> myrecojets);
bool check_dijet_truth(std::vector<struct jet> mytruthjets);
bool passes_time_cut(double calib_lead_time, double calib_delta_time)
{
  double x = calib_lead_time;
  double y = x - calib_delta_time;
  return (fabs(x) < 6 && fabs(y) < 3);
}

void makeJER(const int conesize = 4, const int closure = 0, const int use_herwig = 0, const int pileup = 0, const int ca_sys = 0)
{

  cone_size = conesize*0.1;

  std::string gen[6] = {"PYTHIA","HERWIG","PILEUP","MIX", "CA1","CA2"};

  int index = 0;
  if (use_herwig) index = 1;
  int runnumber_low = 0;
  int runnumber_high = 100000;
  if (ca_sys)
    {      
      if (ca_sys == 1)
	{
	  runnumber_high = 52212;
	  index = 4;
	}
      else
	{
	  runnumber_low = 52212;
	  index = 5;
	}
    }

  int sim_version = (use_herwig? 10 : 11);
  if (pileup == 1)
    {
      index = 2;
      sim_version = 12;
    }
  if (pileup == 2)
    {
      index = 3;
    }
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  float boundary_r4[6];
  boundary_r4[0] = 14;
  boundary_r4[1] = 22;
  boundary_r4[2] = 33;
  boundary_r4[3] = 42;
  boundary_r4[4] = 52;
  boundary_r4[5] = 100;

  TF1 *gjer = nullptr;

  TF1 *fsmear = new TF1("fsmear", "gaus", -1, 1);
  fsmear->SetParameters(1, 0, 0.13);
  TH1D *hjer = nullptr;
  if (closure)
    {
      TFile *finjer = new TFile(Form("r%02d/jer/jer_fits_r%02d_1JES_0_closure_%s.root", conesize, conesize, gen[index].c_str()), "r");
      gjer = (TF1*) finjer->Get("jernom");
      hjer = (TH1D*) finjer->Get("h_nominal_smear");
      hjer->SetName("input");
      fsmear->SetParameter(2, gjer->GetParameter(0));
    }

  etacut_cone = 1.1 - conesize*0.1;; 

  float low_pt = 20;
  float high_pt = 50;

  
  TFile* fin_data[2];
  
  fin_data[0] = new TFile(Form("../trees/TREE_DIJET_SKIM_r%02d_v8_3_ana509_2024p022_v001_gl10-minbias.root", conesize), "r");
  fin_data[1] = new TFile(Form("../trees/TREE_DIJET_SKIM_r%02d_v8_7_ana533_2025p007_v001_gl10-alltime.root", conesize), "r");

  
  int nevents_data[2];
  std::cout << "Getting Data Trees" << std::endl;
  TTree *ttree_data[2];// = (TTree*) fin_data->Get("ttree");
  ULong64_t gl1_scaled_data[2];  
  std::vector<float> *data_jet_pt[2] = {0};
  std::vector<float> *data_jet_t[2] = {0};
  std::vector<float> *data_jet_emcal[2] = {0};
  std::vector<float> *data_jet_e[2] = {0};
  std::vector<float> *data_jet_eta[2] = {0};
  std::vector<float> *data_jet_eta_det[2] = {0};
  std::vector<float> *data_jet_phi[2] = {0};
  int mbd_hit_data[2];
  double calib_lead_time[2];
  double calib_delta_time[2];
  float mbd_vertex_z_data[2];
  int runnumber[2];
  for (int i = 1; i < 2; i++)
    {
      ttree_data[i] = (TTree*) fin_data[i]->Get("ttree");
      ttree_data[i]->SetBranchAddress(Form("jet_pt%s_%d", (use_jes ? "_calib" : ""), conesize), &data_jet_pt[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_t_%d", conesize), &data_jet_t[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_emcal_%d", conesize), &data_jet_emcal[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_e_%d", conesize), &data_jet_e[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_eta_%d", conesize), &data_jet_eta[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_eta_det_%d", conesize), &data_jet_eta_det[i]);
      ttree_data[i]->SetBranchAddress(Form("jet_phi_%d", conesize), &data_jet_phi[i]);
      ttree_data[i]->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z_data[i]);
      ttree_data[i]->SetBranchAddress("mbd_hit", &mbd_hit_data[i]);
      ttree_data[i]->SetBranchAddress("calib_lead_time", &calib_lead_time[i]);
      ttree_data[i]->SetBranchAddress("calib_delta_time", &calib_delta_time[i]);
      ttree_data[i]->SetBranchAddress("mbd_hit", &mbd_hit_data[i]);
      ttree_data[i]->SetBranchAddress("gl1_scaled", &gl1_scaled_data[i]);
      ttree_data[i]->SetBranchAddress("runnumber", &runnumber[i]);
      
      nevents_data[i] = ttree_data[i]->GetEntries();
    }

  std::string j_file[10];
  j_file[0] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v" + std::to_string(sim_version) + "_12_ana509_MDC2-00000028.root";
  j_file[1] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v" + std::to_string(sim_version) + "_20_ana509_MDC2-00000028.root";
  j_file[2] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v" + std::to_string(sim_version) + "_30_ana509_MDC2-00000028.root";
  j_file[3] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v" + std::to_string(sim_version) + "_40_ana509_MDC2-00000028.root";
  j_file[4] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v" + std::to_string(sim_version) + "_50_ana509_MDC2-00000028.root";

  bool use_sample[10] = {1, 1, 1, 1, 1, 0, 0, 0 ,0, 0};

  if (pileup == 2)
    {
      int v = 12;
      j_file[5] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v12_12_ana509_MDC2-00000028.root";
      j_file[6] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v12_20_ana509_MDC2-00000028.root";
      j_file[7] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v12_30_ana509_MDC2-00000028.root";
      j_file[8] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v12_40_ana509_MDC2-00000028.root";
      j_file[9] = "../trees//TREE_JET_SKIM_r0" + std::to_string(conesize) + "_v12_50_ana509_MDC2-00000028.root";
      for (int i = 0; i < 5; i++)
	{
	  use_sample[5+i] = 1;
	}
    }  
  TFile* finsim[10];

  float nevents[10];

  TTree *ttree[10];
  ULong64_t gl1_scaled[10];

  std::vector<float> *truth_jet_pt_ref[10] = {0};
  std::vector<float> *truth_jet_pt[10] = {0};
  std::vector<int> *truth_jet_flavor[10] = {0};
  std::vector<float> *truth_jet_eta[10] = {0};
  std::vector<float> *truth_jet_phi[10] = {0};
  
  std::vector<float> *reco_jet_pt[10] = {0};
  std::vector<float> *reco_jet_emcal[10] = {0};
  std::vector<float> *reco_jet_e[10] = {0};
  std::vector<float> *reco_jet_eta[10] = {0};
  std::vector<float> *reco_jet_eta_det[10] = {0};
  std::vector<float> *reco_jet_phi[10] = {0};

  float truth_vertex_z[10];
  float mbd_vertex_z[10];

  std::cout << "Getting Sim Trees" << std::endl;
  for (int j = 0; j < 10; j++)
    {
      if (!use_sample[j]) continue;
      finsim[j] = new TFile(j_file[j].c_str(),"r");

      ttree[j] = (TTree*) finsim[j]->Get("ttree");
      if (!ttree[j])
	{
	  std::cout << "NO TREE " << j << std::endl;
	  return;
	}

      if (conesize != 4)
	{
	  ttree[j]->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_ref[j]);
	}
      ttree[j]->SetBranchAddress(Form("truth_jet_pt_%d", conesize), &truth_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_eta_%d", conesize), &truth_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("truth_jet_phi_%d", conesize), &truth_jet_phi[j]);
  
      ttree[j]->SetBranchAddress(Form("jet_pt%s_%d", (use_jes ? "_calib" : ""), conesize), &reco_jet_pt[j]);
      //      ttree[j]->SetBranchAddress(Form("jet_pt_%s%d", conesize), &reco_jet_pt[j]);
      ttree[j]->SetBranchAddress(Form("jet_emcal_%d", conesize), &reco_jet_emcal[j]);
      ttree[j]->SetBranchAddress(Form("jet_e_%d", conesize), &reco_jet_e[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_%d", conesize), &reco_jet_eta[j]);
      ttree[j]->SetBranchAddress(Form("jet_eta_det_%d", conesize), &reco_jet_eta_det[j]);
      ttree[j]->SetBranchAddress(Form("jet_phi_%d", conesize), &reco_jet_phi[j]);
      ttree[j]->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z[j]);
      ttree[j]->SetBranchAddress("truth_vertex_z", &truth_vertex_z[j]);
      ttree[j]->SetBranchAddress("gl1_scaled", &gl1_scaled[j]);
      nevents[j] = 10000000.;//ttree[j]->GetEntries();
    }
  
  float cs_12 = (1.49e6);
  float cs_20 = (6.26e4);
  float cs_30 = (2.53e3);
  float cs_40 = (1.355e2);
  float cs_50 = (7.311);


  float csh_10 = 1.5703e-9;
  float csh_30 = 1.473e-12;
  float scale_factor[5];
  scale_factor[0] = cs_12/cs_50;
  scale_factor[1] = cs_20/cs_50;
  scale_factor[2] = cs_30/cs_50;
  scale_factor[3] = cs_40/cs_50; 
  scale_factor[4] = 1;
  float cs_h[5] = {6.7108e+05, 5.2613e4, 2.0694e3, 1.0510e2, 5.2089};
  //float nevents_h[5] = {1063000, 10000000., 10000000., 9671056, 1239000}; 
  float nevents_h[5] = {8057000, 10000000., 10000000., 9671056, 3827000}; 
  if (use_herwig)
    {
      for (int i = 0 ; i < 5; i++)
	{
	  scale_factor[i] = cs_h[i]*nevents_h[4]/(cs_h[4]*nevents_h[i]);
	}
      scale_factor[0]*=2;      
    }

  /* if (use_herwig) */
  /*   { */
  /*     boundary_r4[0] = 14; */
  /*     boundary_r4[1] = 35; */
  /*     boundary_r4[3] = 35; */
  /*     boundary_r4[4] = 70; */
  /*   } */

  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;
  std::vector<struct jet> mytruthjets2;
  std::vector<struct jet> myrecojets2;

  std::vector<std::pair<struct jet, struct jet>> matched_dijets;

  TRandom *rng = new TRandom();  
  std::cout << " Making Histograms"<<std::endl;

  float first_trigger = 23;
  const int nbins_pt = 9;
  double pt_bin_edges[nbins_pt + 1] = {16, 18, 20, 23, 26, 29, 33, 38, 45, 60};

  //const int nbins_pt = 4;
  //double pt_bin_edges[nbins_pt + 1] = {20, 25, 30, 35, 40};

  float dpt = (high_pt - low_pt);
  //float reco_cut = 7.;
  float third_jet_cut_single = 7;
  float dphicut = 6.*TMath::Pi()/8.;

  float third_jet_cuts[8] = {10, 9, 8, 7, 6, 5, 4, 3}; 
  
  TH1D *h_lead_sample[6];
  TH1D *h_lead_combined = new TH1D("h_lead_combined","; Leading Truth Jet #it{p}_{T} [GeV]; lumi scale * counts",100, 0, 100);

  for (int i = 0; i < 6; i++)
    {
      h_lead_sample[i] = new TH1D(Form("h_lead_sample_%d", i), "; Leading Truth Jet #it{p}_{T} [GeV]; lumi scale * counts",100, 0, 100);
    }
  TH1D *h_sigp_sim = new TH1D("h_sigp_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);
  TH1D *h_sigp_data = new TH1D("h_sigp_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);

  TH1D *h_sigv_sim = new TH1D("h_sigv_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);
  TH1D *h_sigv_data = new TH1D("h_sigv_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);

  TH1D *h_jer_sim = new TH1D("h_jer_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
  TH1D *h_jer_data = new TH1D("h_jer_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

  TProfile *hp_avgpt_sim = new TProfile("hp_avgpt_sim",";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);
  TProfile *hp_avgpt_data = new TProfile("hp_avgpt_data",";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);

  TH2D *h_trigger_ptem_pta_sim = new TH2D("h_trigger_ptem_pta_sim",";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
  TH2D *h_trigger_ptem_pta_data = new TH2D("h_trigger_ptem_pta_data",";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);

  TH2D *h_trigger_ptv_pta_sim = new TH2D("h_trigger_ptv_pta_sim",";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);
  TH2D *h_trigger_ptv_pta_data = new TH2D("h_trigger_ptv_pta_data",";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);
  
  TH2D *h_trigger_ptp_pta_sim = new TH2D("h_trigger_ptp_pta_sim",";pta;ptp", nbins_pt, pt_bin_edges, 80, -40, 40);
  TH2D *h_trigger_ptp_pta_data = new TH2D("h_trigger_ptp_pta_data",";pta;ptp", nbins_pt, pt_bin_edges, 80, -40, 40);
  
  TProfile *hp_trigger_ptem_pta_sim = new TProfile("hp_trigger_ptem_pta_sim",";pta;em", nbins_pt, pt_bin_edges);//,"s");
  TProfile *hp_trigger_ptem_pta_data = new TProfile("hp_trigger_ptem_pta_data",";pta;em", nbins_pt, pt_bin_edges);//,"s");

  TProfile *hp_trigger_ptv_pta_sim = new TProfile("hp_trigger_ptv_pta_sim",";pta;ptv", nbins_pt, pt_bin_edges);//,"s");
  TProfile *hp_trigger_ptv_pta_data = new TProfile("hp_trigger_ptv_pta_data",";pta;ptv", nbins_pt, pt_bin_edges);//,"s");

  TProfile *hp_trigger_ptp_pta_sim = new TProfile("hp_trigger_ptp_pta_sim",";pta;ptp", nbins_pt, pt_bin_edges);//,"s");
  TProfile *hp_trigger_ptp_pta_data = new TProfile("hp_trigger_ptp_pta_data",";pta;ptp", nbins_pt, pt_bin_edges);//,"s");

  TProfile *hp_cosdphi_pta_sim = new TProfile("hp_cosdphi_pta_sim",";pta;ptp", nbins_pt, pt_bin_edges,"s");
  TProfile *hp_cosdphi_pta_data = new TProfile("hp_cosdphi_pta_data",";pta;ptp", nbins_pt, pt_bin_edges,"s");


  // 3jet cut
  TH1D *h_sigp_sim_3jet[8]; // new TH1D("h_sigp_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T}, #psi})", 10, low_pt, high_pt);
  TH1D *h_sigp_data_3jet[8]; // new TH1D("h_sigp_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T}, #psi})", 10, low_pt, high_pt);

  TH1D *h_sigv_sim_3jet[8]; // new TH1D("h_sigv_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T}, #lambda})", 10, low_pt, high_pt);
  TH1D *h_sigv_data_3jet[8]; // new TH1D("h_sigv_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T}, #lambda})", 10, low_pt, high_pt);

  TH1D *h_jer_sim_3jet[8]; // new TH1D("h_jer_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_data_3jet[8]; // new TH1D("h_jer_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);

  TH1D *h_jer_truth_im_3jet[8]; // new TH1D("h_jer_truth",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_sim_im_3jet[8]; // new TH1D("h_jer_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_data_im_3jet[8]; // new TH1D("h_jer_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);

  TH1D *h_jer_truth_im_softcorr_3jet[8]; // new TH1D("h_jer_truth",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_sim_im_softcorr_3jet[8]; // new TH1D("h_jer_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_data_im_softcorr_3jet[8]; // new TH1D("h_jer_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);

  TH1D *h_jer_sim_im_softcorr_part_3jet[8]; // new TH1D("h_jer_sim",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);
  TH1D *h_jer_data_im_softcorr_part_3jet[8]; // new TH1D("h_jer_data",";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", 10, low_pt, high_pt);

  TProfile *hp_avgpt_sim_3jet[8]; // new TH1D("h_avgpt_sim",";#LT#it{p}_{T}#GT [GeV];Counts", 10, low_pt, high_pt);
  TProfile *hp_avgpt_data_3jet[8]; // new TH1D("h_avgpt_data",";#LT#it{p}_{T}#GT [GeV];Counts", 10, low_pt, high_pt);
  TProfile *hp_avgpt_truth_3jet[8]; // new TH1D("h_avgpt_truth",";#LT#it{p}_{T}#GT [GeV];Counts", 10, low_pt, high_pt);

  TH2D *h_trigger_ptem_pta_sim_3jet[8]; // new TH2D("h_trigger_ptem_pta_sim",";pta;em", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_trigger_ptem_pta_data_3jet[8]; // new TH2D("h_trigger_ptem_pta_data",";pta;em", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_lead_ptem_pta_sim_3jet[8]; // new TH2D("h_lead_ptem_pta_sim",";pta;em", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_lead_ptem_pta_data_3jet[8]; // new TH2D("h_lead_ptem_pta_data",";pta;em", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_sublead_ptem_pta_sim_3jet[8]; // new TH2D("h_sublead_ptem_pta_sim",";pta;em", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_sublead_ptem_pta_data_3jet[8]; // new TH2D("h_sublead_ptem_pta_data",";pta;em", 10, low_pt, high_pt, 40, -20, 20);

  TH2D *h_trigger_ptv_pta_sim_3jet[8]; // new TH2D("h_trigger_ptv_pta_sim",";pta;ptv", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_trigger_ptv_pta_data_3jet[8]; // new TH2D("h_trigger_ptv_pta_data",";pta;ptv", 10, low_pt, high_pt, 40, -20, 20);

  TH2D *h_trigger_ptp_pta_sim_3jet[8]; // new TH2D("h_trigger_ptp_pta_sim",";pta;ptp", 10, low_pt, high_pt, 40, -20, 20);
  TH2D *h_trigger_ptp_pta_data_3jet[8]; // new TH2D("h_trigger_ptp_pta_data",";pta;ptp", 10, low_pt, high_pt, 40, -20, 20);

  TProfile *hp_trigger_ptem_pta_sim_3jet[8]; // new TProfile("hp_trigger_ptem_pta_sim",";pta;em", 10, low_pt, high_pt);//,"s");
  TProfile *hp_trigger_ptem_pta_data_3jet[8]; // new TProfile("hp_trigger_ptem_pta_data",";pta;em", 10, low_pt, high_pt);//,"s");
  TProfile *hp_lead_ptem_pta_sim_3jet[8]; // new TProfile("hp_lead_ptem_pta_sim",";pta;em", 10, low_pt, high_pt);//,"s");
  TProfile *hp_lead_ptem_pta_data_3jet[8]; // new TProfile("hp_lead_ptem_pta_data",";pta;em", 10, low_pt, high_pt);//,"s");
  TProfile *hp_sublead_ptem_pta_sim_3jet[8]; // new TProfile("hp_sublead_ptem_pta_sim",";pta;em", 10, low_pt, high_pt);//,"s");
  TProfile *hp_sublead_ptem_pta_data_3jet[8]; // new TProfile("hp_sublead_ptem_pta_data",";pta;em", 10, low_pt, high_pt);//,"s");
  
  TProfile *hp_trigger_ptv_pta_sim_3jet[8]; // new TProfile("hp_trigger_ptv_pta_sim",";pta;ptv", 10, low_pt, high_pt);//,"s");
  TProfile *hp_trigger_ptv_pta_data_3jet[8]; // new TProfile("hp_trigger_ptv_pta_data",";pta;ptv", 10, low_pt, high_pt);//,"s");

  TProfile *hp_trigger_ptp_pta_sim_3jet[8]; // new TProfile("hp_trigger_ptp_pta_sim",";pta;ptp", 10, low_pt, high_pt);//,"s");
  TProfile *hp_trigger_ptp_pta_data_3jet[8]; // new TProfile("hp_trigger_ptp_pta_data",";pta;ptp", 10, low_pt, high_pt);//,"s");

  TProfile *hp_cosdphi_pta_sim_3jet[8]; // new TProfile("hp_cosdphi_pta_sim",";pta;ptp", 10, low_pt, high_pt);//,"s");
  TProfile *hp_cosdphi_pta_data_3jet[8]; // new TProfile("hp_cosdphi_pta_data",";pta;ptp", 10, low_pt, high_pt);//,"s");

  for (int i = 0; i < 8; i++)
    {
      h_sigp_sim_3jet[i] = new TH1D(Form("h_sigp_sim_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);
      h_sigp_data_3jet[i] = new TH1D(Form("h_sigp_data_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);

      h_sigv_sim_3jet[i] = new TH1D(Form("h_sigv_sim_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);
      h_sigv_data_3jet[i] = new TH1D(Form("h_sigv_data_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);
      h_jer_sim_3jet[i] = new TH1D(Form("h_jer_sim_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_data_3jet[i] = new TH1D(Form("h_jer_data_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

      h_jer_sim_im_3jet[i] = new TH1D(Form("h_jer_sim_im_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_data_im_3jet[i] = new TH1D(Form("h_jer_data_im_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_truth_im_3jet[i] = new TH1D(Form("h_jer_truth_im_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

      h_jer_sim_im_softcorr_3jet[i] = new TH1D(Form("h_jer_sim_im_softcorr_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_data_im_softcorr_3jet[i] = new TH1D(Form("h_jer_data_im_softcorr_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_truth_im_softcorr_3jet[i] = new TH1D(Form("h_jer_truth_im_softcorr_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

      h_jer_sim_im_softcorr_part_3jet[i] = new TH1D(Form("h_jer_sim_im_softcorr_part_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_jer_data_im_softcorr_part_3jet[i] = new TH1D(Form("h_jer_data_im_softcorr_part_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

      hp_avgpt_sim_3jet[i] = new TProfile(Form("hp_avgpt_sim_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);
      hp_avgpt_data_3jet[i] = new TProfile(Form("hp_avgpt_data_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);
      hp_avgpt_truth_3jet[i] = new TProfile(Form("hp_avgpt_truth_3jet%d", i),";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);

      h_trigger_ptem_pta_sim_3jet[i] = new TH2D(Form("h_trigger_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
      h_trigger_ptem_pta_data_3jet[i] = new TH2D(Form("h_trigger_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
      h_lead_ptem_pta_sim_3jet[i] = new TH2D(Form("h_lead_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
      h_lead_ptem_pta_data_3jet[i] = new TH2D(Form("h_lead_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
      h_sublead_ptem_pta_sim_3jet[i] = new TH2D(Form("h_sublead_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);
      h_sublead_ptem_pta_data_3jet[i] = new TH2D(Form("h_sublead_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges, 14, -0.2, 1.2);

      h_trigger_ptv_pta_sim_3jet[i] = new TH2D(Form("h_trigger_ptv_pta_sim_3jet%d", i),";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);
      h_trigger_ptv_pta_data_3jet[i] = new TH2D(Form("h_trigger_ptv_pta_data_3jet%d", i),";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);
      h_trigger_ptp_pta_sim_3jet[i] = new TH2D(Form("h_trigger_ptp_pta_sim_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges,  80, -40, 40);
      h_trigger_ptp_pta_data_3jet[i] = new TH2D(Form("h_trigger_ptp_pta_data_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges, 80, -40, 40);

      hp_trigger_ptem_pta_sim_3jet[i] = new TProfile(Form("hp_trigger_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");
      hp_trigger_ptem_pta_data_3jet[i] = new TProfile(Form("hp_trigger_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");
      hp_lead_ptem_pta_sim_3jet[i] = new TProfile(Form("hp_lead_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");
      hp_lead_ptem_pta_data_3jet[i] = new TProfile(Form("hp_lead_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");
      hp_sublead_ptem_pta_sim_3jet[i] = new TProfile(Form("hp_sublead_ptem_pta_sim_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");
      hp_sublead_ptem_pta_data_3jet[i] = new TProfile(Form("hp_sublead_ptem_pta_data_3jet%d", i),";pta;em", nbins_pt, pt_bin_edges);//,"s");

      hp_trigger_ptv_pta_sim_3jet[i] = new TProfile(Form("hp_trigger_ptv_pta_sim_3jet%d", i),";pta;ptv", nbins_pt, pt_bin_edges);//,"s");
      hp_trigger_ptv_pta_data_3jet[i] = new TProfile(Form("hp_trigger_ptv_pta_data_3jet%d", i),";pta;ptv", nbins_pt, pt_bin_edges);//,"s");
      hp_trigger_ptp_pta_sim_3jet[i] = new TProfile(Form("hp_trigger_ptp_pta_sim_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges);//,"s");
      hp_trigger_ptp_pta_data_3jet[i] = new TProfile(Form("hp_trigger_ptp_pta_data_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges);//,"s");
      hp_cosdphi_pta_sim_3jet[i] = new TProfile(Form("hp_cosdphi_pta_sim_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges,"s");
      hp_cosdphi_pta_data_3jet[i] = new TProfile(Form("hp_cosdphi_pta_data_3jet%d", i),";pta;ptp", nbins_pt, pt_bin_edges,"s");
    }

  // dijet balance

  TH2D *h2_aj_avgpt_sim_3jet[8];
  TProfile *hp_aj_avgpt_sim_3jet[8];

  TH2D *h2_aj_avgpt_truth_3jet[8];
  TProfile *hp_aj_avgpt_truth_3jet[8];

  TH2D *h2_aj_avgpt_data_3jet[8];
  TProfile *hp_aj_avgpt_data_3jet[8];
  
  for (int i = 0; i < 8; i++)
    {
      h2_aj_avgpt_sim_3jet[i] = new TH2D(Form("h2_aj_avgpt_sim_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, 100, -1, 1);
      hp_aj_avgpt_sim_3jet[i] = new TProfile(Form("hp_aj_avgpt_sim_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, "s");

      h2_aj_avgpt_truth_3jet[i] = new TH2D(Form("h2_aj_avgpt_truth_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, 100, -1, 1);
      hp_aj_avgpt_truth_3jet[i] = new TProfile(Form("hp_aj_avgpt_truth_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, "s");

      h2_aj_avgpt_data_3jet[i] = new TH2D(Form("h2_aj_avgpt_data_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, 100, -1, 1);
      hp_aj_avgpt_data_3jet[i] = new TProfile(Form("hp_aj_avgpt_data_3jet%d", i), ";avg pt; A_{J}",nbins_pt, pt_bin_edges, "s");
    }
  
  // EM fraction

  const int nbins_em = 2;
  float em_bounds[3] = {-0.2, 0.5, 1.2};

  TH1D *h_em_sigp_sim[5];
  TH1D *h_em_sigp_data[5];

  TH1D *h_em_sigv_sim[5];
  TH1D *h_em_sigv_data[5];

  TH1D *h_em_jer_sim[5];
  TH1D *h_em_jer_data[5];

  TProfile *hp_em_avgpt_sim[5];
  TProfile *hp_em_avgpt_data[5];

  TH2D *h_em_trigger_ptv_pta_sim[5];
  TH2D *h_em_trigger_ptv_pta_data[5];

  TH2D *h_em_trigger_ptp_pta_sim[5];
  TH2D *h_em_trigger_ptp_pta_data[5];

  TProfile *hp_em_trigger_ptv_pta_sim[5];
  TProfile *hp_em_trigger_ptv_pta_data[5];

  TProfile *hp_em_trigger_ptp_pta_sim[5];
  TProfile *hp_em_trigger_ptp_pta_data[5];

  TProfile *hp_em_cosdphi_pta_sim[5];
  TProfile *hp_em_cosdphi_pta_data[5];
  for (int i = 0; i < 5; i++)
    {
      h_em_sigp_sim[i] = new TH1D(Form("h_em_sigp_sim_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);
      h_em_sigp_data[i] = new TH1D(Form("h_em_sigp_data_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #psi})", nbins_pt, pt_bin_edges);

      h_em_sigv_sim[i] = new TH1D(Form("h_em_sigv_sim_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);
      h_em_sigv_data[i] = new TH1D(Form("h_em_sigv_data_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T, #lambda})", nbins_pt, pt_bin_edges);

      h_em_jer_sim[i] = new TH1D(Form("h_em_jer_sim_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_em_jer_data[i] = new TH1D(Form("h_em_jer_data_%d", i),";#LT#it{p}_{T}#GT [GeV];#sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);

      hp_em_avgpt_sim[i] = new TProfile(Form("hp_em_avgpt_sim_%d", i),";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);
      hp_em_avgpt_data[i] = new TProfile(Form("hp_em_avgpt_data_%d", i),";#LT#it{p}_{T}#GT [GeV];Counts", nbins_pt, pt_bin_edges);

      h_em_trigger_ptv_pta_sim[i] = new TH2D(Form("h_em_trigger_ptv_pta_sim_%d", i),";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);
      h_em_trigger_ptv_pta_data[i] = new TH2D(Form("h_em_trigger_ptv_pta_data_%d", i),";pta;ptv", nbins_pt, pt_bin_edges, 80, -40, 40);

      h_em_trigger_ptp_pta_sim[i] = new TH2D(Form("h_em_trigger_ptp_pta_sim_%d", i),";pta;ptp", nbins_pt, pt_bin_edges, 80, -40, 40);
      h_em_trigger_ptp_pta_data[i] = new TH2D(Form("h_em_trigger_ptp_pta_data_%d", i),";pta;ptp", nbins_pt, pt_bin_edges, 80, -40, 40);

  
      hp_em_trigger_ptv_pta_sim[i] = new TProfile(Form("hp_em_trigger_ptv_pta_sim_%d", i),";pta;ptv", nbins_pt, pt_bin_edges);//,"s");
      hp_em_trigger_ptv_pta_data[i] = new TProfile(Form("hp_em_trigger_ptv_pta_data_%d", i),";pta;ptv", nbins_pt, pt_bin_edges);//,"s");

      hp_em_trigger_ptp_pta_sim[i] = new TProfile(Form("hp_em_trigger_ptp_pta_sim_%d", i),";pta;ptp", nbins_pt, pt_bin_edges);//,"s");
      hp_em_trigger_ptp_pta_data[i] = new TProfile(Form("hp_em_trigger_ptp_pta_data_%d", i),";pta;ptp", nbins_pt, pt_bin_edges);//,"s");

      hp_em_cosdphi_pta_sim[i] = new TProfile(Form("hp_em_cosdphi_pta_sim_%d", i),";pta;ptp", nbins_pt, pt_bin_edges,"s");
      hp_em_cosdphi_pta_data[i] = new TProfile(Form("hp_em_cosdphi_pta_data_%d", i),";pta;ptp", nbins_pt, pt_bin_edges,"s");

    }
  for (int isample = 1; isample < 2; isample++)
    {
      std::cout << " Starting Data : " << nevents_data[isample] << std::endl;
      for (int i = 0; i < nevents_data[isample]; i++)
	{
	  ttree_data[isample]->GetEntry(i);
	  std::cout << "Event : " << i << " \r" << std::flush;
	  
	  if (runnumber[isample] >= runnumber_high || runnumber[isample] < runnumber_low) continue;

	  if (fabs(mbd_vertex_z_data[isample]) > vertex_cut) continue;


	  // make sure it passes Jet 10 GeV Trigger
	  bool trig = (((gl1_scaled_data[isample] >> 17) & 0x1) == 0x1 );
	  // min bias for lower pt bins < 20
	  bool mbtrig = (((gl1_scaled_data[isample] >> 12) & 0x1) == 0x1 );// || (((gl1_scaled_data[isample] >> 12) & 0x1) == 0x1 ) || (((gl1_scaled_data[isample] >> 14) & 0x1) == 0x1 );
	  if (!trig && !mbtrig) continue;
	  
	  if (NUCLEAR)  std::cout <<  "Looking through jets" << std::endl;
	  int nrecojets = data_jet_pt[isample]->size();
	  if (nrecojets == 0) continue;
	  myrecojets.clear();

	  for (int j = 0; j < nrecojets;j++)
	    {
	      if (data_jet_pt[isample]->at(j) < reco_cut) continue;
	      if (data_jet_e[isample]->at(j) < 0) continue;

	      struct jet tempjet;

	      tempjet.pt = data_jet_pt[isample]->at(j);

	      
	      if (tempjet.pt == 0) continue;

	      tempjet.t = data_jet_t[isample]->at(j);
	      tempjet.eta = data_jet_eta[isample]->at(j);
	      tempjet.eta_det = data_jet_eta_det[isample]->at(j);
	      tempjet.phi = data_jet_phi[isample]->at(j);
	      tempjet.emcal = data_jet_emcal[isample]->at(j);
	      tempjet.id = j;
	      myrecojets.push_back(tempjet);	  
	    }

	  if (myrecojets.size() < 2) continue;
	  if (NUCLEAR)  std::cout <<  "Found Dijet" << std::endl;

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  if (!check_dijet_reco(myrecojets))
	    {
	      continue;
	    }
	  if (!passes_time_cut(calib_lead_time[isample], calib_delta_time[isample])) continue;
	  
	  if (NUCLEAR)  std::cout <<  "Passes Dijet" << std::endl;
	  float pt1_data = myrecojets.at(0).pt;
	  float em1_data = myrecojets.at(0).emcal;
	  float pt2_data = myrecojets.at(1).pt;
	  float em2_data = myrecojets.at(1).emcal;

	  float trigger_jet = myrecojets.at(0).pt;
	  float probe_jet = myrecojets.at(1).pt;
	  float pt1_calib = myrecojets.at(0).pt;
	  float pt2_calib = myrecojets.at(1).pt;

	  float third_jet = 0;
	  float dphi_data = fabs(getDPHI(myrecojets.at(0).phi, myrecojets.at(1).phi));
      
	  if (myrecojets.size() >= 3)
	    {
	      third_jet = myrecojets.at(2).pt;
	    }
      
	  if (rng->Uniform() > 0.5)
	    {
	      float tempy = trigger_jet;
	      trigger_jet = probe_jet;
	      probe_jet = tempy;
	    }

	  bool data_good = check_dijet_reco(myrecojets);
      
	  float ptv = trigger_jet*cos(dphi_data/2) + probe_jet*cos(dphi_data/2);
	  float ptp = trigger_jet*sin(dphi_data/2) - probe_jet*sin(dphi_data/2);

	  float ptem = (em1_data * pt1_calib + em2_data * pt2_calib )/(pt1_calib + pt2_calib);

	  float pta = (trigger_jet + probe_jet)/2.;
	  float aj = (trigger_jet - probe_jet) / (trigger_jet + probe_jet);
	  float xj = (trigger_jet - probe_jet) / (trigger_jet + probe_jet);

	  bool mb_fill = true;
	  bool trig_fill = true;

	  if (mbtrig && (pta < 12 || pta >= first_trigger)) mb_fill = false;//continue;

	  if (trig && (pta < first_trigger)) trig_fill = false;

	  if (!trig_fill && !mb_fill) continue;
	  if (!data_good) continue;
	  
	  if (rng->Uniform() > 0.5)
	    {
	      ptv = -1*ptv;
	    }

	  if (trigger_jet >= reco_cut && probe_jet >= reco_cut && dphi_data >= dphicut)
	    {

	      for (int i3 = 0 ; i3 < 8; i3++)
		{
		  if (third_jet < third_jet_cuts[i3])
		    {
		      h2_aj_avgpt_data_3jet[i3]->Fill(pta, aj);
		      hp_aj_avgpt_data_3jet[i3]->Fill(pta, aj);
		      hp_avgpt_data_3jet[i3]->Fill(pta, pta);

		      h_trigger_ptv_pta_data_3jet[i3]->Fill(pta, ptv);
		      hp_trigger_ptv_pta_data_3jet[i3]->Fill(pta, ptv);

		      h_trigger_ptem_pta_data_3jet[i3]->Fill(pta, ptem);
		      h_lead_ptem_pta_data_3jet[i3]->Fill(pta, em1_data );
		      h_sublead_ptem_pta_data_3jet[i3]->Fill(pta, em2_data );
		      hp_trigger_ptem_pta_data_3jet[i3]->Fill(pta, ptem);
		      hp_lead_ptem_pta_data_3jet[i3]->Fill(pta, em1_data );
		      hp_sublead_ptem_pta_data_3jet[i3]->Fill(pta, em2_data );

		      h_trigger_ptp_pta_data_3jet[i3]->Fill(pta, ptp);
		      hp_trigger_ptp_pta_data_3jet[i3]->Fill(pta, ptp);      
		      hp_cosdphi_pta_data_3jet[i3]->Fill(pta, fabs(cos(dphi_data)));

		    }
		}
	    }



	  hp_avgpt_data->Fill(pta, pta);
	  h_trigger_ptv_pta_data->Fill(pta, ptv);
	  hp_trigger_ptv_pta_data->Fill(pta, ptv);
	  h_trigger_ptp_pta_data->Fill(pta, ptp);
	  hp_trigger_ptp_pta_data->Fill(pta, ptp);      
	  hp_cosdphi_pta_data->Fill(pta, fabs(cos(dphi_data)));

      
	  float ema = (em1_data*pt1_data + em2_data*pt2_data)/(pt1_data+pt2_data);
	  for (int iem = 0; iem < nbins_em; iem++)
	    {
	      if (ema >= em_bounds[iem] && ema < em_bounds[iem+1])
		{
		  hp_em_avgpt_data[iem]->Fill(pta, pta);
		  h_em_trigger_ptv_pta_data[iem]->Fill(pta, ptv);
		  hp_em_trigger_ptv_pta_data[iem]->Fill(pta, ptv);
		  h_em_trigger_ptp_pta_data[iem]->Fill(pta, ptp);
		  hp_em_trigger_ptp_pta_data[iem]->Fill(pta, ptp);      
		  hp_em_cosdphi_pta_data[iem]->Fill(pta, fabs(cos(dphi_data)));	    
		  break;
		}
	    }
	}
    }
  std::cout << "done" << std::endl;
  // Sim
  for (int isample = 0; isample < 10; isample++)
    {

      std::cout << "Sample " << isample << std::endl;
      if (!use_sample[isample]) continue;      
      float nevents = ttree[isample]->GetEntries();
      int iisample = isample%5;
      if (isample >= 5)
	{
	  nevents = ttree[iisample]->GetEntries() * 0.22 * (0.22 + 1.);
	}
      for (int i = 0; i < nevents; i++)
	{
	  ttree[isample]->GetEntry(i);
	  std::cout << "Event : " << i << " \r" << std::flush;	 
	  bool trig = true;//(((gl1_scaled[isample] >> 22) & 0x1) == 0x1 );
	  
	  if (fabs(mbd_vertex_z[isample]) > vertex_cut) continue;
	  float maxpttruth = *std::max_element(truth_jet_pt[isample]->begin(), truth_jet_pt[isample]->end());

	  if (conesize != 4)
	    {
	      maxpttruth = *std::max_element(truth_jet_pt_ref[isample]->begin(), truth_jet_pt_ref[isample]->end());
	    }

	  

	  if (maxpttruth < boundary_r4[iisample] || maxpttruth >= boundary_r4[iisample+1]) continue;
	  h_lead_sample[iisample]->Fill(maxpttruth, scale_factor[iisample]);
	  h_lead_combined->Fill(maxpttruth, scale_factor[iisample]);

	  // fully matched dijets
	  matched_dijets.clear();
	  // all truth jets
	  mytruthjets.clear();
	  // all reco jets
	  myrecojets.clear();
	  // top two of each
	  mytruthjets2.clear();
	  myrecojets2.clear();

	  int ntruthjets = truth_jet_pt[isample]->size();
	  for (int j = 0; j < ntruthjets;j++)
	    {

	      if (truth_jet_pt[isample]->at(j) < truth_cut) continue;

	      struct jet tempjet;
	  
	      tempjet.pt = truth_jet_pt[isample]->at(j);
	      tempjet.eta = truth_jet_eta[isample]->at(j);
	      tempjet.phi = truth_jet_phi[isample]->at(j);
	      tempjet.id = j;

	      mytruthjets.push_back(tempjet);	  	  

	    }

	  std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  bool truth_good = check_dijet_truth(mytruthjets);
	  //truth_good &= (fabs(truth_vertex_z[isample]) < vertex_cut);
	  int nrecojets = reco_jet_pt[isample]->size();
	  for (int j = 0; j < nrecojets;j++)
	    {
	      if (reco_jet_pt[isample]->at(j) < reco_cut) continue;
	      if (reco_jet_e[isample]->at(j) < 0) continue;

	      struct jet tempjet;
	      tempjet.pt = reco_jet_pt[isample]->at(j);

	      tempjet.emcal = reco_jet_emcal[isample]->at(j);
	      tempjet.eta = reco_jet_eta[isample]->at(j);
	      tempjet.eta_det = reco_jet_eta_det[isample]->at(j);
	      tempjet.phi = reco_jet_phi[isample]->at(j);
	      tempjet.t = 0;
	      tempjet.id = j;

	      myrecojets.push_back(tempjet);
	    }

	  std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

	  bool reco_good = check_dijet_reco(myrecojets);

	  // get third jet
	  float third_jet_pt = 0;
	  float third_jet_pt_truth = 0;
	  
	  if (myrecojets.size() > 2)
	    {
	      third_jet_pt = myrecojets.at(2).pt;
	    }
	  if (mytruthjets.size() > 2)
	    {
	      third_jet_pt_truth = mytruthjets.at(2).pt;
	    }

	  float trigger_jet = 0;
	  float probe_jet = 0;
	  float dphi_reco = 0;
	  float em1_reco = 0;
	  float em2_reco = 0;
	  float pt1_reco = 0;
	  float pt2_reco = 0;
	  float ptem = 0;


	  if (reco_good)
	    {
	      dphi_reco = fabs(getDPHI(myrecojets.at(0).phi, myrecojets.at(1).phi));
	      pt1_reco = myrecojets.at(0).pt;
	      em1_reco = myrecojets.at(0).emcal;
	      pt2_reco = myrecojets.at(1).pt;
	      em2_reco = myrecojets.at(1).emcal;
	      ptem = (em1_reco * pt1_reco + em2_reco * pt2_reco )/(pt1_reco + pt2_reco);
	      if (closure)
		{
		  double xmin = 3, xmax = 100;
		  double dx = 0.1;
		  int nbins_smear = (int) ((xmax - xmin) / dx) + 1 ;

		  int bin = floor((pt1_reco - xmin)/dx);//hjer->FindBin(pt1_reco);
		  float par = hjer->GetBinContent(bin);
		  if (par > 0)
		    {
		      fsmear->SetParameter(2, hjer->GetBinContent(bin));
		      float jersmear1 = fsmear->GetRandom();
		      pt1_reco += jersmear1*pt1_reco;
		    }
		  bin = floor((pt2_reco - xmin)/dx) + 1;//hjer->FindBin(pt1_reco);
		  
		  par = hjer->GetBinContent(bin);
		  if (par > 0)
		    {
		      fsmear->SetParameter(2, hjer->GetBinContent(bin));
		      float jersmear2 = fsmear->GetRandom();
		      pt2_reco += jersmear2*pt2_reco;
		    }
		}

	      trigger_jet = pt1_reco;
	      probe_jet = pt2_reco;

	      myrecojets2 = {myrecojets.begin(), myrecojets.begin()+2};
	    }

	  float triggertruth_jet = 0;
	  float probetruth_jet = 0;

	  float dphi_truth = 0;
	  if (truth_good)
	    {
	      dphi_truth = fabs(getDPHI(mytruthjets.at(0).phi, mytruthjets.at(1).phi));
	      triggertruth_jet = mytruthjets.at(0).pt;
	      probetruth_jet = mytruthjets.at(1).pt;
	      mytruthjets2 = {mytruthjets.begin(), mytruthjets.begin()+2};
	    }


	  
	  if (rng->Uniform() > 0.5)
	    {
	      float tempy = trigger_jet;
	      trigger_jet = probe_jet;	
	      probe_jet = tempy;
	    }
	  if (rng->Uniform() > 0.5)
	    {
	      float tempy = triggertruth_jet;
	      triggertruth_jet = probetruth_jet;
	      probetruth_jet = tempy;
	    }

	  matched_dijets = match_dijets(myrecojets2, mytruthjets2);
	  bool matched = (matched_dijets.size() == 2);

	      
	  float ptv = trigger_jet*cos(dphi_reco/2) + probe_jet*cos(dphi_reco/2);
	  float ptp = trigger_jet*sin(dphi_reco/2) - probe_jet*sin(dphi_reco/2);
	  float pta = (trigger_jet + probe_jet)/2.;
	  float aj = (trigger_jet - probe_jet) / (trigger_jet + probe_jet);

	  if (rng->Uniform() > 0.5)
	    {
	      ptv = -1*ptv;
	    }
	  if (pta >= first_trigger && !trig) continue;
	  //	  reco_good &= (pta >= 20);// && !((((gl1_scaled[isample] >> 18) & 0x1) == 0x1 || ((gl1_scaled[isample] >> 22) & 0x1) == 0x1)));
	  if (reco_good)
	    {
	      
	      for (int i3 = 0 ; i3 < 8; i3++)
		{
		  if (third_jet_pt < third_jet_cuts[i3])
		    {
		      h2_aj_avgpt_sim_3jet[i3]->Fill(pta, aj, scale_factor[iisample]);
		      hp_aj_avgpt_sim_3jet[i3]->Fill(pta, aj, scale_factor[iisample]);
		      hp_avgpt_sim_3jet[i3]->Fill(pta, pta, scale_factor[iisample]);

       		      h_trigger_ptem_pta_sim_3jet[i3]->Fill(pta, ptem, scale_factor[iisample]);
		      h_lead_ptem_pta_sim_3jet[i3]->Fill(pta, em1_reco, scale_factor[iisample]);
		      h_sublead_ptem_pta_sim_3jet[i3]->Fill(pta, em2_reco, scale_factor[iisample]);
		      hp_trigger_ptem_pta_sim_3jet[i3]->Fill(pta, ptem, scale_factor[iisample]);
		      hp_lead_ptem_pta_sim_3jet[i3]->Fill(pta, em1_reco, scale_factor[iisample]);
		      hp_sublead_ptem_pta_sim_3jet[i3]->Fill(pta, em2_reco, scale_factor[iisample]);
		      
		      h_trigger_ptv_pta_sim_3jet[i3]->Fill(pta, ptv, scale_factor[iisample]);
		      hp_trigger_ptv_pta_sim_3jet[i3]->Fill(pta, ptv, scale_factor[iisample]);
		      h_trigger_ptp_pta_sim_3jet[i3]->Fill(pta, ptp, scale_factor[iisample]);
		      hp_trigger_ptp_pta_sim_3jet[i3]->Fill(pta, ptp, scale_factor[iisample]);
		      hp_cosdphi_pta_sim_3jet[i3]->Fill(pta, fabs(cos(dphi_reco)), scale_factor[iisample]);
		    }
		}
	      hp_avgpt_sim->Fill(pta,pta, scale_factor[iisample]);
	      h_trigger_ptv_pta_sim->Fill(pta, ptv, scale_factor[iisample]);
	      hp_trigger_ptv_pta_sim->Fill(pta, ptv, scale_factor[iisample]);
	      h_trigger_ptp_pta_sim->Fill(pta, ptp, scale_factor[iisample]);
	      hp_trigger_ptp_pta_sim->Fill(pta, ptp, scale_factor[iisample]);
	      hp_cosdphi_pta_sim->Fill(pta, fabs(cos(dphi_reco)), scale_factor[iisample]);

	      float ema = (em1_reco*pt1_reco + em2_reco*pt2_reco)/(pt1_reco+pt2_reco);

	      for (int iem = 0; iem < 5; iem++)
		{
		  if (ema >= em_bounds[iem] && ema < em_bounds[iem+1])
		    {
		      hp_em_avgpt_sim[iem]->Fill(pta,pta, scale_factor[iisample]);
		      h_em_trigger_ptv_pta_sim[iem]->Fill(pta, ptv, scale_factor[iisample]);
		      hp_em_trigger_ptv_pta_sim[iem]->Fill(pta, ptv, scale_factor[iisample]);
		      h_em_trigger_ptp_pta_sim[iem]->Fill(pta, ptp, scale_factor[iisample]);
		      hp_em_trigger_ptp_pta_sim[iem]->Fill(pta, ptp, scale_factor[iisample]);      
		      hp_em_cosdphi_pta_sim[iem]->Fill(pta, fabs(cos(dphi_reco)), scale_factor[iisample]);	    
		      break;
		    }
		}
	    }
	  if (truth_good)
	    {
	      float ptat = (triggertruth_jet + probetruth_jet)/2.;
	      float ajt = (triggertruth_jet - probetruth_jet) / (triggertruth_jet + probetruth_jet);
	      
	      for (int i3 = 0 ; i3 < 8; i3++)
		{
		  if (third_jet_pt_truth < third_jet_cuts[i3])
		    {
		      h2_aj_avgpt_truth_3jet[i3]->Fill(ptat, ajt, scale_factor[iisample]);
		      hp_aj_avgpt_truth_3jet[i3]->Fill(ptat, ajt, scale_factor[iisample]);
		      hp_avgpt_truth_3jet[i3]->Fill(ptat,ptat, scale_factor[iisample]);
		    }
		}
	    }

	}
    }


  TH1D *h_sim_3jet_width[nbins_pt];
  TH1D *h_data_3jet_width[nbins_pt];
  TH1D *h_truth_3jet_width[nbins_pt];

  TH1D *h_sim_3jet_softcorr[8];
  TH1D *h_data_3jet_softcorr[8];
  TH1D *h_truth_3jet_softcorr[8];

  for (int i = 0; i < 8; i++)
    {
      h_truth_3jet_softcorr[i] = new TH1D(Form("h_truth_3jet_softcorr_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_sim_3jet_softcorr[i] = new TH1D(Form("h_sim_3jet_softcorr_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
      h_data_3jet_softcorr[i] = new TH1D(Form("h_data_3jet_softcorr_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", nbins_pt, pt_bin_edges);
    }
  
  TH1D *h_sim_3jet_sigv[nbins_pt];
  TH1D *h_data_3jet_sigv[nbins_pt];

  TH1D *h_sim_3jet_sigp[nbins_pt];
  TH1D *h_data_3jet_sigp[nbins_pt];

  TH1D *h_sim_3jet_sigma[nbins_pt];
  TH1D *h_data_3jet_sigma[nbins_pt];

  TH1D *h_sim_3jet_bi_width[nbins_pt];
  TH1D *h_data_3jet_bi_width[nbins_pt];

  TF1 *fg = new TF1("fg","gaus", -20, 20);
  TF1 *fg_v = new TF1("fg_v","gaus", -6, 6);
  TF1 *fg_p = new TF1("fg_p","gaus", -20, 20);
  fg->FixParameter(1, 0);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < nbins_pt; i++)
    {
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      h_sim_3jet_sigma[i] = new TH1D(Form("h_sim_3jet_sigma_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #psi})", 10, 0.5, 10.5);
      h_data_3jet_sigma[i] = new TH1D(Form("h_data_3jet_sigma_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #psi})", 10, 0.5, 10.5);
      
      h_sim_3jet_sigp[i] = new TH1D(Form("h_sim_3jet_sigp_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #psi})", 10, 0.5, 10.5);
      h_data_3jet_sigp[i] = new TH1D(Form("h_data_3jet_sigp_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #psi})", 10, 0.5, 10.5);

      h_sim_3jet_sigv[i] = new TH1D(Form("h_sim_3jet_sigv_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #lambda})", 10, 0.5, 10.5);
      h_data_3jet_sigv[i] = new TH1D(Form("h_data_3jet_sigv_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T, #lambda})", 10, 0.5, 10.5);
      
      h_sim_3jet_bi_width[i] = new TH1D(Form("h_sim_3jet_bi_width_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", 10, 0.5, 10.5);
      h_data_3jet_bi_width[i] = new TH1D(Form("h_data_3jet_bi_width_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", 10, 0.5, 10.5);

      for (int j = 0; j < 8; j++)
	{

	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  
	  TH1D *h_simv = (TH1D*) h_trigger_ptv_pta_sim_3jet[j]->ProjectionY("hv2", i+1, i+1);
	  h_simv->Scale(1./h_simv->Integral(), "width");
	  h_simv->Fit("fg", "QR+","", -6, 6);
	  
	  std::cout << " WIDTH: " << fg->GetParameter(2) << " / " << h_simv->GetRMS() << std::endl;

	  //float sigv = fg->GetParameter(2);
	  float sigv = h_simv->GetRMS();
	  float sig2v = h_simv->GetRMSError();


	  TH1D *h_simp = (TH1D*) h_trigger_ptp_pta_sim_3jet[j]->ProjectionY("hp2", i+1, i+1);
	  h_simp->Scale(1./h_simp->Integral(), "width");
	  h_simp->Fit("fg", "QR+","", -20, 20);

	  //float sigp = fg->GetParameter(2);
	  float sigp = h_simp->GetRMS();
	  float sig2p = h_simp->GetRMSError();//float sig2p = fg->GetParError(2);

	  std::cout << " WIDTH (F/R): " << fg->GetParameter(2) << " / " << h_simp->GetRMS() << std::endl;

	  h_sim_3jet_sigv[i]->SetBinContent(int(third_jet_cuts[j]), sigv);
	  h_sim_3jet_sigv[i]->SetBinError(int(third_jet_cuts[j]), sig2v);

	  h_sim_3jet_sigp[i]->SetBinContent(int(third_jet_cuts[j]), sigp);
	  h_sim_3jet_sigp[i]->SetBinError(int(third_jet_cuts[j]), sig2p);

	  float pta = hp_avgpt_sim_3jet[j]->GetBinContent(i+1);
	  
	  float cosdphi = hp_cosdphi_pta_sim_3jet[j]->GetBinContent(i+1);
	  float sigcosdphi = hp_cosdphi_pta_sim_3jet[j]->GetBinError(i+1);
	  float quad_sub = sqrt(sigp*sigp - sigv*sigv);

	  float sig2_quad_sub = TMath::Power(sigp*sig2p/quad_sub, 2) + TMath::Power(sigv*sig2v/quad_sub, 2);
	  float sig_quad_sub = sqrt(sig2_quad_sub);
	  float res = quad_sub /(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig_quad_sub/quad_sub, 2) + TMath::Power(sigcosdphi/cosdphi, 2));
	  
	  if (sigp <= sigv) continue;

	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  std::cout << "sim 3jet : " << j <<  " bin : " << i << " : " << "sigp ( " << sigp << " +- "<<sig2p<<") sige ( " << sigv << " +- " << sig2v << ") cosdphi  ( " << cosdphi << " +- " << sigcosdphi << " ) res ( " <<  res << " +- " << sigres << " ) " << std::endl;
	  
	  h_sim_3jet_sigma[i]->SetBinContent(int(third_jet_cuts[j]), quad_sub);
	  h_sim_3jet_sigma[i]->SetBinError(int(third_jet_cuts[j]), sig_quad_sub);
	  
	  h_sim_3jet_bi_width[i]->SetBinContent(int(third_jet_cuts[j]), res);
	  h_sim_3jet_bi_width[i]->SetBinError(int(third_jet_cuts[j]), sigres);

	  h_jer_sim_3jet[j]->SetBinContent(i+1, res);
	  h_jer_sim_3jet[j]->SetBinError(i+1, sigres);
	}
      for (int j = 0; j < 8; j++)
	{
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  TH1D *h_datav = (TH1D*) h_trigger_ptv_pta_data_3jet[j]->ProjectionY("hv", i+1, i+1);
	  h_datav->Scale(1./h_datav->Integral(), "width");	  
	  h_datav->Fit("fg", "Q+","", -6, 6);
	  //float sigv = fg->GetParameter(2);
	  float sigv = h_datav->GetRMS();
	  float sig2v = h_datav->GetRMSError();
	  //float sig2v = fg->GetParError(2);
	  TH1D *h_datap = (TH1D*) h_trigger_ptp_pta_data_3jet[j]->ProjectionY("hp", i+1, i+1);
	  h_datap->Scale(1./h_datap->Integral(), "width");
	  h_datap->Fit("fg", "Q+","", -20, 20);
	  //float sigp = fg->GetParameter(2);
	  float sigp = h_datap->GetRMS();
	  float sig2p = h_datap->GetRMSError();//	  float sig2p = fg->GetParError(2);

	  /* h_trigger_ptv_pta_data_3jet[j]->FitSlicesY(fg, i+1, i+1); */
	  /* float sigv = fg->GetParameter(2); */
	  /* float sig2v = fg->GetParError(2); */
	  /* h_trigger_ptp_pta_data_3jet[j]->FitSlicesY(fg, i+1, i+1); */
	  /* float sigp = fg->GetParameter(2); */
	  /* float sig2p = fg->GetParError(2); */
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  
	  h_data_3jet_sigv[i]->SetBinContent(int(third_jet_cuts[j]), sigv);//j, third_jet_cuts[j], sigv);
	  h_data_3jet_sigv[i]->SetBinError(int(third_jet_cuts[j]), sig2v);

	  h_data_3jet_sigp[i]->SetBinContent(int(third_jet_cuts[j]), sigp);
	  h_data_3jet_sigp[i]->SetBinError(int(third_jet_cuts[j]), sig2p);
 

	  float pta = hp_avgpt_data_3jet[j]->GetBinContent(i+1);
	  float sigpta = hp_trigger_ptv_pta_data_3jet[j]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_cosdphi_pta_data_3jet[j]->GetBinContent(i+1);
	  float sigcosdphi = hp_cosdphi_pta_data_3jet[j]->GetBinError(i+1);
	  float quad_sub = sqrt(sigp*sigp - sigv*sigv);
	  float sig2_quad_sub = TMath::Power(sigp*sig2p/quad_sub, 2) + TMath::Power(sigv*sig2v/quad_sub, 2);
	  float sig_quad_sub = sqrt(sig2_quad_sub);
	  float res = quad_sub /(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig_quad_sub/quad_sub, 2) + TMath::Power(sigcosdphi/cosdphi, 2));
	  
	  if (sigp <= sigv) continue;
	  std::cout << "data 3jet : " << j <<  " bin : " << i << " : " << "sigp ( " << sigp << " +- "<<sig2p<<") sige ( " << sigv << " +- " << sig2v << ") cosdphi  ( " << cosdphi << " +- " << sigcosdphi << " ) res ( " <<  res << " +- " << sigres << " ) " << std::endl;
	  if (NUCLEAR) std::cout << __LINE__ << std::endl;
	  h_data_3jet_sigma[i]->SetBinContent(int(third_jet_cuts[j]), quad_sub);
	  h_data_3jet_sigma[i]->SetBinError(int(third_jet_cuts[j]), sig_quad_sub);
	  
	  h_data_3jet_bi_width[i]->SetBinContent(int(third_jet_cuts[j]), res);
	  h_data_3jet_bi_width[i]->SetBinError(int(third_jet_cuts[j]), sigres);

	  h_jer_data_3jet[j]->SetBinContent(i+1, res);
	  h_jer_data_3jet[j]->SetBinError(i+1, sigres);
	  
	}

    }

  TH1D *h_aj_data[nbins_pt][8];
  TH1D *h_aj_sim[nbins_pt][8];
  TH1D *h_aj_truth[nbins_pt][8];

  TF1 *fgg = new TF1("fgg", "gaus", -0.2, 0.2);
  TF1 *fgg2 = new TF1("fgg2", "gaus(0) + gaus(3)", -0.2, 0.2);
  fgg->SetParameters(1, 0, 0.2);
  fgg->FixParameter(1, 0);
  fgg->SetParLimits(2, 0, 0.5);
  fgg2->SetParameters(1, 0, 0.2, 0.2, 0, 0.5);
  fgg2->FixParameter(1, 0);
  fgg2->FixParameter(4, 0);

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < nbins_pt; i++)
    {
      

      h_truth_3jet_width[i] = new TH1D(Form("h_truth_3jet_width_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", 10, 0.5, 10.5);
      h_sim_3jet_width[i] = new TH1D(Form("h_sim_3jet_width_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", 10, 0.5, 10.5);
      h_data_3jet_width[i] = new TH1D(Form("h_data_3jet_width_%d", i), ";#it{p}_{T,3}; #sigma(#it{p}_{T})/#it{p}_{T}", 10, 0.5, 10.5);

      
      for (int j = 0; j < 8; j++)
	{
	  h_aj_data[i][j] = (TH1D*) h2_aj_avgpt_data_3jet[j]->ProjectionY(Form("h_aj_data_%d_%d", i, j), i+1, i+1);
	  h_aj_sim[i][j] = (TH1D*) h2_aj_avgpt_sim_3jet[j]->ProjectionY(Form("h_aj_sim_%d_%d", i, j), i+1, i+1);
	  h_aj_truth[i][j] = (TH1D*) h2_aj_avgpt_truth_3jet[j]->ProjectionY(Form("h_aj_truth_%d_%d", i, j), i+1, i+1);

	  h_aj_data[i][j]->Scale(1./h_aj_data[i][j]->Integral(),"width");//;
	  h_aj_sim[i][j]->Scale(1./h_aj_sim[i][j]->Integral(),"width");//;
	  h_aj_truth[i][j]->Scale(1./h_aj_truth[i][j]->Integral(),"width");//;
	  h_aj_sim[i][j]->Fit("fgg","R", "",-0.2, 0.2);
	  /* fgg2->SetParameters(fgg->GetParameter(0), fgg->GetParameter(1), fgg->GetParameter(2), 0.1, 0, 0.5); */
	  /* fgg2->SetParLimits(0, fgg->GetParameter(0) - 0.5,fgg->GetParameter(0) + 0.5); */
	  /* fgg2->FixParameter(1, 0); */
	  /* fgg2->SetParLimits(2, fgg->GetParameter(2) - 0.05,fgg->GetParameter(2) + 0.05); */
	  /* fgg2->SetParLimits(3, 0, 0.5); */
	  /* fgg2->FixParameter(4, 0); */
	  /* fgg2->SetParLimits(5, 0, 1); */

	  /* h_aj_sim[i][j]->Fit("fgg2","R", "",-0.5, 0.5); */
	  float stdev = h_aj_sim[i][j]->GetRMS()*sqrt(2);     //fgg->GetParameter(2)*sqrt(2);
	  float sedev = h_aj_sim[i][j]->GetRMSError()*sqrt(2);// = fgg->GetParError(2)*sqrt(2);
	  //float stdev = fgg->GetParameter(2)*sqrt(2);
	  //float sedev = fgg->GetParError(2)*sqrt(2);

	  h_sim_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_sim_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);

	  h_jer_sim_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_sim_im_3jet[j]->SetBinError(i+1, sedev);
	  h_aj_data[i][j]->Fit("fgg","R", "",-0.2, 0.2);
	  /* fgg2->SetParameters(fgg->GetParameter(0), fgg->GetParameter(1), fgg->GetParameter(2), 0.1, 0, 0.5); */
	  /* fgg2->SetParLimits(0, fgg->GetParameter(0) - 0.5,fgg->GetParameter(0) + 0.5); */
	  /* fgg2->FixParameter(1, 0); */
	  /* fgg2->SetParLimits(2, fgg->GetParameter(2) - 0.05,fgg->GetParameter(2) + 0.05); */
	  /* fgg2->SetParLimits(3, 0, 0.5); */
	  /* fgg2->FixParameter(4, 0); */
	  /* fgg2->SetParLimits(5, 0, 1); */

	  /* h_aj_data[i][j]->Fit("fgg2","R", "",-0.5, 0.5); */

	  stdev = h_aj_data[i][j]->GetRMS()*sqrt(2);//fgg->GetParameter(2)*sqrt(2);
	  sedev = h_aj_data[i][j]->GetRMSError()*sqrt(2);// = fgg->GetParError(2)*sqrt(2);
	  //stdev = fgg->GetParameter(2)*sqrt(2);
	  //sedev = fgg->GetParError(2)*sqrt(2);

	  //stdev = fgg->GetParameter(2)*sqrt(2);
	  //sedev = fgg->GetParError(2)*sqrt(2);

	  h_jer_data_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_data_im_3jet[j]->SetBinError(i+1, sedev);
	  
	  h_data_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_data_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);

	  h_aj_truth[i][j]->Fit("fgg","RQ", "",-0.1, 0.1);

	  fgg2->SetParameters(fgg->GetParameter(0), fgg->GetParameter(1), fgg->GetParameter(2), 0.2, 0, 0.5);
	  fgg2->FixParameter(1, 0);
	  fgg2->FixParameter(4, 0);
	  h_aj_truth[i][j]->Fit("fgg2","RQ", "",-0.4, 0.4);
	  stdev = h_aj_truth[i][j]->GetRMS()*sqrt(2);//fgg->GetParameter(2)*sqrt(2);
	  sedev = h_aj_truth[i][j]->GetRMSError()*sqrt(2);// = fgg->GetParError(2)*sqrt(2);
	  //stdev = fgg2->GetParameter(2)*sqrt(2);
	  //sedev = fgg2->GetParError(2)*sqrt(2);


	  h_jer_truth_im_3jet[j]->SetBinContent(i+1, stdev);
	  h_jer_truth_im_3jet[j]->SetBinError(i+1, sedev);
	  
	  h_truth_3jet_width[i]->SetBinContent(int(third_jet_cuts[j]), stdev);
	  h_truth_3jet_width[i]->SetBinError(int(third_jet_cuts[j]), sedev);
	}
    }
  TF1 *fline = new TF1("fline","pol1", 5, 10);
  // Linear Fit And Extrapolation of soft corrections
  std::cout << "---------------- Estrapolation -------------------" <<std::endl; 
  for (int i = 0; i < nbins_pt ; i++)
    {


      h_sim_3jet_width[i]->Fit("fline","R+","", 5, 10);
      float sim_at_zero = fline->Eval(0);
      for (int j = 0; j < 8; j++)
	{
	  
	  float sim_at_ten = fline->Eval(third_jet_cuts[j]);
	  float ksoft_sim = sim_at_zero/sim_at_ten;
	  h_sim_3jet_softcorr[j]->SetBinContent(i+1, ksoft_sim);
	  float res = h_jer_sim_im_3jet[j]->GetBinContent(i+1)*ksoft_sim;
	  float sigres = h_jer_sim_im_3jet[j]->GetBinError(i+1);
	  	  
	  std::cout << "sim " << i << " , " << j << " --> " << ksoft_sim << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_sim_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_sim_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      std::cout << __LINE__ <<" line " << std::endl; h_data_3jet_width[i]->Fit("fline","R+","", 5, 10);

      float data_at_zero = fline->Eval(0);
      // float data_at_ten = fline->Eval(10);
      // float ksoft_data = data_at_zero/data_at_ten;
      for (int j = 0; j < 8; j++)
	{
	  float data_at_ten = fline->Eval(third_jet_cuts[j]);
	  //	  float data_at_ten = fline->Eval(h_jer_data_im_3jet[j]->GetBinCenter(i+1));
	  float ksoft_data = data_at_zero/data_at_ten;
	  h_data_3jet_softcorr[j]->SetBinContent(i+1, ksoft_data);
	  float res = h_jer_data_im_3jet[j]->GetBinContent(i+1)*ksoft_data;
	  float sigres = h_jer_data_im_3jet[j]->GetBinError(i+1)*ksoft_data;

	  std::cout << "data pt " << i << " , " << third_jet_cuts[j] << " --> " << data_at_zero << " / " << data_at_ten << " = " << ksoft_data << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_data_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_data_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      h_truth_3jet_width[i]->Fit("fline","R+","", 5, 10);
      float truth_at_zero = fline->Eval(0);
      //float truth_at_ten = fline->Eval(10);
      //float ksoft_truth = truth_at_zero/truth_at_ten;
      for (int j = 0; j < 8; j++)
	{
	  float truth_at_ten = fline->Eval(third_jet_cuts[j]);
	  float ksoft_truth = truth_at_zero/truth_at_ten;
	  h_truth_3jet_softcorr[j]->SetBinContent(i+1, ksoft_truth);
	  float res = h_jer_truth_im_3jet[j]->GetBinContent(i+1)*ksoft_truth;
	  float sigres = h_jer_truth_im_3jet[j]->GetBinError(i+1)*ksoft_truth;

	  std::cout << "truth " << i << " , " << j << " --> " << ksoft_truth << " - " << res << " - " << sigres <<  std::endl;
	  h_jer_truth_im_softcorr_3jet[j]->SetBinContent(i+1, res);
	  h_jer_truth_im_softcorr_3jet[j]->SetBinError(i+1, sigres);
	}

      // No begins the subtraction of the truth
  
      for (int j = 0; j < 8; j++)
	{
	  float tres = h_jer_truth_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float tsigres = h_jer_truth_im_softcorr_3jet[j]->GetBinError(i+1);

	  float sres = h_jer_sim_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float ssigres = h_jer_sim_im_softcorr_3jet[j]->GetBinError(i+1);

	  float dres = h_jer_data_im_softcorr_3jet[j]->GetBinContent(i+1);
	  float dsigres = h_jer_data_im_softcorr_3jet[j]->GetBinError(i+1);

	  float sim_new_res = sqrt(sres*sres - tres*tres);
	  float data_new_res = sqrt(dres*dres - tres*tres);
	  float data_new_sigres = sqrt(TMath::Power(tres*tsigres/data_new_res, 2) + TMath::Power(dres*dsigres/data_new_res, 2));
	  float sim_new_sigres = sqrt(TMath::Power(tres*tsigres/sim_new_res, 2) + TMath::Power(sres*ssigres/sim_new_res, 2));

	  std::cout << "data " << i << " , " << j << " --> " << dres << " - " << tres << " --> " << data_new_res << std::endl;
	  std::cout << "sim " << i << " , " << j << " --> " << sres << " - " << tres << " --> " << sim_new_res << std::endl;
	  h_jer_sim_im_softcorr_part_3jet[j]->SetBinContent(i+1, sim_new_res);
	  h_jer_sim_im_softcorr_part_3jet[j]->SetBinError(i+1, sim_new_sigres);
	  h_jer_data_im_softcorr_part_3jet[j]->SetBinContent(i+1, data_new_res);
	  h_jer_data_im_softcorr_part_3jet[j]->SetBinError(i+1, data_new_sigres);
	}

    }

  // 
  int color_sim = kRed;
  int color_truth = kBlack;
  int color_data = kBlue;

  TCanvas *ctt = new TCanvas("ctt","ctt", 500, 500);
  for (int i3 = 0; i3 < 8; i3++)
    {
      dlutility::SetMarkerAtt(h_data_3jet_softcorr[i3], color_data, 1, 8);
      dlutility::SetLineAtt(h_data_3jet_softcorr[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_truth_3jet_softcorr[i3], color_truth, 1, 8);
      dlutility::SetLineAtt(h_truth_3jet_softcorr[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_softcorr[i3], color_sim, 1, 8);
      dlutility::SetLineAtt(h_sim_3jet_softcorr[i3], color_sim, 1, 1);
      h_data_3jet_softcorr[i3]->SetMaximum(1.5);
      h_data_3jet_softcorr[i3]->SetMinimum(0);
      h_data_3jet_softcorr[i3]->Draw("same p");
      h_sim_3jet_softcorr[i3]->Draw("same p");
      h_truth_3jet_softcorr[i3]->Draw("same p");

      dlutility::DrawSPHENIXpp(0.22, 0.82);
      dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
      TLegend *legg = new TLegend(0.6, 0.77, 0.75, 0.87);
      legg->SetLineWidth(0);
      legg->AddEntry(h_data_3jet_softcorr[i3], "Data");
      legg->AddEntry(h_sim_3jet_softcorr[i3], "Sim");
      legg->AddEntry(h_truth_3jet_softcorr[i3], "Truth");

      legg->Draw("same");

      ctt->Print(Form("r%02d/jer//softcorr_3jet%d_r%02d_%d_closure_%s.pdf", conesize , i3, conesize, closure, gen[index].c_str()));
      ctt->Print(Form("r%02d/jer//softcorr_3jet%d_r%02d_%d_closure_%s.png", conesize , i3, conesize, closure, gen[index].c_str()));

    }
  
  
  TH1D *h_blank_res = new TH1D("hh", "; #it{p}_{T3} [GeV]; #sigma(#it{p}_{T})/#it{p}_{T}", 1, 0, 10);
  h_blank_res->SetMinimum(0.0);
  h_blank_res->SetMaximum(0.5);
  
  TCanvas *ctg = new TCanvas("ctg", "ctg", 500, 500);
  for (int i = 0; i < nbins_pt; i++)
    {
      std::cout << __LINE__ << std::endl;
      dlutility::SetLineAtt(h_sim_3jet_bi_width[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_bi_width[i], color_sim, 0.8, 21);
      dlutility::SetLineAtt(h_data_3jet_bi_width[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_bi_width[i], color_data, 0.8, 21);
      std::cout << __LINE__ << std::endl;
      dlutility::SetLineAtt(h_sim_3jet_width[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_width[i], color_sim, 0.8, 20);
      std::cout << __LINE__ << std::endl;
      dlutility::SetLineAtt(h_truth_3jet_width[i], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_truth_3jet_width[i], color_truth, 0.8, 20);
      std::cout << __LINE__ << std::endl;
      dlutility::SetLineAtt(h_data_3jet_width[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_width[i], color_data, 0.8, 20);
      h_blank_res->DrawCopy();
      h_data_3jet_width[i]->Draw("same p");
      h_sim_3jet_width[i]->Draw("same p");
      h_truth_3jet_width[i]->Draw("same p");
      h_data_3jet_bi_width[i]->Draw("same p");
      h_sim_3jet_bi_width[i]->Draw("same p");
      std::cout << __LINE__ << std::endl;
      dlutility::DrawSPHENIXpp(0.2, 0.82);
      dlutility::drawText(Form("%2.0f < #LT#it{p}_{T}#GT < %2.0f GeV", pt_bin_edges[i], pt_bin_edges[i+1]), 0.2, 0.7); 

      TLegend *leg = new TLegend(0.55, 0.67, 0.8, 0.87);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_3jet_width[i], "Data - Imbalance");
      leg->AddEntry(h_sim_3jet_width[i], "Sim - Imbalance");
      leg->AddEntry(h_data_3jet_bi_width[i], "Data - Bisector");
      leg->AddEntry(h_sim_3jet_bi_width[i], "Sim - Bisector");

      leg->Draw("same");

      ctg->Print(Form("r%02d/jer//width3jet_r%02d_pta%d_%d_closure_%s.pdf", conesize , conesize, i, closure, gen[index].c_str()));
      ctg->Print(Form("r%02d/jer//width3jet_r%02d_pta%d_%d_closure_%s.png", conesize , conesize, i, closure, gen[index].c_str()));

    }


  for (int i = 0; i < nbins_pt; i++)
    {
      for (int j = 0; j < 8; j++)
	{
	  ctg->SetLogy(0);
	  dlutility::SetLineAtt(h_aj_sim[i][j], color_sim, 1, 1);
	  dlutility::SetMarkerAtt(h_aj_sim[i][j], color_sim, 0.8, 1);
	  dlutility::SetLineAtt(h_aj_truth[i][j], color_truth, 1, 1);
	  dlutility::SetMarkerAtt(h_aj_truth[i][j], color_truth, 0.8, 1);
	  dlutility::SetLineAtt(h_aj_data[i][j], color_data, 1, 1);
	  dlutility::SetMarkerAtt(h_aj_data[i][j], color_data, 0.8, 8);

	  h_aj_sim[i][j]->SetFillColorAlpha(kRed - 3, 0.3);
	  h_aj_data[i][j]->GetFunction("fgg")->SetLineColor(kBlue + 2);//SetFillColorAlpha(kBlue - 3, 0.3);
	  h_aj_data[i][j]->GetFunction("fgg")->SetLineWidth(2);//SetFillColorAlpha(kBlue - 3, 0.3);
	  h_aj_truth[i][j]->GetFunction("fgg2")->SetLineColor(kBlack);//SetFillColorAlpha(kBlue - 3, 0.3);
	  h_aj_truth[i][j]->GetFunction("fgg2")->SetLineWidth(2);//SetFillColorAlpha(kBlue - 3, 0.3);

	  h_aj_sim[i][j]->GetFunction("fgg")->SetLineColor(kRed + 2);//SetFillColorAlpha(kBlue - 3, 0.3);
	  h_aj_sim[i][j]->GetFunction("fgg")->SetLineWidth(2);//SetFillColorAlpha(kBlue - 3, 0.3);

	  //h_aj_sim[i][j]->SetMaximum(100);
	  h_aj_truth[i][j]->SetTitle("; A_{J} ; #frac{1}{N_{pair}} #frac{dN_{pair}}{dA_{J}}");
	  dlutility::SetFont(h_aj_truth[i][j], 42, 0.05, 0.04, 0.04, 0.04);
	  h_aj_truth[i][j]->GetYaxis()->SetTitleOffset(1.5);
	  h_aj_truth[i][j]->Draw("p");
	  h_aj_sim[i][j]->Draw("same hist");
	  h_aj_sim[i][j]->Draw("same p");
	  h_aj_truth[i][j]->Draw("same p");
	  h_aj_data[i][j]->Draw("same p");
	  dlutility::DrawSPHENIXppInternalsize(0.22, 0.85, 0.04);
	  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.76);
	  dlutility::drawText(Form("%2.1f < #LT#it{p}_{T}#GT #leq %2.1f GeV",  pt_bin_edges[i], pt_bin_edges[i+1]), 0.22, 0.71);
	  dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[j]), 0.22, 0.66);
	  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.87);
	  leg->SetLineWidth(0);
	  leg->AddEntry(h_aj_data[i][j], "Data");
	  leg->AddEntry(h_aj_sim[i][j], Form("Reco. %s", gen[index].c_str()));
	  leg->AddEntry(h_aj_truth[i][j], Form("%s", gen[index].c_str()));
	  leg->Draw("same");

	  ctg->Print(Form("r%02d/jer//aj_r%02d_3jet%d_pta%d_%d_closure_%s.pdf", conesize , conesize, j, i, closure, gen[index].c_str()));
	  ctg->Print(Form("r%02d/jer//aj_r%02d_3jet%d_pta%d_%d_closure_%s.png", conesize , conesize, j, i, closure, gen[index].c_str()));
	}
    }

  h_blank_res->SetMaximum(10);
  TCanvas *ctg2 = new TCanvas("ctg2", "ctg2", 500, 500);
  for (int i = 0; i < nbins_pt; i++)
    {

      dlutility::SetLineAtt(h_sim_3jet_sigp[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigp[i], color_sim, 0.8, 8);

      dlutility::SetLineAtt(h_data_3jet_sigp[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigp[i], color_data, 0.8, 8);

      dlutility::SetLineAtt(h_sim_3jet_sigv[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigv[i], color_sim, 0.8, 21);

      dlutility::SetLineAtt(h_data_3jet_sigv[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigv[i], color_data, 0.8, 21);

      dlutility::SetLineAtt(h_sim_3jet_sigma[i], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_sim_3jet_sigma[i], color_sim, 0.8, 24);

      dlutility::SetLineAtt(h_data_3jet_sigma[i], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_data_3jet_sigma[i], color_data, 0.8, 24);

      h_blank_res->DrawCopy();
      h_data_3jet_sigp[i]->Draw("same p");
      h_sim_3jet_sigp[i]->Draw("same p");
      h_data_3jet_sigv[i]->Draw("same p");
      h_sim_3jet_sigv[i]->Draw("same p");
      h_data_3jet_sigma[i]->Draw("same p");
      h_sim_3jet_sigma[i]->Draw("same p");

      dlutility::DrawSPHENIXpp(0.22, 0.85);
      dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.75);
      dlutility::drawText(Form("%2.1f < #LT#it{p}_{T}#GT #leq %2.1f GeV", pt_bin_edges[i], pt_bin_edges[i+1]), 0.22, 0.7);
      TLegend *legtg = new TLegend(0.6, 0.57, 0.75, 0.877);
      legtg->SetLineWidth(0);
      legtg->SetTextSize(0.04);
      legtg->SetTextFont(42);
      legtg->AddEntry(h_data_3jet_sigp[i], "Data #sigma(#it{p}_{T, #psi})");
      legtg->AddEntry(h_data_3jet_sigv[i],  "Data #sigma(#it{p}_{T, #lambda})");
      legtg->AddEntry(h_data_3jet_sigma[i], "Data #sigma(#it{p}_{T, #psi} - #it{p}_{T, #lambda} ");
      legtg->AddEntry(h_sim_3jet_sigp[i], "Sim #sigma(#it{p}_{T, #psi})");
      legtg->AddEntry(h_sim_3jet_sigv[i],  "Sim #sigma(#it{p}_{T, #lambda})");
      legtg->AddEntry(h_sim_3jet_sigma[i], "Sim #sigma(#it{p}_{T, #psi} - #it{p}_{T, #lambda} ");
      legtg->Draw("same");

      ctg2->Print(Form("r%02d/jer//sigma3jet_r%02d_pta%d_%d_closure_%s.pdf", conesize , conesize, i, closure, gen[index].c_str()));//, i));
      ctg2->Print(Form("r%02d/jer//sigma3jet_r%02d_pta%d_%d_closure_%s.png", conesize , conesize, i, closure, gen[index].c_str()));

    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);

  dlutility::SetLineAtt(hp_trigger_ptv_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_trigger_ptv_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_trigger_ptv_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_trigger_ptv_pta_data, color_data, 1, 8);

  hp_trigger_ptv_pta_sim->Draw("p E1");
  hp_trigger_ptv_pta_data->Draw("same p E1");
  TCanvas *c3 = new TCanvas("c3","c3", 500, 500);

  dlutility::SetLineAtt(hp_trigger_ptp_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_trigger_ptp_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_trigger_ptp_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_trigger_ptp_pta_data, color_data, 1, 8);

  hp_trigger_ptp_pta_sim->Draw("p E1");
  hp_trigger_ptp_pta_data->Draw("same p E1");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c4 = new TCanvas("c4","c4", 500, 500);

  dlutility::SetLineAtt(hp_cosdphi_pta_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(hp_cosdphi_pta_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(hp_cosdphi_pta_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(hp_cosdphi_pta_data, color_data, 1, 8);

  hp_cosdphi_pta_sim->Draw("p E1");
  hp_cosdphi_pta_data->Draw("same p E1");

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int i = 0; i < h_jer_sim->GetNbinsX(); i++)
    {

      TH1D *h_simv = (TH1D*) h_trigger_ptv_pta_sim->ProjectionY("hv", i+1, i+1);
      h_simv->Fit("fg", "R");
      float sigv = fg->GetParameter(2);
      float sig2v = fg->GetParError(2);

      TH1D *h_simp = (TH1D*) h_trigger_ptp_pta_sim->ProjectionY("hp", i+1, i+1);
      h_simp->Fit("fg", "R");
      float sigp = fg->GetParameter(2);
      float sig2p = fg->GetParError(2);

      h_sigv_sim->SetBinContent(i+1, sigv);
      h_sigv_sim->SetBinError(i+1, sig2v);
 
      h_sigp_sim->SetBinContent(i+1, sigp);
      h_sigp_sim->SetBinError(i+1, sig2p);

      float pta = hp_avgpt_sim->GetBinContent(i+1);
      float sigpta = hp_trigger_ptv_pta_sim->GetBinWidth(i+1)/2.;
      float cosdphi = hp_cosdphi_pta_sim->GetBinContent(i+1);
      float sigcosdphi = hp_cosdphi_pta_sim->GetBinError(i+1);
      float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
      float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
      std::cout << "Sim " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
      if (sigp <= sigv) continue;
      h_jer_sim->SetBinContent(i+1, res);
      h_jer_sim->SetBinError(i+1, sigres);

    }


  for (int i = 0; i < h_jer_data->GetNbinsX(); i++)
    {
      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      TH1D *h_datav = (TH1D*) h_trigger_ptv_pta_data->ProjectionY("hv", i+1, i+1);
      h_datav->Fit("fg", "R");
      float sigv = fg->GetParameter(2);
      float sig2v = fg->GetParError(2);
      TH1D *h_datap = (TH1D*) h_trigger_ptp_pta_data->ProjectionY("hp", i+1, i+1);
      h_datap->Fit("fg", "R");
      float sigp = fg->GetParameter(2);
      float sig2p = fg->GetParError(2);
 
      h_sigp_data->SetBinContent(i+1, sigp);
      h_sigp_data->SetBinError(i+1, sig2p);
      h_sigv_data->SetBinContent(i+1, sigv);
      h_sigv_data->SetBinError(i+1, sig2v);

      if (NUCLEAR) std::cout << __LINE__ << std::endl;
      float pta = hp_avgpt_data->GetBinContent(i+1);
      float sigpta = hp_trigger_ptv_pta_data->GetBinWidth(i+1)/2.;
      float cosdphi = hp_cosdphi_pta_data->GetBinContent(i+1);
      float sigcosdphi = hp_cosdphi_pta_data->GetBinError(i+1);
      float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
      float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
      h_jer_data->SetBinContent(i+1, res);
      h_jer_data->SetBinError(i+1, sigres);
      if (sigp <= sigv) continue;
      std::cout << "Data " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
    }

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  for (int iem = 0; iem < 5; iem++)
    {
      for (int i = 0; i < h_em_jer_sim[iem]->GetNbinsX(); i++)
	{
	  h_em_trigger_ptv_pta_sim[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  h_em_trigger_ptp_pta_sim[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);

	  h_em_sigv_sim[iem]->SetBinContent(i+1, sigv);
	  h_em_sigv_sim[iem]->SetBinError(i+1, sig2v);
 
	  h_em_sigp_sim[iem]->SetBinContent(i+1, sigp);
	  h_em_sigp_sim[iem]->SetBinError(i+1, sig2p);
	  float pta = hp_em_avgpt_sim[iem]->GetBinContent(i+1);
	  float sigpta = hp_em_trigger_ptv_pta_sim[iem]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_em_cosdphi_pta_sim[iem]->GetBinContent(i+1);
	  float sigcosdphi = hp_em_cosdphi_pta_sim[iem]->GetBinError(i+1);
	  float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + TMath::Power(sigcosdphi/cosdphi, 2)/4.));
	  h_em_jer_sim[iem]->SetBinContent(i+1, res);
	  h_em_jer_sim[iem]->SetBinError(i+1, sigres);
	  //std::cout << "Sim " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
	}


      for (int i = 0; i < h_em_jer_data[iem]->GetNbinsX(); i++)
	{
	  h_em_trigger_ptv_pta_data[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigv = fg->GetParameter(2);
	  float sig2v = fg->GetParError(2);
	  h_em_trigger_ptp_pta_data[iem]->FitSlicesY(fg, i+1, i+1);
	  float sigp = fg->GetParameter(2);
	  float sig2p = fg->GetParError(2);
	  h_em_sigv_data[iem]->SetBinContent(i+1, sigv);
	  h_em_sigv_data[iem]->SetBinError(i+1, sig2v);
 
	  h_em_sigp_data[iem]->SetBinContent(i+1, sigp);
	  h_em_sigp_data[iem]->SetBinError(i+1, sig2p);

	  float pta = hp_em_avgpt_data[iem]->GetBinContent(i+1);
	  float sigpta = hp_em_trigger_ptv_pta_data[iem]->GetBinWidth(i+1)/2.;
	  float cosdphi = hp_em_cosdphi_pta_data[iem]->GetBinContent(i+1);
	  float sigcosdphi = hp_em_cosdphi_pta_data[iem]->GetBinError(i+1);
	  float res = sqrt(sigp*sigp - sigv*sigv)/(pta*sqrt(2*cosdphi));
	  float sigres = res*sqrt(TMath::Power(sig2p/sigp, 2) + TMath::Power(sig2v/sigv, 2) + 2*(TMath::Power(sigpta/pta, 2) + 0.25*TMath::Power(sigcosdphi/cosdphi, 2)));
	  h_em_jer_data[iem]->SetBinContent(i+1, res);
	  h_em_jer_data[iem]->SetBinError(i+1, sigres);
	  //	  std::cout << "Data " << i << " : " << pta << "/" << sigp << " - " << sigv << " -- " << cosdphi << " --> " << res << " +- " << sigres << std::endl;
	}

    }

  TFile *fjer = new TFile(Form("r%02d/jes/jes_fits_r%02d_1JES.root", conesize, conesize),"r");
  TF1 *f_sim_calib = (TF1*) fjer->Get("fjer");

  TF1 *f = new TF1("f", "sqrt(TMath::Power([2]/x, 2) + TMath::Power([0]/sqrt(x), 2) + TMath::Power([1],2))", low_pt, high_pt);
  f->SetParameters(f_sim_calib->GetParameter(0), f_sim_calib->GetParameter(1), f_sim_calib->GetParameter(2));
  /* f->FixParameter(0, f_sim_calib->GetParameter(0)); */
  /* f->SetParLimits(1, f_sim_calib->GetParameter(1) - 0.02, f_sim_calib->GetParameter(1) + 0.02); */
  /* f->SetParLimits(2, f_sim_calib->GetParameter(2) - 0.4, f_sim_calib->GetParameter(2) + 1); */

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TCanvas *c5 = new TCanvas("c5","c5", 680, 900);
  dlutility::ratioPanelCanvas(c5, 0.4);

  c5->cd(1);
  gPad->SetBottomMargin(0.14);
  dlutility::SetLineAtt(h_jer_sim, color_sim, 1, 1);
  dlutility::SetLineAtt(h_jer_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_jer_sim, color_sim, 1, 8);
  dlutility::SetMarkerAtt(h_jer_data, color_data, 1, 8);
  dlutility::SetFont(h_jer_data, 42, 0.05, 0.06, 0.05, 0.05);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_data->GetYaxis()->SetTitleOffset(1.2);
  h_jer_data->SetMinimum(0);
  h_jer_data->SetMaximum(0.4);
  h_jer_data->Fit("f","R+");
  TF1 *f_data = (TF1*) f->Clone();
  f_data->SetName("f_data");
  f_data->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
  f->SetParameters(0.6313, 0.09508, 2.16);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->Fit("f","R+");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TF1 *f_sim = (TF1*) f->Clone();
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  f_sim->SetName("f_sim");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  f_sim->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  dlutility::SetLineAtt(f_sim_calib, kBlack, 2, 1);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  
  h_jer_data->GetFunction("f")->SetLineColor(kBlue);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->GetFunction("f")->SetLineColor(kRed);
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_data->Draw("p");
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  h_jer_sim->Draw("same p");
  f_sim_calib->Draw("same");
  dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.55, 0.7, 0, kBlack, 0.05);
  TLegend *leg1 = new TLegend(0.20, 0.204, 0.455, 0.354);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.035);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_jer_data,Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data->GetFunction("f")->GetParameter(0), h_jer_data->GetFunction("f")->GetParameter(2), h_jer_data->GetFunction("f")->GetParameter(1)));
  leg1->AddEntry(h_jer_sim,Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim->GetFunction("f")->GetParameter(0), h_jer_sim->GetFunction("f")->GetParameter(2),h_jer_sim->GetFunction("f")->GetParameter(1)));//, h_jer_sim->GetFunction("f")->GetParameter(2)));
  leg1->Draw("same");  
  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  c5->cd(2);
  gPad->SetTopMargin(0.05);

  if (NUCLEAR) std::cout << __LINE__ << std::endl;
  TF1 *fr_compare = new TF1("fr_compare", "(f_data - f_sim)/f_sim", low_pt, high_pt);
  fr_compare->SetLineColor(kRed);
  fr_compare->SetParameters(f_data->GetParameter(0), f_data->GetParameter(1), f_data->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2));

  TF1 *fr_diff = new TF1("fr_diff", "sqrt(fabs(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2)))", low_pt, high_pt);
  fr_diff->SetLineColor(kRed);
  fr_diff->SetParameters(f_data->GetParameter(0), f_data->GetParameter(1), f_data->GetParameter(2), f_sim->GetParameter(0), f_sim->GetParameter(1), f_sim->GetParameter(2));

  TH1D *h_compare = (TH1D*) h_jer_data->Clone();
  h_compare->Reset();

  TH1D *h_dcompare = (TH1D*) h_jer_data->Clone();
  h_dcompare->Reset();
  h_dcompare->SetName("h_dcompare");

  for (int i = 0; i < h_compare->GetNbinsX(); i++)
    {
      float va = h_jer_data->GetBinContent(i+1);
      float ea = h_jer_data->GetBinError(i+1);
      float vb = h_jer_sim->GetBinContent(i+1);
      float eb = h_jer_sim->GetBinError(i+1);
      float rat = (va - vb)/vb;
            
      h_compare->SetBinContent(i+1, rat);
      float err = rat*sqrt(((va*va + vb*vb)*ea*ea + 2*(vb*vb - va*vb)*eb*eb)/(vb*vb*(va - vb)*(va - vb)));

      h_compare->SetBinError(i+1, err);
      
      if (va > vb)
	{
	  float sv = sqrt(va*va - vb*vb);
	  float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	  h_dcompare->SetBinContent(i+1, sv);
	  h_dcompare->SetBinError(i+1, se);
	}
    }
  dlutility::SetLineAtt(h_compare, kRed, 1, 1);
  dlutility::SetMarkerAtt(h_compare, kRed, 1, 8);
  dlutility::SetLineAtt(h_dcompare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_dcompare, kBlack, 1, 8);

  dlutility::SetLineAtt(fr_diff, kBlack, 2, 1);

  dlutility::SetFont(h_compare, 42, 0.07, 0.077, 0.07, 0.07);
  h_compare->SetTitle("; #LT#it{p}_{T}#GT [ GeV ] ; (Data - MC)/MC");
  h_compare->GetYaxis()->SetTitleOffset(0.8);
  h_compare->SetMinimum(-0.3);
  h_compare->SetMaximum(0.6);
  h_compare->Draw("p E0");
  h_dcompare->Draw("p E0 same");
  fr_compare->Draw("same");
  fr_diff->Draw("same");
  TLine *linep = new TLine(20, 0.2, 40, 0.2);
  dlutility::SetLineAtt(linep, kBlack, 3, 4);
  linep->Draw("same");
  TLine *linen = new TLine(20, -0.2, 40, -0.2);
  dlutility::SetLineAtt(linen, kBlack, 3, 4);
  linen->Draw("same");
  TLine *lin0 = new TLine(20, 0, 40, 0);
  dlutility::SetLineAtt(lin0, kBlack, 2, 1);
  lin0->Draw("same");
  TLegend *legl = new TLegend(0.2, 0.78, 0.5, 0.89);
  legl->SetLineWidth(0);
  legl->SetTextFont(42);
  legl->SetTextSize(0.05);
  legl->AddEntry(fr_compare, "Fractional Difference");
  legl->AddEntry(fr_diff, "Quadrature Difference");
  legl->Draw("same");
  c5->Print(Form("r%02d/jer//datasim_jer_r%02d_%d_closure_%s.png", conesize , conesize, closure, gen[index].c_str()));
  c5->Print(Form("r%02d/jer//datasim_jer_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));

  // JER with third jet cuts
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);
      dlutility::SetLineAtt(h_jer_sim_3jet[i3], color_sim, 1, 1);
      dlutility::SetLineAtt(h_jer_data_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_3jet[i3], color_sim, 1, 8);
      dlutility::SetMarkerAtt(h_jer_data_3jet[i3], color_data, 1, 8);
      dlutility::SetFont(h_jer_data_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);
      h_jer_data_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_3jet[i3]->SetMinimum(0);
      h_jer_data_3jet[i3]->SetMaximum(0.4);
      h_jer_data_3jet[i3]->Fit("f","R0");
      TF1 *f_data3 = (TF1*) f->Clone();
      f_data3->SetName("f_data3");
      f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f->SetParameters(0.6313, 0.09508, 2.16);
      h_jer_sim_3jet[i3]->Fit("f","R0");
      TF1 *f_sim3 = (TF1*) f->Clone();
      f_sim3->SetName("f_sim3");
      f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      h_jer_data_3jet[i3]->GetFunction("f")->SetLineColor(kBlue);
      h_jer_sim_3jet[i3]->GetFunction("f")->SetLineColor(kRed);
      h_jer_data_3jet[i3]->Draw("p");
      h_jer_sim_3jet[i3]->Draw("same p");
      f_sim_calib->Draw("same");
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.65, 0, kBlack, 0.05);
      TLegend *leg2 = new TLegend(0.20, 0.204, 0.455, 0.354);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.035);
      leg2->SetTextFont(42);
      leg2->AddEntry(h_jer_data_3jet[i3],Form("Data: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(2), h_jer_data_3jet[i3]->GetFunction("f")->GetParameter(1)));
      leg2->AddEntry(h_jer_sim_3jet[i3],Form("Sim: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(0), h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2),h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(1)));//, h_jer_sim_3jet[i3]->GetFunction("f")->GetParameter(2)));
      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", low_pt, high_pt);
      fr_compare3->SetLineColor(kRed);
      fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", low_pt, high_pt);
      fr_diff3->SetLineColor(kRed);
      fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TH1D *h_compare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");

      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_dcompare3->SetBinContent(i+1, sv);
	      h_dcompare3->SetBinError(i+1, se);
	    }
	}
      dlutility::SetLineAtt(h_compare3, kRed, 1, 1);
      dlutility::SetMarkerAtt(h_compare3, kRed, 1, 8);
      dlutility::SetLineAtt(h_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_dcompare3, kBlack, 1, 8);

      dlutility::SetLineAtt(fr_diff3, kBlack, 2, 1);

      dlutility::SetFont(h_dcompare3, 42, 0.07, 0.077, 0.07, 0.07);
      h_dcompare3->SetTitle("; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference");
      h_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      h_dcompare3->SetMinimum(0);
      h_dcompare3->SetMaximum(0.18);
      h_dcompare3->Draw("p E0");
      fr_diff3->Draw("same");

      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_3jet%d_%d_closure_%s.png", conesize , conesize, i3, closure, gen[index].c_str()));
      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_3jet%d_%d_closure_%s.pdf", conesize , conesize, i3, closure, gen[index].c_str()));
    }
  c5 = new TCanvas("c5","c5", 680, 800);
  dlutility::ratioPanelCanvas(c5, 0.3);

  // JER with third jet cuts
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);

      dlutility::SetLineAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 8);

      dlutility::SetLineAtt(h_jer_truth_im_softcorr_3jet[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_jer_truth_im_softcorr_3jet[i3], color_truth, 1, 8);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 8);
      dlutility::SetFont(h_jer_data_im_softcorr_part_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);

      dlutility::SetLineAtt(h_jer_sim_im_softcorr_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_3jet[i3], color_sim, 1, 24);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_3jet[i3], color_data, 1, 24);

      dlutility::SetLineAtt(h_jer_sim_im_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_3jet[i3], color_sim, 1, 25);

      dlutility::SetLineAtt(h_jer_truth_im_3jet[i3], color_truth, 1, 1);
      dlutility::SetMarkerAtt(h_jer_truth_im_3jet[i3], color_truth, 1, 25);

      dlutility::SetLineAtt(h_jer_data_im_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_3jet[i3], color_data, 1, 25);


      h_jer_data_im_softcorr_part_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMinimum(0);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMaximum(0.5);

      h_jer_data_im_softcorr_part_3jet[i3]->Fit("f","R0");
      
      TF1 *f_data3 = (TF1*) f->Clone();
      f_data3->SetName("f_data3");
      f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f->SetParameters(0.6313, 0.09508, 2.16);
      
      h_jer_sim_im_softcorr_part_3jet[i3]->Fit("f","R0");
      
      TF1 *f_sim3 = (TF1*) f->Clone();
      f_sim3->SetName("f_sim3");
      f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));

      h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kBlue);
      h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kRed);
      h_jer_data_im_softcorr_part_3jet[i3]->Draw("p");
      h_jer_sim_im_softcorr_part_3jet[i3]->Draw("same p");
      h_jer_truth_im_softcorr_3jet[i3]->Draw("same p");
      h_jer_data_im_softcorr_3jet[i3]->Draw("same p"); 
      h_jer_sim_im_softcorr_3jet[i3]->Draw("same p"); 

      /* h_jer_truth_im_3jet[i3]->Draw("same p"); */
      //h_jer_data_im_3jet[i3]->Draw("same p");
      //h_jer_sim_im_3jet[i3]->Draw("same p");
      
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.65, 0, kBlack, 0.05);
      TLegend *leg2 = new TLegend(0.20, 0.55, 0.5, 0.8);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.035);
      leg2->SetTextFont(42);
      leg2->AddEntry(h_jer_truth_im_softcorr_3jet[i3],"Truth after LE");

      leg2->AddEntry(h_jer_data_im_softcorr_3jet[i3],"Data After LE");
      leg2->AddEntry(h_jer_sim_im_softcorr_3jet[i3],"Sim After LE");
      leg2->AddEntry(h_jer_data_im_softcorr_part_3jet[i3],"Data After LE and TS");
      leg2->AddEntry(h_jer_sim_im_softcorr_part_3jet[i3],"Sim After LE and TS");
      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", low_pt, high_pt);
      fr_compare3->SetLineColor(kRed);
      fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", low_pt, high_pt);
      fr_diff3->SetLineColor(kRed);
      fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2));

      TH1D *h_compare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");

      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_im_softcorr_part_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);
      
	  if (va > vb)
	    {
	      float sv = sqrt(va*va - vb*vb);
	      float se = sqrt(TMath::Power(va/sv,2)*ea*ea + TMath::Power(vb/sv, 2)*eb*eb);
	      h_dcompare3->SetBinContent(i+1, sv);
	      h_dcompare3->SetBinError(i+1, se);
	    }
	}
      dlutility::SetLineAtt(h_compare3, kRed, 1, 1);
      dlutility::SetMarkerAtt(h_compare3, kRed, 1, 8);
      dlutility::SetLineAtt(h_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_dcompare3, kBlack, 1, 8);

      dlutility::SetLineAtt(fr_diff3, kBlack, 2, 1);

      dlutility::SetFont(h_dcompare3, 42, 0.09);//, 0.09, 0.07, 0.07);
      h_dcompare3->SetTitle("; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference");
      h_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      h_dcompare3->SetMinimum(0);
      h_dcompare3->SetMaximum(0.2);
      //h_compare3->Draw("p E0");
      h_dcompare3->Draw("p E0");
      //fr_compare3->Draw("same");
      fr_diff3->Draw("same");
      TLegend *legl3 = new TLegend(0.2, 0.78, 0.5, 0.89);
      legl3->SetLineWidth(0);
      legl3->SetTextFont(42);
      legl3->SetTextSize(0.05);
      //      legl3->AddEntry(fr_compare3, "Fractional Difference");
      legl3->AddEntry(fr_diff3, "Quadrature Difference");
      //legl3->Draw("same");
      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_im_3jet%d_%d_closure_%s.png", conesize , conesize, i3, closure, gen[index].c_str()));
      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_im_3jet%d_%d_closure_%s.pdf", conesize , conesize, i3, closure, gen[index].c_str()));
    }


  TF1 *f_sim_im[8];
  TF1 *f_data_im[8];
  TF1 *f_sim_bi[8];
  TF1 *f_data_bi[8];

  TF1 *fr_diff_bi[8];
  TF1 *fr_diff_im[8];
  TF1 *fr_diff_line_bi[8];
  TF1 *fr_diff_line_im[8];

  float err_line_bi[8] = {0};
  float err_line_im[8] = {0};
  float mean_line_bi[8] = {0};
  float mean_line_im[8] = {0};

  // JER with third jet cuts by two methods
  Double_t par_errors[6];

  TF1 *flatline = new TF1("flatline", "[0]", 20, 50);
  
  for (int i3 = 0; i3 < 8 ; i3++)
    {
      c5->cd(1);
      gPad->SetBottomMargin(0.14);

      dlutility::SetLineAtt(h_jer_sim_3jet[i3], color_sim + 2, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_3jet[i3], color_sim + 2, 1, 21);

      dlutility::SetLineAtt(h_jer_data_3jet[i3], color_data + 2, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_3jet[i3], color_data + 2, 1, 21);

      dlutility::SetLineAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 1);
      dlutility::SetMarkerAtt(h_jer_sim_im_softcorr_part_3jet[i3], color_sim, 1, 20);

      dlutility::SetLineAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 1);
      dlutility::SetMarkerAtt(h_jer_data_im_softcorr_part_3jet[i3], color_data, 1, 20);

      dlutility::SetFont(h_jer_data_im_softcorr_part_3jet[i3], 42, 0.05, 0.06, 0.05, 0.05);

      h_jer_data_im_softcorr_part_3jet[i3]->GetYaxis()->SetTitleOffset(1.2);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMinimum(0);
      h_jer_data_im_softcorr_part_3jet[i3]->SetMaximum(0.4);

      f->SetParameters(f_sim_calib->GetParameter(0), f_sim_calib->GetParameter(1), f_sim_calib->GetParameter(2));//      f->SetParameters(0.6313, 0.09508, 2.16);


      f->SetParLimits(0, f_sim_calib->GetParameter(0) - 0.1,f_sim_calib->GetParameter(0) + 0.3);
      f->SetParLimits(1, f_sim_calib->GetParameter(1) - 0.1,f_sim_calib->GetParameter(1) + 1);
      f->SetParLimits(2, f_sim_calib->GetParameter(2) - 0.1,f_sim_calib->GetParameter(2) + 1); 
      /* f->FixParameter(1, f_sim_calib->GetParameter(1)); */
      /* f->SetParLimits(2, std::max(0.0, (float) f_sim_calib->GetParameter(2) - 0.05),f_sim_calib->GetParameter(2) + 0.05); */


      TGraphErrors *tg_sim_im = new TGraphErrors(h_jer_sim_im_softcorr_part_3jet[i3]->GetNbinsX());
      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  tg_sim_im->SetPointX(ib, hp_avgpt_sim_3jet[i3]->GetBinContent(ib+1));
	  tg_sim_im->SetPointY(ib, h_jer_sim_im_softcorr_part_3jet[i3]->GetBinContent(ib+1));
	  tg_sim_im->SetPointError(ib, 0,h_jer_sim_im_softcorr_part_3jet[i3]->GetBinError(ib+1));
	}
      tg_sim_im->Fit("f","R0");

      f_sim_im[i3] = (TF1*) f->Clone();
      f_sim_im[i3]->SetName(Form("f_sim_im_3jet_%d", i3));
      f_sim_im[i3]->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f_sim_im[i3]->SetParErrors(f->GetParErrors());
      f_sim_im[i3]->SetRange(8.3, 50);            


      //f->SetParLimits(0, f_sim_im[i3]->GetParameter(0), f_sim_im[i3]->GetParameter(0) + 0.2);
      f->SetParameter(1, f_sim_im[i3]->GetParameter(1));
      
    
      f->SetParameter(0, f_sim_im[i3]->GetParameter(0));
      f->SetParameter(2, f_sim_im[i3]->GetParameter(2));

      f->SetParLimits(0, f_sim_im[i3]->GetParameter(0) - 0.1, f_sim_im[i3]->GetParameter(0) + 0.3);
      f->SetParLimits(1, f_sim_im[i3]->GetParameter(1) - 0.1, f_sim_im[i3]->GetParameter(1) + 1);
      f->SetParLimits(2, f_sim_im[i3]->GetParameter(2) - 0.1, f_sim_im[i3]->GetParameter(2) + 1);
      //f->FixParameter(0, f->GetParameter(0));
      //f->FixParameter(2, f->GetParameter(2));

      TGraphErrors *tg_data_im = new TGraphErrors(h_jer_data_im_softcorr_part_3jet[i3]->GetNbinsX());
      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  tg_data_im->SetPointX(ib, hp_avgpt_data_3jet[i3]->GetBinContent(ib+1));
	  tg_data_im->SetPointY(ib, h_jer_data_im_softcorr_part_3jet[i3]->GetBinContent(ib+1));
	  tg_data_im->SetPointError(ib, 0,h_jer_data_im_softcorr_part_3jet[i3]->GetBinError(ib+1));
	}
      tg_data_im->Fit("f","R0");

      f_data_im[i3] = (TF1*) f->Clone();
      f_data_im[i3]->SetName(Form("f_data_im_3jet_%d", i3));
      f_data_im[i3]->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f_data_im[i3]->SetParErrors(f->GetParErrors());
      f_data_im[i3]->SetRange(8.3, 50);            
      //f->SetParameters(0.6313, 0.09508, 2.16);
      f->SetParameters(f_sim_calib->GetParameter(0), f_sim_calib->GetParameter(1), f_sim_calib->GetParameter(2));
      f->SetParLimits(0, f_sim_calib->GetParameter(0) - 0.1,f_sim_calib->GetParameter(0) + 0.3);
      f->SetParLimits(1, f_sim_calib->GetParameter(1) - 0.1,f_sim_calib->GetParameter(1) + 1);
      f->SetParLimits(2, f_sim_calib->GetParameter(2) - 0.1,f_sim_calib->GetParameter(2) + 1); 

      /* f->FixParameter(0, f_sim_calib->GetParameter(0)); */
      /* f->SetParLimits(1, f_sim_calib->GetParameter(1) - 0.02, f_sim_calib->GetParameter(1) + 0.02); */
      /* f->SetParLimits(2, f_sim_calib->GetParameter(2) - 0.4, f_sim_calib->GetParameter(2) + 1); */
      /* f->SetParLimits(0, std::max(0.0, (float) f_sim_calib->GetParameter(0) - 0.3),f_sim_calib->GetParameter(0) + 0.3); */
      /* f->FixParameter(1, f_sim_calib->GetParameter(1)); */
      /* f->SetParLimits(2, std::max(0.0, (float) f_sim_calib->GetParameter(2) - 0.05),f_sim_calib->GetParameter(2) + 0.05); */

      
      TGraphErrors *tg_sim_bi = new TGraphErrors(h_jer_sim_3jet[i3]->GetNbinsX());
      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  tg_sim_bi->SetPointX(ib, hp_avgpt_sim_3jet[i3]->GetBinContent(ib+1));
	  tg_sim_bi->SetPointY(ib, h_jer_sim_3jet[i3]->GetBinContent(ib+1));
	  tg_sim_bi->SetPointError(ib, 0, h_jer_sim_3jet[i3]->GetBinError(ib+1));
	}

      tg_sim_bi->Fit("f","R0");

      f_sim_bi[i3] = (TF1*) f->Clone();
      f_sim_bi[i3]->SetName(Form("f_sim_bi_3jet_%d", i3));
      f_sim_bi[i3]->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f_sim_bi[i3]->SetParErrors(f->GetParErrors());
      f_sim_bi[i3]->SetRange(8.3, 50);            
      f->SetParameter(1, f_sim_bi[i3]->GetParameter(1));
      f->SetParameter(0, f_sim_bi[i3]->GetParameter(0));
      f->SetParameter(2, f_sim_bi[i3]->GetParameter(2));

      TGraphErrors *tg_data_bi = new TGraphErrors(h_jer_data_3jet[i3]->GetNbinsX());
      for (int ib = 0; ib < nbins_pt; ib++)
	{
	  tg_data_bi->SetPointX(ib, hp_avgpt_data_3jet[i3]->GetBinContent(ib+1));
	  tg_data_bi->SetPointY(ib, h_jer_data_3jet[i3]->GetBinContent(ib+1));
	  tg_data_bi->SetPointError(ib, 0, h_jer_data_3jet[i3]->GetBinError(ib+1));
	}
      tg_data_bi->Fit("f","R0");

      f_data_bi[i3] = (TF1*) f->Clone();
      f_data_bi[i3]->SetName(Form("f_data_bi_3jet_%d", i3));
      f_data_bi[i3]->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
      f_data_bi[i3]->SetParErrors(f->GetParErrors());
      f_data_bi[i3]->SetRange(8.3, 50);
      
      dlutility::SetLineAtt(f_data_im[i3], color_data, 2, 1);
      dlutility::SetLineAtt(f_sim_im[i3], color_sim, 2, 1);
      dlutility::SetLineAtt(f_data_bi[i3], color_data + 2, 2, 1);
      dlutility::SetLineAtt(f_sim_bi[i3], color_sim + 2, 2, 1);

      dlutility::SetLineAtt(tg_sim_bi, color_sim + 2, 1, 1);
      dlutility::SetMarkerAtt(tg_sim_bi, color_sim + 2, 1, 21);

      dlutility::SetLineAtt(tg_sim_im, color_sim, 1, 1);
      dlutility::SetMarkerAtt(tg_sim_im, color_sim, 1, 20);

      dlutility::SetLineAtt(tg_data_bi, color_data + 2, 1, 1);
      dlutility::SetMarkerAtt(tg_data_bi, color_data + 2, 1, 21);

      dlutility::SetLineAtt(tg_data_im, color_data, 1, 1);
      dlutility::SetMarkerAtt(tg_data_im, color_data, 1, 20);

      dlutility::SetFont(tg_data_im, 42, 0.05, 0.06, 0.05, 0.05);
      TH1D *hb1 = new TH1D("hb1", "; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference", 1, 8, 50);


      hb1->SetTitle(";#LT#it{p}_{T}#GT [GeV]; #sigma(#it{p}_{T})/#it{p}_{T}");
      hb1->GetYaxis()->SetTitleOffset(1.2);
      hb1->SetMinimum(0);
      hb1->SetMaximum(0.4);
      hb1->Draw();
      //h_jer_data_im_softcorr_part_3jet[i3]->Fit("f","R0");
      /* TF1 *f_data3 = (TF1*) f->Clone(); */
      /* f_data3->SetName("f_data3"); */
      /* f_data3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2)); */
      /* f->SetParameters(0.6313, 0.09508, 2.16); */
      /* h_jer_sim_im_softcorr_part_3jet[i3]->Fit("f","R"); */
      /* TF1 *f_sim3 = (TF1*) f->Clone(); */
      /* f_sim3->SetName("f_sim3"); */
      /* f_sim3->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2)); */
      /* h_jer_data_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kBlue); */
      /* h_jer_sim_im_softcorr_part_3jet[i3]->GetFunction("f")->SetLineColor(kRed); */

      tg_sim_bi->Draw("same p");
      tg_sim_im->Draw("same p");
      tg_data_bi->Draw("same p");
      tg_data_im->Draw("same p");
      

      f_sim_calib->Draw("same");

      f_sim_im[i3]->Draw("same ");
      f_data_im[i3]->Draw("same ");
      f_sim_bi[i3]->Draw("same ");
      f_data_bi[i3]->Draw("same ");

      TLine *linesample = new TLine(23, 0, 23, 0.4);
      linesample->SetLineColor(kBlack);
      linesample->SetLineStyle(2);
      linesample->SetLineWidth(2);
      linesample->Draw("same");


      fr_diff_im[i3] = new TF1(Form("fr_diff_im_3jet_%d", i3), "sqrt(fabs(TMath::Power(f, 2) - TMath::Power(f, 2)))", low_pt, high_pt);
      fr_diff_im[i3]->SetLineColor(kRed);
      fr_diff_im[i3]->SetParameters(f_data_im[i3]->GetParameter(0), f_data_im[i3]->GetParameter(1), f_data_im[i3]->GetParameter(2), f_sim_im[i3]->GetParameter(0), f_sim_im[i3]->GetParameter(1), f_sim_im[i3]->GetParameter(2));
      fr_diff_im[i3]->SetRange(8.3, 50);
      
      fr_diff_bi[i3] = new TF1(Form("fr_diff_bi_3jet_%d", i3), "sqrt(fabs(TMath::Power(f, 2) - TMath::Power(f, 2)))", low_pt, high_pt);
      fr_diff_bi[i3]->SetLineColor(kRed);
      fr_diff_bi[i3]->SetParameters(f_data_bi[i3]->GetParameter(0), f_data_bi[i3]->GetParameter(1), f_data_bi[i3]->GetParameter(2), f_sim_bi[i3]->GetParameter(0), f_sim_bi[i3]->GetParameter(1), f_sim_bi[i3]->GetParameter(2));
      fr_diff_bi[i3]->SetRange(8.3, 50);

      dlutility::SetLineAtt(fr_diff_bi[i3], color_data + 2, 2, 1);
      dlutility::SetLineAtt(fr_diff_im[i3], color_data, 2, 1);

      
      dlutility::DrawSPHENIXppInternalsize(0.55, 0.8, 0.05);
      dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.55, 0.7, 0, kBlack, 0.05);
      dlutility::drawText("Isolation #Delta#it{R} < 1.0", 0.55, 0.65, 0, kBlack, 0.05);
      dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.55, 0.6, 0, kBlack, 0.05);
      dlutility::drawText("Min. Bias. #leftarrow  ",0.445, 0.37, 1);
      dlutility::drawText("   #rightarrow Jet Trigger",0.445, 0.37, 0);
      TLegend *leg2 = new TLegend(0.20, 0.184, 0.455, 0.334);
      leg2->SetLineWidth(0);
      leg2->SetTextSize(0.03);
      leg2->SetTextFont(42);
      leg2->AddEntry(f_sim_calib,"MC JER Analysis");
      leg2->AddEntry(tg_data_im,Form("Data Im.: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", f_data_im[i3]->GetParameter(0), f_data_im[i3]->GetParameter(2), f_data_im[i3]->GetParameter(1)));
      leg2->AddEntry(tg_sim_im,Form("%s Im.: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f",  gen[index].c_str(), f_sim_im[i3]->GetParameter(0), f_sim_im[i3]->GetParameter(2), f_sim_im[i3]->GetParameter(1)));
      leg2->AddEntry(tg_data_bi,Form("Data Bi.: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", f_data_bi[i3]->GetParameter(0), f_data_bi[i3]->GetParameter(2), f_data_bi[i3]->GetParameter(1)));
      leg2->AddEntry(tg_sim_bi,Form("%s Bi.: %2.2f / #sqrt{E} #oplus %2.2f / E #oplus %2.4f", gen[index].c_str(), f_sim_bi[i3]->GetParameter(0), f_sim_bi[i3]->GetParameter(2), f_sim_bi[i3]->GetParameter(1)));

      leg2->Draw("same");  

      c5->cd(2);
      gPad->SetTopMargin(0.05);

  
      /* TF1 *fr_compare3 = new TF1("fr_compare3", "(f_data - f_sim)/f_sim", low_pt, high_pt); */
      /* fr_compare3->SetLineColor(kRed); */
      /* fr_compare3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2)); */

      /* TF1 *fr_diff3 = new TF1("fr_diff", "sqrt(TMath::Power(f_data, 2) - TMath::Power(f_sim, 2))", low_pt, high_pt); */
      /* fr_diff3->SetLineColor(kRed); */
      /* fr_diff3->SetParameters(f_data3->GetParameter(0), f_data3->GetParameter(1), f_data3->GetParameter(2), f_sim3->GetParameter(0), f_sim3->GetParameter(1), f_sim3->GetParameter(2)); */

      TH1D *h_im_compare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_im_compare3->Reset();
      h_im_compare3->SetName("h_im_compare3");
      TH1D *h_im_dcompare3 = (TH1D*) h_jer_data_im_softcorr_part_3jet[i3]->Clone();
      h_im_dcompare3->Reset();
      h_im_dcompare3->SetName("h_im_dcompare3");
      TH1D *h_compare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_compare3->Reset();
      h_compare3->SetName("h_compare3");
      TH1D *h_dcompare3 = (TH1D*) h_jer_data_3jet[i3]->Clone();
      h_dcompare3->Reset();
      h_dcompare3->SetName("h_dcompare3");      

      TGraphAsymmErrors *g_dcompare3 = new TGraphAsymmErrors(nbins_pt);
      TGraphAsymmErrors *g_im_dcompare3 = new TGraphAsymmErrors(nbins_pt);

      double min_ax = 0;
      double max_ax = 0;

      
      for (int i = 0; i < h_compare3->GetNbinsX(); i++)
	{
	  float va = h_jer_data_3jet[i3]->GetBinContent(i+1);
	  float ea = h_jer_data_3jet[i3]->GetBinError(i+1);
	  float vb = h_jer_sim_3jet[i3]->GetBinContent(i+1);
	  float eb = h_jer_sim_3jet[i3]->GetBinError(i+1);

	  float rat = (va - vb)/vb;
            
	  h_compare3->SetBinContent(i+1, rat);
	  float err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_compare3->SetBinError(i+1, err);

	  int a_or_b = 1;
	  float nom_high = va;
	  float nom_low = vb;
	  if (va < vb)
	    {
	      a_or_b = -1;
	      float temp = nom_high;
	      nom_high = nom_low;
	      nom_low = temp;
	    }

	  float nom_diff = a_or_b * sqrt(nom_high * nom_high - nom_low * nom_low);

	  float min_err = 0;
	  float max_err = 0;

	  for (int i4 = 0; i4 < 4; i4++)
	    {
	      int a_or_b_err = 1;
	      float err_high = va -  ea + 2 * (i4/2) * ea;
	      float err_low = vb -  eb + 2 * (i4%2) * eb;
	      if (err_high < err_low)
		{
		  a_or_b_err = -1;
		  float temp = err_high;
		  err_high = err_low;
		  err_low = temp;
		}

	      float err_diff = a_or_b_err * sqrt(err_high * err_high - err_low * err_low);

	      float derr = err_diff - nom_diff;
	     
	      if (derr > 0 && derr > max_err)
		{
		  max_err = derr;
		}
	      if (derr < 0 && derr < min_err)
		{
		  min_err = derr;
		}

	      std::cout << "NOM_ERR_ERR = " << i4 <<  " - " << i << " : " << nom_diff << " // " << err_diff << " = " << derr << std::endl;
	    }

	  if (nom_diff - fabs(min_err) < min_ax)
	    {
	      min_ax = nom_diff  - fabs(min_err);
	    }
	  if (nom_diff + fabs(max_err) > max_ax)
	    {
	      max_ax = nom_diff  + fabs(max_err);
	    }
	  
	  g_dcompare3->SetPointX(i, hp_avgpt_sim_3jet[i3]->GetBinContent(i+1));
	  g_dcompare3->SetPointY(i, nom_diff);
	  g_dcompare3->SetPointError(i, 0, 0, fabs(min_err), fabs(max_err));

	  va = h_jer_data_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  ea = h_jer_data_im_softcorr_part_3jet[i3]->GetBinError(i+1);
	  vb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinContent(i+1);
	  eb = h_jer_sim_im_softcorr_part_3jet[i3]->GetBinError(i+1);

	  rat = (va - vb)/vb;
            
	  h_im_compare3->SetBinContent(i+1, rat);
	  err = rat*sqrt(TMath::Power(sqrt(ea*ea + eb*eb)/(va - vb) , 2) + TMath::Power(eb/vb , 2));

	  h_im_compare3->SetBinError(i+1, err);

	  a_or_b = 1;
	  nom_high = va;
	  nom_low = vb;

	  if (va < vb)
	    {
	      a_or_b = -1;
	      float temp = nom_high;
	      nom_high = nom_low;
	      nom_low = temp;
	    }

	  nom_diff = a_or_b * sqrt(nom_high * nom_high - nom_low * nom_low);

	  min_err = 0;
	  max_err = 0;

	  for (int i4 = 0; i4 < 4; i4++)
	    {
	      int a_or_b_err = 1;
	      float err_high = va -  ea + 2 * (i4/2) * ea;
	      float err_low = vb -  eb + 2 * (i4%2) * eb;
	      if (err_high < err_low)
		{
		  a_or_b_err = -1;
		  float temp = err_high;
		  err_high = err_low;
		  err_low = temp;
		}

	      float err_diff = a_or_b_err * sqrt(err_high * err_high - err_low * err_low);

	      float derr = err_diff - nom_diff;

	      if (derr > 0 && derr > max_err)
		{
		  max_err = derr;
		}
	      if (derr < 0 && derr < min_err)
		{
		  min_err = derr;
		}

	    }

	  if (nom_diff - fabs(min_err) < min_ax)
	    {
	      min_ax = nom_diff  - fabs(min_err);
	    }
	  if (nom_diff + fabs(max_err) > max_ax)
	    {
	      max_ax = nom_diff  + fabs(max_err);
	    }
	  
	  g_im_dcompare3->SetPointX(i, hp_avgpt_sim_3jet[i3]->GetBinContent(i+1));
	  g_im_dcompare3->SetPointY(i, nom_diff);
	  g_im_dcompare3->SetPointError(i, 0, 0, fabs(min_err), fabs(max_err));

	}

      dlutility::SetLineAtt(g_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(g_dcompare3, kBlack, 1, 21);

      dlutility::SetLineAtt(g_im_dcompare3, kBlack, 1, 1);
      dlutility::SetMarkerAtt(g_im_dcompare3, kBlack, 1, 20);

      dlutility::SetFont(g_dcompare3, 42, 0.12, 0.08, 0.08, 0.08);
      g_dcompare3->SetTitle("; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference");
      g_dcompare3->GetYaxis()->SetTitleOffset(0.8);
      
      flatline->SetParameter(0, 0.1);
      g_dcompare3->Fit("flatline","NDRQ0","goff", 20, 50);
      fr_diff_line_bi[i3] = (TF1*) flatline->Clone();
      fr_diff_line_bi[i3]->SetName(Form("fr_diff_line_bi_%d", i3));
      fr_diff_line_bi[i3]->SetParameter(0, flatline->GetParameter(0));
      fr_diff_line_bi[i3]->SetRange(8, 50);
      float mean_bi = flatline->GetParameter(0);
      float chi2_bi = 0;
      float sigma_bi = 0;
      float wmean_bi = 0;
      float wsigma_bi = 0;
      float wsum_bi = 0;
      float count_bi;
      for (int i = 0; i < nbins_pt - 1; i++)
	{

	  if (g_dcompare3->GetPointY(i) < 0) continue;
	  count_bi += 1;
	  wmean_bi += g_dcompare3->GetPointY(i) / g_dcompare3->GetErrorY(i);
	  wsum_bi += 1/(g_dcompare3->GetErrorY(i));
	  std::cout << " POINT " << i << " : " << g_dcompare3->GetPointY(i) << " +- " << g_dcompare3->GetErrorY(i) <<std::endl;
	}

      wmean_bi /= wsum_bi;
      wsum_bi = 0;      

      for (int i = 0; i < nbins_pt - 1; i++)
	{
	  if (g_dcompare3->GetPointY(i) < 0) continue;
	  wsum_bi += 1/(g_dcompare3->GetErrorY(i));
	  wsigma_bi += (g_dcompare3->GetPointY(i) - wmean_bi)*(g_dcompare3->GetPointY(i) - wmean_bi)/(g_dcompare3->GetErrorY(i));
	}

      // error on the wesighted mean
      wsigma_bi /= wsum_bi;
      wsigma_bi *= ((float)(count_bi))/((float)(count_bi - 1));
      wsigma_bi = sqrt(wsigma_bi);

      err_line_bi[i3] = wsigma_bi;
      mean_line_bi[i3] = wmean_bi;

      std::cout << "CALC : " << mean_line_bi[i3] << " / " << err_line_bi[i3] << std::endl;
      for (int i = 0; i < nbins_pt; i++)
	{
	  float d = g_dcompare3->GetPointY(i) - mean_bi;
	  float w = (d > 0? g_dcompare3->GetErrorYlow(i) : g_dcompare3->GetErrorYhigh(i));
	  w = 1./(w*w);
	  chi2_bi += w * d * d;
	  sigma_bi += w;
	}
      chi2_bi /= (float) (nbins_pt - 1);
      sigma_bi = 1./sqrt(sigma_bi);
      
      float sigma_int_bi = sigma_bi * sqrt(chi2_bi - 1);

      std::cout << i3 << " : SIGMA INT BI = " << chi2_bi <<" = "  << sigma_int_bi << " / " << sigma_bi << std::endl;
      flatline->SetParameter(0, 0.1);
      g_im_dcompare3->Fit("flatline","NDRQ0","goff", 20, 50);

      fr_diff_line_im[i3] = (TF1*) flatline->Clone();
      fr_diff_line_im[i3]->SetName(Form("fr_diff_line_im_%d", i3));
      fr_diff_line_im[i3]->SetParameter(0, flatline->GetParameter(0));
      fr_diff_line_im[i3]->SetRange(8, 50);
      
      float mean_im = flatline->GetParameter(0);
      float chi2_im = 0;
      float sigma_im = 0;
      float wmean_im = 0;
      float wsigma_im = 0;
      float wsum_im = 0;
      float count_im;
      for (int i = 0; i < nbins_pt - 1; i++)
	{
	  if (g_im_dcompare3->GetPointY(i) < 0) continue;
	  count_im += 1;
	  wmean_im += g_im_dcompare3->GetPointY(i) / g_im_dcompare3->GetErrorY(i);
	  wsum_im += 1/(g_im_dcompare3->GetErrorY(i));
	  std::cout << " POINT " << i << " : " << g_im_dcompare3->GetPointY(i) << " +- " << g_im_dcompare3->GetErrorY(i) <<std::endl;
	}

      wmean_im /= wsum_im;
      wsum_im = 0;      

      for (int i = 0; i < nbins_pt - 1; i++)
	{
	  if (g_im_dcompare3->GetPointY(i) < 0) continue;
	  wsum_im += 1/(g_im_dcompare3->GetErrorY(i));
	  wsigma_im += (g_im_dcompare3->GetPointY(i) - wmean_im)*(g_im_dcompare3->GetPointY(i) - wmean_im)/(g_im_dcompare3->GetErrorY(i));
	}
      
      wsigma_im /= wsum_im;
      wsigma_im *= ((float)(count_im))/((float)(count_im - 1));
      wsigma_im = sqrt(wsigma_im);

      err_line_im[i3] = wsigma_im;
      mean_line_im[i3] = wmean_im;

      std::cout << "CALC : " << mean_line_im[i3] << " / " << err_line_im[i3] << std::endl;
      for (int i = 0; i < nbins_pt; i++)
	{
	  float d = g_im_dcompare3->GetPointY(i) - mean_im;
	  float w = (d > 0? g_im_dcompare3->GetErrorYlow(i) : g_im_dcompare3->GetErrorYhigh(i));
	  w = 1./(w*w);
	  chi2_im += w * d * d;
	  sigma_im += w;
	}
      chi2_im /= (float) (nbins_pt - 1);
      sigma_im = 1./sqrt(sigma_im);
      
      float sigma_int_im = sigma_im * sqrt(chi2_im - 1);

      std::cout << i3 << " : SIGMA INT IM = " << chi2_im << " = " << sigma_int_im << " / " << sigma_im << std::endl;

      dlutility::SetLineAtt(fr_diff_line_bi[i3], color_data + 2, 2, 3);
      dlutility::SetLineAtt(fr_diff_line_im[i3], color_data, 2, 4);

      TH1D *hb = new TH1D("hb", "; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference", 1, 8, 50);
      dlutility::SetFont(hb, 42, 0.12, 0.08, 0.08, 0.08);
      hb->SetTitle("; #LT#it{p}_{T}#GT [ GeV ] ; Quadrature Difference");
      hb->GetYaxis()->SetTitleOffset(0.8);
      hb->SetMinimum(0);
      hb->SetMaximum(0.2);
      hb->Draw();
      g_im_dcompare3->Draw("same p E1");      
      g_dcompare3->Draw("p E1 same");

      fr_diff_line_bi[i3]->Draw("same");
      fr_diff_line_im[i3]->Draw("same");

      //fr_diff_im[i3]->Draw("same");
      //fr_diff_bi[i3]->Draw("same");

      TLegend *legl3 = new TLegend(0.25, 0.7, 0.45, 0.9);
      legl3->SetLineWidth(0);
      legl3->SetTextFont(42);
      legl3->SetTextSize(0.07);
      legl3->AddEntry(g_dcompare3, "Bisector");
      legl3->AddEntry(g_im_dcompare3, "Imbalance");
      legl3->Draw("same");
      TLegend *legl4 = new TLegend(0.45, 0.7, 0.65, 0.9);
      legl4->SetLineWidth(0);
      legl4->SetTextFont(42);
      legl4->SetTextSize(0.07);
      legl4->AddEntry(fr_diff_line_bi[i3]," Constant Fit: Bisector");
      legl4->AddEntry(fr_diff_line_im[i3]," Constant Fit: Imbalance");
      legl4->Draw("same");

      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_im_v_bi_3jet%d_%d_closure_%s.png", conesize , conesize, i3, closure, gen[index].c_str()));
      c5->Print(Form("r%02d/jer//datasim_jer_r%02d_im_v_bi_3jet%d_%d_closure_%s.pdf", conesize , conesize, i3, closure, gen[index].c_str()));
    }


  TFile *fint = new TFile("/Users/daniel/pythia8/plotting/truth_sig.root", "r");
  TH1D *h_sigp_truth = (TH1D*) fint->Get("h_sigp");
  TH1D *h_sigv_truth = (TH1D*) fint->Get("h_sigv");
  TFile *fintf = new TFile("/Users/daniel/pythia8/plotting/truth_sig_noFSR.root", "r");
  TH1D *h_sigp_truthfsr = (TH1D*) fintf->Get("h_sigp");
  TH1D *h_sigv_truthfsr = (TH1D*) fintf->Get("h_sigv");

  h_sigp_truth->SetTitle(";#LT#it{p}_{T}#GT [GeV]; #sigma(#it{p}_{T,#psi}) or #sigma(#it{p}_{T,#lambda}) [GeV]");

  TCanvas *c10 = new TCanvas("c10","c10", 500, 500);

  dlutility::SetLineAtt(h_sigp_data, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_data, kGreen + 2, 1, 20);
  dlutility::SetLineAtt(h_sigv_data, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_data, kMagenta + 2, 1, 20);
  dlutility::SetLineAtt(h_sigp_sim, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_sim, kGreen + 2, 1, 24);
  dlutility::SetLineAtt(h_sigv_sim, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_sim, kMagenta + 2, 1, 24);
  dlutility::SetLineAtt(h_sigp_truth, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_truth, kGreen + 2, 1, 25);
  dlutility::SetLineAtt(h_sigv_truth, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_truth, kMagenta + 2, 1, 25);

  dlutility::SetLineAtt(h_sigp_truthfsr, kGreen + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigp_truthfsr, kGreen + 2, 1, 21);
  dlutility::SetLineAtt(h_sigv_truthfsr, kMagenta + 2, 1, 1);
  dlutility::SetMarkerAtt(h_sigv_truthfsr, kMagenta + 2, 1, 21);

  h_sigp_truth->GetXaxis()->SetRangeUser(10, 50);
  h_sigp_truth->SetMaximum(10);
  h_sigp_truth->SetMinimum(0);
  dlutility::SetFont(h_sigp_truth, 42, 0.05, 0.05, 0.05, 0.05);
  h_sigp_truth->Draw("p");
  h_sigv_truth->Draw("same p");
  dlutility::DrawSPHENIXppInternalsize(0.22, 0.85, 0.05);
  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.22, 0.70);
  dlutility::drawText("#it{p}_{T} > 10 GeV", 0.22, 0.65);
  TLegend *leg = new TLegend(0.67, 0.65, 0.83, 0.877);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_sigp_truth, "#sigma(#it{p}_{T, #psi}^{Truth})");
  leg->AddEntry(h_sigv_truth, "#sigma(#it{p}_{T, #lambda}^{Truth})");

  leg->Draw("same");

  c10->Print(Form("r%02d/jer//sigpsigv_truth_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));
  c10->Print(Form("r%02d/jer//sigpsigv_truth_r%02d_%d_closure_%s.png", conesize , conesize, closure, gen[index].c_str()));



  h_sigp_sim->SetMaximum(20);
  h_sigp_sim->SetMinimum(0);
  //ydlutility::SetFont(h_sigp_truth, 42, 0.05);
  h_sigp_sim->Draw("p");
  h_sigv_sim->Draw("same p");
  h_sigp_truth->Draw("same p");
  h_sigv_truth->Draw("same p");

  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("#it{p}_{T} > 10 GeV", 0.22, 0.75);
  leg = new TLegend(0.67, 0.65, 0.83, 0.877);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_sigp_sim, "#sigma(#it{p}_{T, #psi}^{Reco})");
  leg->AddEntry(h_sigp_truth, "#sigma(#it{p}_{T, #psi}^{Truth})");
  leg->AddEntry(h_sigv_sim, "#sigma(#it{p}_{T, #lambda}^{Reco})");

  leg->Draw("same");

  c10->Print(Form("r%02d/jer//sigpsigv_sim_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));
  c10->Print(Form("r%02d/jer//sigpsigv_sim_r%02d_%d_closure_%s.png", conesize , conesize, closure, gen[index].c_str()));

  h_sigp_sim->SetMaximum(20);
  h_sigp_sim->SetMinimum(0);
  h_sigp_sim->Draw("p");
  h_sigv_sim->Draw("same p");
  h_sigp_data->Draw("same p");
  h_sigv_data->Draw("same p");
  h_sigp_truth->Draw("same p");
  h_sigv_truth->Draw("same p");

  dlutility::DrawSPHENIXppInternalsize(0.22, 0.85, 0.04);
  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = 0.%d", conesize), 0.22, 0.75);
  dlutility::drawText("#it{p}_{T} > 10 GeV", 0.22, 0.7);
  //  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg = new TLegend(0.55, 0.64, 0.7, 0.89);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_sigp_data, "#sigma(#it{p}_{T, #psi}^{Data})","p");
  leg->AddEntry(h_sigp_sim, "#sigma(#it{p}_{T, #psi}^{Reco})","p");
  leg->AddEntry(h_sigp_truth, "#sigma(#it{p}_{T, #psi}^{Truth})","p");
  leg->Draw("same");
  leg = new TLegend(0.72, 0.64, 0.91, 0.89);
  leg->SetLineWidth(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_sigv_data, "#sigma(#it{p}_{T, #lambda}^{Data})","p");
  leg->AddEntry(h_sigv_sim, "#sigma(#it{p}_{T, #lambda}^{Reco})","p");
  leg->AddEntry(h_sigv_truth, "#sigma(#it{p}_{T, #lambda}^{Truth})","p");

  leg->Draw("same");

  c10->Print(Form("r%02d/jer//sigpsigv_data_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));
  c10->Print(Form("r%02d/jer//sigpsigv_data_r%02d_%d_closure_%s.png", conesize , conesize, closure, gen[index].c_str()));

  // Here is the average EM fraction result

  TCanvas *c_em = new TCanvas("c_em", "c_em", 800, 400);
  int color_em[5] = {kRed, kGreen, kRed, kBlack, kBlack};

  h_jer_data->SetMinimum(0);
  h_jer_data->SetMaximum(0.4);

  h_jer_sim->SetMinimum(0);
  h_jer_sim->SetMaximum(0.4);

  dlutility::SetLineAtt(h_jer_sim, color_data, 1, 1);
  dlutility::SetLineAtt(h_jer_data, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_jer_sim, color_data, 1, 20);
  dlutility::SetMarkerAtt(h_jer_data, color_data, 1, 20);

  dlutility::SetFont(h_jer_sim, 42, 0.04);
  dlutility::SetFont(h_jer_data, 42, 0.04);
  
  for (int iem = 0; iem < nbins_em; iem++)
    {
      dlutility::SetLineAtt(h_em_jer_data[iem], color_em[iem], 1, 1);
      dlutility::SetMarkerAtt(h_em_jer_data[iem], color_em[iem], 1, 20);
      dlutility::SetLineAtt(h_em_jer_sim[iem], color_em[iem], 1, 1);
      dlutility::SetMarkerAtt(h_em_jer_sim[iem], color_em[iem], 1, 20);
    }  

  c_em->Divide(2, 1);
  c_em->cd(1);
  h_jer_sim->Draw("p");
  h_em_jer_sim[0]->Draw("p same");
  h_em_jer_sim[1]->Draw("p same");
  dlutility::DrawSPHENIXppInternal(0.22, 0.85);
  dlutility::drawText("#it{p}_{T} > 10 GeV", 0.22, 0.75);
  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg->SetLineWidth(0);  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_jer_sim, "ALL SIM","p");
  leg->AddEntry(h_em_jer_sim[0], "Sim <E_{T}^{EM}/E_{T}> < 0.5","p");
  leg->AddEntry(h_em_jer_sim[1], "Sim <E_{T}^{EM}/E_{T}> #geq 0.5","p");
  leg->Draw("same");
  c_em->cd(2);
  h_jer_data->Draw("p");
  h_em_jer_data[0]->Draw("p same");
  h_em_jer_data[1]->Draw("p same");
  leg = new TLegend(0.61, 0.65, 0.77, 0.87);
  leg->SetLineWidth(0);  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h_jer_data, "ALL DATA","p");
  leg->AddEntry(h_em_jer_data[0], "Data <E_{T}^{EM}/E_{T}> < 0.5","p");
  leg->AddEntry(h_em_jer_data[1], "Data <E_{T}^{EM}/E_{T}> #geq 0.5","p");
  leg->Draw("same");

  c_em->Print(Form("r%02d/jer//em_fraction_jer_r%02d_%d_closure_%s.png", conesize , conesize, closure, gen[index].c_str()));
  c_em->Print(Form("r%02d/jer//em_fraction_jer_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));  


  /// Widths fo sigeta and sigpsi
  TCanvas *cj3 = new TCanvas("cj3","cj3", 500, 500);
  cj3->SetLogy(0);
  for (int i3 = 0 ; i3 < 8 ; i3++)
    {
      for (int i = 0; i < nbins_pt; i++)
	{
	  
	  TH1D *h_sigp_sim = (TH1D*) h_trigger_ptp_pta_sim_3jet[i3]->ProjectionY("h_sigp_sim", i+1, i+1);
	  TH1D *h_sigp_data = (TH1D*) h_trigger_ptp_pta_data_3jet[i3]->ProjectionY("h_sigp_data", i+1, i+1);

	  dlutility::SetLineAtt(h_sigp_sim, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_sigp_sim, kRed, 1, 8);

	  dlutility::SetLineAtt(h_sigp_data, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_sigp_data, kBlue, 1, 8);

	  h_sigp_data->Scale(1./h_sigp_data->Integral(), "width");
	  h_sigp_sim->Scale(1./h_sigp_sim->Integral(), "width");

	  /* h_sigp_data->Fit("fg", "QR+", "",-20, 20); */
	  /* h_sigp_sim->Fit("fg", "QR+", "",-20, 20); */

	  /* h_sigp_data->GetFunction("fg")->SetLineColor(kBlue + 2); */
	  /* h_sigp_sim->GetFunction("fg")->SetLineColor(kRed + 2); */
	  
	  h_sigp_data->SetTitle(";#sigma(#it{p}_{T, #psi}) [GeV]; Arb.");
	  h_sigp_data->SetFillColorAlpha(kBlue - 3, .3);
	  h_sigp_data->SetMaximum(0.3);
	  h_sigp_data->Draw("hist");
	  h_sigp_data->Draw("same p");
	  h_sigp_sim->Draw("same");
	  
	  dlutility::DrawSPHENIXpp(0.22, 0.85);
	  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.75);
	  dlutility::drawText(Form("%2.1f < #LT#it{p}_{T}#GT #leq %2.1f GeV",pt_bin_edges[i], pt_bin_edges[i+1]) , 0.22, 0.7);
	  dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
	  //dlutility::drawText(Form("Sim Width = %2.2f (%2.2f)", h_sigp_sim->GetFunction("fg")->GetParameter(2), h_sigp_sim->GetRMS()), 0.22, 0.6);
	  TLegend *legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
	  legj3->SetLineWidth(0);
	  legj3->SetTextFont(42);
	  legj3->SetTextSize(0.04);
	  legj3->AddEntry(h_sigp_sim, Form("%s Reco", gen[index].c_str()));
	  legj3->AddEntry(h_sigp_data, "Data");

	  legj3->Draw("same");

	  cj3->Print(Form("r%02d/jer//sigp_r%02d_3jet%d_pta%d_%d_closure_%s.pdf", conesize , conesize, i3, i, closure, gen[index].c_str()));
	  cj3->Print(Form("r%02d/jer//sigp_r%02d_3jet%d_pta%d_%d_closure_%s.png", conesize , conesize, i3, i, closure, gen[index].c_str()));

	  TH1D *h_sige_sim = (TH1D*) h_trigger_ptv_pta_sim_3jet[i3]->ProjectionY("h_sige_sim", i+1, i+1);
	  TH1D *h_sige_data = (TH1D*) h_trigger_ptv_pta_data_3jet[i3]->ProjectionY("h_sige_data", i+1, i+1);


	  dlutility::SetLineAtt(h_sige_sim, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_sige_sim, kRed, 1, 8);

	  dlutility::SetLineAtt(h_sige_data, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_sige_data, kBlue, 1, 8);

	  h_sige_data->Scale(1./h_sige_data->Integral(), "width");
	  h_sige_sim->Scale(1./h_sige_sim->Integral(), "width");

	  /* h_sige_data->Fit("fg", "QR+","", -10, 10); */
	  /* h_sige_sim->Fit("fg", "Q+R","", -10, 10); */

	  /* h_sige_data->GetFunction("fg")->SetLineColor(kBlue + 2); */
	  /* h_sige_sim->GetFunction("fg")->SetLineColor(kRed + 2); */

	  h_sige_data->SetTitle(";#sigma(#it{p}_{T, #lambda}) [GeV]; Arb.");
	  h_sige_data->SetFillColorAlpha(kBlue - 3, .3);
	  h_sige_data->SetMaximum(0.3);

	  h_sige_data->Draw("hist");
	  h_sige_data->Draw("same p");
	  h_sige_sim->Draw("same");
	  
	  dlutility::DrawSPHENIXpp(0.22, 0.85);
	  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.75);
	  dlutility::drawText(Form("%2.1f < #LT#it{p}_{T}#GT #leq %2.1f GeV", pt_bin_edges[i], pt_bin_edges[i+1]), 0.22, 0.7);
	  dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);
	  //dlutility::drawText(Form("Sim Width = %2.2f (%2.2f)", h_sige_sim->GetFunction("fg")->GetParameter(2), h_sige_sim->GetRMS()), 0.22, 0.6);
	  legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
	  legj3->SetLineWidth(0);
	  legj3->SetTextSize(0.04);
	  legj3->SetTextFont(42);
	  legj3->AddEntry(h_sige_sim, Form("%s Reco", gen[index].c_str()));
	  legj3->AddEntry(h_sige_data, "Data");

	  legj3->Draw("same");

	  cj3->Print(Form("r%02d/jer//sige_r%02d_3jet%d_pta%d_%d_closure_%s.pdf", conesize , conesize, i3, i, closure, gen[index].c_str()));
	  cj3->Print(Form("r%02d/jer//sige_r%02d_3jet%d_pta%d_%d_closure_%s.png", conesize , conesize, i3, i, closure, gen[index].c_str()));

	  
	}
    }

  /// Widths fo sigeta and sigpsi
  TCanvas *cem4 = new TCanvas("cem4","cem4", 1000, 300);
  cem4->Divide(3, 1);
  TCanvas *cem3 = new TCanvas("cem3","cem3", 500, 500);
  cem3->SetLogy(0);
  for (int i3 = 0 ; i3 < 8 ; i3++)
    {
      TH2D *h_divide_trigger = (TH2D*) h_trigger_ptem_pta_data_3jet[i3]->Clone();
      TH2D *h_divide_lead = (TH2D*) h_lead_ptem_pta_data_3jet[i3]->Clone();
      TH2D *h_divide_sublead = (TH2D*) h_sublead_ptem_pta_data_3jet[i3]->Clone();

      TH2D *h_denominator_trigger = (TH2D*) h_trigger_ptem_pta_sim_3jet[i3]->Clone();
      TH2D *h_denominator_lead = (TH2D*) h_lead_ptem_pta_sim_3jet[i3]->Clone();
      TH2D *h_denominator_sublead = (TH2D*) h_sublead_ptem_pta_sim_3jet[i3]->Clone();

      h_divide_trigger->Scale(1./h_divide_trigger->Integral());
      h_divide_lead->Scale(1./h_divide_lead->Integral());
      h_divide_sublead->Scale(1./h_divide_sublead->Integral());

      h_denominator_trigger->Scale(1./h_denominator_trigger->Integral());
      h_denominator_lead->Scale(1./h_denominator_lead->Integral());
      h_denominator_sublead->Scale(1./h_denominator_sublead->Integral());

      h_divide_trigger->Divide(h_denominator_trigger);
      h_divide_lead->Divide(h_denominator_lead);
      h_divide_sublead->Divide(h_denominator_sublead);

      cem4->cd(1);
      gPad->SetRightMargin(0.2);
      h_divide_trigger->SetMaximum(2.0);
      h_divide_trigger->Draw("colz");
      
      cem4->cd(2);
      gPad->SetRightMargin(0.2);
      h_divide_lead->SetMaximum(2.0);
      h_divide_lead->Draw("colz");
      cem4->cd(3);
      gPad->SetRightMargin(0.2);
      h_divide_sublead->SetMaximum(2.0);
      h_divide_sublead->Draw("colz");
 
      cem4->Print(Form("r%02d/jer//ptem2d_r%02d_3jet%d_%d_closure_%s.pdf", conesize , conesize, i3, closure, gen[index].c_str()));
      cem4->Print(Form("r%02d/jer//ptem2d_r%02d_3jet%d_%d_closure_%s.png", conesize , conesize, i3, closure, gen[index].c_str()));

      cem3->cd();
      dlutility::SetLineAtt(hp_trigger_ptem_pta_sim_3jet[i3], kRed, 2, 1);
      dlutility::SetMarkerAtt(hp_trigger_ptem_pta_sim_3jet[i3], kRed, 1, 20);

      dlutility::SetLineAtt(hp_trigger_ptem_pta_data_3jet[i3], kBlue, 2, 1);
      dlutility::SetMarkerAtt(hp_trigger_ptem_pta_data_3jet[i3], kBlue, 1, 20);

      dlutility::SetLineAtt(hp_lead_ptem_pta_sim_3jet[i3], kRed, 2, 1);
      dlutility::SetMarkerAtt(hp_lead_ptem_pta_sim_3jet[i3], kRed, 1, 24);

      dlutility::SetLineAtt(hp_lead_ptem_pta_data_3jet[i3], kBlue, 2, 1);
      dlutility::SetMarkerAtt(hp_lead_ptem_pta_data_3jet[i3], kBlue, 1, 24);

      dlutility::SetLineAtt(hp_sublead_ptem_pta_sim_3jet[i3], kRed, 2, 1);
      dlutility::SetMarkerAtt(hp_sublead_ptem_pta_sim_3jet[i3], kRed, 1, 25);

      dlutility::SetLineAtt(hp_sublead_ptem_pta_data_3jet[i3], kBlue, 2, 1);
      dlutility::SetMarkerAtt(hp_sublead_ptem_pta_data_3jet[i3], kBlue, 1, 25);
      
      
      hp_trigger_ptem_pta_data_3jet[i3]->SetTitle(";#sigma(#it{p}_{T, #psi}) [GeV]; Arb.");
      hp_trigger_ptem_pta_data_3jet[i3]->SetMinimum(0.5);
      hp_trigger_ptem_pta_data_3jet[i3]->SetMaximum(1.0);
      hp_trigger_ptem_pta_data_3jet[i3]->Draw("p");
      hp_trigger_ptem_pta_sim_3jet[i3]->Draw("p same");
      hp_sublead_ptem_pta_data_3jet[i3]->Draw("p same");
      hp_sublead_ptem_pta_sim_3jet[i3]->Draw("p same");
      hp_lead_ptem_pta_data_3jet[i3]->Draw("p same");
      hp_lead_ptem_pta_sim_3jet[i3]->Draw("p same");
      
      dlutility::DrawSPHENIXpp(0.22, 0.85);
      dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.75);
      TLegend *legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
      legj3->SetLineWidth(0);
      legj3->SetTextFont(42);
      legj3->SetTextSize(0.04);
      //legj3->AddEntry(h_sige_sim, 
      legj3->AddEntry(hp_trigger_ptem_pta_sim_3jet[i3], Form("%s Reco", gen[index].c_str()));
      legj3->AddEntry(hp_trigger_ptem_pta_data_3jet[i3], "Data");
      
      legj3->Draw("same");
      
      cem3->Print(Form("r%02d/jer//ptem_r%02d_3jet%d_%d_closure_%s.pdf", conesize , conesize, i3, closure, gen[index].c_str()));
      cem3->Print(Form("r%02d/jer//ptem_r%02d_3jet%d_%d_closure_%s.png", conesize , conesize, i3, closure, gen[index].c_str()));

      for (int i = 0; i < nbins_pt; i++)
	{
	  TH1D *h_ptem_sim = (TH1D*) h_trigger_ptem_pta_sim_3jet[i3]->ProjectionY("h_ptem_sim", i+1, i+1);
	  TH1D *h_ptem_data = (TH1D*) h_trigger_ptem_pta_data_3jet[i3]->ProjectionY("h_ptem_data", i+1, i+1);

	  TH1D *h_ptem_sim_lead = (TH1D*) h_lead_ptem_pta_sim_3jet[i3]->ProjectionY("h_ptem_sim_lead", i+1, i+1);
	  TH1D *h_ptem_data_lead = (TH1D*) h_lead_ptem_pta_data_3jet[i3]->ProjectionY("h_ptem_data_lead", i+1, i+1);

	  TH1D *h_ptem_sim_sublead = (TH1D*) h_sublead_ptem_pta_sim_3jet[i3]->ProjectionY("h_ptem_sim_sublead", i+1, i+1);
	  TH1D *h_ptem_data_sublead = (TH1D*) h_sublead_ptem_pta_data_3jet[i3]->ProjectionY("h_ptem_data_sublead", i+1, i+1);

	  dlutility::SetLineAtt(h_ptem_sim, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_sim, kRed, 1, 8);

	  dlutility::SetLineAtt(h_ptem_data, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_data, kBlue, 1, 8);

	  dlutility::SetLineAtt(h_ptem_sim_lead, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_sim_lead, kRed, 1, 24);

	  dlutility::SetLineAtt(h_ptem_data_lead, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_data_lead, kBlue, 1, 24);

	  dlutility::SetLineAtt(h_ptem_sim_sublead, kRed, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_sim_sublead, kRed, 1, 25);

	  dlutility::SetLineAtt(h_ptem_data_sublead, kBlue, 2, 1);
	  dlutility::SetMarkerAtt(h_ptem_data_sublead, kBlue, 1, 25);

	  h_ptem_data->Scale(1./h_ptem_data->Integral(), "width");
	  h_ptem_sim->Scale(1./h_ptem_sim->Integral(), "width");
	
	  h_ptem_data_lead->Scale(1./h_ptem_data_lead->Integral(), "width");
	  h_ptem_sim_lead->Scale(1./h_ptem_sim_lead->Integral(), "width");

	  h_ptem_data_sublead->Scale(1./h_ptem_data_sublead->Integral(), "width");
	  h_ptem_sim_sublead->Scale(1./h_ptem_sim_sublead->Integral(), "width");

	  h_ptem_data->SetTitle(";p_{T}^{EM}; Arb.");
	  h_ptem_data->SetFillColorAlpha(kBlue - 3, .3);
	  h_ptem_data->SetMaximum(4.0);;
	  h_ptem_data->Draw("");
	  h_ptem_sim->Draw("same");
	  h_ptem_data_lead->Draw("same");
	  h_ptem_sim_lead->Draw("same");
	  h_ptem_data_sublead->Draw("same");
	  h_ptem_sim_sublead->Draw("same");
	  
	  dlutility::DrawSPHENIXpp(0.22, 0.85);
	  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.22, 0.75);
	  dlutility::drawText(Form("%2.1f < #LT#it{p}_{T}#GT #leq %2.1f GeV", pt_bin_edges[i], pt_bin_edges[i+1]), 0.22, 0.7);
	  dlutility::drawText(Form("#it{p}_{T,3} < %2.0f GeV", third_jet_cuts[i3]), 0.22, 0.65);

	  TLegend *legj3 = new TLegend(0.6, 0.77, 0.75, 0.877);
	  legj3->SetLineWidth(0);
	  legj3->SetTextFont(42);
	  legj3->SetTextSize(0.04);
	  legj3->AddEntry(h_ptem_sim, Form("%s Reco", gen[index].c_str()));//"Pythia-8 Reco");
	  legj3->AddEntry(h_ptem_data, "Data");
	  legj3->AddEntry(h_ptem_sim_lead, Form("Lead %s Reco", gen[index].c_str()));//"Lead Pythia-8 Reco");
	  legj3->AddEntry(h_ptem_data_lead, "Lead Data");
	  legj3->AddEntry(h_ptem_sim_sublead, Form("Sub %s Reco", gen[index].c_str()));//"Sub Pythia-8 Reco");
	  legj3->AddEntry(h_ptem_data_sublead, "Sub Data");

	  legj3->Draw("same");

	  cem3->Print(Form("r%02d/jer//ptem_r%02d_3jet%d_pta%d_%d_closure_%s.pdf", conesize , conesize, i3, i, closure, gen[index].c_str()));
	  cem3->Print(Form("r%02d/jer//ptem_r%02d_3jet%d_pta%d_%d_closure_%s.png", conesize , conesize, i3, i, closure, gen[index].c_str()));
	}
    }

  int color_sample[5] = {kBlue, kGreen, kRed, kOrange, kViolet};
  int color_combined = kBlack;
  for (int i = 0; i < 5; i++)
    {
      dlutility::SetLineAtt(h_lead_sample[i], color_sample[i], 2, 1);
      dlutility::SetMarkerAtt(h_lead_sample[i], color_sample[i], 1.1, 20);
    }
  dlutility::SetLineAtt(h_lead_combined, color_combined, 2, 1);
  dlutility::SetMarkerAtt(h_lead_combined, color_combined, 1.1, 24);

  TCanvas *ccom = new TCanvas("ccom","ccom", 500, 500);
  ccom->SetLogy();
  h_lead_combined->Draw("p");
  h_lead_sample[0]->Draw("p same");
  h_lead_sample[1]->Draw("p same");
  h_lead_sample[2]->Draw("p same");
  h_lead_sample[3]->Draw("p same");
  h_lead_sample[4]->Draw("p same");
  h_lead_combined->Draw("p same");

  ccom->Print(Form("r%02d/jer//sample_combine_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));


  double xmin = 3, xmax = 100;
  double dx = 0.1;
  int nbins_smear = (int) ((xmax - xmin) / dx) ;
  TCanvas *cjj = new TCanvas("cjj","cjj", 500, 500);
  
  TH1D *hhh = new TH1D("hhh","; #it{p}_{T}^{truth} ; Quad. Diff. #frac{#sigma(#it{p}_{T})}{#it{p}_{T}}", 1 , xmin, xmax);
  hhh->SetMaximum(0.4);
  hhh->Draw();

  TF1 *band = (TF1*) flatline->Clone();
  band->SetName("band");
  band->SetRange(xmin, xmax);
  TF1 *band_low = (TF1*) flatline->Clone();
  band_low->SetName("band_low");
  band_low->SetRange(xmin, xmax);
  TF1 *band_high = (TF1*) flatline->Clone();
  band_high->SetName("band_high");
  band_high->SetRange(xmin, xmax);
  TF1 *band_im = (TF1*) flatline->Clone();
  band_im->SetName("band_im");
  band_im->SetRange(xmin, xmax);
  TF1 *band_bi = (TF1*) flatline->Clone();
  band_bi->SetName("band_bi");
  band_bi->SetRange(xmin, xmax);

  
  float average_bi = 0;
  float average_im = 0;
  float sum_bi = 0;
  float sum_im = 0;
  float sigma_bi = 0;
  float sigma_im = 0;
  std::cout << "CALCULATE" << std::endl;
  for (int i3 = 0; i3 < 5; i3++)
    { 
      std::cout << i3 << " : " << mean_line_bi[i3] << " / " << err_line_bi[i3] << " / " << mean_line_im[i3] << " / " << err_line_im[i3] << std::endl;
      float wb = 1./TMath::Power(err_line_bi[i3], 2);
      float wi = 1./TMath::Power(err_line_im[i3], 2);
      average_bi += mean_line_bi[i3]*wb;
      sum_bi += wb;//1./(err_line_bi[i3]);
      average_im += mean_line_im[i3]*wi;///err_line_im[3];
      sum_im += wi;///(err_line_im[i3]);
    }
  average_bi /= sum_bi;
  average_im /= sum_im;
  sigma_bi = err_line_bi[3];
  sigma_im = err_line_im[3];//1./sqrt(sum_im);
  std::cout << average_bi << " / " << average_im << std::endl;
  std::cout << sigma_bi << " / " << sigma_im << std::endl;
  
  
  band_im->SetParameter(0, average_im);
  band_bi->SetParameter(0, average_bi);

  float total_average = (average_im + average_bi)/2.;
  float total_diff = (average_im - average_bi)/2.;
  
  float total_sys = sqrt(sigma_bi*sigma_bi + sigma_im*sigma_im + total_diff*total_diff);
  band->SetParameter(0, total_average);
  band_high->SetParameter(0, total_average + total_sys);
  band_low->SetParameter(0, total_average - total_sys);

  // The pt dependent
  TProfile *h_pt_dependent_smear_im = new TProfile("h_pt_dependent_smear_im", ";p_{T};smear", nbins_smear, xmin, xmax, "s");
  TProfile *h_pt_dependent_smear_bi = new TProfile("h_pt_dependent_smear_bi", ";p_{T};smear", nbins_smear, xmin, xmax, "s");
  for (int i = 0; i < 4; i++)
    {

      for (int ib = 0; ib < nbins_smear; ib++)
	{
	  float v = h_pt_dependent_smear_im->GetBinCenter(ib+1);

	  float im_diff = 0;
	  float data_im_v = f_data_im[i]->Eval(v);
	  float sim_im_v = f_sim_im[i]->Eval(v);
	  if (data_im_v > sim_im_v)
	    {
	      im_diff = sqrt(data_im_v*data_im_v - sim_im_v*sim_im_v);
	    }
	  float bi_diff = 0;
	  float data_bi_v = f_data_bi[i]->Eval(v);
	  float sim_bi_v = f_sim_bi[i]->Eval(v);
	  if (data_bi_v > sim_bi_v)
	    {
	      bi_diff = sqrt(data_bi_v*data_bi_v - sim_bi_v*sim_bi_v);
	    }
	  h_pt_dependent_smear_im->Fill(v, im_diff);
	  h_pt_dependent_smear_bi->Fill(v, bi_diff);
	}
    }


  
  TH1D *h_nominal_smear =  new TH1D("h_nominal_smear", ";p_{T};smear", nbins_smear, xmin, xmax);
  TH1D *h_sys_up_smear =  new TH1D("h_sys_up_smear", ";p_{T};smear", nbins_smear, xmin, xmax);
  TH1D *h_sys_down_smear =  new TH1D("h_sys_down_smear", ";p_{T};smear", nbins_smear, xmin, xmax);

  for (int ib = 0; ib < nbins_smear; ib++)
    {
      float imv = h_pt_dependent_smear_im->GetBinContent(ib+1);
      float biv = h_pt_dependent_smear_bi->GetBinContent(ib+1);
      float average_nom = (biv + imv)/2.;

      float diff_sys = (biv - imv)/2.;
      float err_bi = h_pt_dependent_smear_bi->GetBinError(ib+1);
      float err_im = h_pt_dependent_smear_im->GetBinError(ib+1);

      h_nominal_smear->SetBinContent(ib+1, average_nom);
      float sys = sqrt(diff_sys*diff_sys + err_bi*err_bi + err_im*err_im + sigma_im*sigma_im + sigma_bi*sigma_bi);
      h_sys_up_smear->SetBinContent(ib+1, average_nom + sys);
      h_sys_down_smear->SetBinContent(ib+1, (average_nom - sys > 0 ? average_nom - sys : 0 ));
    }

 
  dlutility::SetLineAtt(band, kRed, 3, 1);
  dlutility::SetLineAtt(band_low, kRed + 2, 3, 1);
  dlutility::SetLineAtt(band_high, kRed + 2, 3, 1);

  dlutility::SetLineAtt(h_nominal_smear, kGreen, 3, 1);
  dlutility::SetLineAtt(h_sys_down_smear, kGreen + 2, 3, 1);
  dlutility::SetLineAtt(h_sys_up_smear, kGreen + 2, 3, 1);

  
  dlutility::SetLineAtt(band_im, kBlue, 3, 1);
  dlutility::SetLineAtt(band_bi, kBlue + 2, 3, 1);
  TBox *bsys = new TBox(xmin, band_low->GetParameter(0), xmax, band_high->GetParameter(0));
  bsys->SetLineWidth(0);
  bsys->SetFillStyle(1001);
  bsys->SetFillColorAlpha(kRed - 2, 0.5);
  bsys->Draw("same");
  band->Draw("same");
  band_low->Draw("same");
  band_high->Draw("same");
  band_im->Draw("same");
  band_bi->Draw("same");
  
  dlutility::DrawSPHENIXppInternalsize(0.6, 0.85, 0.04);
  dlutility::drawText(Form("anti-#it{k_{t}} #it{R} = %0.1f", conesize*0.1), 0.6, 0.75);

  TLegend *lejj = new TLegend(0.6, 0.55, 0.84, 0.7);
  lejj->SetLineWidth(0);
  lejj->SetTextSize(0.03);
  lejj->SetTextFont(42);

  lejj->AddEntry(band_im, "Imbalance Method","l");
  lejj->AddEntry(band_bi, "Bisector Method","l");
  lejj->AddEntry(band, "Nominal +/- Sys","l");

  lejj->Draw("same");
  
  
  cjj->Print(Form("r%02d/jer//jer_fits_r%02d_%d_closure_%s.pdf", conesize , conesize, closure, gen[index].c_str()));

  TFile *fout_jer = new TFile(Form("r%02d/jer/jer_fits_r%02d_1JES_%d_closure_%s.root", conesize, conesize, closure, gen[index].c_str()),"recreate");
  band->SetName("jernom");
  band_low->SetName("jerneg");
  band_high->SetName("jerpos");
  h_nominal_smear->Write();
  h_sys_up_smear->Write();
  h_sys_down_smear->Write();
  band->Write();
  band_low->Write();
  band_high->Write();
  
  fout_jer->Close();
  return;
}

std::vector<std::pair<struct jet, struct jet>>  match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto jet : myrecojets)
    {
      for (auto tjet : mytruthjets)
	{
	  if (std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;}) != matched_dijets.end()) continue;
	  float dR = getDR(jet, tjet);
	      
	  if (dR < dRcut*cone_size)
	    {
	      if (std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;}) != matched_dijets.end()) continue;


	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;

	      matched_dijets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }

  return matched_dijets;
}

bool isIsolated(struct jet myjet, std::vector<struct jet> myrecojets)
{
  for (auto jet : myrecojets)
    {
      if (myjet.id == jet.id) continue;
      if (jet.pt < 5) continue;
      float dR = getDR(jet, myjet);
	      
      if (dR < (isocut))//*cone_size))
	{
	  return false;
	}
    }

  return true;
}



float getDPHI(float phi1, float phi2)
{
  float dphitr = phi1 - phi2;
  if (dphitr > TMath::Pi())
    {
      dphitr -= 2*TMath::Pi();
    }
  if (dphitr < -1*TMath::Pi())
    {
      dphitr += 2*TMath::Pi();
    }
  return dphitr;
  
}
bool check_dijet_reco(std::vector<struct jet> myrecojets)
{

  if (myrecojets.size() < 2) return false;
    
  auto leading_iter = myrecojets.begin();
  
  auto subleading_iter = myrecojets.begin() + 1;
  
  if (fabs(leading_iter->eta) > etacut_cone || fabs(subleading_iter->eta) > etacut_cone) return false;
  if (fabs(leading_iter->eta_det) > etacut_cone || fabs(subleading_iter->eta_det) > etacut_cone) return false;
	  
  bool isolated = true;
  isolated &= isIsolated(*leading_iter, myrecojets);
  isolated &= isIsolated(*subleading_iter, myrecojets);
  if (!isolated) return false;

  float dphir = fabs(getDPHI(leading_iter->phi, subleading_iter->phi));

  if (!(leading_iter->pt >= 7 && subleading_iter->pt >= 3 && dphir >= dphicut)) return false;

  /* double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t); */
  /* double jetleadtime = 17.6*(leading_iter->t); */
  /* bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < 6.0 ); */
  /* bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);	   */
  
  /* bool passbothtime = (passdijettime) && (passleadtime); */
  
  /* if (!passbothtime) return false; */

  return true;
  
}
bool check_dijet_truth(std::vector<struct jet> mytruthjets)
{

  if (mytruthjets.size() < 2) return false;
    
  auto leading_iter = mytruthjets.begin();
  
  auto subleading_iter = mytruthjets.begin() + 1;
  
  if (fabs(leading_iter->eta) > etacut_cone || fabs(subleading_iter->eta) > etacut_cone) return false;
	  
  bool isolated = true;
  isolated &= isIsolated(*leading_iter, mytruthjets);
  isolated &= isIsolated(*subleading_iter, mytruthjets);

  if (!isolated) return false;

  float dphir = fabs(getDPHI(leading_iter->phi, subleading_iter->phi));

  if (!(leading_iter->pt >= 7 && subleading_iter->pt >= 3 && dphir >= dphicut)) return false;

  return true;

  

}