#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"

#include "read_binning.h"


using std::cout;
using std::endl;

struct jet
{
  int id;
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
};

const bool Debug = false;

const float dRcut = 1.3;
float cone_size = 4;
const float truth_cut = 3;
const float reco_cut = 8;
const float etacut = 0.8;

const float dphicut = 0;//3*TMath::Pi()/4.;
const float dphicutloose = 0;//3*TMath::Pi()/4.;

const float vertex_cut = 60;

void makeDataTreeAuAu(const int cone_size_int = 3, const int isAuAu = 1,
                      const std::string infile = "/home/tmengel/PPG14/rootfiles/v001_20260720/run2auau_rho_jet.root",
                      const std::string outfile = "",
                      const std::string configfile = "binning_AA.config")
{

  read_binning rb(configfile.c_str());

  cone_size = (float) cone_size_int;
  std::cout << cone_size << std::endl;
  std::string mycopy = infile;
  mycopy = mycopy.substr(mycopy.rfind("/")+1);
  std::cout << mycopy << std::endl;
  std::string newname = "TNTUPLE_DIJET_r0" + std::to_string(cone_size_int);
  std::string oldname = "run2auau";
  const std::size_t oldpos = mycopy.find(oldname);
  if (oldpos != std::string::npos)
    mycopy.replace(oldpos, oldname.length(), newname);
  else
    mycopy = newname + "_" + mycopy;
  std::string newfile = outfile.empty() ? rb.get_tntuple_location() + "/" + mycopy : outfile;
  std::cout << newfile << std::endl;

  TFile *f = new TFile(infile.c_str(), "r");
  TTree *t = (TTree*) f->Get("T");
  if (!t)
    {
      std::cerr << "Could not find tree T in " << infile << std::endl;
      return;
    }

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("cent", 1);
  t->SetBranchStatus("sumeT", 1);
  t->SetBranchStatus("jet_pT", 1);
  t->SetBranchStatus("jet_unsub_pT", 1);
  t->SetBranchStatus("jet_E", 1);
  t->SetBranchStatus("jet_unsub_E", 1);
  t->SetBranchStatus("jet_eta", 1);
  t->SetBranchStatus("jet_phi", 1);
  t->SetBranchStatus("zvrtx", 1);
  t->SetCacheSize(256 * 1024 * 1024);
  t->AddBranchToCache("*", kTRUE);

  

  int centrality;
  int minbias;
  float mbd_vertex_z;
  ULong64_t gl1_scaled;
  float sumeT;
  float mbd_time_zero;
    std::vector<float> *reco_jet_pt = 0;
  std::vector<float> *reco_jet_pt_unsub = 0;
  std::vector<float> *reco_jet_e = 0;
  std::vector<float> *reco_jet_e_unsub = 0;
  std::vector<float> *reco_jet_eta = 0;
  std::vector<float> *reco_jet_phi = 0;


  t->SetBranchAddress("cent", &centrality);
  t->SetBranchAddress("sumeT", &sumeT);
  t->SetBranchAddress("jet_pT", &reco_jet_pt);
  t->SetBranchAddress("jet_unsub_pT", &reco_jet_pt_unsub);
  t->SetBranchAddress("jet_E", &reco_jet_e);
  t->SetBranchAddress("jet_unsub_E", &reco_jet_e_unsub);
  t->SetBranchAddress("jet_eta", &reco_jet_eta);
  t->SetBranchAddress("jet_phi", &reco_jet_phi);
  t->SetBranchAddress("zvrtx", &mbd_vertex_z);

  // t->SetBranchAddress("gl1_scaled", &gl1_scaled);
  // t->SetBranchAddress("centrality", &centrality);
  // t->SetBranchAddress("minbias", &minbias);
  // t->SetBranchAddress("mbd_time_zero", &mbd_time_zero);
  // t->SetBranchAddress(Form("jet_pt_%d_sub", cone_size_int), &reco_jet_pt);
  // t->SetBranchAddress(Form("jet_pt_unsub_%d_sub", cone_size_int), &reco_jet_pt_unsub);
  // t->SetBranchAddress(Form("jet_e_%d_sub", cone_size_int), &reco_jet_e);
  // t->SetBranchAddress(Form("jet_e_unsub_%d_sub", cone_size_int), &reco_jet_e_unsub);
  // t->SetBranchAddress(Form("jet_eta_%d_sub", cone_size_int), &reco_jet_eta);
  // t->SetBranchAddress(Form("jet_phi_%d_sub", cone_size_int), &reco_jet_phi);
  // t->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);

  int total_entries = t->GetEntries();
  int first_entry = 0;
  if (const char *first_entry_env = std::getenv("DIJET_DATA_FIRST_ENTRY"))
    {
      const int requested_first_entry = std::atoi(first_entry_env);
      if (requested_first_entry > 0 && requested_first_entry < total_entries)
        first_entry = requested_first_entry;
    }
  int entries = total_entries - first_entry;
  if (const char *max_events_env = std::getenv("DIJET_MAX_DATA_EVENTS"))
    {
      const int max_events = std::atoi(max_events_env);
      if (max_events > 0 && max_events < entries)
        entries = max_events;
    }
  if (first_entry > 0 || entries < total_entries)
    {
      std::cout << "Processing entries [" << first_entry << ", " << first_entry + entries
                << ") from " << total_entries << " total entries" << std::endl;
    }

  TFile *fout = new TFile(newfile.c_str(), "recreate");
  TNtuple *tn_dijet = new TNtuple("tn_dijet","matched truth and reco","pt1_reco:pt2_reco:dphi_reco:deta_reco:trigger:njets:centrality:mbd_vertex:sumeT");

  
  std::pair<int, float> id_leaders[2];
  TF1 *fcut = new TF1("fcut","[0]+[1]*TMath::Exp(-[2]*x)",0.0,100.0);
  //  fcut->SetParameters(2.5,36.2,0.035);
  fcut->SetParameters(0.0,40,0.038);
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(first_entry + i);

      if (i % 100000 == 0)
        std::cout << "Event: " << i << " / " << entries << "\r" << std::flush;

      // if (fabs(mbd_vertex_z) > vertex_cut) continue;
      // if (!minbias) continue;

      int trigger_fired = 1;

      // if ((gl1_scaled >> 10) & 0x1U) trigger_fired = 10;

      // if (trigger_fired == 0) continue;

      int nrecojets = reco_jet_pt->size();
      // this is the reco index for the leading and subleading jet
      if (Debug) std::cout << __LINE__ << std::endl;
      id_leaders[0] = std::make_pair(0, 0);
      id_leaders[1] = std::make_pair(0, 0);
      int njet_good = 0;
      int bad_event = 0;
      float cut_value = fcut->Eval(centrality);
      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt->at(j) < reco_cut) continue;

	  if (reco_jet_e->at(j) < 0) continue;
	  if (reco_jet_e_unsub->at(j) < 0) continue;
	  if (fabs(reco_jet_eta->at(j)) > etacut) continue;
	  
	  float pt_unsub = reco_jet_pt_unsub->at(j) - reco_jet_pt->at(j);

	  if (pt_unsub > cut_value)
	    {
	      bad_event = 1;
	      continue;
	    }


	  njet_good++;
	  if (reco_jet_pt->at(j) > id_leaders[0].second)
	    {
	      id_leaders[1] = id_leaders[0];
	      id_leaders[0] = std::make_pair(j, reco_jet_pt->at(j));
	    }
	  else if (reco_jet_pt->at(j) > id_leaders[1].second)
	    {
	      id_leaders[1] = std::make_pair(j, reco_jet_pt->at(j));
	    } 	  

	}
      if (bad_event) continue;
      if (njet_good < 2) continue;

      const int lead = id_leaders[0].first;
      const int sub = id_leaders[1].first;
      const float lead_pt = reco_jet_pt->at(lead);
      const float sub_pt = reco_jet_pt->at(sub);
      double dphir = fabs(reco_jet_phi->at(lead) - reco_jet_phi->at(sub));
      if (dphir > TMath::Pi())
        dphir = 2*TMath::Pi() - dphir;
      const double detar = fabs(reco_jet_eta->at(lead) - reco_jet_eta->at(sub));

      if (dphir >= dphicut && lead_pt >= 15 && sub_pt >= 5)
        tn_dijet->Fill(lead_pt, sub_pt, dphir, detar, trigger_fired, nrecojets, centrality, mbd_vertex_z, sumeT);
    
    }
  std::cout << std::endl;

  tn_dijet->Write();

  fout->Close();

}
