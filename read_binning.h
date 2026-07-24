#ifndef READ_BINNING_H
#define READ_BINNING_H
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

class read_binning
{
public:
  read_binning(const std::string configfile)
    {
      fjer = new TFile(jerfile.c_str(), "r");
      
      penv = new TEnv(configfile.c_str());
      tntuple_location = std::getenv("DIJET_TNTUPLE_PATH");
      code_location = std::getenv("DIJET_UNFOLDING_PATH");
      sim_location = std::getenv("AUAU_SIM_PATH");
      auau_config = std::getenv("AUAU_CONFIG");
    }
  std::string get_tntuple_location() 
  {
    return tntuple_location;      
  }
  std::string get_code_location() 
  {
    return code_location;      
  }
  std::string get_sim_location() 
  {
    return sim_location;      
  }
  std::string get_jesr_location() 
  {
    return jesr_location;      
  }

  std::vector<int> getRunnumbers()
  {
    // Define vectors to store the data                                                                                                                                                                                                    
    std::vector<int> v;
    // Open the CSV file                                                                                                                                                                                                                   
    std::ifstream file("grl24.txt");
    
    if (!file.is_open()) {
      std::cerr << "Error opening file!" << std::endl;
      return v;
    }

    std::string line;
    while (std::getline(file, line)) {
      std::istringstream ss(line);
      std::string runNumberStr;

      int runNumber = std::stoi(line);

      v.push_back(runNumber);
    }

    file.close();

    return v;
  }

  Int_t get_nbins(){ return penv->GetValue("nbins", 1); }
  Int_t get_bbins(){ return penv->GetValue("bbins", 1); }
  Int_t get_minentries(){ return penv->GetValue("minentries", 1); }
  Int_t get_measure_bins(){ return penv->GetValue("measure_bins", 1); }
  Double_t get_minimum(){ return penv->GetValue("minimum", 1.0); }
  Double_t get_maximum_reco(){ return penv->GetValue("maximum", 68.0); }
  Int_t get_maximum_reco_bin() {return max_reco_bin; }
  
  Double_t get_fixed(){ return penv->GetValue("fixed", 1.0); }
  Int_t get_primer(){ return penv->GetValue("primer", 0); }
  Int_t get_njet_sys(){ return penv->GetValue("NJET", 0); }
  Int_t get_zyam_sys(){ return penv->GetValue("ZYAM", 0); }
  Double_t get_flow_sys(){ return penv->GetValue("FLOW_V22_SCALE", 1.0); }
  Double_t get_flow_v33_sys(){ return penv->GetValue("FLOW_V33_SCALE", 1.0); }
  std::string get_flow_systematic_name()
  {
    return penv->GetValue("FLOW_SYSTEMATIC_NAME", "COMB");
  }
  Int_t get_prior_sys(){ return penv->GetValue("PRIOR", 0); }
  Double_t get_pileup_sys(){ return penv->GetValue("PILEUP", 0.0); }
  Int_t get_full_sys(){ return penv->GetValue("FULL", 0); }
  Int_t get_emfrac_sys(){ return penv->GetValue("EMFRAC", 0); }
  Int_t get_trigger_sys(){ return penv->GetValue("TRIGEFF", 0); }
  Int_t get_crossingangle_sys(){ return penv->GetValue("CA", 0); }
  Int_t get_philoc_sys(){ return penv->GetValue("PHILOC", 0); }
  Int_t get_vtx_sys(){ return penv->GetValue("VTX", 0); }
  Int_t get_herwig(){ return penv->GetValue("HERWIG", 0); }
  
  Int_t get_inclusive_sys(){ return penv->GetValue("INCLUSIVE", 0); }
  Double_t get_jes_sys(){ return penv->GetValue("JES", 0.0); }
  Double_t get_jer_sys(){ return penv->GetValue("JER", 0.0); }
  Double_t get_first_xj(){ return penv->GetValue("first_xj", 0.0); }
  
  Double_t get_pt_regions(int reg){ return penv->GetValue(Form("pt1_%d", reg), 1.0); }
  Double_t get_pt2_regions(int reg){ return penv->GetValue(Form("pt2_%d", reg), 1.0); }
  
  Double_t get_low_trigger(int reg){ return low_trigger[reg];}
  Double_t get_low_trigger_goal() { return penv->GetValue("low_trigger", 1.0);}
  
  Int_t get_measure_region(int reg){ return measure_region[reg]; }
  Int_t get_subleading_measure_region(int reg){ return subleading_measure_region[reg]; }

  Double_t get_sample_boundary(int ib){ return sample_boundary[ib]; }

  Double_t get_vtx_cut() { return 60; }
  Double_t get_njet_cut() { return penv->GetValue("njet_cut", 3); }
  std::string get_dphi_string(){ return Form("%d#pi/%d", (int)penv->GetValue("dphi_top", 1.0), (int) penv->GetValue("dphi_bottom", 2.0)); }

  Double_t get_dphicut() { return penv->GetValue("dphi_top", 1.0) * TMath::Pi() / penv->GetValue("dphi_bottom", 2.0); }
  Double_t get_truth_leading_goal(){ return penv->GetValue("truth_leading_goal", 14.0); }
  Double_t get_truth_subleading_goal(){ return penv->GetValue("truth_subleading_goal", 14.0); }

  Int_t get_reco_leading_bin_buffer(){ return penv->GetValue("reco_leading_bin_buffer", 14.0); }
  Int_t get_reco_subleading_bin_buffer(){ return penv->GetValue("reco_subleading_bin_buffer", 14.0); }

  Int_t get_measure_leading_bin_buffer(){ return penv->GetValue("measure_leading_bin_buffer", 14.0); }
  Int_t get_measure_subleading_bin_buffer(){ return penv->GetValue("measure_subleading_bin_buffer", 14.0); }


  int get_truth_leading_bin(){ return truth_leading_bin;}
  int get_truth_subleading_bin(){ return truth_subleading_bin;}

  int get_reco_leading_bin(){ return reco_leading_bin;}
  int get_reco_subleading_bin(){ return reco_subleading_bin;}

  int get_measure_leading_bin(){ return measure_leading_bin;}
  int get_measure_subleading_bin(){ return measure_subleading_bin;}

  
  float get_truth_leading_cut(){ return truth_leading_cut;}
  float get_truth_subleading_cut(){ return truth_subleading_cut;}

  float get_reco_leading_cut(){ return reco_leading_cut;}
  float get_reco_subleading_cut(){ return reco_subleading_cut;}

  float get_measure_leading_cut(){ return measure_leading_cut;}
  float get_measure_subleading_cut(){ return measure_subleading_cut;}

  float get_zyam_low() {return 10.*TMath::Pi()/32.;};
  float get_zyam_high() {return 14.*TMath::Pi()/32.;};

  int get_number_centrality_bins() { return centrality_bins;}

  TF1* get_smear_function(int centrality_bin)
  {
    Double_t jer_sys = get_jer_sys();

    std::string sys_name = "f_nominal";
    if (jer_sys < 0)
      {
	sys_name = "f_negative";
      }
    else if (jer_sys > 0)
      {
	sys_name = "f_positive";
      }
    f_smear_function = (TF1*) fjer->Get(Form("%s_%d", sys_name.c_str(), centrality_bin));

    return f_smear_function;
  }
  void get_centrality_bins(float icentrality_bins[])
  {
    icentrality_bins[0] = 0;
    icentrality_bins[1] = 10;
    icentrality_bins[2] = 30;
    icentrality_bins[3] = 50;
    icentrality_bins[4] = 90;
  }
  void get_pt_bins(float ipt_bins[])
  {
    Float_t minimum = get_minimum();
    Float_t fixed = get_fixed();
    Int_t beforebin = get_bbins();
    Float_t low_trigger_goal = get_low_trigger_goal();

    Float_t reco_max_goal = get_maximum_reco();
    float alpha = TMath::Power(fixed/minimum, 1/(float)beforebin);

    Float_t truth_leading_goal = get_truth_leading_goal();

    Int_t nbins = get_nbins();
    int ireg = 0;
    int ireg2 = 0;
    for (int i = 0; i < nbins+1; i++)
      {
	float ipt = minimum*TMath::Power(alpha, (float)i);
	ipt_bins[i] = ipt;
	for (int ib = 0; ib < 4; ib++)
	  {
	    if (sample_boundary[ib] == 0 && ipt >= sample_boundary_goal[ib])
	      {
		sample_boundary[ib] = ipt;
	      }
	  }
	if (low_trigger_goal <= ipt && low_trigger_bin == 0) low_trigger_bin = i;
	if (ireg2 <= get_measure_bins() && subleading_measure_region[ireg2] == 0 && ipt >= get_pt2_regions(ireg2)) subleading_measure_region[ireg2++] = i;
	if (ireg <= get_measure_bins() && measure_region[ireg] == 0 && ipt >= get_pt_regions(ireg)) measure_region[ireg++] = i;
	if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
	if (ipt >= reco_max_goal && max_reco_bin == 0) max_reco_bin = i;
      }
    if (max_reco_bin == 0) max_reco_bin = nbins - 1;
    low_trigger[0] = ipt_bins[low_trigger_bin];
    low_trigger[1] = ipt_bins[low_trigger_bin+1];
    low_trigger[2] = ipt_bins[low_trigger_bin+2];

    truth_subleading_bin = 0;
    
    reco_leading_bin = truth_leading_bin + get_reco_leading_bin_buffer();
    measure_leading_bin = truth_leading_bin + get_measure_leading_bin_buffer();
    reco_subleading_bin = truth_subleading_bin + get_reco_subleading_bin_buffer();
    measure_subleading_bin = truth_subleading_bin + get_measure_subleading_bin_buffer();

    truth_leading_cut = ipt_bins[truth_leading_bin];
    reco_leading_cut = ipt_bins[reco_leading_bin];
    measure_leading_cut = ipt_bins[measure_leading_bin];
    truth_subleading_cut = ipt_bins[truth_subleading_bin];
    reco_subleading_cut = ipt_bins[reco_subleading_bin];
    measure_subleading_cut = ipt_bins[measure_subleading_bin];


    return;
  }

  void get_xj_bins(float ixj_bins[])
  {
    Int_t nbins = get_nbins();
    Float_t minimum = get_minimum();
    Float_t fixed = get_fixed();
    Int_t beforebin = get_bbins();

    float alpha = TMath::Power(fixed/minimum, 1/(float)beforebin);
    float maximum = minimum*TMath::Power(alpha, nbins );
    //float bin_xj10 = 1.0;
    float bin_xj1 = 1.0*(minimum/maximum);

    for (int i = 0; i < nbins+1; i++)
      {
	float ixj = bin_xj1*TMath::Power(alpha, (float)i);
	ixj_bins[i] = ixj;
      }
    return;
  }

  
 private:

  std::string jerfile = "jer/jer_smear_functions.root";
  TFile *fjer{nullptr};
  TEnv *penv{nullptr};
  TF1 *f_smear_function{nullptr};

  std::string tntuple_location = ".";
  std::string sim_location = ".";
  std::string jesr_location = ".";
  std::string code_location = ".";
  std::string auau_config = ".";
  
  float sample_boundary_goal[4] = {9.99, 19.99, 29.99, 100};

  float sample_boundary[4] = {0, 0, 0, 100};

  int max_reco_bin = 0;
  float low_trigger[3] = {0};
  int low_trigger_bin = 0;
  int measure_region[10] = {0};
  int subleading_measure_region[10] = {0};

  float truth_leading_cut = 0;
  float truth_subleading_cut = 0;

  float reco_leading_cut = 0;
  float reco_subleading_cut = 0;

  float measure_leading_cut = 0;
  float measure_subleading_cut = 0;

  
  int truth_leading_bin = 0;
  int truth_subleading_bin = 0;

  int reco_leading_bin = 0;
  int reco_subleading_bin = 0;

  int measure_leading_bin = 0;
  int measure_subleading_bin = 0;

  int centrality_bins = 4;  
};
#endif
