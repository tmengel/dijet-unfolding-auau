#include <array>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TStyle.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
#include "PlotUtils.h"
#include "priorReweightQA.h"

// v2 variant of createResponse_exclusive_AA.cxx: reads the pair-level TTree
// "T" written by dijet_pair_matching.C -- one row per EVENT (not per leg),
// already carrying a final pair-level category -- instead of the per-leg
// TTree "T" written by dijet_matching.C. dijet_pair_matching.C classifies
// each event directly as kFill/kMiss/kFake/kSkip/kUESub (DijetPair::Category
// in dijet_pair_matching.C, same numeric values as the enum below), building
// the reco dijet the way the data analysis does -- leading/subleading
// ACCEPTED reco jets by pT, then the dphi cut -- rather than only ever
// looking at each truth leg's own dR match, and it separates out UE
// fluctuations that have taken over a leg (kUESub) instead of silently
// folding them into Fill or Miss. Since the classification is already
// pair-level, there is no per-leg combination step here: category is read
// straight off the row and mapped onto real/miss/fake/skip below, with
// kUESub skipped alongside kSkip (see the file header comment on kUESub in
// dijet_pair_matching.C -- those pairs are covered by the subleading-jet
// efficiency and the inclusive cross-check, so counting them here too would
// double count them). Everything downstream of that -- reweighting,
// closure-test split, smearing, the RooUnfoldResponse/histogram fills, the
// noempty skim, and the closure-test unfold -- is unchanged from
// createResponse_exclusive_AA.cxx and operates purely on the resulting
// histograms/response, agnostic to where they came from.
enum DijetPairCategory { kFill = 0, kMiss = 1, kFake = 2, kSkip = 3, kUESub = 4 , kFakeMiss = 5 };

// Angle of a jet to the 2nd-order event plane, folded into [0, pi/2]. Same
// definition as AnaUtils::get_dpsi2 -- which is what wrote the tree's dpsi2
// branch, but only for the LEADING truth leg -- reimplemented here because
// this macro runs in the ROOT interpreter and does not link libmyana. The
// subleading leg's angle is rebuilt from the psi2 and truth_phi2 branches
// with this, for the jet-v2 cross-check weight below. Folding by pi is
// harmless there: cos(2*dpsi2) is pi-periodic in the underlying dphi.
static double jetv2_dpsi2( const double psi2, const double phi_jet )
{
	auto dphi_wrap = [](double a, double b) -> double
	{
		double dphi = b - a;
		while ( dphi >  TMath::Pi() ) { dphi -= 2.0 * TMath::Pi(); }
		while ( dphi < -TMath::Pi() ) { dphi += 2.0 * TMath::Pi(); }
		return std::fabs(dphi);
	};
	const double psi2_mod = ( psi2 < 0 ) ? psi2 + TMath::Pi() : psi2 - TMath::Pi();
	return std::fabs( std::min( dphi_wrap(phi_jet, psi2), dphi_wrap(phi_jet, psi2_mod) ) );
}

int createResponse_exclusive_v2_AA (
	const std::string configfile = "binning.config",
	const int full_or_half = 0,
	const int niterations = 10,
	const int cone_size = 4,
	const int centrality_bin = 0,
	const int primer = 0,
	const std::string exclusive_dir = "/home/tmengel/PPG14/rootfiles/dijet_match_08_31_2026/exclusive"
)
{

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
  	gStyle->SetOptTitle(0);
 	PlotUtils::set_sphenix_style();

	const bool OVERRIDE_EVENT_WEIGHT = false;
	const bool DO_CENT_EVENT_WEIGHT = false;
	const bool DO_CENT_CUT = false;
	
	const bool DO_FAKES = true;
	const bool VETO_NULL_WEIGHTS = false;
	
	const double prior_var = 0.5;
	const double prior_norm = 1.0;
	// Space in which the prior fraction is applied (see priorReweightQA.h):
	//   kBlendXJ     - blend the xJ shapes and map the ratio back onto the
	//                  (pt1, pt2) bins. One weight per pT-bin separation d, so
	//                  the prior moves in xJ without chasing the per-bin
	//                  statistics of the unfolded data, and the pt1 spectrum is
	//                  left alone.
	//   kBlendPerBin - the original per-(pt1, pt2)-bin blend. Hits the flat
	//                  target exactly, but at fraction 1 it puts weights of up
	//                  to 23x on corner bins with a few hundred raw counts.
	const int prior_blend_space = prior_qa::kBlendPerBin;
	// Nominal Bayesian-iteration index (0-indexed) -> N_iter = 2. The iteration
	// scan minimises the quadrature sum of sigma_sim, sigma_data and
	// sigma_bin-by-bin (see makeIterationPlot_AA.C); using absolute changes, as
	// the AA note specifies, that minimum is N_iter = 2 for every centrality.
	const int prior_iteration = 1;

	const bool not_a_closure_test = (full_or_half == 0);
	const bool half_closure  = (full_or_half == 1);
	const bool full_closure  = (full_or_half == 2);
	std::cout << "not_a_closure_test = " << not_a_closure_test << std::endl;
	std::cout << "half_closure = " << half_closure << std::endl;
	std::cout << "full_closure = " << full_closure << std::endl;
	const bool NO_PRIMER = (primer == 0);
	const bool PRIMER1 = (primer == 1);
	const bool PRIMER2 = (primer == 2);
	std::cout << "PRIMER1 = " << PRIMER1 << std::endl;
	std::cout << "PRIMER2 = " << PRIMER2 << std::endl;
	std::cout << "NO_PRIMER = " << NO_PRIMER << std::endl;
	read_binning rb(configfile.c_str());

	std::string system_string = rb.get_system_string(centrality_bin);
	std::cout << "System string: " << system_string << std::endl;
	
	std::cout << "Using exclusive dijet-matching files in: " << exclusive_dir << std::endl;

	// Flavor-tagged cross-check (qq vs qg/gg leading-dijet parton origin,
	// from dijet_matching_flavor.C) -- 0 selects the plain jetNN_scaled.root
	// files, 1/2 select the qq/qg_gg-suffixed ones written alongside them, 3
	// mixes both flavor files into one response at a controlled QQ share
	// (see get_flavor_qq_fraction() / the normalization pre-pass below).
	// Read once here (used for both the filenames below and sys_name).
	const int flavor_sys = rb.get_flavor_sys();
	const bool flavor_mix = (flavor_sys == 3);
	const double flavor_qq_fraction = rb.get_flavor_qq_fraction();
	const std::string flavor_suffix = (flavor_sys == 1) ? "_qq" : (flavor_sys == 2) ? "_qg_gg" : "";


	
	// Every other mode reads one flavor stream per pT sample (3 streams,
	// stream index == pT-sample index, matching every array below
	// historically). Mix mode instead reads BOTH flavor-tagged files for
	// each of the 3 pT samples (6 streams): streams [0..2] are the qq side
	// and [3..5] the qg/gg side of the same 3 pT samples. stream_sample[]
	// maps a stream back to its pT sample (for n_events/scale_factor and
	// the QA histograms below, which stay indexed by pT sample) and
	// stream_is_qq[] to its flavor (for the normalization pre-pass).
	const int n_streams = flavor_mix ? 6 : 3;
	const int stream_sample[6] = {0, 1, 2, 0, 1, 2};
	const bool stream_is_qq[6]  = {true, true, true, false, false, false};

	std::array<std::string, 6> exclusive_names;
	if (flavor_mix)
	{
		exclusive_names = {
			"jet10_scaled_qq.root", "jet20_scaled_qq.root", "jet30_scaled_qq.root",
			"jet10_scaled_qg_gg.root", "jet20_scaled_qg_gg.root", "jet30_scaled_qg_gg.root"
		};
	}
	else
	{
		exclusive_names = {
			"jet10_scaled" + flavor_suffix + ".root",
			"jet20_scaled" + flavor_suffix + ".root",
			"jet30_scaled" + flavor_suffix + ".root",
			"", "", ""
		};
	}

	// No per-event "weight" branch in this schema -- event weighting here
	// comes entirely from the mbd/sumeT/centrality reweighting below, same
	// as the tn_match chain when it has no "weight" branch either.
	float event_weight[3] = {1, 1, 1};
	bool has_event_weight[3] = {false, false, false};

	// Total processed events per sample (match_standalone.C's raw merged
	// jetNN_scaled.root entry counts) -- the cross-section scale_factor
	// denominator below. Hardcoded rather than read from a file: it's
	// pinned to this specific simulation production, cross-checked against
	// both the raw merged files directly and by counting event_id
	// transitions in the derived exclusive tree, and hardcoding it avoids
	// depending on the (sizeable, and at one point actually corrupted) raw
	// files just for this one number.
	float n_events[3] = { 2892123, 2871303, 2854291 };

	// One row per EVENT (pair-level, not per-leg): truth_pt1/truth_pt2 are
	// the truth-leading/truth-subleading legs (truth_idx1=0, truth_idx2=1
	// in dijet_pair_matching.C -- pT-sorted by construction, so
	// truth_pt1 >= truth_pt2 always), reco_pt1/reco_pt2 the corresponding
	// selected reco legs, valid only when reco_pair is set. truth_lead_pt
	// (dijet_pair_matching_v3.C) is the event's hardest truth jet's pT
	// regardless of category/acceptance -- unlike truth_pt1, it is not
	// gated on a truth dijet candidate existing at all, so it is what the
	// pT-hat sample-boundary veto below keys on. Sized for the mix-mode
	// 6-stream case; non-mix modes only ever touch indices [0..2].
	int   b_cent[6];
	float b_zvrtx[6];
	float b_sumeT[6];
	int   b_category[6];
	int   b_reco_pair[6];
	float b_truth_pt1[6];
	float b_truth_pt2[6];
	float b_truth_lead_pt[6];
	float b_reco_pt1[6];
	float b_reco_pt2[6];

	float b_psi2[6];
	float b_truth_phi1[6];
	float b_truth_phi2[6];
	float b_dpsi2[6];
	int b_truth_match_idx1[6];
	int b_truth_match_idx2[6];

	TFile * fin[6];
	TTree * texcl[6];
	for (int i = 0 ; i < n_streams; i++)
	{
		const std::string path = exclusive_dir + "/" + exclusive_names[i];
		fin[i] = new TFile(path.c_str(), "READ");
		if (!fin[i] || fin[i]->IsZombie())
		{
			std::cerr << "Cannot open exclusive dijet-matching file " << path << std::endl;
			return 1;
		}

		texcl[i] = (TTree*) fin[i]->Get("T");
		if (!texcl[i])
		{
			std::cerr << "Missing tree T in " << path << std::endl;
			return 1;
		}

		texcl[i]->SetBranchAddress("cent", &b_cent[i]);
		texcl[i]->SetBranchAddress("zvrtx", &b_zvrtx[i]);
		texcl[i]->SetBranchAddress("sumeT", &b_sumeT[i]);
		texcl[i]->SetBranchAddress("dpsi2", &b_dpsi2[i]);
		texcl[i]->SetBranchAddress("category", &b_category[i]);
		texcl[i]->SetBranchAddress("reco_pair", &b_reco_pair[i]);
		texcl[i]->SetBranchAddress("truth_pt1", &b_truth_pt1[i]);
		texcl[i]->SetBranchAddress("truth_pt2", &b_truth_pt2[i]);
		texcl[i]->SetBranchAddress("truth_lead_pt", &b_truth_lead_pt[i]);
		texcl[i]->SetBranchAddress("reco_pt1", &b_reco_pt1[i]);
		texcl[i]->SetBranchAddress("reco_pt2", &b_reco_pt2[i]);
		texcl[i]->SetBranchAddress("psi2", &b_psi2[i]);
		texcl[i]->SetBranchAddress("truth_match_idx1", &b_truth_match_idx1[i]);
		texcl[i]->SetBranchAddress("truth_match_idx2", &b_truth_match_idx2[i]);
		texcl[i]->SetBranchAddress("truth_phi1", &b_truth_phi1[i]);
		texcl[i]->SetBranchAddress("truth_phi2", &b_truth_phi2[i]);
	}
	std::cout << "n_events (hardcoded): " << n_events[0] << " " << n_events[1] << " " << n_events[2] << std::endl;
	std::cout << "has_event_weight: " << has_event_weight[0] << " " << has_event_weight[1] << " " << has_event_weight[2] << std::endl;

	const float cs_10 = 0.000003997;
  	const float cs_20 = 6.218e-8;
  	const float cs_30 = 2.505e-9;
  
	float scale_factor[3];
	scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
	scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30; 
	scale_factor[2] = 1.0;
	std::cout << "Scale factors: " << scale_factor[0] << " " << scale_factor[1] << " " << scale_factor[2] << std::endl;

	const int minentries = rb.get_minentries();
	const int minentries_truth = rb.get_minentries_truth();
	const int minentries_link = rb.get_minentries_link();
	std::cout << "Minimum entries: " << minentries
			  << "  (truth-side: " << minentries_truth
			  << ", truth-reco link: " << minentries_link << ")" << std::endl;

	const int nbins = rb.get_nbins();
	float ipt_bins[nbins+1];
	float ixj_bins[nbins+1];
	rb.get_pt_bins(ipt_bins);
	rb.get_xj_bins(ixj_bins);
	std::cout << "Pt bins -- Xj bins: " << std::endl;
	for (int i = 0 ; i < nbins + 1; i++)
	{
		std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
	}

	const int centrality_bins = rb.get_number_centrality_bins();
	float icentrality_bins[centrality_bins + 1];
	rb.get_centrality_bins(icentrality_bins);
	std::cout << "Centrality bins: " << std::endl;
	for (int i = 0 ; i < centrality_bins + 1; i++)
	{
		std::cout << icentrality_bins[i] << std::endl;
	}

	const int max_reco_bin = rb.get_maximum_reco_bin();
	std::cout << "Max reco bin: " << max_reco_bin << std::endl;

	const double dphicut = rb.get_dphicut();
	const double dphicuttruth = dphicut;
	std::cout << "dphicut: " << dphicut << std::endl;
	std::cout << "dphicuttruth: " << dphicuttruth << std::endl;

	const float truth_leading_cut = rb.get_truth_leading_cut();
	const float truth_subleading_cut = rb.get_truth_subleading_cut();
	std::cout << "Truth leading cut: " << truth_leading_cut << std::endl;
	std::cout << "Truth subleading cut: " << truth_subleading_cut << std::endl;

	const float reco_leading_cut = rb.get_reco_leading_cut();
	const float reco_subleading_cut = rb.get_reco_subleading_cut();
	std::cout << "Reco leading cut: " << reco_leading_cut << std::endl;
	std::cout << "Reco subleading cut: " << reco_subleading_cut << std::endl;

	const float measure_leading_cut = rb.get_measure_leading_cut();
	const float measure_subleading_cut = rb.get_measure_subleading_cut();
	std::cout << "Measure leading cut: " << measure_leading_cut << std::endl;
	std::cout << "Measure subleading cut: " << measure_subleading_cut << std::endl;

	const int measure_leading_bin = rb.get_measure_leading_bin();
	const int measure_subleading_bin = rb.get_measure_subleading_bin();
	std::cout << "Measure leading bin: " << measure_leading_bin << std::endl;
	std::cout << "Measure subleading bin: " << measure_subleading_bin << std::endl;

	float sample_boundary[4] = {0};
	std::cout << "Sample boundaries: " << std::endl;
	for (int ib = 0; ib < 4; ib++)
	{
		sample_boundary[ib] = rb.get_sample_boundary(ib);
		std::cout <<  sample_boundary[ib] << std::endl;
	}


	const double flow_v22_scale = rb.get_flow_sys();
	const double flow_v33_scale = rb.get_flow_v33_sys();
	const bool flow_sys = std::fabs(flow_v22_scale - 1.0) > 1e-6 || std::fabs(flow_v33_scale - 1.0) > 1e-6;
	
	
	// Jet-v2 cross-check: HIJING+Pythia embeds the signal jets with no azimuthal
	// correlation to the event plane (the tree's dpsi2 is flat, mean pi/4), so a
	// hypothetical jet v2 is injected here by reweighting each MC pair by
	// 1 + amp*(cos(2*dpsi2_leg1) + cos(2*dpsi2_leg2)), with amp fixed by the
	// normalization pre-pass below so that each leg's realized v2 is exactly
	// JETV2_SCALE. That changes which part of the modulated UE each jet sits on, so it moves
	// the reco legs' background subtraction and therefore the response itself --
	// which is the whole point of running it here rather than as a flat scaling of
	// the final spectrum. Like the flavor cross-check this is a "what if", not a
	// measured up/down uncertainty, so it is deliberately NOT folded into
	// drawSys_AA.C's total band; it only needs a distinct output name.
	//
	// JETV2_SCALE is the jet v2 itself (0.03 for the 3% cross-check), NOT a
	// multiplier -- the 1.0 default is just an out-of-range sentinel meaning
	// "off", matching how FLOW_V22_SCALE flags itself.
	const double jetv2_scale = rb.get_jetv2_scale();
	const bool jetv2_sys = std::fabs(jetv2_scale - 1.0) > 1e-6;
	if ( jetv2_sys )
	{
		std::cout << "JetV2 cross-check ACTIVE: injecting jet v2 = " << jetv2_scale
		          << " (weight 1 + 2*v2*cos(2*dpsi2) per truth leg)" << std::endl;
		if ( !(jetv2_scale > 0.0) || jetv2_scale >= 0.5 )
		{
			std::cerr << "JETV2_SCALE must be a jet v2 in (0, 0.5) -- got " << jetv2_scale
			          << ". Outside that the weight 1 + 2*v2*cos(2*dpsi2) goes negative." << std::endl;
			return 1;
		}
	}
	else
	{
		std::cout << "JetV2 cross-check off (JETV2_SCALE = " << jetv2_scale << ")" << std::endl;
	}
	
	const int inclusive_sys = rb.get_inclusive_sys();


	const double JES_sys = rb.get_jes_sys();
	const double JER_sys = rb.get_jer_sys();
	const bool prior_sys = rb.get_prior_sys();

	std::string sys_name = "nominal";
	std::cout << "JES = " << JES_sys << std::endl;
	std::cout << "JER = " << JER_sys << std::endl;
	
	int using_sys = 0;

	auto * f_smear = (TF1*) rb.get_smear_function(centrality_bin);
	if ( !f_smear )
	{
		std::cout << " NO SMEAR " << std::endl;
		exit(-1);
	}
	
	// Fractional-pT smearing: a Gaussian of width f_smear->Eval(pt), truncated to
	// [-0.5, 0.5]. This used to be a TF1("gaus") with SetParameter(2, sigma)
	// followed by GetRandom(). SetParameter calls TF1::Update(), which throws away
	// the cumulative-distribution table, so GetRandom rebuilt a 100-point numerical
	// CDF on EVERY call: ~8 us per smear, ~6M smears per invocation, ~50 s of the
	// event loop. Below is the closed-form inverse CDF of the same truncated
	// Gaussian -- ~0.07 us per smear, and exact rather than TF1's 100-bin
	// interpolation. It draws exactly one uniform per call, as GetRandom did.
	//
	// JER systematics enter only through f_smear (get_smear_function returns
	// f_positive/f_negative for JER_sys != 0). The old `width = 0.8 + JER_sys`
	// was dead: SetParameter(2, ...) overwrote it before every single draw.
	const double smear_trunc = 0.5;
	auto smear_random = [&](double pt) -> double
	{
		const double sigma = f_smear->Eval(pt);
		if ( !(sigma > 0) ) { return 0.0; }
		// P(|x| < trunc) for N(0, sigma), via Phi(trunc/sigma).
		const double phi_hi = TMath::Freq(smear_trunc / sigma);
		double u = (1.0 - phi_hi) + gRandom->Rndm() * (2.0 * phi_hi - 1.0);
		// NormQuantile is only defined on the open interval; sigma << trunc makes
		// phi_hi round to exactly 1 and u can then reach the endpoints.
		if ( u <= 0.0 ) { u = std::numeric_limits<double>::epsilon(); }
		if ( u >= 1.0 ) { u = 1.0 - std::numeric_limits<double>::epsilon(); }
		return sigma * TMath::NormQuantile(u);
	};

	if ( prior_sys )
	{
		using_sys = 1;
		sys_name = "PRIOR";
	}
	if ( flow_sys )
	{
		using_sys = 1;
		sys_name = rb.get_flow_systematic_name();
	}
	if (inclusive_sys)
	{
		using_sys = 1;
		sys_name = "INCLUSIVE";
	}
	if ( jetv2_sys )
	{
		using_sys = 1;
		sys_name = rb.get_jetv2_systematic_name();
	}

	// Flavor-tagged cross-check (qq vs qg/gg leading-dijet parton origin) --
	// see dijet_matching_flavor.C. Not a JES/JER-style up/down uncertainty,
	// so it's deliberately not folded into drawSys_AA.C's total systematic
	// band; this only needs to pick the right input files (below) and give
	// the response/unfold chain a distinct, non-colliding output name.
	if (flavor_sys == 1)
	{
		using_sys = 1;
		sys_name = "QQ";
	}
	else if (flavor_sys == 2)
	{
		using_sys = 1;
		sys_name = "QGGG";
	}
	else if (flavor_sys == 3)
	{
		using_sys = 1;
		// Percent-QQ tag so e.g. 0.5 and 0.2 land on distinct, non-colliding
		// response/unfold/QA output names (MIX50, MIX20, ...).
		sys_name = Form("MIX%02d", (int) std::lround(flavor_qq_fraction * 100.0));
	}
	if (JER_sys != 0)
	{
		using_sys = 1;
		if (JER_sys < 0) sys_name = "negJER";
		else if (JER_sys > 0) sys_name = "posJER";
		std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
	}
	if (JES_sys != 0)
	{
		using_sys = 1;
		if (JES_sys < 0) sys_name = "negJES";
		else if (JES_sys > 0) sys_name = "posJES";
		std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
	}

  	const std::string diagnostic_name = half_closure ? "HALF_" + sys_name : (full_closure ? "FULL_" + sys_name : sys_name);

	TH1D * h_centrality_reweight = nullptr;
	TH1D * h_mbd_reweight 		 = nullptr;
	TH1D * h_sumeT_reweight 	 = nullptr;
	
	if ( !PRIMER1 || OVERRIDE_EVENT_WEIGHT )
	{
		std::cout << "Loading centrality, vertex, and sumeT reweighting histograms for " << system_string << " " << sys_name << std::endl;
		
		auto * fcent = new TFile(Form("centrality/centrality_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
		if ( !fcent || fcent->IsZombie() )
		{
			std::cerr << "Missing required centrality reweight for " << system_string << " " << sys_name << std::endl;
			if ( !PRIMER1 ) { return 1; }
		}
		else
		{
			h_centrality_reweight = (TH1D*) fcent->Get("h_centrality_reweight") -> Clone(Form("h_centrality_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
			h_centrality_reweight -> SetDirectory(0);
			fcent -> Close();
		}
		
		auto * fvtx = new TFile(Form("vertex/vertex_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
		if ( !fvtx || fvtx->IsZombie()  && !PRIMER1 )
		{
			std::cerr << "Missing required vertex reweight for " << system_string << " "<< sys_name << std::endl;
			if ( !PRIMER1 ) { return 1; }
		}
		else
		{
			h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight") -> Clone(Form("h_mbd_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
			h_mbd_reweight -> SetDirectory(0);
			fvtx -> Close();
		}
		
		auto * fsumeT = new TFile(Form("sumeT/sumeT_reweight_%s_r%02d_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
		if ( !fsumeT || fsumeT->IsZombie() && !PRIMER1 )
		{
			std::cerr << "Missing required sumeT reweight for " << system_string << " "<< sys_name << std::endl;
			if ( !PRIMER1 ) { return 1; }
		}
		else
		{
			h_sumeT_reweight = (TH1D*) fsumeT->Get("h_sumeT_reweight") -> Clone(Form("h_sumeT_reweight_%s_r%02d_%s", system_string.c_str(), cone_size, sys_name.c_str()));
			h_sumeT_reweight -> SetDirectory(0);
			fsumeT -> Close();
		}
	}
	
  	TH1D * h_truth_lead_sample[3];
	for (int i = 0; i < 3; i++)
	{
		h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	}
  
	TH1D *h_truth_lead = new TH1D("h_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_truth_sublead = new TH1D("h_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_reco_lead = new TH1D("h_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_reco_sublead = new TH1D("h_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);

	TH1D *h_match_truth_lead = new TH1D("h_match_truth_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_truth_sublead = new TH1D("h_match_truth_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_reco_lead = new TH1D("h_match_reco_lead", " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
	TH1D *h_match_reco_sublead = new TH1D("h_match_reco_sublead", " ; Subleading Jet p_{T} [GeV]; counts", 100, 0, 100);
	
	TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
	TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
	TH1D *h_sumeT = new TH1D("h_sumeT", ";#Sigma E_{T}; counts", 200, 0, 2000);

	// Debug: distribution of the jet-v2 cross-check weight
	// (1 + jetv2_amp*(c1 + c2)) actually applied per pair, unweighted --
	// only filled when the jetv2 cross-check is active. Mean should sit at
	// 1 and the spread should be consistent with jetv2_amp*sqrt(2), the
	// same closure the console printout at the end of the event loop checks
	// numerically.
	TH1D *h_jetv2_weight = new TH1D("h_jetv2_weight", ";jet-v2 cross-check weight;counts", 200, 0.5, 1.5);

	// pure fills
	TH1D *h_truth_xj = new TH1D("h_truth_xj",";A_{J};1/N", nbins, ixj_bins);
	TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

	TH1D *h_linear_truth_xj = new TH1D("h_lineartruth_xj",";A_{J};1/N", 20, 0, 1.0);
	TH1D *h_linear_reco_xj = new TH1D("h_linearreco_xj",";A_{J};1/N", 20, 0, 1.0);

	TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);
	TH2D *h_e1e2 = new TH2D("h_e1e2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

	TH1D *h_flat_truth_pt1pt2 = new TH1D("h_truth_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_count_flat_truth_pt1pt2 = new TH1D("h_truth_count_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_flat_truth_to_response_pt1pt2 = new TH1D("h_truth_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

	TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_count_flat_reco_pt1pt2 = new TH1D("h_count_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
	TH1D *h_flat_reco_to_response_pt1pt2 = new TH1D("h_reco_flat_to_response_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

	TH2D *h_flat_response_pt1pt2 = new TH2D("h_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);
	// Unweighted response fills, used by the minentries_link truth-bin criterion.
	TH2D *h_count_flat_response_pt1pt2 = new TH2D("h_count_flat_response_pt1pt2",";p_{T,1, reco} + p_{T,2, reco};p_{T,1, truth} + p_{T,2, truth}", nbins*nbins, 0, nbins*nbins, nbins*nbins, 0, nbins*nbins);

	// For the prior sensitivity. The weights themselves are built (and QA'd) by
	// prior_qa so that testPriorReweight_AA.C exercises the identical code path.
	prior_qa::weights prior_weights;
	TH1D * h_flatreweight_pt1pt2 = nullptr;
	auto * h2_flatreweight_pt1pt2 = new TH2D("h2_flatreweight_pt1pt2",";p_{T,1}^{truth} [GeV]; p_{T,2}^{truth} [GeV] ; Reweight Factor", nbins, ipt_bins, nbins, ipt_bins);

	// Prior-sensitivity configuration.
	// fraction = 0.0 : keep the Pythia(+HIJING) prior (all weights exactly 1)
	// fraction = 0.5 : reweight the prior halfway toward the unfolded data
	// fraction = 1.0 : full data-driven prior reweighting (preliminary p+p style)
	//
	// Nominal uses 0.0: every result produced and presented for the Run-24
	// to 1 (this matches the committed Jul 30 value prior_fraction = 0.0).
	// Turning the full reweighting on with the quenched Au+Au data produces a
	// runaway feedback (weights 0.08-1.3) and a nominal far outside the band
	// of everything reviewed; enabling it is a physics decision, not a bugfix.
	//
	// The PRIOR systematic is then |unfold(fraction 0.5) - unfold(fraction 0)|,
	// i.e. half of the preliminary full-vs-none prior excursion, implementing
	// the intended "reduce the prior-reweighting systematic by 50%" without
	// touching the nominal.
	//
	// Ternary order matters: the first value applies to the PRIOR-systematic
	// config (PRIOR: 1), the second to the nominal. It must be 0.5 : 0.0 —
	// writing 0.0 : 0.5 makes the NOMINAL half-reweighted (central value moves
	// halfway to the steep solution) and the PRIOR variation unreweighted.
	const double prior_fraction = prior_sys ? prior_var : prior_norm;
	std::cout << "Prior fraction: " << prior_fraction << std::endl;

	
	const int nbin_response = nbins*nbins;
  	RooUnfoldResponse rooResponse(nbin_response, 0, nbin_response);

	if ( not_a_closure_test && NO_PRIMER )
	{

		std::cout << "doing prior" << std::endl;

		// Flavor cross-check (QQ/QGGG): data itself is never flavor-tagged, so
		// "the data" the flavor-tagged truth should be pulled toward is the
		// nominal unfolded pt1pt2, not the QQ/QGGG-response-unfolded version of
		// the same data -- that would just reintroduce the flavor-tagged
		// response's own migration pattern into the target. The truth below
		// stays flavor-specific; only the unfold/data target changes.
		const std::string prior_data_sys_name = (flavor_sys != 0) ? "nominal" : sys_name;

		auto * fun = new TFile(Form("unfolding_hists/unfolding_hists_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, prior_data_sys_name.c_str()), "READ");
		if ( !fun || fun->IsZombie() )
		{
			std::cerr << "Missing required unfolding histograms for " << system_string << " "<< prior_data_sys_name << std::endl;
			return 1;
		}
		auto * h_unfold_flat_in = (TH1D*) fun->Get(Form("h_flat_unfold_pt1pt2_%d", prior_iteration)) -> Clone(Form("h_flat_unfold_pt1pt2_%d", prior_iteration));
		h_unfold_flat_in -> SetDirectory(0);
		fun -> Close();

		auto * ftr = new TFile(Form("response_matrices/response_matrix_%s_r%02d_PRIMER2_%s.root", system_string.c_str(), cone_size, sys_name.c_str()), "READ");
		if ( !ftr || ftr->IsZombie() )
		{
			std::cerr << "Missing required response matrix for " << system_string << " "<< sys_name << std::endl;
			return 1;
		}
		auto * h_truth_flat_in = (TH1D*) ftr->Get("h_truth_flat_pt1pt2") -> Clone(Form("h_truth_flat_pt1pt2_%d", prior_iteration));
		h_truth_flat_in -> SetDirectory(0);
		ftr -> Close();

		// Normalizes both inputs, blends them by prior_fraction and forms the
		// per-truth-bin weights in ROOT bin k+1. Projection onto xJ is linear, so
		// prior_weights.target is the same blended target either way -- only the
		// weights differ.
		std::cout << "Prior blend space: " << prior_qa::blend_name(prior_blend_space) << std::endl;
		// One xJ ratio per leading-pT measurement range, so the reweighted prior
		// reproduces each range's xJ target and not just the inclusive one.
		std::vector<prior_qa::lead_group> prior_lead_groups;
		for (int irange = 0; irange < rb.get_measure_bins(); irange++)
		{
			prior_lead_groups.push_back({rb.get_measure_region(irange), rb.get_measure_region(irange + 1)});
		}

		prior_weights = ( prior_blend_space == prior_qa::kBlendXJ )
			? prior_qa::build_xj(h_truth_flat_in, h_unfold_flat_in, prior_fraction, nbins,
			                     ipt_bins, ixj_bins, prior_lead_groups,
			                     measure_subleading_bin, nbins - 2)
			: prior_qa::build(h_truth_flat_in, h_unfold_flat_in, prior_fraction, nbins);
		h_flatreweight_pt1pt2 = prior_weights.weight;
		if ( !h_flatreweight_pt1pt2 )
		{
			std::cerr << "Failed to build the prior reweighting." << std::endl;
			return 1;
		}

		TH1D * h_truth_flat = prior_weights.truth;
		TH1D * h_unfold_flat = prior_weights.unfold;
		TH1D * h_fractional_unfold_flat = prior_weights.target;

		// 2D view of the weights. Flat index k = ibin-1 maps to (k/nbins, k%nbins);
		// using ibin itself here shifted the whole matrix by one truth bin.
		for (int ibin = 1; ibin <= nbin_response; ibin++)
		{
			const int flat_index = ibin - 1;
			const int gbin = h2_flatreweight_pt1pt2 -> GetBin( flat_index/nbins + 1, flat_index%nbins + 1 );
			h2_flatreweight_pt1pt2->SetBinContent( gbin, h_flatreweight_pt1pt2->GetBinContent(ibin) );
		}

		TCanvas *cre = new TCanvas("cre","cre", 500, 500);
		cre->SetLeftMargin(0.1);
		cre->SetRightMargin(0.19);
		// h2_flatreweight_pt1pt2->GetYaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
		// h2_flatreweight_pt1pt2->GetXaxis()->SetRangeUser(ipt_bins[0], ipt_bins[13]);
		h2_flatreweight_pt1pt2->Draw("colz");
		dlutility::DrawSPHENIX(0.2, 0.87);
		dlutility::drawText("Prior Reweighting Matrix", 0.2, 0.77);
		dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.2, 0.72);
		// cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
		// need a primer and sys and diagnostic name to save the prior matrix
		cre->Print(Form("%s/unfolding_plots/prior_matrix_%s_r%02d_%s_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str(), diagnostic_name.c_str()));
		// draw 1D truth, unfold, and fractional unfold
		TCanvas *c1d = new TCanvas("c1d","c1d", 500, 500);
		c1d->SetLeftMargin(0.1);
		c1d->SetRightMargin(0.19);
		gPad->SetLogy(0);
		h_truth_flat->SetLineColor(kBlack);
		h_truth_flat->SetLineWidth(2);
		h_unfold_flat->SetLineColor(kRed);
		h_unfold_flat->SetLineWidth(2);
		h_fractional_unfold_flat->SetLineColor(kBlue);
		h_fractional_unfold_flat->SetLineWidth(2);
		h_truth_flat->Draw("hist");
		h_unfold_flat->Draw("hist same");
		h_fractional_unfold_flat->Draw("hist same");
		dlutility::DrawSPHENIX(0.2, 0.87);
		dlutility::drawText("Prior Reweighting", 0.2, 0.77);
		dlutility::drawText(Form("%d - %d %%", (int) icentrality_bins[centrality_bin], (int) icentrality_bins[centrality_bin+1]), 0.2, 0.72);
		// c1d->Print(Form("%s/unfolding_plots/prior_1D_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str()));
		c1d->Print(Form("%s/unfolding_plots/prior_1D_%s_r%02d_%s_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, sys_name.c_str(), diagnostic_name.c_str()));
		gPad->SetLogy(0);

		// ---------------------------------------------------------------------
		// Prior QA: apply the weights to the MC truth exactly as the event loop
		// below will, and compare the result with the target the reweighting is
		// supposed to reach. "reweighted / target" must be 1 in every bin; any
		// structure there is a bug in the reweighting, not physics. The same
		// three distributions are also projected onto xJ so the effect is visible
		// in the observable the systematic is actually quoted on.
		// ---------------------------------------------------------------------
		{
			TH1D * h_prior_reweighted_truth = prior_qa::apply(prior_weights.truth, h_flatreweight_pt1pt2,
			                                                  nbins, "h_prior_reweighted_truth");
			const std::string qa_label = Form("%s r%02d %s %s blend=%s", system_string.c_str(),
			                                  cone_size, sys_name.c_str(), diagnostic_name.c_str(),
			                                  prior_qa::blend_name(prior_blend_space));
			prior_qa::print_summary(prior_weights, h_prior_reweighted_truth, nbins, qa_label);
			// The number that matters for kBlendXJ: flat-space closure is not
			// expected there, xJ closure is.
			prior_qa::print_xj_summary(prior_weights, h_prior_reweighted_truth, nbins,
			                           ipt_bins, ixj_bins, measure_leading_bin, nbins - 2,
			                           measure_subleading_bin, nbins - 2, qa_label);

			const std::string cent_label = Form("%d - %d %%", (int) icentrality_bins[centrality_bin],
			                                                  (int) icentrality_bins[centrality_bin+1]);
			const std::vector<std::string> captions = {
				Form("Prior fraction %.2f", prior_fraction),
				Form("blended in %s", prior_blend_space == prior_qa::kBlendXJ
				                      ? "#it{x}_{J}" : "(#it{p}_{T,1}, #it{p}_{T,2})"),
				cent_label
			};

			prior_qa::draw_flat(prior_weights, h_prior_reweighted_truth, captions,
				Form("%s/unfolding_plots/priorQA_flat_%s_r%02d_%s_%s.pdf", rb.get_code_location().c_str(),
				     system_string.c_str(), cone_size, sys_name.c_str(), diagnostic_name.c_str()));

			const int mbins = rb.get_measure_bins();
			for (int irange = 0; irange < mbins; irange++)
			{
				const int lead_lo = rb.get_measure_region(irange);
				const int lead_hi = rb.get_measure_region(irange + 1);
				std::vector<std::string> xj_captions = captions;
				xj_captions.push_back(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV",
				                           ipt_bins[lead_lo], ipt_bins[lead_hi]));
				prior_qa::draw_xj(prior_weights, h_prior_reweighted_truth, nbins, ipt_bins, ixj_bins,
				                  lead_lo, lead_hi, measure_subleading_bin, nbins - 2, xj_captions,
					Form("%s/unfolding_plots/priorQA_xj_%s_r%02d_%s_%s_range%d.pdf", rb.get_code_location().c_str(),
					     system_string.c_str(), cone_size, sys_name.c_str(), diagnostic_name.c_str(), irange));
			}
		}
	}


	// Mix-mode normalization: the two flavor streams for the same pT sample
	// are disjoint subsets of the same underlying MC and share scale_factor[]
	// -- but Pythia+HIJING doesn't hand out QQ vs QG+GG leading-dijet-flavor
	// events in any controlled ratio, so hitting an exact target QQ share of
	// the response's cross-section-weighted yield needs one extra factor per
	// flavor. This is a fast pre-pass over the same trees (same in-range /
	// non-Skip/UESub cuts the main loop below applies, no smearing or
	// histogram fills) that sums each flavor's raw weighted yield across all
	// 3 pT samples, then solves for the two factors that (a) hit
	// flavor_qq_fraction exactly and (b) together reproduce the same total
	// yield as summing both flavors unweighted -- so e.g. MIX50 sits at the
	// same absolute scale as combining the QQ-only and QGGG-only samples,
	// not an arbitrarily rescaled one.
	double extra_scale[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	if (flavor_mix)
	{
		double yield_qq = 0.0;
		double yield_qg_gg = 0.0;
		for (int istream = 0; istream < n_streams; istream++)
		{
			const int isample = stream_sample[istream];
			const Long64_t entries_norm = texcl[istream]->GetEntries();
			double stream_yield = 0.0;
			for (Long64_t i = 0; i < entries_norm; i++)
			{
				texcl[istream] -> GetEntry(i);
				const int cat0 = b_category[istream];
				if (cat0 == kSkip || cat0 == kUESub) { continue; }
				const float truth_lead_pt0 = b_truth_lead_pt[istream];
				if ( !(truth_lead_pt0 >= sample_boundary[isample] && truth_lead_pt0 < sample_boundary[isample+1]) ) { continue; }
				stream_yield += scale_factor[isample];
			}
			if (stream_is_qq[istream]) { yield_qq += stream_yield; } else { yield_qg_gg += stream_yield; }
		}
		std::cout << "Flavor mix normalization -- raw weighted yield: QQ = " << yield_qq
		          << ", QG+GG = " << yield_qg_gg << std::endl;
		if ( !(yield_qq > 0.0) || !(yield_qg_gg > 0.0) )
		{
			std::cerr << "Flavor mix requested (QQ fraction " << flavor_qq_fraction
			          << ") but one flavor has zero weighted yield in this centrality/cone-size bin -- cannot mix."
			          << std::endl;
			return 1;
		}
		const double yield_total = yield_qq + yield_qg_gg;
		const double scale_qq = flavor_qq_fraction * yield_total / yield_qq;
		const double scale_qg_gg = (1.0 - flavor_qq_fraction) * yield_total / yield_qg_gg;
		std::cout << "Flavor mix target QQ fraction " << flavor_qq_fraction
		          << " -> extra per-event scale: QQ x" << scale_qq << ", QG+GG x" << scale_qg_gg << std::endl;
		for (int istream = 0; istream < n_streams; istream++)
		{
			extra_scale[istream] = stream_is_qq[istream] ? scale_qq : scale_qg_gg;
		}
	}

	// Jet-v2 normalization pre-pass (same idea as the flavor-mix pre-pass above:
	// one cheap extra sweep over the same trees under the same in-range /
	// non-Skip/UESub cuts, no smearing or histogram fills, solving for a scale the
	// main loop then just applies).
	//
	// The pair weight is w = 1 + amp*(c1 + c2), with ci = cos(2*dpsi2_i) for each
	// truth leg. amp is NOT 2*v2: the two legs of a back-to-back dijet are one
	// azimuthal degree of freedom, not two independent ones, and cos(2*dphi) is
	// pi-periodic, so c2 ~ c1 and each leg picks up the OTHER leg's modulation as
	// well as its own. Reweighting with the textbook single-particle form
	// (1 + 2*v2*c1)(1 + 2*v2*c2) therefore lands at a realized per-jet v2 of
	// ~1.8*v2 -- 5.3% for a requested 3% in the jet10 sample -- which would put the
	// wrong number on the cross-check. Solving for amp instead makes the realized
	// v2 equal JETV2_SCALE by construction:
	//
	//   <c1>_w = amp * ( <c1^2> + <c1*c2> ) = v2   ->   amp = v2 / (<c1^2> + <c1*c2>)
	//
	// and by the c1 <-> c2 symmetry of both the sample and the weight, <c2>_w = v2
	// as well, so both legs come out at the requested v2. <w> = 1 exactly, since
	// <c1> = <c2> = 0 in the unweighted MC (HIJING+Pythia embeds the signal jets
	// flat in azimuth relative to Psi_2). The expectation values are taken
	// cross-section weighted, matching how the response is actually filled.
	double jetv2_amp = 0.0;
	if ( jetv2_sys )
	{
		double sum_w = 0.0, sum_c11 = 0.0, sum_c12 = 0.0, sum_c1 = 0.0;
		for (int istream = 0; istream < n_streams; istream++)
		{
			const int isample = stream_sample[istream];
			const Long64_t entries_norm = texcl[istream]->GetEntries();
			for (Long64_t i = 0; i < entries_norm; i++)
			{
				texcl[istream] -> GetEntry(i);
				const int cat0 = b_category[istream];
				if (cat0 == kSkip || cat0 == kUESub) { continue; }
				const float truth_lead_pt0 = b_truth_lead_pt[istream];
				if ( !(truth_lead_pt0 >= sample_boundary[isample] && truth_lead_pt0 < sample_boundary[isample+1]) ) { continue; }
				const double w = scale_factor[isample] * extra_scale[istream];
				const double c1 = std::cos( 2.0 * b_dpsi2[istream] );
				const double c2 = std::cos( 2.0 * jetv2_dpsi2( b_psi2[istream], b_truth_phi2[istream] ) );
				sum_w += w;
				sum_c1 += w * c1;
				sum_c11 += w * c1 * c1;
				sum_c12 += w * c1 * c2;
			}
		}
		if ( !(sum_w > 0.0) )
		{
			std::cerr << "JetV2 cross-check requested but no pairs survive the in-range cuts "
			          << "in this centrality/cone-size bin -- cannot set the v2 amplitude." << std::endl;
			return 1;
		}
		const double mean_c1 = sum_c1 / sum_w;
		const double response_c = ( sum_c11 + sum_c12 ) / sum_w;
		if ( !(response_c > 0.0) )
		{
			std::cerr << "JetV2 cross-check: <c1^2> + <c1*c2> = " << response_c
			          << " is not positive -- cannot solve for the v2 amplitude." << std::endl;
			return 1;
		}
		jetv2_amp = jetv2_scale / response_c;
		std::cout << "JetV2 normalization -- unweighted <cos(2 dpsi2)> = " << mean_c1
		          << " (should be ~0: no jet v2 in the input MC)" << std::endl;
		std::cout << "JetV2 normalization -- <c1^2> + <c1 c2> = " << response_c
		          << " -> weight amplitude = " << jetv2_amp
		          << " (naive 2*v2 would be " << 2.0 * jetv2_scale << ")" << std::endl;
		if ( std::fabs(mean_c1) > 0.01 )
		{
			std::cout << "  WARNING: the input MC already carries a jet v2 of ~" << mean_c1
			          << "; the injected v2 adds to it rather than replacing it." << std::endl;
		}
		if ( 2.0 * std::fabs(jetv2_amp) >= 1.0 )
		{
			std::cerr << "JetV2 weight amplitude " << jetv2_amp
			          << " makes 1 + amp*(c1 + c2) go negative -- v2 = " << jetv2_scale
			          << " is too large for this sample." << std::endl;
			return 1;
		}
	}

	// Fixed seeds: gRandom drives smear_random() — 4357 is ROOT's
	// TRandom3 default, pinned here so nominal/systematic chains stay event-aligned
	// (common random numbers) even if ROOT's default changes. rng drives only the
	// closure-test split; a fixed seed makes closure tests reproducible.
	gRandom->SetSeed(4357);
	auto * rng = new TRandom(12345);
	// Unweighted pair-level counts (post sample-boundary/centrality/reweight
	// cuts), for a quick sanity check against dijet_pair_matching.C's own
	// Fill/Miss/Fake/Skip/UESub summary.
	long n_real[3] = {0, 0, 0};
	long n_miss[3] = {0, 0, 0};
	long n_fake[3] = {0, 0, 0};
	long n_skip[3] = {0, 0, 0};
	long n_uesub[3] = {0, 0, 0};
	// Realized jet-v2 weight, reported after the loop as a closure check on the
	// normalization pre-pass: <w> should be 1 and the weighted <cos(2 dpsi2)> of
	// EACH leg should come back at JETV2_SCALE.
	double jetv2_weight_sum = 0.0;
	double jetv2_c1_sum = 0.0;
	double jetv2_c2_sum = 0.0;
	long   jetv2_weight_n = 0;
  	for (int istream = 0; istream < n_streams; istream++)
    {
		const int isample = stream_sample[istream];
		const Long64_t entries2 = texcl[istream]->GetEntries();
		const Long64_t print_every = (entries2/10 > 0) ? entries2/10 : 1;

		for (Long64_t i = 0; i < entries2; i++)
		{
			texcl[istream] -> GetEntry(i);

			const int cat0 = b_category[istream];
			
			const float truth_pt0 = b_truth_pt1[istream];
			const float truth_pt1 = b_truth_pt2[istream];
			const float truth_lead_pt_val = b_truth_lead_pt[istream];

			const bool has_reco_pair = ( b_reco_pair[istream] != 0 );
			const float reco_pt0 = has_reco_pair ? b_reco_pt1[istream] : 0.0F;
			const float reco_pt1 = has_reco_pair ? b_reco_pt2[istream] : 0.0F;
			
			const int cent_val = b_cent[istream];
			
			const float zvrtx_val = b_zvrtx[istream];
			const float sumeT_val = b_sumeT[istream];

			const float psi2_val = b_psi2[istream];
			const int truth_match_idx1_val = b_truth_match_idx1[istream];
			const int truth_match_idx2_val = b_truth_match_idx2[istream];
			const float truth_phi1_val = b_truth_phi1[istream];
			const float truth_phi2_val = b_truth_phi2[istream];

			if ( i % print_every == 0 )
			{
				std::cout << "Sample " << isample << " (stream " << istream << ") : " << i << " / " << entries2 << "\r" << std::flush;
			}

			// Jet-v2 cross-check weight, w = 1 + amp*(c1 + c2), with amp solved in
			// the normalization pre-pass above so each leg realizes exactly
			// JETV2_SCALE. Both legs enter, so the subleading jet gets the same v2
			// as the leading one. Leg 1's angle to the event plane is the tree's
			// own dpsi2 branch; leg 2's is rebuilt from psi2 and truth_phi2, which
			// the tree stores no dpsi2 for. Both use the TRUTH phi, so the weight
			// is a property of the generated event and does not depend on whether
			// the pair reconstructed -- Fill, Miss and Fake pairs are all weighted
			// the same way, which is what keeps the miss/fake fractions (and hence
			// the efficiency folded into the response) internally consistent.
			double psi2_weight = 1.0;
			if ( jetv2_sys )
			{
				const double c1 = std::cos( 2.0 * b_dpsi2[istream] );
				const double c2 = std::cos( 2.0 * jetv2_dpsi2( psi2_val, truth_phi2_val ) );
				psi2_weight = 1.0 + jetv2_amp * ( c1 + c2 );
				jetv2_weight_sum += psi2_weight;
				jetv2_c1_sum += psi2_weight * c1;
				jetv2_c2_sum += psi2_weight * c2;
				jetv2_weight_n++;
				h_jetv2_weight->Fill(psi2_weight);
			}

			double event_scale 			= scale_factor[isample] * extra_scale[istream] * psi2_weight;
			
			double mbd_vertex_scale 	= 1.0;
			double sumeT_scale 			= 1.0;
			double centrality_scale 	= 1.0;

			// The event's highest-pT truth jet, no acceptance cut applied --
			// the same quantity the old ntuple called maxpttruth. Read
			// straight off dijet_pair_matching_v3.C's truth_lead_pt branch
			// rather than std::max(truth_pt0, truth_pt1): those two are only
			// the truth dijet CANDIDATE's legs (truth_idx1=0, truth_idx2=1),
			// which sit at -999 whenever the event has fewer than two truth
			// jets or the leading leg was dropped by the slimming -- both
			// cases where the real leading truth jet still exists and this
			// veto should still see it.
			const float maxpttruth_val = truth_lead_pt_val;

			const bool in_maxpttruth_range = ( maxpttruth_val >= sample_boundary[isample] && maxpttruth_val < sample_boundary[isample+1] );
			if ( !in_maxpttruth_range )
			{
				continue;
			}

			const bool in_centrality_bin = ( cent_val >= icentrality_bins[centrality_bin] && cent_val < icentrality_bins[centrality_bin+1] );
			if ( !in_centrality_bin && DO_CENT_CUT )
			{
				continue;
			}

			// No hard centrality cut on the MC: the response is built from the full
			// 0-90% embedded sample and the final (non-PRIMER) pass weights each
			// event by the data/MC z-vertex and SumET ratios for this centrality
			// bin (AA note Fig. 14/15 procedure). Restricting the MC to the target
			// window starves the response (~8.5x fewer pairs in 0-10%), which the
			// minentries trimming then turns into a different unfolding space.

			// Closure semantics: FULL closure trains on ALL events and unfolds the
			// same sample (response = test = everything); HALF closure trains on a
			// random half and unfolds the other half. The previous logic split
			// 50/50 for both modes, making FULL bit-identical to HALF.
			const bool accept_prob			= ( rng->Uniform() >= 0.5 );
			const bool fill_response 		= full_closure || ( half_closure && accept_prob );
			const bool use_for_response 	= not_a_closure_test || fill_response;
			const bool use_for_test 		= not_a_closure_test || full_closure || !fill_response;

			if ( !PRIMER1 || OVERRIDE_EVENT_WEIGHT )
			{
				const int imbd_bin = h_mbd_reweight -> FindBin( zvrtx_val );
				const int isumeT_bin = h_sumeT_reweight -> FindBin( sumeT_val );
				const int icent_bin = h_centrality_reweight -> FindBin( cent_val );

				if ( imbd_bin > 0 && imbd_bin <= h_mbd_reweight->GetNbinsX() )
				{
					mbd_vertex_scale = h_mbd_reweight -> GetBinContent(imbd_bin);
				}
				else
				{
					continue;
				}

				if ( isumeT_bin > 0 && isumeT_bin <= h_sumeT_reweight->GetNbinsX() )
				{
					sumeT_scale = h_sumeT_reweight -> GetBinContent(isumeT_bin);
				}
				else
				{
					if ( !DO_CENT_EVENT_WEIGHT ) continue;
				}

				if ( icent_bin > 0 && icent_bin <= h_centrality_reweight->GetNbinsX() )
				{
					centrality_scale = h_centrality_reweight -> GetBinContent(icent_bin);
				}
				else
				{
					if ( DO_CENT_EVENT_WEIGHT ) continue;
				}

				if ( mbd_vertex_scale < 0 ) { mbd_vertex_scale = 0; }
				if ( sumeT_scale < 0 ) { sumeT_scale = 0; }
				if ( centrality_scale < 0 ) { centrality_scale = 0; }

				double this_w = mbd_vertex_scale * sumeT_scale;
				if ( DO_CENT_EVENT_WEIGHT )
				{
					this_w *= centrality_scale;
					// this_w = centrality_scale * mbd_vertex_scale;
				}

				if ( not_a_closure_test )
				{
					event_scale *= this_w;
				}
				
			}

			if ( event_scale <= 0 && VETO_NULL_WEIGHTS ) { continue; }
			
			// leg 1 (index 0 in dijet_matching.C) is always the
			// truth-leading jet and leg 2 the truth-subleading one --
			// idxT = {0, 1} indexes directly into the pT-sorted
			// truth_jet_pT vector, so truth_pt0 >= truth_pt1 is guaranteed
			// by construction, not just the common case.
			
			const float e1 = truth_pt0;
			const float e2 = truth_pt1;
			float es1 = reco_pt0;
			float es2 = reco_pt1;

			float pt1_truth_bin = nbins;
			float pt2_truth_bin = nbins;
			float pt1_reco_bin = nbins;
			float pt2_reco_bin = nbins;

			h_truth_lead_sample[isample]->Fill(e1, event_scale);

			double smear1 = 0;
			double smear2 = 0;
			if ( not_a_closure_test )
			{
				smear1 = smear_random(e1);
				smear2 = smear_random(e2);
			}
			if (JES_sys != 0)
			{
				es1 = es1 + (JES_sys + smear1)*e1;
				es2 = es2 + (JES_sys + smear2)*e2;
			}
			else
			{
				es1 = es1 + smear1*e1;
				es2 = es2 + smear2*e2;
			}

			float maxi = std::max(es1, es2);
			float mini = std::min(es1, es2);
			float maxit = std::max(e1, e2);
			float minit = std::min(e1, e2);

			if (maxi >= ipt_bins[max_reco_bin]){ continue; }

			if (maxit > truth_leading_cut){ h_truth_lead->Fill(e1, event_scale);}
			if (minit >  truth_subleading_cut){ h_truth_sublead->Fill(e2, event_scale);}

			if (maxi >  reco_leading_cut){ h_reco_lead->Fill(maxi, event_scale);}
			if (mini >  reco_subleading_cut){ h_reco_sublead->Fill(mini, event_scale);}


			for (int ib = 0; ib < nbins; ib++)
			{
				if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
				{
					pt1_truth_bin = ib;
				}
				if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
				{
					pt2_truth_bin = ib;
				}
				if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
				{
					pt1_reco_bin = ib;
				}
				if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
				{
					pt2_reco_bin = ib;
				}
			}

			// Pair-level status read straight off dijet_pair_matching.C's
			// `category` -- no cuts are re-derived here, and no per-leg
			// combination is needed since the classification is already
			// pair-level (see the file header comment). kUESub pairs are
			// dropped alongside kSkip: dijet_pair_matching.C's header
			// explains they are neither a clean Fill, Miss, nor Fake, and
			// are already covered by the subleading-jet efficiency and the
			// inclusive cross-check.
			const bool skip_pair = (cat0 == kSkip) || (cat0 == kUESub);
			if ( skip_pair )
			{
				if (cat0 == kUESub) { ++n_uesub[isample]; }
				else { ++n_skip[isample]; }
				continue;
			}

			const bool miss_pair  = (cat0 == kMiss) || ( cat0 == kFakeMiss );
			const bool fake_pair  = (cat0 == kFake) || ( cat0 == kFakeMiss );
			const bool real_pair  = (cat0 == kFill);

			if (miss_pair) { ++n_miss[isample]; }
			if (fake_pair) { ++n_fake[isample]; }
			if (real_pair) { ++n_real[isample]; }
			
			// prior_qa::weight_bin() owns the flat-index-to-ROOT-bin convention and is
			// the same call the QA above uses, so the plots cannot silently disagree
			// with what is applied here. It returns k+1 for flat index
			// k = pt1_bin*nbins + pt2_bin: reading GetBinContent(k) — as version0 also
			// did — returns the weight of the NEIGHBORING truth bin (pt2-1), an
			// off-by-one that leaves the fraction-0 nominal untouched (all weights are
			// 1) but distorts every reweighted prior it is used for.

			if ( not_a_closure_test && NO_PRIMER )
			{
				const int prior_weight_bin = prior_qa::weight_bin(pt1_truth_bin, pt2_truth_bin, nbins);
				double flat_scale = h_flatreweight_pt1pt2->GetBinContent( prior_weight_bin );
				if ( !(flat_scale >= 0) )
				{
					std::cerr << "Warning: prior reweighting factor is negative or NaN for truth bins ("
					          << pt1_truth_bin << ", " << pt2_truth_bin << ")" << std::endl;
					flat_scale = 1.0;
				}
				event_scale *= flat_scale;

				// rewight jet
			}

			if ( miss_pair )
			{
				if (use_for_response)
				{
					h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					
					h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
					h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);

					rooResponse.Miss(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					rooResponse.Miss(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				}

				if (use_for_test)
				{
					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
					{
						h_truth_xj->Fill(minit/maxit, event_scale);
						h_linear_truth_xj->Fill(minit/maxit, event_scale);
					}
					
					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					
					h_e1e2->Fill(e1, e2, event_scale);
					h_e1e2->Fill(e2, e1, event_scale);
				}
			}
			
			if ( fake_pair  && DO_FAKES )
			{
			
				if (use_for_response)
				{
					h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
					
					h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
					h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
					
					rooResponse.Fake(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					rooResponse.Fake(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
				}

				if (use_for_test)
				{

					h_pt1pt2->Fill(es1, es2, event_scale);
					h_pt1pt2->Fill(es2, es1, event_scale);

					if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
					{
						h_reco_xj->Fill(mini/maxi, event_scale);
						h_linear_reco_xj->Fill(mini/maxi, event_scale);
					}
					h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
					
					// h_e1e2->Fill(es1, es2, event_scale);
					// h_e1e2->Fill(es2, es1, event_scale);
				}
			}
			
			if ( real_pair )
			{
			
				if (maxit > truth_leading_cut) { h_match_truth_lead->Fill(maxit, event_scale); }
				if (minit > truth_subleading_cut) { h_match_truth_sublead->Fill(minit, event_scale);}
				if (maxi >  reco_leading_cut){ h_match_reco_lead->Fill(maxi, event_scale);}
				if (mini >  reco_subleading_cut) { h_match_reco_sublead->Fill(mini, event_scale);}
				
			
				if (use_for_response)
				{
					
					h_flat_reco_to_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_to_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
					h_flat_truth_to_response_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_to_response_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					rooResponse.Fill(pt1_reco_bin + nbins*pt2_reco_bin,pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					rooResponse.Fill(pt2_reco_bin + nbins*pt1_reco_bin,pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
					h_count_flat_response_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, pt1_truth_bin + nbins*pt2_truth_bin);
					h_count_flat_response_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, pt2_truth_bin + nbins*pt1_truth_bin);
					h_count_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin);
					h_count_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin);
					h_count_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
					h_count_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);

				}
				// Independent of the response fill (version0 convention): the test/prior
				// truth distribution must contain ALL truth-good pairs, not only misses.
				// With an else-if here the nominal (both flags true) never adds matched
				// pairs to h_flat_truth_pt1pt2, so the prior reweighting divides the
				// unfolded data by a misses-only truth and the FULL closure compares
				// against an incomplete truth reference.
				if (use_for_test)
				{
					if (minit >= measure_subleading_cut && maxit >= measure_leading_cut && maxit < ipt_bins[measure_leading_bin + 2])
					{
						h_truth_xj->Fill(minit/maxit, event_scale);
						h_linear_truth_xj->Fill(minit/maxit, event_scale);
					}
					h_pt1pt2->Fill(es1, es2, event_scale);
					h_pt1pt2->Fill(es2, es1, event_scale);
					h_e1e2->Fill(e1, e2, event_scale);
					h_e1e2->Fill(e2, e1, event_scale);

					h_flat_truth_pt1pt2->Fill(pt1_truth_bin + nbins*pt2_truth_bin, event_scale);
					h_flat_truth_pt1pt2->Fill(pt2_truth_bin + nbins*pt1_truth_bin, event_scale);
				}

				// Reco spectrum of matched pairs: filled whenever the pair is usable as a
				// "measured" distribution (nominal: every real pair; closure: test half only).
				// Keeping this outside the use_for_response/use_for_test else-if matches the
				// version0 reference (h_reco_flat entries = 2 per real pair).
				if (use_for_test)
				{
					h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, event_scale);
					h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, event_scale);
				}

				if (maxi >= measure_leading_cut && mini>= measure_subleading_cut)
				{
					h_reco_xj->Fill(mini/maxi, event_scale);
					h_linear_reco_xj->Fill(mini/maxi, event_scale);
				}

				h_mbd_vertex->Fill(zvrtx_val, event_scale);
				h_sumeT->Fill(sumeT_val, event_scale);
				h_centrality->Fill(static_cast<float>(cent_val), event_scale);

			}
		}
		std::cout << std::endl;
    }

	for (int i = 0; i < 3; i++)
	{
		std::cout << "Sample " << i << " pairs -- Fill: " << n_real[i] << ", Miss: " << n_miss[i]
		          << ", Fake: " << n_fake[i] << ", Skip: " << n_skip[i] << ", UESub: " << n_uesub[i] << std::endl;
	}

	if ( jetv2_sys && jetv2_weight_n > 0 )
	{
		// Closure on the pre-pass: this is the unweighted mean of w over the
		// pairs the main loop actually kept (a slightly tighter set than the
		// pre-pass, which does not apply the centrality/reweight cuts), and the
		// realized per-leg v2 that the response was built with.
		const double mean_w = jetv2_weight_sum / (double) jetv2_weight_n;
		const double v2_leg1 = jetv2_c1_sum / jetv2_weight_sum;
		const double v2_leg2 = jetv2_c2_sum / jetv2_weight_sum;
		std::cout << "JetV2 closure -- requested v2 = " << jetv2_scale
		          << ", realized leg1 = " << v2_leg1 << ", leg2 = " << v2_leg2
		          << ", <weight> = " << mean_w
		          << " over " << jetv2_weight_n << " pairs" << std::endl;
		if ( std::fabs(v2_leg1 - jetv2_scale) > 0.1 * jetv2_scale
		  || std::fabs(v2_leg2 - jetv2_scale) > 0.1 * jetv2_scale )
		{
			std::cout << "  WARNING: realized jet v2 is more than 10% off the requested value."
			          << std::endl;
		}
	}

	h_flat_reco_pt1pt2->Scale(.5);
	h_flat_truth_pt1pt2->Scale(.5);
	h_flat_reco_to_response_pt1pt2->Scale(.5);
	h_flat_truth_to_response_pt1pt2->Scale(.5);
	h_flat_response_pt1pt2->Scale(.5);

	int number_of_mins = 1;
	if (minentries == 0){ number_of_mins = 20;}
  	for (int imin = minentries; imin < minentries + number_of_mins; imin++)       
    {
		std::cout << "Minimum Entries == " << imin << std::endl;
		int nbinsx = h_flat_response_pt1pt2->GetXaxis()->GetNbins();
		int nbinsy = h_flat_response_pt1pt2->GetYaxis()->GetNbins();

		TH1D *h_flat_truth_mapping = (TH1D*) h_flat_truth_to_response_pt1pt2->Clone();
		TH1D *h_flat_reco_mapping = (TH1D*) h_flat_reco_to_response_pt1pt2->Clone();
		if (minentries)
		{
			h_flat_truth_mapping->SetName("h_flat_truth_mapping");
			h_flat_reco_mapping->SetName("h_flat_reco_mapping");
		}
		else
		{
			h_flat_truth_mapping->SetName(Form("h_flat_truth_mapping_min%d", imin));
			h_flat_reco_mapping->SetName(Form("h_flat_reco_mapping_min%d", imin));
		}
      
		h_flat_reco_mapping->Reset();
		h_flat_truth_mapping->Reset();
	
		int nempty_reco = 0;
		std::vector<int> binnumbers_reco{};
		int nempty_truth = 0;
		std::vector<int> binnumbers_truth{};
		int nrecobins = 0;
		int ntruthbins = 0;
		// Pass 1: reco bins, unweighted-count threshold (unchanged behavior).
		for (int ib = 0; ib < nbins*nbins; ib++)
		{
			if (h_count_flat_reco_pt1pt2->GetBinContent(ib+1) < imin)
			{
				nempty_reco++;
				binnumbers_reco.push_back(0);
				h_flat_reco_mapping->SetBinContent(ib+1, 0);
			}
	  		else
			{
				nrecobins++;
				binnumbers_reco.push_back(nrecobins);
				h_flat_reco_mapping->SetBinContent(ib+1, nrecobins);
			}
		}
		// Pass 2: truth bins. Two criteria on top of the historical count cut:
		//   (a) minentries_truth: raise the truth-side count threshold until the
		//       number of kept truth bins is comparable to the kept reco bins
		//       (117 truth vs 95 reco leaves ~22 null-space directions in which
		//       Bayes iterations wander, i.e. non-converging prior dependence);
		//   (b) minentries_link: require the truth bin to have at least this many
		//       unweighted response fills into KEPT reco bins, so no retained
		//       truth bin is (nearly) unconstrained by the measurement.
		// With minentries_truth == minentries and minentries_link == 0 this is
		// identical to the historical trimming.
		const int imin_truth = (minentries == 0) ? imin : std::max(imin, minentries_truth);
		int ndropped_by_link = 0;
		for (int ib = 0; ib < nbins*nbins; ib++)
		{
			bool keep = h_count_flat_truth_pt1pt2->GetBinContent(ib+1) >= imin_truth;
			if (keep && minentries_link > 0)
			{
				double nlink = 0;
				for (int ibr = 0; ibr < nbins*nbins; ibr++)
				{
					if (binnumbers_reco.at(ibr) == 0) continue;
					nlink += h_count_flat_response_pt1pt2->GetBinContent(h_count_flat_response_pt1pt2->GetBin(ibr+1, ib+1));
				}
				if (nlink < minentries_link)
				{
					keep = false;
					ndropped_by_link++;
				}
			}
	  		if (!keep)
	    	{
				nempty_truth++;
				binnumbers_truth.push_back(0);
				h_flat_truth_mapping->SetBinContent(ib+1, 0);
			}
			else
			{
				ntruthbins++;
				binnumbers_truth.push_back(ntruthbins);
				h_flat_truth_mapping->SetBinContent(ib+1, ntruthbins);
			}
		}
		std::cout << "Trimming: kept " << ntruthbins << " truth / " << nrecobins
				  << " reco bins (truth threshold " << imin_truth
				  << ", dropped by link criterion: " << ndropped_by_link << ")" << std::endl;
		TH1D *h_flat_truth_skim = new TH1D(Form("h_flat_truth_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
		TH1D *h_flat_reco_skim = new TH1D(Form("h_flat_reco_skim_min%d", imin),"", nrecobins, 0, nrecobins);
		TH1D *h_flat_truth_to_unfold_skim = new TH1D(Form("h_flat_truth_to_unfold_skim_min%d", imin),"", ntruthbins, 0, ntruthbins);
		TH1D *h_flat_reco_to_unfold_skim = new TH1D(Form("h_flat_reco_to_unfold_skim_min%d", imin),"", nrecobins, 0, nrecobins);

		TH2D *h_flat_response_skim = new TH2D(Form("h_flat_response_skim_min%d", imin),"", nrecobins, 0, nrecobins, ntruthbins, 0, ntruthbins);

		if (minentries)
		{
			h_flat_response_skim->SetName("h_flat_response_skim");
			h_flat_truth_skim->SetName("h_flat_truth_skim");
			h_flat_reco_skim->SetName("h_flat_reco_skim");
			h_flat_truth_to_unfold_skim->SetName("h_flat_truth_to_unfold_skim");
			h_flat_reco_to_unfold_skim->SetName("h_flat_reco_to_unfold_skim");
		}

      	for (int ib = 0; ib < nbins*nbins; ib++)
		{
			int bintruth = binnumbers_truth.at(ib);      
			int binreco = binnumbers_reco.at(ib);
			if (binreco)
			{
				h_flat_reco_skim->SetBinContent(binreco, h_flat_reco_to_response_pt1pt2->GetBinContent(ib+1));
				h_flat_reco_skim->SetBinError(binreco, h_flat_reco_to_response_pt1pt2->GetBinError(ib+1));
				h_flat_reco_to_unfold_skim->SetBinContent(binreco, h_flat_reco_pt1pt2->GetBinContent(ib+1));
				h_flat_reco_to_unfold_skim->SetBinError(binreco, h_flat_reco_pt1pt2->GetBinError(ib+1));

			}

			if (bintruth)
			{
			
				h_flat_truth_skim->SetBinContent(bintruth, h_flat_truth_to_response_pt1pt2->GetBinContent(ib+1));
				h_flat_truth_skim->SetBinError(bintruth, h_flat_truth_to_response_pt1pt2->GetBinError(ib+1));
				h_flat_truth_to_unfold_skim->SetBinContent(bintruth, h_flat_truth_pt1pt2->GetBinContent(ib+1));
				h_flat_truth_to_unfold_skim->SetBinError(bintruth, h_flat_truth_pt1pt2->GetBinError(ib+1));

				for (int ibr = 0; ibr < nbins*nbins; ibr++)
				{
					int binreco2 = binnumbers_reco.at(ibr);
					if (binreco2)
					{
						int rbin = h_flat_response_skim->GetBin(binreco2, bintruth);
						int tbin = h_flat_response_pt1pt2->GetBin(ibr+1, ib+1);
						h_flat_response_skim->SetBinContent(rbin, h_flat_response_pt1pt2->GetBinContent(tbin));
						h_flat_response_skim->SetBinError(rbin, h_flat_response_pt1pt2->GetBinError(tbin));
					}
				}
	    	}
		}

		std::cout <<" Nbins skim reco = "<<h_flat_reco_skim->GetNbinsX()<<std::endl;
		std::cout <<" Nbins skim reco = "<<h_flat_truth_skim->GetNbinsX()<<std::endl;
		std::cout <<" Nbins skim reco = "<<h_flat_response_skim->GetXaxis()->GetNbins()<<"  --  " << h_flat_response_skim->GetYaxis()->GetNbins()<<std::endl; 
		RooUnfoldResponse rooResponsehist(h_flat_reco_skim, h_flat_truth_skim, h_flat_response_skim);
		if (minentries)
		{
			rooResponsehist.SetName("response_noempty");
		}
		else
		{
			rooResponsehist.SetName(Form("response_noempty_min%d", imin));
		}

		if (minentries)
		{
			TH1D* h_flat_unfold_skim[niterations];
			TH1D* h_flat_unfold_pt1pt2[niterations];
			TH1D* h_flat_refold_skim[niterations];
			TH1D* h_flat_refold_pt1pt2[niterations];
			TH1D* h_flat_truth_fold_skim;
			TH1D* h_flat_truth_fold_pt1pt2;
			int niter = 3;
			
			h_flat_truth_fold_skim = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_truth_skim, "hTruth_Folded");

			h_flat_truth_fold_pt1pt2 = (TH1D*) h_flat_truth_pt1pt2->Clone();
			h_flat_truth_fold_pt1pt2->Reset();
			h_flat_truth_fold_pt1pt2->SetName("h_flat_truth_fold_pt1pt2");
			for (int ib = 0; ib < nbins*nbins; ib++)
			{
				int bin = binnumbers_truth.at(ib);
				if (bin)
				{
					h_flat_truth_fold_pt1pt2->SetBinContent(ib+1, h_flat_truth_fold_skim->GetBinContent(bin));
					h_flat_truth_fold_pt1pt2->SetBinError(ib+1, h_flat_truth_fold_skim->GetBinError(bin));
				}
			}

			for (int iter = 0; iter < niterations; iter++ )
			{
		
				// handleFakes=true: h_flat_reco_skim (the "measured" side of
				// rooResponsehist) includes fake_pair content that
				// h_flat_response_skim (built from real_pair only) doesn't,
				// so RooUnfoldResponse::Setup already computes a non-empty
				// Fakes histogram (measured - response.ProjectionX()) -- but
				// RooUnfoldBayes defaults handleFakes to false and silently
				// unfolds that fake content as if it were signal instead of
				// excluding it, which grows with iteration count instead of
				// converging. This was invisible while fake_pair was
				// permanently disabled (no fakes to mishandle); it isn't now.
				RooUnfoldBayes   unfold (&rooResponsehist, h_flat_reco_to_unfold_skim, iter + 1, false, true);    // OR

				h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
				std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
				h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
				h_flat_unfold_pt1pt2[iter]->Reset();
				h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
				for (int ib = 0; ib < nbins*nbins; ib++)
				{
					int bin = binnumbers_truth.at(ib);
					if (bin)
					{
						h_flat_unfold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_unfold_skim[iter]->GetBinContent(bin));
						h_flat_unfold_pt1pt2[iter]->SetBinError(ib+1, h_flat_unfold_skim[iter]->GetBinError(bin));
					}
				}

				h_flat_refold_skim[iter] = (TH1D*) rooResponsehist.ApplyToTruth(h_flat_unfold_skim[iter], "hRefolded");

				std::cout <<" Nbins skim reco = "<<h_flat_refold_skim[iter]->GetNbinsX()<<std::endl;
				h_flat_refold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
				h_flat_refold_pt1pt2[iter]->Reset();
				h_flat_refold_pt1pt2[iter]->SetName(Form("h_flat_refold_pt1pt2_%d",iter));
				for (int ib = 0; ib < nbins*nbins; ib++)
				{
					int bin = binnumbers_truth.at(ib);
					if (bin)
					{
						h_flat_refold_pt1pt2[iter]->SetBinContent(ib+1, h_flat_refold_skim[iter]->GetBinContent(bin));
						h_flat_refold_pt1pt2[iter]->SetBinError(ib+1, h_flat_refold_skim[iter]->GetBinError(bin));
					}
				}

	    	}


			TCanvas *c = new TCanvas("c","c", 500, 500);

			for (int iter = 0; iter < niterations; iter++)
			{
				dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
				dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
			}
		
			dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
			dlutility::SetLineAtt(h_flat_reco_pt1pt2, kBlue, 2, 1);

			h_flat_truth_pt1pt2->Draw("hist");
			h_flat_reco_pt1pt2->Draw("hist same");
			h_flat_unfold_pt1pt2[niter]->Draw("same p");

			TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_truth_fold = new TH2D("h_pt1pt2_truth_fold", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
			TH2D *h_pt1pt2_unfold[niterations];
			TH2D *h_pt1pt2_refold[niterations];
			for (int iter = 0; iter < niterations; iter++)
			{
				h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
				h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
				h_pt1pt2_refold[iter] = new TH2D("h_pt1pt2_refold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
				h_pt1pt2_refold[iter]->SetName(Form("h_pt1pt2_refold_iter%d", iter));
			}
			TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
			TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_truth_fold = new TH1D("h_xj_truth_fold", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_truth_direct = new TH1D("h_xj_truth_direct", ";x_{J};",nbins, ixj_bins);
			TH1D *h_xj_unfold[niterations];
			TH1D *h_xj_refold[niterations];
			for (int iter = 0; iter < niterations; iter++)
			{
				h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
				h_xj_refold[iter] = new TH1D(Form("h_xj_refold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
			}

			h_pt1pt2_reco->SetTitle(";Reco p_{T, 1} [GeV]; Reco p_{T, 2} [GeV]; Counts * lumi scale ");
			h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
			h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");

			histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
			histo_opps::make_sym_pt1pt2(h_flat_truth_fold_pt1pt2, h_pt1pt2_truth_fold, nbins);
			histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins);
			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
				histo_opps::make_sym_pt1pt2(h_flat_refold_pt1pt2[iter], h_pt1pt2_refold[iter], nbins);
			}

			TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
			cpt1pt2->Divide(3, 1);
			cpt1pt2->cd(1);
			gPad->SetLogz();
			gPad->SetRightMargin(0.2);
			h_pt1pt2_reco->Draw("colz");
			dlutility::DrawSPHENIXpp(0.22, 0.85);
		
			cpt1pt2->cd(2);
			gPad->SetRightMargin(0.2);
			gPad->SetLogz();
			h_pt1pt2_truth->Draw("colz");
			cpt1pt2->cd(3);
			gPad->SetRightMargin(0.2);
			gPad->SetLogz();
			h_pt1pt2_unfold[niter]->Draw("colz");


			cpt1pt2->Print(Form("%s/unfolding_plots/pt1pt2_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(), cone_size, diagnostic_name.c_str()));

			histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_pt1pt2_truth_fold, h_xj_truth_fold, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			histo_opps::project_xj(h_e1e2, h_xj_truth_direct, nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
				histo_opps::project_xj(h_pt1pt2_refold[iter], h_xj_refold[iter], nbins, measure_leading_bin, measure_leading_bin + 2, measure_subleading_bin, nbins - 2);
			}
		
			TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
			dlutility::ratioPanelCanvas(cxj);
			cxj->cd(1);
			dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
			dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

			dlutility::SetLineAtt(h_truth_xj, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_truth_xj, kRed, 2, 24);

			dlutility::SetLineAtt(h_linear_truth_xj, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_linear_truth_xj, kRed, 2, 25);

			dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_reco_xj, kBlue, 2, 24);

			dlutility::SetLineAtt(h_linear_reco_xj, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_linear_reco_xj, kBlue, 2, 25);

		
			dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
			dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

			dlutility::SetLineAtt(h_xj_reco, kBlue, 2, 1);
			dlutility::SetMarkerAtt(h_xj_reco, kBlue, 1, 8);
		
			histo_opps::normalize_histo(h_xj_truth, nbins);
			histo_opps::normalize_histo(h_xj_reco, nbins);
			histo_opps::normalize_histo(h_truth_xj, nbins);
			histo_opps::normalize_histo(h_reco_xj, nbins);

			//histo_opps::normalize_histo(h_linear_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
			h_truth_xj->Scale(1./h_linear_truth_xj->Integral(), "width");
			//h_linear_reco_xj->Scale(1./h_linear_reco_xj->Integral(), "width");

			for (int iter = 0; iter < niterations; iter++)
			{
				histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
			}

			dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
			h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
			//h_xj_truth->SetMaximum(5);
			h_xj_truth->Draw("p");
			h_xj_unfold[niter]->Draw("same p");
			h_xj_reco->Draw("hist same");
			h_xj_truth->Draw("same hist");
			h_xj_reco->Draw("p same");
			h_xj_unfold[niter]->Draw("same hist");
			h_xj_unfold[niter]->Draw("same p");
			dlutility::DrawSPHENIXpp(0.22, 0.84);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.22, 0.74);
			dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
			dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
			TLegend *leg = new TLegend(0.22, 0.3, 0.4, 0.55);
			leg->SetLineWidth(0);
			leg->AddEntry(h_xj_reco, "Pythia Reco");
			leg->AddEntry(h_xj_truth, "Pythia8 Truth");
			leg->AddEntry(h_xj_unfold[niter], "Unfolded");
			leg->Draw("same");
			dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
			
			cxj->cd(2);

			TH1D *h_reco_compare = (TH1D*) h_xj_unfold[niter]->Clone();
			h_reco_compare->Divide(h_xj_truth);
			h_reco_compare->SetTitle(";x_{J}; Unfold / Truth");
			dlutility::SetFont(h_reco_compare, 42, 0.06);
			dlutility::SetLineAtt(h_reco_compare, kBlack, 1,1);
			dlutility::SetMarkerAtt(h_reco_compare, kBlack, 1,8);

		
			h_reco_compare->Draw("p");
			TLine *line = new TLine(0.1, 1, 1, 1);
			line->SetLineStyle(4);
			line->SetLineColor(kRed + 3);
			line->SetLineWidth(2);
			line->Draw("same");

			TCanvas *c_iter = new TCanvas("c_iter", "c_iter");
			TH1D *h_closure[niterations];

			for (int iter = 0; iter < niterations; iter++)
			{
				h_closure[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
				h_closure[iter]->SetName(Form("h_closure_%d", iter));
				
				h_closure[iter]->Divide(h_xj_truth);
				h_closure[iter]->SetTitle(";x_{J}; Unfold / Truth");
				dlutility::SetFont(h_closure[iter], 42, 0.06);
			}
			int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 

			for (int iter = 0; iter < 5; iter++)
			{
				dlutility::SetLineAtt(h_closure[1 + 2*iter], colors[iter], 1,1);
				dlutility::SetMarkerAtt(h_closure[1 + 2*iter], colors[iter], 1,8);
			}
			h_closure[1]->Draw("p");
			h_closure[3]->Draw("p same");
			h_closure[5]->Draw("p same");
			h_closure[7]->Draw("p same");
			h_closure[9]->Draw("p same");
			TLine *line2 = new TLine(0.1, 1, 1, 1);
			line2->SetLineStyle(4);
			line2->SetLineColor(kRed + 3);
			line2->SetLineWidth(2);
			line2->Draw("same");


			TCanvas *cjet = new TCanvas("cjet","cjet", 700, 500);
			dlutility::createCutCanvas(cjet);
			cjet->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_truth_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead, kRed, 1, 24);

			dlutility::SetLineAtt(h_truth_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_truth_sublead, kBlue, 1, 24);

			dlutility::SetLineAtt(h_reco_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_reco_lead, kRed, 1, 20);

			dlutility::SetLineAtt(h_reco_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_reco_sublead, kBlue, 1, 20);

			h_truth_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
			h_truth_lead->Draw();
			h_truth_sublead->Draw("same");
			h_reco_lead->Draw("same");
			h_reco_sublead->Draw("same");
			cjet->cd(2);
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", 0.1*cone_size), 0.1, 0.74, 0, kBlack, 0.1);
			dlutility::drawText("All Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);

			TLegend *leg33 = new TLegend(0.01, 0.13, 0.7, 0.62);
			leg33->SetLineWidth(0);
			leg33->AddEntry(h_reco_lead, "Leading Reco");
			leg33->AddEntry(h_reco_sublead, "Subleading Reco");
			leg33->AddEntry(h_truth_lead, "Leading Truth");
			leg33->AddEntry(h_truth_sublead, "Subleading Truth");
			leg33->SetTextSize(0.1);
			leg33->Draw("same");

			TCanvas *cjetdiv = new TCanvas("cjetdiv","cjetdiv", 700, 500);
			dlutility::createCutCanvas(cjetdiv);
			cjetdiv->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_truth_lead, kBlack, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead, kBlack, 0.8, 24);

			dlutility::SetLineAtt(h_truth_lead_sample[0], kRed, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[0], kRed, 1, 20);
			dlutility::SetLineAtt(h_truth_lead_sample[1], kGreen, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[1], kGreen, 1, 20);
			dlutility::SetLineAtt(h_truth_lead_sample[2], kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_truth_lead_sample[2], kBlue, 1, 20);

			h_truth_lead->SetTitle(";Leading Jet p_{T} [GeV];counts * lumiscale");
			h_truth_lead_sample[0]->SetMinimum(1);
			h_truth_lead_sample[0]->Draw("p");
			h_truth_lead->Draw("p same");
			h_truth_lead_sample[1]->Draw("same p");
			h_truth_lead_sample[2]->Draw("same p");
			h_truth_lead->Draw(" same p");


			TLine  *line1 = new TLine(14, 1, 14, h_truth_lead->GetBinContent(15));
			dlutility::SetLineAtt(line1, kRed, 1.5, 1);
			line1->Draw("same");
			TLine  *lin2 = new TLine(20, 1, 20, h_truth_lead->GetBinContent(21));
			dlutility::SetLineAtt(lin2, kGreen, 1.5, 1);
			lin2->Draw("same");
			TLine  *lin3 = new TLine(30, 1, 30, h_truth_lead->GetBinContent(31));
			dlutility::SetLineAtt(lin3, kBlue, 1.5, 1);
			lin3->Draw("same");
			cjetdiv->cd(2);

			dlutility::DrawSPHENIXppsize(0.05, 0.84, 0.08, 1, 0, 1, "Pythia8");
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.05, 0.74, 0, kBlack, 0.08);
			dlutility::drawText("All Dijet Pairs", 0.05, 0.69, 0, kBlack, 0.08);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.05, 0.64, 0, kBlack, 0.08);

			TLegend *leg3 = new TLegend(0.01, 0.13, 0.7, 0.62);
			leg3->SetLineWidth(0);
			leg3->AddEntry(h_truth_lead, "Combined Sample");
			leg3->AddEntry(h_truth_lead_sample[0], "10 GeV Sample");
			leg3->AddEntry(h_truth_lead_sample[1], "20 GeV Sample");
			leg3->AddEntry(h_truth_lead_sample[2], "30 GeV Sample");
			leg3->SetTextSize(0.08);
			leg3->Draw("same");

			cjetdiv->Print(Form("%s/unfolding_plots/combined_sample_%s_r%02d_%s.pdf", rb.get_code_location().c_str(), system_string.c_str(),  cone_size, diagnostic_name.c_str()));

			TCanvas *cmjet = new TCanvas("cmjet","cmjet", 700, 500);
			dlutility::createCutCanvas(cmjet);
			cmjet->cd(1);
			gPad->SetLogy();

			dlutility::SetLineAtt(h_match_truth_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_match_truth_lead, kRed, 1, 24);

			dlutility::SetLineAtt(h_match_truth_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_match_truth_sublead, kBlue, 1, 24);

			dlutility::SetLineAtt(h_match_reco_lead, kRed, 1, 1);
			dlutility::SetMarkerAtt(h_match_reco_lead, kRed, 1, 20);

			dlutility::SetLineAtt(h_match_reco_sublead, kBlue, 1, 1);
			dlutility::SetMarkerAtt(h_match_reco_sublead, kBlue, 1, 20);

			h_match_reco_lead->SetTitle(";Jet p_{T} [GeV];counts * lumiscale");
			
			h_match_reco_lead->Draw();
			h_match_truth_sublead->Draw("same");
			h_match_truth_lead->Draw("same");
			h_match_reco_sublead->Draw("same");
			cmjet->cd(2);
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);
			dlutility::drawText(Form("anti-k_{T} R = %0.1f", cone_size*0.1), 0.1, 0.74, 0, kBlack, 0.1);
			dlutility::drawText("Matched Dijet Pairs", 0.1, 0.69, 0, kBlack, 0.1);
			dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.1, 0.64, 0, kBlack, 0.1);
			TLegend *lef3 = new TLegend(0.01, 0.13, 0.7, 0.62);
			lef3->SetLineWidth(0);
			lef3->AddEntry(h_match_reco_lead, "Leading Reco");
			lef3->AddEntry(h_match_reco_sublead, "Subleading Reco");
			lef3->AddEntry(h_match_truth_lead, "Leading Truth");
			lef3->AddEntry(h_match_truth_sublead, "Subleading Truth");
			lef3->SetTextSize(0.1);
			lef3->Draw("same");


			TCanvas *cresponse = new TCanvas("fd","fd", 500, 500);
			h_flat_response_pt1pt2->Draw("colz");
			TCanvas *cresponseskim = new TCanvas("fds","fds", 500, 500);
			cresponseskim->SetLogz();
			h_flat_response_skim->Draw("colz");
			dlutility::DrawSPHENIXppsize(0.1, 0.84, 0.1);


			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size);
			
			if (PRIMER1)
			{
				responsepath += "_PRIMER1";
			}
			if (PRIMER2)
			{
				responsepath += "_PRIMER2";
			}
			

			if (half_closure)
			{
				responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_HALF";
			}
			else if (full_closure)
			{
				responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_FULL";
			}
			responsepath += "_" + sys_name;

			responsepath += ".root";
			TFile *fr = new TFile(responsepath.Data(),"recreate");
			rooResponsehist.Write();
			h_flat_reco_pt1pt2->Write();
			h_flat_truth_pt1pt2->Write();
			h_flat_truth_to_response_pt1pt2->Write();
			h_flat_truth_fold_pt1pt2->Write();
			h_flat_response_pt1pt2->Write();
			h_flat_truth_mapping->Write();
			h_flat_reco_mapping->Write();
			h_flat_truth_skim->Write();
			h_flat_reco_skim->Write();
			h_flat_response_skim->Write();
			h_mbd_vertex->Write();
			h_centrality->Write();
			h_sumeT->Write();
			h_pt1pt2->Write();
			h_e1e2->Write();
			h_jetv2_weight->Write();
			// Unweighted trimming inputs, saved so minentries_truth/minentries_link
			// can be tuned offline without re-running the event loop.
			h_count_flat_truth_pt1pt2->Write();
			h_count_flat_reco_pt1pt2->Write();
			h_count_flat_response_pt1pt2->Write();
			// Full closure drawing
			h_linear_truth_xj->Write();
			h_xj_truth->Write();
			h_xj_truth_fold->Write();
			h_xj_truth_direct->Write();
			h_truth_xj->Write();
			for (int iter = 0 ; iter < niterations; iter++)
			{
				h_flat_unfold_pt1pt2[iter]->Write();
				h_xj_unfold[iter]->Write();
				h_flat_refold_pt1pt2[iter]->Write();
				h_xj_refold[iter]->Write();
			}
			h_xj_reco->Write();
			
				
			fr->Close();

		}
		else
		{
			TString responsepath = "response_matrices/response_matrix_" + system_string + "_r0" + std::to_string(cone_size) + "_min_" + std::to_string(imin) + ".root";

			TFile *fr = new TFile(responsepath.Data(),"recreate");
			rooResponsehist.Write();
			
			TH1D *h1 = (TH1D*)h_flat_reco_pt1pt2->Clone();
			TH1D *h2 = (TH1D*)h_flat_truth_pt1pt2->Clone();
			TH2D *h3 = (TH2D*)h_flat_response_pt1pt2->Clone();
			TH1D *h4 = (TH1D*)h_flat_truth_mapping->Clone();
			TH1D *h5 = (TH1D*)h_flat_reco_mapping->Clone();
			TH1D *h6 = (TH1D*)h_flat_truth_skim->Clone();
			TH1D *h7 = (TH1D*)h_flat_reco_skim->Clone();
			TH2D *h8 = (TH2D*)h_flat_response_skim->Clone();
			h1->Write();
			h2->Write();
			h3->Write();
			h4->Write();
			h5->Write();
			h6->Write();
			h7->Write();
			h8->Write();
						
			fr->Close();
		}
		std::cout << "Reco empty: " << nempty_reco << std::endl;
		std::cout << "Truth empty: " << nempty_truth << std::endl;
    }
  	
	return 0;
}
