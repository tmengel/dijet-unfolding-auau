#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "TFile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "PlotUtils.h"

enum SubleadFlag {
  kSubleadMatched = -1,
  kNoRecoCandidate = 0,
  kRealElsewhere = 1,
  kTruthOutsideFiducial = 2
};

// QA for the pair categories used by createResponse_noempty_AA.cxx.
//
// Runs the identical selection/classification pipeline (same cuts, same
// smearing, same sublead_flag logic) but only tallies counts instead of
// filling a response matrix, so the category fractions reported here match
// exactly what the response builder sees. Does NOT apply the mbd/sumET/
// centrality event reweighting the production macro applies in its final
// pass -- that reweighting only rescales events to match data, it does not
// change which category a pair falls into, so skipping it keeps this macro
// self-contained (no reweighting histogram files required).
int checkSubleadPopulations_AA (
	const std::string configfile = "/home/tmengel/PPG14/version1/dijet-unfolding-final/configs/binning_AA.config",
	const int cone_size = 3,
	const int centrality_bin = 0
)
{
	gStyle->SetOptStat(0);
	PlotUtils::set_sphenix_style();

	read_binning rb(configfile.c_str());
	std::string system_string = rb.get_system_string(centrality_bin);
	std::cout << "System string: " << system_string << std::endl;

	std::string j10_file = std::getenv("TNUPLE_SIM_FILE_JET10");
	std::string j20_file = std::getenv("TNUPLE_SIM_FILE_JET20");
	std::string j30_file = std::getenv("TNUPLE_SIM_FILE_JET30");

	float maxpttruth[3];
	float pt1_truth[3];
	float pt2_truth[3];
	float dphi_truth[3];
	float pt1_reco[3];
	float pt2_reco[3];
	float dphi_reco[3];
	float match[3];
	float centrality[3];
	float sublead_flag[3] = { -1, -1, -1 };

	float n_events[3];
	float b_n_events = 0;

	float event_weight[3] = {1, 1, 1};
	bool has_event_weight[3] = {false, false, false};
	bool has_sublead_flag[3] = {false, false, false};

	TFile * fin[3];
	fin[0] = new TFile(j10_file.c_str(), "r");
	fin[1] = new TFile(j20_file.c_str(), "r");
	fin[2] = new TFile(j30_file.c_str(), "r");
	TNtuple *tn[3];
	for (int i = 0; i < 3; i++)
	{
		if (!fin[i] || fin[i]->IsZombie())
		{
			std::cerr << "Cannot open matched simulation file " << i << std::endl;
			return 1;
		}

		tn[i] = (TNtuple*) fin[i]->Get("tn_match");
		if (!tn[i])
		{
			std::cerr << "Missing tn_match in " << fin[i]->GetName() << std::endl;
			return 1;
		}

		tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
		tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
		tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
		tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
		tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
		tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
		tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
		tn[i]->SetBranchAddress("matched", &match[i]);
		tn[i]->SetBranchAddress("centrality", &centrality[i]);

		has_sublead_flag[i] = tn[i]->GetBranch("sublead_flag") != nullptr;
		if (has_sublead_flag[i]) tn[i]->SetBranchAddress("sublead_flag", &sublead_flag[i]);
		else std::cerr << "Warning: no sublead_flag branch in sample " << i
		                << "; those entries default to kSubleadMatched (no miss/fake reclassification, no UE-flux skip)." << std::endl;

		has_event_weight[i] = tn[i]->GetBranch("weight") != nullptr;
		if (has_event_weight[i]) tn[i]->SetBranchAddress("weight", &event_weight[i]);

		TNtuple *tn_stats = (TNtuple*) fin[i]->Get("tn_stats");
		if (!tn_stats)
		{
			std::cerr << "Missing tn_stats in " << fin[i]->GetName() << std::endl;
			return 1;
		}
		tn_stats->SetBranchAddress("nevents", &b_n_events);
		n_events[i] = 0;
		for (int j = 0; j < tn_stats->GetEntries(); j++)
		{
			tn_stats->GetEntry(j);
			n_events[i] += b_n_events;
		}
	}
	std::cout << "n_events: " << n_events[0] << " " << n_events[1] << " " << n_events[2] << std::endl;

	const float cs_10 = 0.000003997;
	const float cs_20 = 6.218e-8;
	const float cs_30 = 2.505e-9;

	float scale_factor[3];
	scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
	scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30;
	scale_factor[2] = 1.0;
	std::cout << "Scale factors: " << scale_factor[0] << " " << scale_factor[1] << " " << scale_factor[2] << std::endl;

	const int nbins = rb.get_nbins();
	float ipt_bins[nbins+1];
	rb.get_pt_bins(ipt_bins);
	const int max_reco_bin = rb.get_maximum_reco_bin();

	const double dphicut = rb.get_dphicut();

	const float truth_leading_cut = rb.get_truth_leading_cut();
	const float truth_subleading_cut = rb.get_truth_subleading_cut();
	const float reco_leading_cut = rb.get_reco_leading_cut();
	const float reco_subleading_cut = rb.get_reco_subleading_cut();

	float sample_boundary[4] = {0};
	for (int ib = 0; ib < 4; ib++) { sample_boundary[ib] = rb.get_sample_boundary(ib); }

	const double JES_sys = rb.get_jes_sys();

	auto * f_smear = (TF1*) rb.get_smear_function(centrality_bin);
	if (!f_smear)
	{
		std::cerr << "NO SMEAR FUNCTION" << std::endl;
		return 1;
	}

	// Same closed-form truncated-Gaussian smear as createResponse_noempty_AA.cxx.
	const double smear_trunc = 0.5;
	auto smear_random = [&](double pt) -> double
	{
		const double sigma = f_smear->Eval(pt);
		if ( !(sigma > 0) ) { return 0.0; }
		const double phi_hi = TMath::Freq(smear_trunc / sigma);
		double u = (1.0 - phi_hi) + gRandom->Rndm() * (2.0 * phi_hi - 1.0);
		if ( u <= 0.0 ) { u = std::numeric_limits<double>::epsilon(); }
		if ( u >= 1.0 ) { u = 1.0 - std::numeric_limits<double>::epsilon(); }
		return sigma * TMath::NormQuantile(u);
	};

	gRandom->SetSeed(4357);

	// Category counters. "considered" = passed the maxpttruth sample-boundary
	// cut and the reco-overflow cut, i.e. every pair that reaches the
	// classification block in createResponse_noempty_AA.cxx.
	double n_considered = 0, w_considered = 0;
	double n_skip_pair = 0, w_skip_pair = 0;
	double n_ue_flux = 0, w_ue_flux = 0;           // sublead_flag == kRealElsewhere, dropped regardless of truth/reco_good
	double n_miss_pair = 0, w_miss_pair = 0;
	double n_missed_from_fid = 0, w_missed_from_fid = 0;
	double n_fake_pair = 0, w_fake_pair = 0;
	double n_real_pair = 0, w_real_pair = 0;
	double n_other = 0, w_other = 0;               // truth_good && reco_good && !match -- not filled anywhere in production

	// Raw sublead_flag composition among all considered pairs, independent of
	// the truth_good/reco_good classification above.
	double n_flag_matched = 0, n_flag_no_reco = 0, n_flag_real_elsewhere = 0, n_flag_outside_fid = 0;

	for (int isample = 0; isample < 3; isample++)
	{
		const int entries = tn[isample]->GetEntries();
		for (int i = 0; i < entries; i++)
		{
			tn[isample]->GetEntry(i);

			const bool in_maxpttruth_range = ( maxpttruth[isample] >= sample_boundary[isample] && maxpttruth[isample] < sample_boundary[isample+1] );
			if ( !in_maxpttruth_range ) { continue; }

			const double event_scale = scale_factor[isample] * (has_event_weight[isample] ? event_weight[isample] : 1.0);

			float max_truth, min_truth, max_reco, min_reco;
			if (pt1_truth[isample] >= pt2_truth[isample])
			{
				max_truth = pt1_truth[isample]; max_reco = pt1_reco[isample];
				min_truth = pt2_truth[isample]; min_reco = pt2_reco[isample];
			}
			else
			{
				max_truth = pt2_truth[isample]; max_reco = pt2_reco[isample];
				min_truth = pt1_truth[isample]; min_reco = pt1_reco[isample];
			}

			const float e1 = max_truth, e2 = min_truth;
			float es1 = max_reco, es2 = min_reco;

			const double smear1 = smear_random(e1);
			const double smear2 = smear_random(e2);
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

			const float maxi = std::max(es1, es2);
			const float mini = std::min(es1, es2);
			const float maxit = std::max(e1, e2);
			const float minit = std::min(e1, e2);

			if (maxi >= ipt_bins[max_reco_bin]) { continue; }

			const bool truth_good = (maxit >= truth_leading_cut && minit >= truth_subleading_cut && dphi_truth[isample] >= dphicut);
			const bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] >= dphicut);

			const bool sublead_is_ue_flux = ( sublead_flag[isample] == 1 );
			const bool sublead_truth_outside_range = ( sublead_flag[isample] == 2 );
			const bool skip_pair = !truth_good && !reco_good;

			n_considered += 1; w_considered += event_scale;

			if      (sublead_flag[isample] == kSubleadMatched)      { n_flag_matched++; }
			else if (sublead_flag[isample] == kNoRecoCandidate)     { n_flag_no_reco++; }
			else if (sublead_flag[isample] == kRealElsewhere)       { n_flag_real_elsewhere++; }
			else if (sublead_flag[isample] == kTruthOutsideFiducial){ n_flag_outside_fid++; }

			// Mutually exclusive: an entry that is both skip_pair and UE-flux
			// (common -- a missing/failed reco pair often also carries the
			// UE-flux flag) is booked once, as skip_pair, so the categories
			// below partition n_considered exactly instead of double-counting.
			if (skip_pair)                            { n_skip_pair++; w_skip_pair += event_scale; }
			else if (sublead_is_ue_flux)               { n_ue_flux++;   w_ue_flux   += event_scale; }
			if (skip_pair || sublead_is_ue_flux) { continue; }

			const bool miss_pair       = truth_good && !reco_good;
			const bool missed_from_fid = !truth_good && reco_good && sublead_truth_outside_range;
			const bool fake_pair       = !truth_good && reco_good && !sublead_truth_outside_range;
			const bool real_pair       = truth_good && reco_good && match[isample];
			const bool other           = truth_good && reco_good && !match[isample];

			if (miss_pair)       { n_miss_pair++;       w_miss_pair       += event_scale; }
			if (missed_from_fid) { n_missed_from_fid++; w_missed_from_fid += event_scale; }
			if (fake_pair)       { n_fake_pair++;       w_fake_pair       += event_scale; }
			if (real_pair)       { n_real_pair++;       w_real_pair       += event_scale; }
			if (other)           { n_other++;           w_other           += event_scale; }
		}
	}

	const auto pct = [](double x, double total) { return total > 0 ? 100.0 * x / total : 0.0; };

	std::cout << "\n=== Pair populations (" << system_string << ", r0" << cone_size << ") ===\n";
	std::cout << "Considered pairs: " << n_considered << " raw, " << w_considered << " lumi-weighted\n\n";
	printf("%-18s %12s %8s   %14s %8s\n", "category", "raw", "raw %", "weighted", "wtd %");
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "skip_pair",       n_skip_pair,       pct(n_skip_pair, n_considered),       w_skip_pair,       pct(w_skip_pair, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "ue_flux (skip)",  n_ue_flux,         pct(n_ue_flux, n_considered),         w_ue_flux,         pct(w_ue_flux, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "miss_pair",       n_miss_pair,       pct(n_miss_pair, n_considered),       w_miss_pair,       pct(w_miss_pair, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "missed_from_fid", n_missed_from_fid, pct(n_missed_from_fid, n_considered), w_missed_from_fid, pct(w_missed_from_fid, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "fake_pair",       n_fake_pair,       pct(n_fake_pair, n_considered),       w_fake_pair,       pct(w_fake_pair, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "real_pair",       n_real_pair,       pct(n_real_pair, n_considered),       w_real_pair,       pct(w_real_pair, w_considered));
	printf("%-18s %12.0f %7.2f%%   %14.4g %7.2f%%\n", "other (dropped)", n_other,           pct(n_other, n_considered),           w_other,           pct(w_other, w_considered));

	std::cout << "\n=== sublead_flag composition (all considered pairs, unweighted) ===\n";
	printf("%-22s %12.0f %7.2f%%\n", "kSubleadMatched",       n_flag_matched,        pct(n_flag_matched, n_considered));
	printf("%-22s %12.0f %7.2f%%\n", "kNoRecoCandidate",      n_flag_no_reco,        pct(n_flag_no_reco, n_considered));
	printf("%-22s %12.0f %7.2f%%\n", "kRealElsewhere",        n_flag_real_elsewhere, pct(n_flag_real_elsewhere, n_considered));
	printf("%-22s %12.0f %7.2f%%\n", "kTruthOutsideFiducial", n_flag_outside_fid,    pct(n_flag_outside_fid, n_considered));

	// Bar-chart summary, lumi-weighted fractions.
	const char* labels[7] = {"skip_pair", "ue_flux", "miss", "missed_fid", "fake", "real", "other"};
	const double vals[7] = {
		pct(w_skip_pair, w_considered), pct(w_ue_flux, w_considered), pct(w_miss_pair, w_considered),
		pct(w_missed_from_fid, w_considered), pct(w_fake_pair, w_considered), pct(w_real_pair, w_considered),
		pct(w_other, w_considered)
	};
	TH1D *h_pop = new TH1D("h_pop", ";;Fraction of considered pairs [%]", 7, 0, 7);
	for (int ib = 0; ib < 7; ib++)
	{
		h_pop->SetBinContent(ib+1, vals[ib]);
		h_pop->GetXaxis()->SetBinLabel(ib+1, labels[ib]);
	}
	h_pop->SetFillColor(kAzure-4);
	h_pop->GetXaxis()->SetLabelSize(0.05);
	h_pop->SetMinimum(0);

	TCanvas *c = new TCanvas("c_pop", "c_pop", 700, 500);
	c->SetBottomMargin(0.12);
	h_pop->Draw("hist text0");
	dlutility::drawText(Form("%s, anti-k_{T} R = %0.1f", system_string.c_str(), 0.1*cone_size), 0.4, 0.85);
	c->Print(Form("%s/unfolding_plots/sublead_populations_%s_r%02d_cent%d.pdf",
	              rb.get_code_location().c_str(), system_string.c_str(), cone_size, centrality_bin));

	return 0;
}
