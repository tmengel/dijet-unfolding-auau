#ifndef _DIJET_PAIR_MATCHING_INCLUSIVE_V3_C_
#define _DIJET_PAIR_MATCHING_INCLUSIVE_V3_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Inclusive-pairing variant of dijet_pair_matching_v3.C, run on the
// SLIMMED tree written by slim_matched_jets.C.
//
// The output schema is branch-for-branch the one
// dijet_pair_matching_inclusive_v2.C writes (same names, same category
// numbering, same pair_partner_rank convention), so
// createResponse_exclusive_v2_AA_inclusive.cxx reads it unchanged. The
// v3-only branches are appended at the end.
//
// THE ALGORITHM -- RECO-DRIVEN PAIRING (one row per (event, partner rank))
//
// The leading reco jet is paired with every other reco jet:
// (pT1,pT2), (pT1,pT3), (pT1,pT4), ... and each pairing is categorized on
// its own. The pairing is built on the RECO side and the truth pair is
// then read off the two legs' truth partners:
//
//   reco pair  = ( accepted[0], accepted[k] ), k = 1, 2, 3, ...
//                where `accepted` is the event's matched reco jets that
//                are inside the calorimeter acceptance with E > 0, pT
//                descending. NOTE: no pT threshold enters the list, so a
//                pairing still gets a row (and can be a Miss) when its
//                legs are too soft to be a reco dijet candidate.
//                accepted[0] is the leading such jet.
//       accepted when pT_1(2) > 25(8) GeV on the two legs and
//       |dphi_12| >= 7pi/8.
//
//   truth pair = the two legs' OWN truth partners, ordered by truth pT --
//                truth leg 1 is the harder of the two. Every reco jet in
//                the slim tree is truth-matched by construction, so this
//                pair always exists.
//       accepted when pT_1(2) > 14(7) GeV, truth accept_eta on both legs,
//       and |dphi_12| > 7pi/8.
//
// THREE populations are written (kSkip is dropped, never filled):
//
//     truth cand | reco cand | category
//     -----------+-----------+-----------
//         yes    |    yes    | Fill      (0)
//         yes    |    no     | Miss      (1)
//         no     |    yes    | Fake      (2)
//         no     |    no     | (not written)
//
// There is no FakeMiss here, and it is not a category that was dropped:
// it cannot arise. FakeMiss is the exclusive statement "the reco dijet
// candidate is not the event's real leading truth pair", and inclusively
// there is no leading-pair requirement to violate -- the truth pair IS
// the reco pair's own two matches by construction. A resolution-driven
// pT-ordering flip between the two legs is still a Fill, flagged by
// legs_swapped; reco_prov1/reco_prov2 are correspondingly always
// { kProvLead, kProvSub } or { kProvSub, kProvLead }, and
// kProvOtherTruth / kProvUE / kUESub are never assigned.
//
// COVERAGE, stated plainly: the enumeration is reco-driven, so every
// accepted reco pair built on the leading reco jet gets a row, but a
// truth pair whose harder leg is not the leading reco jet's own truth
// partner never gets one. Truth pairs are therefore complete only for
// events where the leading accepted reco jet belongs to the truth jet
// heading that pair. dijet_pair_matching_v3.C (exclusive) is the macro
// that accounts for that population, as its FakeMiss category.
//
// The acc_truth_* collection of the slim tree (truth jets in acceptance
// with NO reco match) is deliberately not read here: an unmatched truth
// jet has no reco leg, so it cannot enter a reco-driven pairing. Its
// effect on the leading pair is what dijet_pair_matching_v3.C reports.
//
// Input:  the tree written by slim_matched_jets.C.
// Output: one row per (event, pair_partner_rank) that is a Fill, Miss or
//         Fake. pair_partner_rank is k+1, so rank 2 is the leading pair
//         (pT1,pT2), rank 3 is (pT1,pT3), and so on -- the same
//         convention as dijet_pair_matching_inclusive_v2.C.
namespace DijetPairInclusiveV3
{
    // Same numbering as DijetPairV2 / DijetPairV3 --
    // createResponse_exclusive_v2_AA_inclusive.cxx hard-codes it.
    enum Category
    {
        kFill     = 0,
        kMiss     = 1,
        kFake     = 2,
        kSkip     = 3, // never written out
        kUESub    = 4, // never assigned
        kFakeMiss = 5  // never assigned -- see the header
    };

    // What a reco leg is, relative to THIS pairing's truth pair.
    enum Prov
    {
        kProvNone       = -1, // no reco leg
        kProvLead       =  0, // holds the harder truth leg
        kProvSub        =  1, // holds the softer truth leg
        kProvOtherTruth =  2, // never assigned here
        kProvUE         =  3  // never assigned here
    };

    // Why a truth leg is out of acceptance; zero means it is in.
    // kFailDphi is a pair-level condition, set on both legs together.
    enum TruthFail
    {
        kFailNone = 0,
        kFailPt   = 1,
        kFailEta  = 2,
        kFailDphi = 4
    };

    struct Config
    {
        float truth_pt_thresh[2] = { 14.0f, 7.0f };
        float reco_pt_thresh[2]  = { 25.0f, 8.0f };
        float min_dphi           = 7.0f * TMath::Pi() / 8.0f;
    };

    struct Result
    {
        int   category      = kSkip;
        int   truth_in_acc  = 0;
        int   truth_fail[2] = { kFailNone, kFailNone };
        int   truth_idx[2]  = { -1, -1 }; // original truth indices
        float truth_dphi    = -999.0f;
        int   truth_match_idx[2] = { -1, -1 };
        float truth_match_pt[2]  = { -999.0f, -999.0f };

        int   reco_pair         = 0;
        int   reco_slot[2]      = { -1, -1 }; // slots in match_*
        int   reco_idx[2]       = { -1, -1 }; // original reco indices
        float reco_dphi         = -999.0f;
        int   reco_in_acc       = 0;
        int   reco_prov[2]      = { kProvNone, kProvNone };
        int   reco_truth_idx[2] = { -1, -1 };
        int   legs_swapped      = 0;
        int   n_accepted_reco   = 0;
        // which of the two reco legs (0 or 1) holds truth leg 1 -- the
        // bookkeeping legs_swapped is derived from
        int   truth_lead_leg    = 0;
    };

    // The event's reco jets available for pairing, as slots in the match_*
    // collection, pT descending. match_reco_pT is pT-descending on input,
    // so accepted[0] is the leading matched reco jet in the calorimeter
    // acceptance.
    //
    // Unlike dijet_pair_matching_v3.C's list, NO pT threshold is applied
    // here: the thresholds decide whether the formed pair is a reco dijet
    // CANDIDATE (reco_pair / reco_in_acc), not whether the pairing exists.
    // Dropping soft jets from the list would silently delete the Miss
    // population in which the subleading leg fell below 8 GeV.
    inline void build_accepted_reco(
        const std::vector< float > & match_reco_E,
        const std::vector< int >   & match_reco_accept_eta,
        std::vector< int > & accepted
    )
    {
        accepted.clear();
        for ( size_t s = 0; s < match_reco_accept_eta.size(); ++s )
        {
            if ( match_reco_accept_eta.at( s ) && match_reco_E.at( s ) > 0.0f )
            {
                accepted.push_back( static_cast< int >( s ) );
            }
        }
    }

    // Classifies the pairing ( accepted[0], accepted[partner] ). partner
    // == 1 is the leading pair (pT1,pT2), partner == 2 is (pT1,pT3), and
    // so on.
    //
    // `accepted` is the event-level list from build_accepted_reco(),
    // computed once by the caller and reused across every pairing.
    //
    // Returns false when the pairing cannot be formed at all (fewer than
    // partner+1 entries in `accepted`) or when it is a kSkip, in which
    // case nothing should be written. Pure otherwise: no ROOT I/O state.
    inline bool classify(
        const int partner,
        const std::vector< int >   & accepted,
        const std::vector< int >   & match_reco_idx,
        const std::vector< float > & match_reco_pT,
        const std::vector< float > & match_reco_phi,
        const std::vector< int >   & match_truth_idx,
        const std::vector< float > & match_truth_pT,
        const std::vector< float > & match_truth_phi,
        const std::vector< int >   & match_truth_accept_eta,
        const Config & cfg,
        Result & r
    )
    {
        if ( partner < 1 || partner >= static_cast< int >( accepted.size() ) ) return false;

        r.n_accepted_reco = static_cast< int >( accepted.size() );

        //------------------------------------------------------------
        // reco pair, in reco pT order. The pair is always FORMED; the pT
        // thresholds set reco_pair and the dphi requirement then sets
        // reco_in_acc -- whether it is a reco dijet CANDIDATE.
        //------------------------------------------------------------
        const int s1 = accepted[0];
        const int s2 = accepted[partner];

        r.reco_slot[0] = s1;
        r.reco_slot[1] = s2;
        r.reco_idx[0]  = match_reco_idx.at( s1 );
        r.reco_idx[1]  = match_reco_idx.at( s2 );
        r.reco_dphi    = AnaUtils::dphi_wrap( match_reco_phi.at( s1 ), match_reco_phi.at( s2 ) );

        r.reco_pair = ( match_reco_pT.at( s1 ) > cfg.reco_pt_thresh[0]
                     && match_reco_pT.at( s2 ) > cfg.reco_pt_thresh[1] ) ? 1 : 0;
        r.reco_in_acc = ( r.reco_pair && r.reco_dphi >= cfg.min_dphi ) ? 1 : 0;

        //------------------------------------------------------------
        // truth pair: the two legs' own truth partners, ordered by truth
        // pT so truth leg 1 is the harder one and the 14/7 GeV thresholds
        // land on the right legs. The reco legs stay in reco pT order, so
        // a crossing between the two orderings is legs_swapped -- a
        // resolution effect between two real legs, not a failure.
        //------------------------------------------------------------
        const bool swapped = ( match_truth_pT.at( s2 ) > match_truth_pT.at( s1 ) );
        r.truth_lead_leg = swapped ? 1 : 0;
        r.legs_swapped   = swapped ? 1 : 0;

        const int ts[2] = { swapped ? s2 : s1, swapped ? s1 : s2 };

        r.truth_idx[0] = match_truth_idx.at( ts[0] );
        r.truth_idx[1] = match_truth_idx.at( ts[1] );

        r.truth_dphi = AnaUtils::dphi_wrap( match_truth_phi.at( ts[0] ),
                                            match_truth_phi.at( ts[1] ) );
        const bool truth_dphi_ok = ( r.truth_dphi > cfg.min_dphi );

        for ( int leg = 0; leg < 2; ++leg )
        {
            int fail = kFailNone;
            if ( !( match_truth_pT.at( ts[leg] ) > cfg.truth_pt_thresh[leg] ) ) fail |= kFailPt;
            if ( !match_truth_accept_eta.at( ts[leg] ) )                        fail |= kFailEta;
            if ( !truth_dphi_ok )                                               fail |= kFailDphi;
            r.truth_fail[leg] = fail;

            // each truth leg's own match is the reco leg it came from --
            // the matching is 1-to-1, so there is nothing else it could be.
            r.truth_match_idx[leg] = match_reco_idx.at( ts[leg] );
            r.truth_match_pt[leg]  = match_reco_pT.at( ts[leg] );
        }

        r.truth_in_acc = ( r.truth_fail[0] == kFailNone && r.truth_fail[1] == kFailNone ) ? 1 : 0;

        // reco leg -> which truth leg it holds. Always {Lead,Sub} or
        // {Sub,Lead}: the truth pair IS this reco pair's own matches.
        r.reco_prov[0] = swapped ? kProvSub  : kProvLead;
        r.reco_prov[1] = swapped ? kProvLead : kProvSub;
        r.reco_truth_idx[0] = match_truth_idx.at( s1 );
        r.reco_truth_idx[1] = match_truth_idx.at( s2 );

        //------------------------------------------------------------
        // pair-level category -- three populations, kFakeMiss impossible
        //------------------------------------------------------------
        const bool has_truth_dijet_candidate = ( r.truth_in_acc == 1 );
        const bool has_reco_dijet_candidate  = ( r.reco_in_acc == 1 );

        if ( !has_truth_dijet_candidate && !has_reco_dijet_candidate ) return false;
        if      ( has_truth_dijet_candidate && has_reco_dijet_candidate ) r.category = kFill;
        else if ( has_truth_dijet_candidate )                             r.category = kMiss;
        else                                                              r.category = kFake;

        return true;
    }
}

int dijet_pair_matching_inclusive_v3(
    const std::string & infile  = "slim_matched_jets.root",
    const std::string & outfile = "dijet_pair_matching_inclusive_v3.root",
    const int max_partner_rank  = -1,   // -1 = every reco jet; else cap the
                                        // pairing at (1, max_partner_rank)
    const float reco_pt1_thresh  = 25.0f,
    const float reco_pt2_thresh  = 8.0f,
    const float truth_pt1_thresh = 14.0f,
    const float truth_pt2_thresh = 7.0f,
    const long  max_entries      = -1   // -1 = the whole chain; else stop early
)
{
    auto * t = new TChain( "T" );
    // infile can be either a single merged .root file (as written by
    // slim_matched_jets.C) or a .list file listing several of them --
    // AnaUtils::getFilelist always reads its argument as a text file, so a
    // .root file must be added to the chain directly rather than routed
    // through it.
    if ( infile.size() >= 5 && infile.compare( infile.size() - 5, 5, ".root" ) == 0 )
    {
        std::cout << "Adding file: " << infile << std::endl;
        t -> Add( infile.c_str() );
    }
    else
    {
        auto files = AnaUtils::getFilelist( infile, ".root" );
        for ( const auto & file : files )
        {
            std::cout << "Adding file: " << file << std::endl;
            t -> Add( file.c_str() );
        }
    }

    int nentries = t -> GetEntries();
    if ( max_entries > 0 && max_entries < nentries ) nentries = static_cast< int >( max_entries );
    if ( nentries < 1 )
    {
        std::cerr << "Error: No entries in TChain." << std::endl;
        return -1;
    }
    std::cout << "Total entries in TChain: " << t -> GetEntries();
    if ( max_entries > 0 && max_entries < t -> GetEntries() )
    {
        std::cout << "  (capped at " << nentries << ")";
    }
    std::cout << std::endl;

    t -> SetBranchStatus( "*", false );

    auto has_branch = [ t ]( const std::string & name )
    {
        return t -> GetBranch( name.c_str() ) != nullptr;
    };
    auto enable = [ t ]( const std::string & name, void * addr )
    {
        t -> SetBranchStatus ( name.c_str(), true );
        t -> SetBranchAddress( name.c_str(), addr );
    };

    //----------------------------------------------------------------
    // event-level passthrough
    //----------------------------------------------------------------
    int   event_id = -1, cent = -1, is_minbias = 0;
    float zvrtx = 0.0, mbd_q = -999.0, sumeT = -999.0;
    float psi2 = -999.0, psi3 = -999.0;
    int   n_truth_jets = 0, n_reco_jets = 0;

    const bool has_event_id = has_branch( "event_id" );
    const bool has_cent     = has_branch( "cent" );
    const bool has_zvrtx    = has_branch( "zvrtx" );
    const bool has_minbias  = has_branch( "is_minbias" );
    const bool has_mbd_q    = has_branch( "mbd_q" );
    const bool has_sumeT    = has_branch( "sumeT" );
    const bool has_psi2     = has_branch( "psi2" );
    const bool has_psi3     = has_branch( "psi3" );

    if ( has_event_id ) enable( "event_id", &event_id );
    if ( has_cent )     enable( "cent", &cent );
    if ( has_zvrtx )    enable( "zvrtx", &zvrtx );
    if ( has_minbias )  enable( "is_minbias", &is_minbias );
    if ( has_mbd_q )    enable( "mbd_q", &mbd_q );
    if ( has_sumeT )    enable( "sumeT", &sumeT );
    if ( has_psi2 )     enable( "psi2", &psi2 );
    if ( has_psi3 )     enable( "psi3", &psi3 );

    //----------------------------------------------------------------
    // the match_* collection of the slim tree -- the whole reco-side and,
    // through the truth partner carried alongside each reco jet, the whole
    // truth side of a reco-driven pairing.
    //----------------------------------------------------------------
    std::vector< int >   * match_reco_idx         = nullptr;
    std::vector< int >   * match_truth_idx        = nullptr;
    std::vector< float > * match_reco_pT          = nullptr;
    std::vector< float > * match_reco_eta         = nullptr;
    std::vector< float > * match_reco_phi         = nullptr;
    std::vector< float > * match_reco_E           = nullptr;
    std::vector< float > * match_reco_unsub_pT    = nullptr;
    std::vector< int >   * match_reco_accept_eta  = nullptr;
    std::vector< float > * match_truth_pT         = nullptr;
    std::vector< float > * match_truth_eta        = nullptr;
    std::vector< float > * match_truth_phi        = nullptr;
    std::vector< float > * match_truth_E          = nullptr;
    std::vector< int >   * match_truth_accept_eta = nullptr;
    std::vector< int >   * match_truth_flavor     = nullptr;
    std::vector< float > * match_truth_parton_pT  = nullptr;
    std::vector< float > * match_dr               = nullptr;

    const bool has_all = has_branch( "match_reco_idx" )
                      && has_branch( "match_truth_idx" )
                      && has_branch( "match_reco_pT" )
                      && has_branch( "match_reco_accept_eta" )
                      && has_branch( "match_truth_accept_eta" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the slim match_* "
                  << "collection written by slim_matched_jets.C." << std::endl;
        return -1;
    }

    const bool has_flavor = has_branch( "match_truth_flavor" );
    const bool has_unsub  = has_branch( "match_reco_unsub_pT" );
    const bool has_dr     = has_branch( "match_dr" );

    if ( has_branch( "n_truth_jets" ) ) enable( "n_truth_jets", &n_truth_jets );
    if ( has_branch( "n_reco_jets" ) )  enable( "n_reco_jets", &n_reco_jets );

    enable( "match_reco_idx", &match_reco_idx );
    enable( "match_truth_idx", &match_truth_idx );
    enable( "match_reco_pT", &match_reco_pT );
    enable( "match_reco_eta", &match_reco_eta );
    enable( "match_reco_phi", &match_reco_phi );
    enable( "match_reco_E", &match_reco_E );
    enable( "match_reco_accept_eta", &match_reco_accept_eta );
    enable( "match_truth_pT", &match_truth_pT );
    enable( "match_truth_eta", &match_truth_eta );
    enable( "match_truth_phi", &match_truth_phi );
    enable( "match_truth_E", &match_truth_E );
    enable( "match_truth_accept_eta", &match_truth_accept_eta );
    if ( has_unsub ) enable( "match_reco_unsub_pT", &match_reco_unsub_pT );
    if ( has_dr )    enable( "match_dr", &match_dr );
    if ( has_flavor )
    {
        enable( "match_truth_flavor", &match_truth_flavor );
        enable( "match_truth_parton_pT", &match_truth_parton_pT );
    }

    //----------------------------------------------------------------
    // output -- branch-for-branch the dijet_pair_matching_inclusive_v2.C
    // schema, with the v3-only branches appended at the end.
    //----------------------------------------------------------------
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "inclusive-pairing dijet matching / response (slim input)" );

    int   o_event_id = -1, o_cent = -1, o_is_minbias = 0;
    float o_zvrtx = 0.0, o_mbd_q = -999.0, o_sumeT = -999.0;
    float o_psi2 = -999.0, o_psi3 = -999.0, o_dpsi2 = -999.0;

    int   o_pair_partner_rank = -1, o_n_truth_jets = 0;
    int   o_category = -1, o_truth_in_acc = 0, o_reco_pair = 0, o_reco_in_acc = 0;
    int   o_truth_fail1 = 0, o_truth_fail2 = 0;
    int   o_truth_idx1 = -1, o_truth_idx2 = -1;
    float o_truth_pt1 = -999.0, o_truth_eta1 = -999.0, o_truth_phi1 = -999.0, o_truth_E1 = -999.0;
    float o_truth_pt2 = -999.0, o_truth_eta2 = -999.0, o_truth_phi2 = -999.0, o_truth_E2 = -999.0;
    float o_truth_dphi = -999.0, o_truth_xj = -999.0;
    int   o_truth_match_idx1 = -1, o_truth_match_idx2 = -1;
    float o_truth_match_pt1 = -999.0, o_truth_match_pt2 = -999.0;
    int   o_reco_idx1 = -1, o_reco_idx2 = -1;
    float o_reco_pt1 = -999.0, o_reco_eta1 = -999.0, o_reco_phi1 = -999.0, o_reco_E1 = -999.0;
    float o_reco_pt2 = -999.0, o_reco_eta2 = -999.0, o_reco_phi2 = -999.0, o_reco_E2 = -999.0;
    float o_reco_dphi = -999.0, o_reco_xj = -999.0;
    int   o_reco_prov1 = -1, o_reco_prov2 = -1;
    int   o_reco_truth_idx1 = -1, o_reco_truth_idx2 = -1;
    int   o_legs_swapped = 0, o_n_accepted_reco = 0;

    // v3 additions
    int   o_n_reco_jets = 0, o_n_match = 0;
    int   o_truth_flavor1 = -999, o_truth_flavor2 = -999;
    float o_truth_parton_pt1 = -999.0, o_truth_parton_pt2 = -999.0;
    float o_reco_unsub_pt1 = -999.0, o_reco_unsub_pt2 = -999.0;
    float o_reco_dr1 = -999.0, o_reco_dr2 = -999.0;

    if ( has_event_id ) tout -> Branch( "event_id", &o_event_id, "event_id/I" );
    if ( has_cent )     tout -> Branch( "cent", &o_cent, "cent/I" );
    if ( has_zvrtx )    tout -> Branch( "zvrtx", &o_zvrtx, "zvrtx/F" );
    if ( has_minbias )  tout -> Branch( "is_minbias", &o_is_minbias, "is_minbias/I" );
    if ( has_mbd_q )    tout -> Branch( "mbd_q", &o_mbd_q, "mbd_q/F" );
    if ( has_sumeT )    tout -> Branch( "sumeT", &o_sumeT, "sumeT/F" );
    if ( has_psi2 )
    {
        tout -> Branch( "psi2", &o_psi2, "psi2/F" );
        // dpsi2 = wrapped angle between truth leg 1 and the 2nd-order
        // event plane; see AnaUtils::get_dpsi2.
        tout -> Branch( "dpsi2", &o_dpsi2, "dpsi2/F" );
    }
    if ( has_psi3 )     tout -> Branch( "psi3", &o_psi3, "psi3/F" );

    // k+1, so rank 2 is the leading pair (pT1,pT2), rank 3 is (pT1,pT3).
    // Here k is the RECO rank of the subleading leg -- the pairing is
    // reco-driven.
    tout -> Branch( "pair_partner_rank", &o_pair_partner_rank, "pair_partner_rank/I" );
    tout -> Branch( "n_truth_jets", &o_n_truth_jets, "n_truth_jets/I" );
    tout -> Branch( "category", &o_category, "category/I" );
    tout -> Branch( "truth_in_acc", &o_truth_in_acc, "truth_in_acc/I" );
    tout -> Branch( "reco_pair", &o_reco_pair, "reco_pair/I" );
    // reco_pair alone is two legs above threshold; reco_in_acc is the
    // reco dijet CANDIDATE, i.e. those two legs plus the dphi
    // requirement. The category keys on reco_in_acc.
    tout -> Branch( "reco_in_acc", &o_reco_in_acc, "reco_in_acc/I" );
    // bitmask: 1 = pT below threshold, 2 = outside truth accept_eta,
    // 4 = pair fails the dphi requirement (set on both legs together).
    tout -> Branch( "truth_fail1", &o_truth_fail1, "truth_fail1/I" );
    tout -> Branch( "truth_fail2", &o_truth_fail2, "truth_fail2/I" );
    tout -> Branch( "truth_idx1", &o_truth_idx1, "truth_idx1/I" );
    tout -> Branch( "truth_idx2", &o_truth_idx2, "truth_idx2/I" );
    tout -> Branch( "truth_pt1", &o_truth_pt1, "truth_pt1/F" );
    tout -> Branch( "truth_eta1", &o_truth_eta1, "truth_eta1/F" );
    tout -> Branch( "truth_phi1", &o_truth_phi1, "truth_phi1/F" );
    tout -> Branch( "truth_E1", &o_truth_E1, "truth_E1/F" );
    tout -> Branch( "truth_pt2", &o_truth_pt2, "truth_pt2/F" );
    tout -> Branch( "truth_eta2", &o_truth_eta2, "truth_eta2/F" );
    tout -> Branch( "truth_phi2", &o_truth_phi2, "truth_phi2/F" );
    tout -> Branch( "truth_E2", &o_truth_E2, "truth_E2/F" );
    tout -> Branch( "truth_dphi", &o_truth_dphi, "truth_dphi/F" );
    tout -> Branch( "truth_xj", &o_truth_xj, "truth_xj/F" );
    tout -> Branch( "truth_match_idx1", &o_truth_match_idx1, "truth_match_idx1/I" );
    tout -> Branch( "truth_match_idx2", &o_truth_match_idx2, "truth_match_idx2/I" );
    tout -> Branch( "truth_match_pt1", &o_truth_match_pt1, "truth_match_pt1/F" );
    tout -> Branch( "truth_match_pt2", &o_truth_match_pt2, "truth_match_pt2/F" );
    tout -> Branch( "reco_idx1", &o_reco_idx1, "reco_idx1/I" );
    tout -> Branch( "reco_idx2", &o_reco_idx2, "reco_idx2/I" );
    tout -> Branch( "reco_pt1", &o_reco_pt1, "reco_pt1/F" );
    tout -> Branch( "reco_eta1", &o_reco_eta1, "reco_eta1/F" );
    tout -> Branch( "reco_phi1", &o_reco_phi1, "reco_phi1/F" );
    tout -> Branch( "reco_E1", &o_reco_E1, "reco_E1/F" );
    tout -> Branch( "reco_pt2", &o_reco_pt2, "reco_pt2/F" );
    tout -> Branch( "reco_eta2", &o_reco_eta2, "reco_eta2/F" );
    tout -> Branch( "reco_phi2", &o_reco_phi2, "reco_phi2/F" );
    tout -> Branch( "reco_E2", &o_reco_E2, "reco_E2/F" );
    tout -> Branch( "reco_dphi", &o_reco_dphi, "reco_dphi/F" );
    tout -> Branch( "reco_xj", &o_reco_xj, "reco_xj/F" );
    // 0 = this reco leg holds truth leg 1, 1 = it holds truth leg 2.
    // Always {0,1} or {1,0} here -- the truth pair is this reco pair's
    // own two matches.
    tout -> Branch( "reco_prov1", &o_reco_prov1, "reco_prov1/I" );
    tout -> Branch( "reco_prov2", &o_reco_prov2, "reco_prov2/I" );
    tout -> Branch( "reco_truth_idx1", &o_reco_truth_idx1, "reco_truth_idx1/I" );
    tout -> Branch( "reco_truth_idx2", &o_reco_truth_idx2, "reco_truth_idx2/I" );
    tout -> Branch( "legs_swapped", &o_legs_swapped, "legs_swapped/I" );
    tout -> Branch( "n_accepted_reco", &o_n_accepted_reco, "n_accepted_reco/I" );

    // ---- v3 additions, appended so the v2 schema above is untouched ----
    tout -> Branch( "n_reco_jets", &o_n_reco_jets, "n_reco_jets/I" );
    tout -> Branch( "n_match", &o_n_match, "n_match/I" );
    if ( has_flavor )
    {
        tout -> Branch( "truth_flavor1", &o_truth_flavor1, "truth_flavor1/I" );
        tout -> Branch( "truth_flavor2", &o_truth_flavor2, "truth_flavor2/I" );
        tout -> Branch( "truth_parton_pt1", &o_truth_parton_pt1, "truth_parton_pt1/F" );
        tout -> Branch( "truth_parton_pt2", &o_truth_parton_pt2, "truth_parton_pt2/F" );
    }
    if ( has_unsub )
    {
        tout -> Branch( "reco_unsub_pt1", &o_reco_unsub_pt1, "reco_unsub_pt1/F" );
        tout -> Branch( "reco_unsub_pt2", &o_reco_unsub_pt2, "reco_unsub_pt2/F" );
    }
    if ( has_dr )
    {
        tout -> Branch( "reco_dr1", &o_reco_dr1, "reco_dr1/F" );
        tout -> Branch( "reco_dr2", &o_reco_dr2, "reco_dr2/F" );
    }

    DijetPairInclusiveV3::Config cfg;
    cfg.truth_pt_thresh[0] = truth_pt1_thresh;
    cfg.truth_pt_thresh[1] = truth_pt2_thresh;
    cfg.reco_pt_thresh[0]  = reco_pt1_thresh;
    cfg.reco_pt_thresh[1]  = reco_pt2_thresh;

    std::cout << "Truth acceptance: pT_1(2) > " << cfg.truth_pt_thresh[0]
              << "(" << cfg.truth_pt_thresh[1] << ") GeV, accept_eta, dphi > "
              << cfg.min_dphi << std::endl;
    std::cout << "Reco  acceptance: pT_1(2) > " << cfg.reco_pt_thresh[0]
              << "(" << cfg.reco_pt_thresh[1] << ") GeV, accept_eta, E > 0, dphi >= "
              << cfg.min_dphi << std::endl;

    long n_events = 0, n_rows = 0;
    // one slot per Category -- ++n_cat[r.category] indexes it with the raw
    // enum value, so it stays sized to the full enum even though only
    // 0/1/2 can be filled here.
    long n_cat[6] = { 0, 0, 0, 0, 0, 0 };
    // same, restricted to the leading pairing (rank 2) so it can be lined
    // up against dijet_pair_matching_v3.C on the same input
    long n_cat_r2[6] = { 0, 0, 0, 0, 0, 0 };
    long n_swapped = 0;
    // kMiss breakdown: why the reco pair was not a candidate
    long n_miss_reco_pt = 0, n_miss_reco_dphi = 0;
    // kFake breakdown: why the truth pair was out of acceptance
    long n_fake_truth_eta = 0, n_fake_truth_pt = 0, n_fake_truth_dphi = 0;

    // reused across events so the per-event rebuild does not reallocate
    std::vector< int > accepted;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        //------------------------------------------------------------
        // event-level reco jet list, computed once and reused across every
        // pairing in this event. The pairing rank indexes straight into
        // it, accepted[0] being the leading matched reco jet in the
        // calorimeter acceptance.
        //------------------------------------------------------------
        DijetPairInclusiveV3::build_accepted_reco(
            *match_reco_E, *match_reco_accept_eta, accepted );

        if ( accepted.size() < 2 ) continue; // no pairing possible
        ++n_events;

        int last_partner = static_cast< int >( accepted.size() ) - 1;
        if ( max_partner_rank > 0 && ( max_partner_rank - 1 ) < last_partner )
        {
            last_partner = max_partner_rank - 1;
        }

        for ( int partner = 1; partner <= last_partner; ++partner )
        {
            DijetPairInclusiveV3::Result r;
            if ( !DijetPairInclusiveV3::classify(
                     partner, accepted,
                     *match_reco_idx, *match_reco_pT, *match_reco_phi,
                     *match_truth_idx, *match_truth_pT, *match_truth_phi,
                     *match_truth_accept_eta,
                     cfg, r ) )
            {
                continue; // kSkip -- neither candidate in acceptance
            }

            ++n_rows;
            ++n_cat[ r.category ];
            if ( partner == 1 ) ++n_cat_r2[ r.category ];
            if ( r.category == DijetPairInclusiveV3::kFill ) n_swapped += r.legs_swapped;

            if ( r.category == DijetPairInclusiveV3::kMiss )
            {
                if ( !r.reco_pair )                       ++n_miss_reco_pt;
                else if ( r.reco_dphi < cfg.min_dphi )    ++n_miss_reco_dphi;
            }
            else if ( r.category == DijetPairInclusiveV3::kFake )
            {
                const int f = r.truth_fail[0] | r.truth_fail[1];
                if ( f & DijetPairInclusiveV3::kFailEta )  ++n_fake_truth_eta;
                if ( f & DijetPairInclusiveV3::kFailPt )   ++n_fake_truth_pt;
                if ( f & DijetPairInclusiveV3::kFailDphi ) ++n_fake_truth_dphi;
            }

            const int s1 = r.reco_slot[0];
            const int s2 = r.reco_slot[1];
            // the two match_* slots in TRUTH pT order -- truth leg 1 first
            const int t1 = ( r.truth_lead_leg == 0 ) ? s1 : s2;
            const int t2 = ( r.truth_lead_leg == 0 ) ? s2 : s1;

            o_event_id   = event_id;
            o_cent       = cent;
            o_zvrtx      = zvrtx;
            o_is_minbias = is_minbias;
            o_mbd_q      = mbd_q;
            o_sumeT      = sumeT;
            o_psi2       = psi2;
            o_psi3       = psi3;
            o_dpsi2      = has_psi2
                         ? AnaUtils::get_dpsi2( psi2, match_truth_phi -> at( t1 ) )
                         : -999.0f;

            o_pair_partner_rank = partner + 1;
            o_n_truth_jets = n_truth_jets;
            o_category     = r.category;
            o_truth_in_acc = r.truth_in_acc;
            o_reco_pair    = r.reco_pair;
            o_reco_in_acc  = r.reco_in_acc;
            o_truth_fail1  = r.truth_fail[0];
            o_truth_fail2  = r.truth_fail[1];
            o_truth_idx1   = r.truth_idx[0];
            o_truth_idx2   = r.truth_idx[1];
            o_truth_dphi   = r.truth_dphi;

            o_truth_pt1  = match_truth_pT  -> at( t1 );
            o_truth_eta1 = match_truth_eta -> at( t1 );
            o_truth_phi1 = match_truth_phi -> at( t1 );
            o_truth_E1   = match_truth_E   -> at( t1 );
            o_truth_pt2  = match_truth_pT  -> at( t2 );
            o_truth_eta2 = match_truth_eta -> at( t2 );
            o_truth_phi2 = match_truth_phi -> at( t2 );
            o_truth_E2   = match_truth_E   -> at( t2 );
            o_truth_xj   = ( o_truth_pt1 > 0.0f ) ? o_truth_pt2 / o_truth_pt1 : -999.0f;

            o_truth_match_idx1 = r.truth_match_idx[0];
            o_truth_match_idx2 = r.truth_match_idx[1];
            o_truth_match_pt1  = r.truth_match_pt[0];
            o_truth_match_pt2  = r.truth_match_pt[1];

            o_reco_idx1       = r.reco_idx[0];
            o_reco_idx2       = r.reco_idx[1];
            o_reco_prov1      = r.reco_prov[0];
            o_reco_prov2      = r.reco_prov[1];
            o_reco_truth_idx1 = r.reco_truth_idx[0];
            o_reco_truth_idx2 = r.reco_truth_idx[1];
            o_reco_dphi       = r.reco_dphi;
            o_legs_swapped    = r.legs_swapped;
            o_n_accepted_reco = r.n_accepted_reco;

            o_reco_pt1  = match_reco_pT  -> at( s1 );
            o_reco_eta1 = match_reco_eta -> at( s1 );
            o_reco_phi1 = match_reco_phi -> at( s1 );
            o_reco_E1   = match_reco_E   -> at( s1 );
            o_reco_pt2  = match_reco_pT  -> at( s2 );
            o_reco_eta2 = match_reco_eta -> at( s2 );
            o_reco_phi2 = match_reco_phi -> at( s2 );
            o_reco_E2   = match_reco_E   -> at( s2 );
            o_reco_xj   = ( o_reco_pt1 > 0.0f ) ? o_reco_pt2 / o_reco_pt1 : -999.0f;

            if ( has_flavor )
            {
                o_truth_flavor1    = match_truth_flavor    -> at( t1 );
                o_truth_flavor2    = match_truth_flavor    -> at( t2 );
                o_truth_parton_pt1 = match_truth_parton_pT -> at( t1 );
                o_truth_parton_pt2 = match_truth_parton_pT -> at( t2 );
            }
            if ( has_unsub )
            {
                o_reco_unsub_pt1 = match_reco_unsub_pT -> at( s1 );
                o_reco_unsub_pt2 = match_reco_unsub_pT -> at( s2 );
            }
            if ( has_dr )
            {
                o_reco_dr1 = match_dr -> at( s1 );
                o_reco_dr2 = match_dr -> at( s2 );
            }

            o_n_reco_jets = n_reco_jets;
            o_n_match     = static_cast< int >( match_reco_pT -> size() );

            tout -> Fill();
        }
    }

    std::cout << "\n---- pair-level dijet categories (inclusive, reco-driven) ----" << std::endl;
    std::cout << "Events with at least 2 accepted matched reco jets: " << n_events << std::endl;
    std::cout << "Rows written: " << n_rows
              << "  (pairings that are neither candidate are dropped as Skip)" << std::endl;
    std::cout << "  Fill (0): " << n_cat[0] << "   <- the only category that fills the response" << std::endl;
    std::cout << "            of which reco legs pT-swapped w.r.t. truth: " << n_swapped << std::endl;
    std::cout << "  Miss (1): " << n_cat[1] << std::endl;
    std::cout << "            reco pair below threshold: " << n_miss_reco_pt
              << ", reco pair not back-to-back: " << n_miss_reco_dphi << std::endl;
    std::cout << "  Fake (2): " << n_cat[2] << std::endl;
    std::cout << "            truth pair failed eta: " << n_fake_truth_eta
              << ", pT: " << n_fake_truth_pt
              << ", dphi: " << n_fake_truth_dphi << std::endl;
    std::cout << "  FakeMiss (5): " << n_cat[5] << "   <- impossible by construction: the truth"
              << " pair IS the reco pair's own matches" << std::endl;

    std::cout << "  ---- leading pairing only (pair_partner_rank == 2) ----" << std::endl;
    std::cout << "  Fill: " << n_cat_r2[0]
              << ", Miss: " << n_cat_r2[1]
              << ", Fake: " << n_cat_r2[2] << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " rows to " << outfile << std::endl;

    return 0;
}

#endif
