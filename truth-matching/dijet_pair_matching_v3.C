#ifndef _DIJET_PAIR_MATCHING_V3_C_
#define _DIJET_PAIR_MATCHING_V3_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Pair-level dijet matching / response bookkeeping, EXCLUSIVE (leading
// pair only), run on the SLIMMED tree written by slim_matched_jets.C
// rather than on the full match_standalone.C tree.
//
// The output schema is branch-for-branch the one dijet_pair_matching_v2.C
// writes (same names, same category numbering), so createResponse_*.cxx
// and everything else downstream reads it unchanged. A handful of extra
// branches are appended at the end; nothing existing was renamed or
// dropped.
//
// THE ALGORITHM (one row per event)
//
//   truth dijet candidate = truth jets 0 and 1 -- the two hardest truth
//     jets in the event, NO eta or pT selection at selection time.
//       accepted when pT_1(2) > 14(7) GeV, truth_jet_accept_eta on both
//       legs, and |dphi_12| > 7pi/8.
//
//   reco dijet candidate = the two hardest MATCHED reco jets that are in
//     the calorimeter acceptance with E > 0 and pT > 8 GeV, the leading
//     one additionally above 25 GeV.
//       accepted when both legs exist and |dphi_12| >= 7pi/8.
//     Every reco jet in the slim tree is truth-matched by construction,
//     so -- exactly as in v2 -- an unmatched UE fluctuation can never
//     take a leg and kProvUE / kUESub are never assigned.
//
//   FOUR populations are written (kSkip is dropped, never filled):
//
//     truth cand | reco cand | reco IS the truth pair | category
//     -----------+-----------+------------------------+-----------
//         yes    |    yes    |          yes           | Fill      (0)
//         yes    |    no     |          --            | Miss      (1)
//         no     |    yes    |          --            | Fake      (2)
//         yes    |    yes    |          no            | FakeMiss  (5)
//         no     |    no     |          --            | (not written)
//
//   "reco IS the truth pair" means the two reco legs are truth jet 0's
//   and truth jet 1's own dR matches, in either order -- a resolution
//   driven pT-ordering flip between two real legs is a Fill, flagged by
//   legs_swapped, not a failure. FakeMiss is the slide-2 population:
//   both candidates are in acceptance but the reco candidate was built
//   from the wrong jets, so it is a fake (that reco dijet) AND a miss
//   (the real truth dijet) at once.
//
// READING THE SLIM TREE
//
// slim_matched_jets.C keeps two disjoint collections and drops
// everything else:
//
//   match_*     one entry per truth-matched reco jet, carrying the reco
//               jet AND its truth partner side by side, pT(reco)
//               descending. This is the complete reco-side input --
//               unmatched reco jets are gone, which is exactly the
//               selection v2 applied anyway.
//   acc_truth_* one entry per UNMATCHED truth jet inside the eta
//               acceptance.
//
// A truth jet is therefore recoverable from the slim tree unless it is
// BOTH unmatched AND outside the eta acceptance. n_truth_jets records the
// original multiplicity, so those jets are still *countable*: any index
// in [0,n_truth_jets) missing from both collections is one of them, and
// its acceptance flag is known to be 0 even though its kinematics are
// not. build_truth_view() below rebuilds the truth jet list this way and
// marks such entries known = 0; truth_known1 / truth_known2 carry the
// flag into the output, and a leg with known = 0 gets kFailEta (which is
// certain) but not kFailPt or kFailDphi (which are not measurable). Such
// a leg can never be part of a truth dijet candidate, so it only ever
// moves an event into Fake or out of the tree entirely.
//
// Input:  the tree written by slim_matched_jets.C.
// Output: one row per event that is a Fill, Miss, Fake or FakeMiss.
namespace DijetPairV3
{
    // Same numbering as DijetPairV2 -- createResponse_exclusive_v2_AA.cxx
    // hard-codes it.
    enum Category
    {
        kFill     = 0,
        kMiss     = 1,
        kFake     = 2,
        kSkip     = 3, // never written out
        kUESub    = 4, // never assigned (every reco jet is truth-matched)
        kFakeMiss = 5
    };

    // What a selected reco dijet leg actually is.
    enum Prov
    {
        kProvNone       = -1, // no reco leg
        kProvLead       =  0, // truth jet 0's match
        kProvSub        =  1, // truth jet 1's match
        kProvOtherTruth =  2, // matched to a softer truth jet
        kProvUE         =  3  // no truth match -- impossible on slim input
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

    // One truth jet, addressed by its ORIGINAL index in the
    // match_standalone.C collection. known = 0 means the jet was dropped
    // by the slimming (unmatched AND outside the eta acceptance): its
    // acceptance flag is known to be 0, its kinematics are not known.
    struct TruthJet
    {
        int   known      = 0;
        float pT         = -999.0f;
        float eta        = -999.0f;
        float phi        = -999.0f;
        float E          = -999.0f;
        int   accept_eta = 0;
        int   flavor     = -999;
        float parton_pT  = -999.0f;
        // slot in the match_* collection, -1 when this truth jet has no
        // reco match
        int   match_slot = -1;
    };

    // Rebuilds the event's truth jet list from the two slim collections.
    // n_truth is the ORIGINAL truth multiplicity (the n_truth_jets
    // branch), so indices dropped by the slimming stay addressable.
    inline void build_truth_view(
        const int n_truth,
        const std::vector< int >   & match_truth_idx,
        const std::vector< float > & match_truth_pT,
        const std::vector< float > & match_truth_eta,
        const std::vector< float > & match_truth_phi,
        const std::vector< float > & match_truth_E,
        const std::vector< int >   & match_truth_accept_eta,
        const std::vector< int >   * match_truth_flavor,
        const std::vector< float > * match_truth_parton_pT,
        const std::vector< int >   & acc_truth_idx,
        const std::vector< float > & acc_truth_pT,
        const std::vector< float > & acc_truth_eta,
        const std::vector< float > & acc_truth_phi,
        const std::vector< float > & acc_truth_E,
        const std::vector< int >   * acc_truth_flavor,
        const std::vector< float > * acc_truth_parton_pT,
        std::vector< TruthJet > & truth
    )
    {
        truth.assign( n_truth < 0 ? 0 : n_truth, TruthJet() );

        // matched truth jets, from collection (1)
        for ( size_t s = 0; s < match_truth_idx.size(); ++s )
        {
            const int ti = match_truth_idx.at( s );
            if ( ti < 0 || ti >= static_cast< int >( truth.size() ) ) continue;

            TruthJet & tj = truth[ti];
            tj.known      = 1;
            tj.pT         = match_truth_pT.at( s );
            tj.eta        = match_truth_eta.at( s );
            tj.phi        = match_truth_phi.at( s );
            tj.E          = match_truth_E.at( s );
            tj.accept_eta = match_truth_accept_eta.at( s );
            tj.match_slot = static_cast< int >( s );
            if ( match_truth_flavor )    tj.flavor    = match_truth_flavor    -> at( s );
            if ( match_truth_parton_pT ) tj.parton_pT = match_truth_parton_pT -> at( s );
        }

        // unmatched truth jets inside the eta acceptance, from collection
        // (2). Disjoint from (1) by construction, so nothing is
        // overwritten; the guard is defensive only.
        for ( size_t k = 0; k < acc_truth_idx.size(); ++k )
        {
            const int ti = acc_truth_idx.at( k );
            if ( ti < 0 || ti >= static_cast< int >( truth.size() ) ) continue;
            if ( truth[ti].known ) continue;

            TruthJet & tj = truth[ti];
            tj.known      = 1;
            tj.pT         = acc_truth_pT.at( k );
            tj.eta        = acc_truth_eta.at( k );
            tj.phi        = acc_truth_phi.at( k );
            tj.E          = acc_truth_E.at( k );
            tj.accept_eta = 1; // the defining cut of this collection
            tj.match_slot = -1;
            if ( acc_truth_flavor )    tj.flavor    = acc_truth_flavor    -> at( k );
            if ( acc_truth_parton_pT ) tj.parton_pT = acc_truth_parton_pT -> at( k );
        }
    }

    // The event's accepted reco jets, as slots in the match_* collection,
    // pT descending. match_reco_pT is pT-descending on input, so the
    // result is too: accepted[0] is the leading accepted matched reco jet.
    //
    // Every entry of match_* is truth-matched by construction, so this is
    // exactly the "accepted AND truth-matched" list v2 built by hand.
    inline void build_accepted_reco(
        const std::vector< float > & match_reco_pT,
        const std::vector< float > & match_reco_E,
        const std::vector< int >   & match_reco_accept_eta,
        const Config & cfg,
        std::vector< int > & accepted
    )
    {
        accepted.clear();
        for ( size_t s = 0; s < match_reco_pT.size(); ++s )
        {
            if (
                match_reco_accept_eta.at( s )
                && match_reco_E.at( s ) > 0.0f
                && match_reco_pT.at( s ) > cfg.reco_pt_thresh[1]
            )
            {
                accepted.push_back( static_cast< int >( s ) );
            }
        }
    }

    struct Result
    {
        int   category      = kSkip;
        int   truth_in_acc  = 0;
        int   truth_fail[2] = { kFailNone, kFailNone };
        int   truth_idx[2]  = { -1, -1 };
        int   truth_known[2] = { 0, 0 };
        float truth_dphi    = -999.0f;
        // each truth leg's own dR match, kept even when it is not one of
        // the selected reco dijet legs -- this is what shows "leg 2's real
        // match was at 6 GeV and a 12 GeV jet took its place".
        int   truth_match_idx[2]  = { -1, -1 }; // reco index of that match
        int   truth_match_slot[2] = { -1, -1 }; // its slot in match_*
        float truth_match_pt[2]   = { -999.0f, -999.0f };

        int   reco_pair         = 0;
        int   reco_slot[2]      = { -1, -1 }; // slots in match_*
        int   reco_idx[2]       = { -1, -1 }; // original reco indices
        float reco_dphi         = -999.0f;
        int   reco_in_acc       = 0;
        int   reco_prov[2]      = { kProvNone, kProvNone };
        int   reco_truth_idx[2] = { -1, -1 };
        int   legs_swapped      = 0;
        int   n_accepted_reco   = 0;
    };

    // Pure classification -- no ROOT I/O state. Returns false when the
    // event is a kSkip (neither candidate in acceptance), in which case
    // nothing should be written for it.
    inline bool classify(
        const std::vector< TruthJet > & truth,
        const std::vector< int >      & accepted,
        const std::vector< int >      & match_reco_idx,
        const std::vector< float >    & match_reco_pT,
        const std::vector< float >    & match_reco_phi,
        const std::vector< int >      & match_truth_idx,
        const Config & cfg,
        Result & r
    )
    {
        r.n_accepted_reco = static_cast< int >( accepted.size() );

        //------------------------------------------------------------
        // truth pair: always truth jets 0 and 1, NO eta cut at selection
        // time -- the acceptance decision is recorded, not applied.
        //------------------------------------------------------------
        if ( truth.size() >= 2 )
        {
            r.truth_idx[0] = 0;
            r.truth_idx[1] = 1;

            const bool both_known = ( truth[0].known && truth[1].known );
            if ( both_known )
            {
                r.truth_dphi = AnaUtils::dphi_wrap( truth[0].phi, truth[1].phi );
            }
            const bool truth_dphi_ok = both_known && ( r.truth_dphi > cfg.min_dphi );

            for ( int leg = 0; leg < 2; ++leg )
            {
                const TruthJet & tj = truth[ r.truth_idx[leg] ];
                r.truth_known[leg] = tj.known;

                int fail = kFailNone;
                if ( !tj.known )
                {
                    // dropped by the slimming => unmatched AND outside the
                    // eta acceptance. The eta failure is certain; the pT
                    // and dphi bits are not measurable, so they are left
                    // unset rather than guessed.
                    fail |= kFailEta;
                }
                else
                {
                    if ( !( tj.pT > cfg.truth_pt_thresh[leg] ) ) fail |= kFailPt;
                    if ( !tj.accept_eta )                        fail |= kFailEta;
                    if ( !truth_dphi_ok )                        fail |= kFailDphi;
                }
                r.truth_fail[leg] = fail;

                const int s = tj.match_slot;
                r.truth_match_slot[leg] = s;
                r.truth_match_idx[leg]  = ( s >= 0 ) ? match_reco_idx.at( s ) : -1;
                r.truth_match_pt[leg]   = ( s >= 0 ) ? match_reco_pT.at( s )  : -999.0f;
            }

            r.truth_in_acc = ( r.truth_fail[0] == kFailNone && r.truth_fail[1] == kFailNone ) ? 1 : 0;
        }

        const bool has_truth_dijet_candidate = ( r.truth_in_acc == 1 );

        //------------------------------------------------------------
        // reco pair: the two hardest accepted matched reco jets, built
        // the way the data analysis builds it -- by pT rank, NOT from the
        // truth legs' matches. This is what lets the wrong jet take a leg
        // and produce a FakeMiss.
        //------------------------------------------------------------
        int s1 = -1, s2 = -1;
        if ( !accepted.empty() && match_reco_pT.at( accepted[0] ) > cfg.reco_pt_thresh[0] )
        {
            s1 = accepted[0];
        }
        if ( accepted.size() > 1 && match_reco_pT.at( accepted[1] ) > cfg.reco_pt_thresh[1] )
        {
            s2 = accepted[1];
        }
        if ( s1 >= 0 && s2 >= 0 )
        {
            r.reco_pair = 1;
            r.reco_slot[0] = s1;
            r.reco_slot[1] = s2;
            r.reco_idx[0]  = match_reco_idx.at( s1 );
            r.reco_idx[1]  = match_reco_idx.at( s2 );
            r.reco_dphi    = AnaUtils::dphi_wrap( match_reco_phi.at( s1 ), match_reco_phi.at( s2 ) );

            for ( int leg = 0; leg < 2; ++leg )
            {
                const int ti = match_truth_idx.at( r.reco_slot[leg] );
                r.reco_truth_idx[leg] = ti;
                if      ( ti < 0 )  r.reco_prov[leg] = kProvUE;  // impossible here
                else if ( ti == 0 ) r.reco_prov[leg] = kProvLead;
                else if ( ti == 1 ) r.reco_prov[leg] = kProvSub;
                else                r.reco_prov[leg] = kProvOtherTruth;
            }
            r.reco_in_acc = ( r.reco_dphi >= cfg.min_dphi );
        }

        const bool has_reco_dijet_candidate = ( r.reco_in_acc == 1 );
        const bool reco_is_matched_to_truth =
               ( r.reco_prov[0] == kProvLead && r.reco_prov[1] == kProvSub )
            || ( r.reco_prov[0] == kProvSub  && r.reco_prov[1] == kProvLead );

        //------------------------------------------------------------
        // pair-level category
        //------------------------------------------------------------
        if ( !has_truth_dijet_candidate && !has_reco_dijet_candidate ) return false;
        if ( !has_truth_dijet_candidate && has_reco_dijet_candidate ) r.category = kFake;
        if ( has_truth_dijet_candidate )
        {
            if ( !has_reco_dijet_candidate )
            {
                r.category = kMiss;
            }
            else
            {
                r.category = reco_is_matched_to_truth ? kFill : kFakeMiss;
                r.legs_swapped = ( r.reco_prov[0] == kProvSub && r.reco_prov[1] == kProvLead ) ? 1 : 0;
            }
        }
        return true;
    }
}

int dijet_pair_matching_v3(
    const std::string & infile  = "slim_matched_jets.root",
    const std::string & outfile = "dijet_pair_matching_v3.root",
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
    // the two slim collections
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

    std::vector< int >   * acc_truth_idx          = nullptr;
    std::vector< float > * acc_truth_pT           = nullptr;
    std::vector< float > * acc_truth_eta          = nullptr;
    std::vector< float > * acc_truth_phi          = nullptr;
    std::vector< float > * acc_truth_E            = nullptr;
    std::vector< int >   * acc_truth_flavor       = nullptr;
    std::vector< float > * acc_truth_parton_pT    = nullptr;

    const bool has_all = has_branch( "n_truth_jets" )
                      && has_branch( "match_reco_idx" )
                      && has_branch( "match_truth_idx" )
                      && has_branch( "match_reco_pT" )
                      && has_branch( "match_reco_accept_eta" )
                      && has_branch( "match_truth_accept_eta" )
                      && has_branch( "acc_truth_idx" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the slim collections "
                  << "written by slim_matched_jets.C (match_* / acc_truth_* "
                  << "/ n_truth_jets)." << std::endl;
        return -1;
    }

    const bool has_flavor = has_branch( "match_truth_flavor" );
    const bool has_unsub  = has_branch( "match_reco_unsub_pT" );
    const bool has_dr     = has_branch( "match_dr" );

    enable( "n_truth_jets", &n_truth_jets );
    if ( has_branch( "n_reco_jets" ) ) enable( "n_reco_jets", &n_reco_jets );

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
    if ( has_unsub )  enable( "match_reco_unsub_pT", &match_reco_unsub_pT );
    if ( has_dr )     enable( "match_dr", &match_dr );
    if ( has_flavor )
    {
        enable( "match_truth_flavor", &match_truth_flavor );
        enable( "match_truth_parton_pT", &match_truth_parton_pT );
    }

    enable( "acc_truth_idx", &acc_truth_idx );
    enable( "acc_truth_pT", &acc_truth_pT );
    enable( "acc_truth_eta", &acc_truth_eta );
    enable( "acc_truth_phi", &acc_truth_phi );
    enable( "acc_truth_E", &acc_truth_E );
    if ( has_flavor )
    {
        enable( "acc_truth_flavor", &acc_truth_flavor );
        enable( "acc_truth_parton_pT", &acc_truth_parton_pT );
    }

    //----------------------------------------------------------------
    // output -- branch-for-branch the dijet_pair_matching_v2.C schema,
    // with the v3-only branches appended at the end.
    //----------------------------------------------------------------
    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "pair-level dijet matching / response (slim input)" );

    int   o_event_id = -1, o_cent = -1, o_is_minbias = 0;
    float o_zvrtx = 0.0, o_mbd_q = -999.0, o_sumeT = -999.0;
    float o_psi2 = -999.0, o_psi3 = -999.0, o_dpsi2 = -999.0;

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
    int   o_truth_known1 = 0, o_truth_known2 = 0;
    int   o_n_truth_jets = 0, o_n_reco_jets = 0, o_n_match = 0;
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
        // dpsi2 = wrapped angle between the truth leading jet and the
        // 2nd-order event plane; see AnaUtils::get_dpsi2.
        tout -> Branch( "dpsi2", &o_dpsi2, "dpsi2/F" );
    }
    if ( has_psi3 )     tout -> Branch( "psi3", &o_psi3, "psi3/F" );

    tout -> Branch( "category", &o_category, "category/I" );
    tout -> Branch( "truth_in_acc", &o_truth_in_acc, "truth_in_acc/I" );
    tout -> Branch( "reco_pair", &o_reco_pair, "reco_pair/I" );
    // reco_pair alone is two legs above threshold; reco_in_acc is the
    // reco dijet CANDIDATE, i.e. those two legs plus the dphi
    // requirement. The category keys on reco_in_acc.
    tout -> Branch( "reco_in_acc", &o_reco_in_acc, "reco_in_acc/I" );
    // bitmask: 1 = pT below threshold, 2 = outside truth_jet_accept_eta,
    // 4 = pair fails the dphi requirement (set on both legs together).
    // A leg with truth_known == 0 carries only bit 2 -- see the header.
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
    // each truth leg's own dR match, whether or not it ended up being one
    // of the selected reco dijet legs.
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
    // -1 no leg, 0 truth leg 1's match, 1 truth leg 2's match,
    //  2 a softer truth jet outside the leading pair, 3 no truth match
    //  (never on slim input).
    tout -> Branch( "reco_prov1", &o_reco_prov1, "reco_prov1/I" );
    tout -> Branch( "reco_prov2", &o_reco_prov2, "reco_prov2/I" );
    tout -> Branch( "reco_truth_idx1", &o_reco_truth_idx1, "reco_truth_idx1/I" );
    tout -> Branch( "reco_truth_idx2", &o_reco_truth_idx2, "reco_truth_idx2/I" );
    tout -> Branch( "legs_swapped", &o_legs_swapped, "legs_swapped/I" );
    tout -> Branch( "n_accepted_reco", &o_n_accepted_reco, "n_accepted_reco/I" );

    // ---- v3 additions, appended so the v2 schema above is untouched ----
    // 0 when the truth leg was dropped by the slimming (unmatched AND
    // outside the eta acceptance): its kinematics are -999 and only the
    // eta bit of truth_fail is set.
    tout -> Branch( "truth_known1", &o_truth_known1, "truth_known1/I" );
    tout -> Branch( "truth_known2", &o_truth_known2, "truth_known2/I" );
    tout -> Branch( "n_truth_jets", &o_n_truth_jets, "n_truth_jets/I" );
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

    DijetPairV3::Config cfg;
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

    long n_events = 0, n_read = 0;
    // one slot per Category, kFakeMiss (5) included -- sizing this to the
    // number of categories is load-bearing, ++n_cat[r.category] indexes it
    // with the raw enum value.
    long n_cat[6] = { 0, 0, 0, 0, 0, 0 };
    long n_swapped = 0;
    // events whose truth leg 1 or 2 was dropped by the slimming
    long n_truth_leg_unknown = 0;
    // kFake breakdown by why the truth pair was out of acceptance
    long n_fake_truth_eta = 0, n_fake_truth_pt = 0, n_fake_truth_dphi = 0;
    long n_fake_legs_are_truth = 0;
    // kFakeMiss breakdown: which reco leg is not the truth leg it should
    // be. Both legs are truth-matched by construction, so the culprit is
    // always a softer truth jet (kProvOtherTruth).
    long n_fm_lead_wrong = 0, n_fm_sub_wrong = 0, n_fm_both_wrong = 0;
    long n_fm_truth_lead_unmatched = 0, n_fm_truth_sub_unmatched = 0;

    std::vector< DijetPairV3::TruthJet > truth;
    std::vector< int > accepted;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        ++n_read;

        DijetPairV3::build_truth_view(
            n_truth_jets,
            *match_truth_idx, *match_truth_pT, *match_truth_eta, *match_truth_phi,
            *match_truth_E, *match_truth_accept_eta,
            match_truth_flavor, match_truth_parton_pT,
            *acc_truth_idx, *acc_truth_pT, *acc_truth_eta, *acc_truth_phi,
            *acc_truth_E, acc_truth_flavor, acc_truth_parton_pT,
            truth );

        DijetPairV3::build_accepted_reco(
            *match_reco_pT, *match_reco_E, *match_reco_accept_eta, cfg, accepted );

        DijetPairV3::Result r;
        if ( !DijetPairV3::classify(
            truth, accepted,
            *match_reco_idx, *match_reco_pT, *match_reco_phi, *match_truth_idx,
            cfg, r ) )
        {
            continue; // kSkip -- neither candidate in acceptance
        }

        ++n_events;
        ++n_cat[ r.category ];
        if ( r.category == DijetPairV3::kFill ) n_swapped += r.legs_swapped;
        if ( truth.size() >= 2 && ( !r.truth_known[0] || !r.truth_known[1] ) ) ++n_truth_leg_unknown;

        if ( r.category == DijetPairV3::kFakeMiss )
        {
            // slide-2 population: the truth pair IS in acceptance and a
            // reco dijet candidate IS in acceptance, but the candidate is
            // not the truth pair's own two matches -- so the event is a
            // fake (that reco dijet) AND a miss (the truth dijet) at once.
            const bool lead_ok = ( r.reco_prov[0] == DijetPairV3::kProvLead
                                || r.reco_prov[0] == DijetPairV3::kProvSub );
            const bool sub_ok  = ( r.reco_prov[1] == DijetPairV3::kProvLead
                                || r.reco_prov[1] == DijetPairV3::kProvSub );
            if      ( !lead_ok && !sub_ok ) ++n_fm_both_wrong;
            else if ( !lead_ok )            ++n_fm_lead_wrong;
            else                            ++n_fm_sub_wrong;

            // did the truth leg that lost its slot have a reco match at
            // all? "no" means the leg was never reconstructed; "yes"
            // means it was reconstructed but outranked.
            if ( r.truth_match_idx[0] < 0 ) ++n_fm_truth_lead_unmatched;
            if ( r.truth_match_idx[1] < 0 ) ++n_fm_truth_sub_unmatched;
        }
        else if ( r.category == DijetPairV3::kFake )
        {
            const int f = r.truth_fail[0] | r.truth_fail[1];
            if ( f & DijetPairV3::kFailEta )  ++n_fake_truth_eta;
            if ( f & DijetPairV3::kFailPt )   ++n_fake_truth_pt;
            if ( f & DijetPairV3::kFailDphi ) ++n_fake_truth_dphi;

            const bool both_are_truth_legs =
                   ( r.reco_prov[0] == DijetPairV3::kProvLead || r.reco_prov[0] == DijetPairV3::kProvSub )
                && ( r.reco_prov[1] == DijetPairV3::kProvLead || r.reco_prov[1] == DijetPairV3::kProvSub );
            if ( both_are_truth_legs ) ++n_fake_legs_are_truth;
        }

        o_event_id   = event_id;
        o_cent       = cent;
        o_zvrtx      = zvrtx;
        o_is_minbias = is_minbias;
        o_mbd_q      = mbd_q;
        o_sumeT      = sumeT;
        o_psi2       = psi2;
        o_psi3       = psi3;

        o_category     = r.category;
        o_truth_in_acc = r.truth_in_acc;
        o_reco_pair    = r.reco_pair;
        o_reco_in_acc  = r.reco_in_acc;
        o_truth_fail1  = r.truth_fail[0];
        o_truth_fail2  = r.truth_fail[1];
        o_truth_idx1   = r.truth_idx[0];
        o_truth_idx2   = r.truth_idx[1];
        o_truth_dphi   = r.truth_dphi;
        o_truth_known1 = r.truth_known[0];
        o_truth_known2 = r.truth_known[1];

        if ( r.truth_idx[0] >= 0 )
        {
            const DijetPairV3::TruthJet & t1 = truth[ r.truth_idx[0] ];
            const DijetPairV3::TruthJet & t2 = truth[ r.truth_idx[1] ];
            o_truth_pt1  = t1.pT;  o_truth_eta1 = t1.eta;
            o_truth_phi1 = t1.phi; o_truth_E1   = t1.E;
            o_truth_pt2  = t2.pT;  o_truth_eta2 = t2.eta;
            o_truth_phi2 = t2.phi; o_truth_E2   = t2.E;
            o_truth_xj   = ( o_truth_pt1 > 0.0f && o_truth_pt2 > 0.0f )
                         ? o_truth_pt2 / o_truth_pt1 : -999.0f;
            o_truth_flavor1    = t1.flavor;    o_truth_flavor2    = t2.flavor;
            o_truth_parton_pt1 = t1.parton_pT; o_truth_parton_pt2 = t2.parton_pT;
            o_dpsi2 = ( has_psi2 && t1.known ) ? AnaUtils::get_dpsi2( psi2, t1.phi ) : -999.0f;
        }
        else
        {
            o_truth_pt1 = o_truth_eta1 = o_truth_phi1 = o_truth_E1 = -999.0;
            o_truth_pt2 = o_truth_eta2 = o_truth_phi2 = o_truth_E2 = -999.0;
            o_truth_xj  = -999.0;
            o_truth_flavor1 = o_truth_flavor2 = -999;
            o_truth_parton_pt1 = o_truth_parton_pt2 = -999.0;
            o_dpsi2 = -999.0;
        }

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

        if ( r.reco_pair )
        {
            const int s1 = r.reco_slot[0];
            const int s2 = r.reco_slot[1];
            o_reco_pt1  = match_reco_pT  -> at( s1 );
            o_reco_eta1 = match_reco_eta -> at( s1 );
            o_reco_phi1 = match_reco_phi -> at( s1 );
            o_reco_E1   = match_reco_E   -> at( s1 );
            o_reco_pt2  = match_reco_pT  -> at( s2 );
            o_reco_eta2 = match_reco_eta -> at( s2 );
            o_reco_phi2 = match_reco_phi -> at( s2 );
            o_reco_E2   = match_reco_E   -> at( s2 );
            o_reco_xj   = ( o_reco_pt1 > 0.0f ) ? o_reco_pt2 / o_reco_pt1 : -999.0f;
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
        }
        else
        {
            o_reco_pt1 = o_reco_eta1 = o_reco_phi1 = o_reco_E1 = -999.0;
            o_reco_pt2 = o_reco_eta2 = o_reco_phi2 = o_reco_E2 = -999.0;
            o_reco_xj  = -999.0;
            o_reco_unsub_pt1 = o_reco_unsub_pt2 = -999.0;
            o_reco_dr1 = o_reco_dr2 = -999.0;
        }

        o_n_truth_jets = n_truth_jets;
        o_n_reco_jets  = n_reco_jets;
        o_n_match      = static_cast< int >( match_reco_pT -> size() );

        tout -> Fill();
    }

    std::cout << "\n---- pair-level dijet categories (exclusive) ----" << std::endl;
    std::cout << "Events read: " << n_read << ", written: " << n_events
              << "  (the rest are Skip: neither candidate in acceptance)" << std::endl;
    std::cout << "  Fill  (0): " << n_cat[0] << "   <- the only category that fills the response" << std::endl;
    std::cout << "             of which reco legs pT-swapped: " << n_swapped << std::endl;
    std::cout << "  Miss  (1): " << n_cat[1] << std::endl;
    std::cout << "  Fake  (2): " << n_cat[2] << std::endl;
    std::cout << "             truth pair failed eta: " << n_fake_truth_eta
              << ", pT: " << n_fake_truth_pt
              << ", dphi: " << n_fake_truth_dphi << std::endl;
    std::cout << "             both reco legs are the truth pair: " << n_fake_legs_are_truth << std::endl;
    std::cout << "  FakeMiss (5): " << n_cat[5] << "   <- reco dijet candidate exists and is in"
              << " acceptance, but is not the truth pair: counts as BOTH a fake and a miss"
              << std::endl;
    std::cout << "             wrong leg -- leading: " << n_fm_lead_wrong
              << ", subleading: " << n_fm_sub_wrong
              << ", both: " << n_fm_both_wrong << std::endl;
    std::cout << "             truth leg with no reco match at all -- leading: " << n_fm_truth_lead_unmatched
              << ", subleading: " << n_fm_truth_sub_unmatched << std::endl;
    std::cout << "  events with a truth leg dropped by the slimming (truth_known == 0): "
              << n_truth_leg_unknown << std::endl;
    std::cout << "  UESub (4): " << n_cat[4] << "   <- never assigned; every reco jet on slim"
              << " input is truth-matched" << std::endl;

    std::cout << "  ---- cross-checks ----" << std::endl;
    std::cout << "  categorized: "
              << ( n_cat[0] + n_cat[1] + n_cat[2] + n_cat[3] + n_cat[4] + n_cat[5] )
              << " / " << n_events << std::endl;
    // response denominators, stated the way the unfolding uses them
    std::cout << "  misses (Miss + FakeMiss): " << ( n_cat[1] + n_cat[5] )
              << " | fakes (Fake + FakeMiss): " << ( n_cat[2] + n_cat[5] ) << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " rows to " << outfile << std::endl;

    return 0;
}

#endif
