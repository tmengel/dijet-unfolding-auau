#ifndef _SLIM_MATCHED_JETS_C_
#define _SLIM_MATCHED_JETS_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Slims the tree written by match_standalone.C down to the two jet
// populations a matching / response study actually needs, one row per
// event, with every event-level branch passed through unchanged.
//
// TWO COLLECTIONS, deliberately kept separate:
//
//  (1) match_*  -- one entry per RECO jet that HAS a truth match
//                  (jet_truth_match_idx >= 0), carrying the reco jet AND
//                  its corresponding truth jet side by side. No eta or pT
//                  requirement is applied on either side, so the reco
//                  jet's and truth jet's acceptance flags are carried
//                  along (match_reco_accept_eta / match_truth_accept_eta)
//                  for cutting on downstream.
//
//  (2) acc_truth_*  -- one entry per TRUTH jet inside the eta acceptance
//                  (truth_jet_accept_eta != 0), REGARDLESS of match
//                  status. acc_truth_is_matched says whether it was
//                  matched; acc_truth_match_reco_pT gives the matched
//                  reco pT (-999 when unmatched), so a matching
//                  efficiency and a pT response can both be read off this
//                  one collection with a bare TTree::Draw.
//
// NEITHER COLLECTION CONTAINS THE OTHER, which is the whole point of
// keeping them apart:
//   - a matched reco jet's truth partner can be OUTSIDE the eta
//     acceptance, so it appears in (1) but not in (2);
//   - a truth jet in acceptance can be unmatched, so it appears in (2)
//     but not in (1).
// acc_truth_match_slot cross-references the two: for a matched accepted
// truth jet it is that pair's index in the match_* collection, else -1.
//
// Indices back into the ORIGINAL match_standalone.C collections are kept
// too (match_reco_idx, match_truth_idx, acc_truth_idx) so a row can
// always be traced back to the unslimmed tree. n_truth_jets / n_reco_jets
// record the original multiplicities, so what was dropped is recoverable
// without reopening the input.
//
// Definitions are inherited from match_standalone.C and NOT re-derived
// here: the match is its strict 1-to-1 highest-pT-within-dR < 0.75*0.3
// assignment, and the eta acceptance is its
// AnaUtils::accept_jet_eta( eta, zvrtx, jet_R ) flag. match_dr is
// recomputed only so the pair's dR is available without a join.
//
// Input: the tree written by match_standalone.C (needs truth_jet_*,
// jet_*, truth_jet_accept_eta, truth_jet_reco_match_idx and
// jet_truth_match_idx).
int slim_matched_jets(
    const std::string & infile  = "output.root",
    const std::string & outfile = "slim_matched_jets.root"
)
{
    auto * t = new TChain( "T" );
    // infile can be either a single merged .root file (as written by
    // match_standalone.C) or a .list file listing several of them --
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
    if ( nentries < 1 )
    {
        std::cerr << "Error: No entries in TChain." << std::endl;
        return -1;
    }
    std::cout << "Total entries in TChain: " << nentries << std::endl;

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
    float sumeT_cemc = -999.0, sumeT_ihcal = -999.0, sumeT_ohcal = -999.0;
    float psi2 = -999.0, psi3 = -999.0;
    float truth_zvtx = -999.0, truth_jet_R = -999.0, jet_R = -999.0;
    float truth_jet_maxpT_r04 = -999.0;

    const bool has_event_id    = has_branch( "event_id" );
    const bool has_cent        = has_branch( "cent" );
    const bool has_zvrtx       = has_branch( "zvrtx" );
    const bool has_minbias     = has_branch( "is_minbias" );
    const bool has_mbd_q       = has_branch( "mbd_q" );
    const bool has_sumeT       = has_branch( "sumeT" );
    const bool has_sumeT_parts = has_branch( "sumeT_cemc" );
    const bool has_psi2        = has_branch( "psi2" );
    const bool has_psi3        = has_branch( "psi3" );
    const bool has_truth_zvtx  = has_branch( "truth_zvtx" );
    const bool has_truth_R     = has_branch( "truth_jet_R" );
    const bool has_jet_R       = has_branch( "jet_R" );
    const bool has_maxpT_r04   = has_branch( "truth_jet_maxpT_r04" );

    if ( has_event_id )    enable( "event_id", &event_id );
    if ( has_cent )        enable( "cent", &cent );
    if ( has_zvrtx )       enable( "zvrtx", &zvrtx );
    if ( has_minbias )     enable( "is_minbias", &is_minbias );
    if ( has_mbd_q )       enable( "mbd_q", &mbd_q );
    if ( has_sumeT )       enable( "sumeT", &sumeT );
    if ( has_sumeT_parts )
    {
        enable( "sumeT_cemc", &sumeT_cemc );
        enable( "sumeT_ihcal", &sumeT_ihcal );
        enable( "sumeT_ohcal", &sumeT_ohcal );
    }
    if ( has_psi2 )       enable( "psi2", &psi2 );
    if ( has_psi3 )       enable( "psi3", &psi3 );
    if ( has_truth_zvtx ) enable( "truth_zvtx", &truth_zvtx );
    if ( has_truth_R )    enable( "truth_jet_R", &truth_jet_R );
    if ( has_jet_R )      enable( "jet_R", &jet_R );
    if ( has_maxpT_r04 )  enable( "truth_jet_maxpT_r04", &truth_jet_maxpT_r04 );

    //----------------------------------------------------------------
    // truth / reco jets and the match_standalone.C matching branches
    //----------------------------------------------------------------
    std::vector< float > * truth_jet_pT  = nullptr;
    std::vector< float > * truth_jet_eta = nullptr;
    std::vector< float > * truth_jet_phi = nullptr;
    std::vector< float > * truth_jet_E   = nullptr;
    std::vector< int >   * truth_jet_accept_eta = nullptr;
    std::vector< int >   * truth_jet_reco_match_idx = nullptr;

    std::vector< int >   * truth_jet_flavor    = nullptr;
    std::vector< float > * truth_jet_parton_pT = nullptr;
    std::vector< float > * truth_jet_parton_dr = nullptr;

    std::vector< float > * jet_pT  = nullptr;
    std::vector< float > * jet_eta = nullptr;
    std::vector< float > * jet_phi = nullptr;
    std::vector< float > * jet_E   = nullptr;
    std::vector< float > * jet_unsub_pT = nullptr;
    std::vector< float > * jet_unsub_E  = nullptr;
    std::vector< int >   * jet_accept_eta = nullptr;
    std::vector< int >   * jet_truth_match_idx = nullptr;

    std::vector< float > * jet_cemcfrac  = nullptr;
    std::vector< float > * jet_ihcalfrac = nullptr;
    std::vector< float > * jet_ohcalfrac = nullptr;

    const bool has_all = has_branch( "truth_jet_pT" )
                      && has_branch( "truth_jet_accept_eta" )
                      && has_branch( "truth_jet_reco_match_idx" )
                      && has_branch( "jet_pT" )
                      && has_branch( "jet_truth_match_idx" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the truth/reco jet or "
                  << "matching branches written by match_standalone.C "
                  << "(truth_jet_accept_eta / truth_jet_reco_match_idx / "
                  << "jet_truth_match_idx)." << std::endl;
        return -1;
    }

    const bool has_partons    = has_branch( "truth_jet_flavor" );
    const bool has_unsub      = has_branch( "jet_unsub_pT" );
    const bool has_calofracs  = has_branch( "jet_cemcfrac" );
    const bool has_jet_acc    = has_branch( "jet_accept_eta" );

    enable( "truth_jet_pT", &truth_jet_pT );
    enable( "truth_jet_eta", &truth_jet_eta );
    enable( "truth_jet_phi", &truth_jet_phi );
    enable( "truth_jet_E", &truth_jet_E );
    enable( "truth_jet_accept_eta", &truth_jet_accept_eta );
    enable( "truth_jet_reco_match_idx", &truth_jet_reco_match_idx );
    enable( "jet_pT", &jet_pT );
    enable( "jet_eta", &jet_eta );
    enable( "jet_phi", &jet_phi );
    enable( "jet_E", &jet_E );
    enable( "jet_truth_match_idx", &jet_truth_match_idx );
    if ( has_jet_acc ) enable( "jet_accept_eta", &jet_accept_eta );
    if ( has_partons )
    {
        enable( "truth_jet_flavor", &truth_jet_flavor );
        enable( "truth_jet_parton_pT", &truth_jet_parton_pT );
        enable( "truth_jet_parton_dr", &truth_jet_parton_dr );
    }
    if ( has_unsub )
    {
        enable( "jet_unsub_pT", &jet_unsub_pT );
        enable( "jet_unsub_E", &jet_unsub_E );
    }
    if ( has_calofracs )
    {
        enable( "jet_cemcfrac", &jet_cemcfrac );
        enable( "jet_ihcalfrac", &jet_ihcalfrac );
        enable( "jet_ohcalfrac", &jet_ohcalfrac );
    }

    //----------------------------------------------------------------
    // generic per-node rho values (rho_val_<node>, rho_sigma_<node>) --
    // branch names are dynamic, so wire these up in a loop instead of
    // one variable per node.
    //----------------------------------------------------------------
    std::vector< std::string > rho_names;
    for ( auto * obj : *t -> GetListOfBranches() )
    {
        const std::string bname = obj -> GetName();
        if ( bname.rfind( "rho_val_", 0 ) == 0 )
        {
            rho_names.push_back( bname.substr( std::string( "rho_val_" ).size() ) );
        }
    }
    std::vector< float > rho_vals  ( rho_names.size(), -999.0f );
    std::vector< float > rho_sigmas( rho_names.size(), -999.0f );

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );
    auto * tout = new TTree( "T", "matched reco/truth jet pairs + accepted truth jets" );

    if ( has_event_id )    tout -> Branch( "event_id", &event_id, "event_id/I" );
    if ( has_cent )        tout -> Branch( "cent", &cent, "cent/I" );
    if ( has_zvrtx )       tout -> Branch( "zvrtx", &zvrtx, "zvrtx/F" );
    if ( has_minbias )     tout -> Branch( "is_minbias", &is_minbias, "is_minbias/I" );
    if ( has_mbd_q )       tout -> Branch( "mbd_q", &mbd_q, "mbd_q/F" );
    if ( has_sumeT )       tout -> Branch( "sumeT", &sumeT, "sumeT/F" );
    if ( has_sumeT_parts )
    {
        tout -> Branch( "sumeT_cemc", &sumeT_cemc, "sumeT_cemc/F" );
        tout -> Branch( "sumeT_ihcal", &sumeT_ihcal, "sumeT_ihcal/F" );
        tout -> Branch( "sumeT_ohcal", &sumeT_ohcal, "sumeT_ohcal/F" );
    }
    if ( has_psi2 )       tout -> Branch( "psi2", &psi2, "psi2/F" );
    if ( has_psi3 )       tout -> Branch( "psi3", &psi3, "psi3/F" );
    if ( has_truth_zvtx ) tout -> Branch( "truth_zvtx", &truth_zvtx, "truth_zvtx/F" );
    if ( has_truth_R )    tout -> Branch( "truth_jet_R", &truth_jet_R, "truth_jet_R/F" );
    if ( has_jet_R )      tout -> Branch( "jet_R", &jet_R, "jet_R/F" );
    if ( has_maxpT_r04 )  tout -> Branch( "truth_jet_maxpT_r04", &truth_jet_maxpT_r04, "truth_jet_maxpT_r04/F" );

    for ( size_t i = 0; i < rho_names.size(); ++i )
    {
        const std::string val_name   = "rho_val_"   + rho_names[i];
        const std::string sigma_name = "rho_sigma_" + rho_names[i];
        enable( val_name.c_str(), &rho_vals[i] );
        enable( sigma_name.c_str(), &rho_sigmas[i] );
        tout -> Branch( val_name.c_str(), &rho_vals[i], ( val_name + "/F" ).c_str() );
        tout -> Branch( sigma_name.c_str(), &rho_sigmas[i], ( sigma_name + "/F" ).c_str() );
    }

    //----------------------------------------------------------------
    // original multiplicities, so what was slimmed away stays knowable
    //----------------------------------------------------------------
    int o_n_truth_jets = 0, o_n_reco_jets = 0;
    int o_n_match = 0, o_n_acc_truth = 0;

    tout -> Branch( "n_truth_jets", &o_n_truth_jets, "n_truth_jets/I" );
    tout -> Branch( "n_reco_jets", &o_n_reco_jets, "n_reco_jets/I" );
    tout -> Branch( "n_match", &o_n_match, "n_match/I" );
    tout -> Branch( "n_acc_truth", &o_n_acc_truth, "n_acc_truth/I" );

    //----------------------------------------------------------------
    // (1) matched reco jets and their corresponding truth jets
    //----------------------------------------------------------------
    std::vector< int >   m_reco_idx, m_truth_idx;
    std::vector< float > m_reco_pT, m_reco_eta, m_reco_phi, m_reco_E;
    std::vector< float > m_reco_unsub_pT, m_reco_unsub_E;
    std::vector< int >   m_reco_accept_eta;
    std::vector< float > m_reco_cemcfrac, m_reco_ihcalfrac, m_reco_ohcalfrac;
    std::vector< float > m_truth_pT, m_truth_eta, m_truth_phi, m_truth_E;
    std::vector< int >   m_truth_accept_eta;
    std::vector< int >   m_truth_flavor;
    std::vector< float > m_truth_parton_pT, m_truth_parton_dr;
    std::vector< float > m_dr, m_response;

    // index back into the unslimmed match_standalone.C collections
    tout -> Branch( "match_reco_idx", &m_reco_idx );
    tout -> Branch( "match_truth_idx", &m_truth_idx );
    tout -> Branch( "match_reco_pT", &m_reco_pT );
    tout -> Branch( "match_reco_eta", &m_reco_eta );
    tout -> Branch( "match_reco_phi", &m_reco_phi );
    tout -> Branch( "match_reco_E", &m_reco_E );
    if ( has_unsub )
    {
        tout -> Branch( "match_reco_unsub_pT", &m_reco_unsub_pT );
        tout -> Branch( "match_reco_unsub_E", &m_reco_unsub_E );
    }
    if ( has_jet_acc ) tout -> Branch( "match_reco_accept_eta", &m_reco_accept_eta );
    if ( has_calofracs )
    {
        tout -> Branch( "match_reco_cemcfrac", &m_reco_cemcfrac );
        tout -> Branch( "match_reco_ihcalfrac", &m_reco_ihcalfrac );
        tout -> Branch( "match_reco_ohcalfrac", &m_reco_ohcalfrac );
    }
    tout -> Branch( "match_truth_pT", &m_truth_pT );
    tout -> Branch( "match_truth_eta", &m_truth_eta );
    tout -> Branch( "match_truth_phi", &m_truth_phi );
    tout -> Branch( "match_truth_E", &m_truth_E );
    // the truth partner of a matched reco jet is NOT required to be in the
    // eta acceptance -- carry the flag so it can be cut on downstream
    tout -> Branch( "match_truth_accept_eta", &m_truth_accept_eta );
    if ( has_partons )
    {
        tout -> Branch( "match_truth_flavor", &m_truth_flavor );
        tout -> Branch( "match_truth_parton_pT", &m_truth_parton_pT );
        tout -> Branch( "match_truth_parton_dr", &m_truth_parton_dr );
    }
    tout -> Branch( "match_dr", &m_dr );
    // reco pT / truth pT -- the pT response, the quantity this collection
    // exists to make plottable in one Draw
    tout -> Branch( "match_response", &m_response );

    //----------------------------------------------------------------
    // (2) truth jets in the eta acceptance, which are not matched
    //----------------------------------------------------------------
    std::vector< int >   a_truth_idx, a_is_matched, a_match_slot;
    std::vector< float > a_pT, a_eta, a_phi, a_E;
    std::vector< int >   a_flavor;
    std::vector< float > a_parton_pT, a_parton_dr;
    std::vector< float > a_match_reco_pT;

    tout -> Branch( "acc_truth_idx", &a_truth_idx );
    tout -> Branch( "acc_truth_pT", &a_pT );
    tout -> Branch( "acc_truth_eta", &a_eta );
    tout -> Branch( "acc_truth_phi", &a_phi );
    tout -> Branch( "acc_truth_E", &a_E );
    if ( has_partons )
    {
        tout -> Branch( "acc_truth_flavor", &a_flavor );
        tout -> Branch( "acc_truth_parton_pT", &a_parton_pT );
        tout -> Branch( "acc_truth_parton_dr", &a_parton_dr );
    }
    // 1 when this truth jet has a reco match -- the matching efficiency
    // numerator, with the collection itself as the denominator
    tout -> Branch( "acc_truth_is_matched", &a_is_matched );
    // this pair's index in the match_* collection above, -1 if unmatched
    tout -> Branch( "acc_truth_match_slot", &a_match_slot );
    // matched reco pT, -999 when unmatched, so efficiency AND response are
    // both readable off this one collection without a join
    tout -> Branch( "acc_truth_match_reco_pT", &a_match_reco_pT );

    long n_events = 0;
    long tot_truth = 0, tot_reco = 0, tot_match = 0, tot_acc_truth = 0;
    long tot_acc_truth_matched = 0, tot_match_truth_outside_acc = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );

        const int n_truth = static_cast< int >( truth_jet_pT -> size() );
        const int n_reco  = static_cast< int >( jet_pT -> size() );

        m_reco_idx.clear();       m_truth_idx.clear();
        m_reco_pT.clear();        m_reco_eta.clear();
        m_reco_phi.clear();       m_reco_E.clear();
        m_reco_unsub_pT.clear();  m_reco_unsub_E.clear();
        m_reco_accept_eta.clear();
        m_reco_cemcfrac.clear();  m_reco_ihcalfrac.clear();
        m_reco_ohcalfrac.clear();
        m_truth_pT.clear();       m_truth_eta.clear();
        m_truth_phi.clear();      m_truth_E.clear();
        m_truth_accept_eta.clear();
        m_truth_flavor.clear();
        m_truth_parton_pT.clear(); m_truth_parton_dr.clear();
        m_dr.clear();             m_response.clear();

        a_truth_idx.clear();  a_is_matched.clear();  a_match_slot.clear();
        a_pT.clear();  a_eta.clear();  a_phi.clear();  a_E.clear();
        a_flavor.clear();  a_parton_pT.clear();  a_parton_dr.clear();
        a_match_reco_pT.clear();

        //------------------------------------------------------------
        // (1) every reco jet in that has a truth match, with its match.
        // jet_pT is pT-descending on input (sort_reco_jets_by_pt in
        // match_standalone.C), so this collection stays pT-ordered.
        //------------------------------------------------------------
        // truth index -> slot in the match_* collection, so collection (2)
        // can point at the same pair rather than duplicate it
        std::vector< int > truth_to_slot( n_truth, -1 );

        for ( int ir = 0; ir < n_reco; ++ir )
        {
            const int it = jet_truth_match_idx -> at( ir );
            if ( it < 0 ) continue;             // unmatched reco jet -- dropped
            if ( it >= n_truth ) continue;      // defensive: stale index

            truth_to_slot[it] = static_cast< int >( m_reco_idx.size() );

            m_reco_idx.push_back( ir );
            m_truth_idx.push_back( it );

            m_reco_pT.push_back ( jet_pT  -> at( ir ) );
            m_reco_eta.push_back( jet_eta -> at( ir ) );
            m_reco_phi.push_back( jet_phi -> at( ir ) );
            m_reco_E.push_back  ( jet_E   -> at( ir ) );
            if ( has_unsub )
            {
                m_reco_unsub_pT.push_back( jet_unsub_pT -> at( ir ) );
                m_reco_unsub_E.push_back ( jet_unsub_E  -> at( ir ) );
            }
            if ( has_jet_acc ) m_reco_accept_eta.push_back( jet_accept_eta -> at( ir ) );
            if ( has_calofracs )
            {
                m_reco_cemcfrac.push_back ( jet_cemcfrac  -> at( ir ) );
                m_reco_ihcalfrac.push_back( jet_ihcalfrac -> at( ir ) );
                m_reco_ohcalfrac.push_back( jet_ohcalfrac -> at( ir ) );
            }

            m_truth_pT.push_back ( truth_jet_pT  -> at( it ) );
            m_truth_eta.push_back( truth_jet_eta -> at( it ) );
            m_truth_phi.push_back( truth_jet_phi -> at( it ) );
            m_truth_E.push_back  ( truth_jet_E   -> at( it ) );
            m_truth_accept_eta.push_back( truth_jet_accept_eta -> at( it ) );
            if ( has_partons )
            {
                m_truth_flavor.push_back   ( truth_jet_flavor    -> at( it ) );
                m_truth_parton_pT.push_back( truth_jet_parton_pT -> at( it ) );
                m_truth_parton_dr.push_back( truth_jet_parton_dr -> at( it ) );
            }

            m_dr.push_back( AnaUtils::calc_dr(
                truth_jet_eta -> at( it ), truth_jet_phi -> at( it ),
                jet_eta       -> at( ir ), jet_phi       -> at( ir ) ) 
            );
            const float tpt = truth_jet_pT -> at( it );
            m_response.push_back( ( tpt > 0.0f ) ? jet_pT -> at( ir ) / tpt : -999.0f );

            if ( !truth_jet_accept_eta -> at( it ) ) ++tot_match_truth_outside_acc;
        }

        //------------------------------------------------------------
        // (2) every truth jet inside the eta acceptance, without a given reco match
        // truth_jet_pT is pT-descending on input too.
        //------------------------------------------------------------
        for ( int it = 0; it < n_truth; ++it )
        {
            if ( !truth_jet_accept_eta -> at( it ) ) continue;

            const int mi = truth_jet_reco_match_idx -> at( it );
            const bool matched = ( mi >= 0 && mi < n_reco );
            if ( matched ) continue; // skip truth jets that already have a reco match

            a_truth_idx.push_back( it );
            a_pT.push_back ( truth_jet_pT  -> at( it ) );
            a_eta.push_back( truth_jet_eta -> at( it ) );
            a_phi.push_back( truth_jet_phi -> at( it ) );
            a_E.push_back  ( truth_jet_E   -> at( it ) );
            if ( has_partons )
            {
                a_flavor.push_back   ( truth_jet_flavor    -> at( it ) );
                a_parton_pT.push_back( truth_jet_parton_pT -> at( it ) );
                a_parton_dr.push_back( truth_jet_parton_dr -> at( it ) );
            }

            a_is_matched.push_back( matched ? 1 : 0 );
            a_match_slot.push_back( truth_to_slot[it] );
            a_match_reco_pT.push_back( matched ? jet_pT -> at( mi ) : -999.0f );

            if ( matched ) ++tot_acc_truth_matched;
        }

        o_n_truth_jets = n_truth;
        o_n_reco_jets  = n_reco;
        o_n_match      = static_cast< int >( m_reco_idx.size() );
        o_n_acc_truth  = static_cast< int >( a_truth_idx.size() );

        ++n_events;
        tot_truth     += n_truth;
        tot_reco      += n_reco;
        tot_match     += o_n_match;
        tot_acc_truth += o_n_acc_truth;

        tout -> Fill();
    }

    std::cout << "\n---- slimmed jet collections ----" << std::endl;
    std::cout << "Events: " << n_events << std::endl;
    std::cout << "  input jets            -- truth: " << tot_truth
              << ", reco: " << tot_reco << std::endl;
    std::cout << "  (1) matched reco jets : " << tot_match
              << "  (" << ( tot_reco ? 100.0 * tot_match / tot_reco : 0.0 )
              << "% of reco jets kept)" << std::endl;
    std::cout << "        of which the truth partner is OUTSIDE the eta acceptance: "
              << tot_match_truth_outside_acc << std::endl;
    std::cout << "  (2) accepted truth jets: " << tot_acc_truth
              << "  (" << ( tot_truth ? 100.0 * tot_acc_truth / tot_truth : 0.0 )
              << "% of truth jets kept)" << std::endl;
    std::cout << "        of which matched: " << tot_acc_truth_matched
              << "  -> matching efficiency "
              << ( tot_acc_truth ? 100.0 * tot_acc_truth_matched / tot_acc_truth : 0.0 )
              << "%" << std::endl;

    fout -> cd();
    tout -> Write();
    const long n_written = tout -> GetEntries();
    fout -> Close();

    std::cout << "\nWrote " << n_written << " events to " << outfile << std::endl;

    return 0;
}

#endif
