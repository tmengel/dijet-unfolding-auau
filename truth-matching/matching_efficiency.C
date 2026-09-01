#ifndef _MATCHING_EFFICIENCY_C_
#define _MATCHING_EFFICIENCY_C_

#include <myana/AnaUtils.h>

#include <TChain.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY( libmyana.so )

// Standalone QA of the truth<->reco jet matching done by
// match_standalone.C. Nothing here knows about dijets, categories or
// response -- it only answers "how often does a truth jet get a reco
// match, and how often does a reco jet have a truth partner".
//
// The matcher under test (match_standalone.C, "strict 1-to-1"): truth
// jets are walked in descending pT and each takes the first not-yet-used
// reco jet within dR < 0.75 * 0.3 = 0.225. No pT or eta selection is
// applied on either side, and no reco jet can be claimed twice.
//
// WHAT IS MEASURED
//
//   Efficiency (truth side) = N(truth jets with truth_jet_reco_match_idx
//   >= 0) / N(truth jets), binned in truth pT, |eta| and centrality. This
//   is the quantity that turns into the Miss rate downstream.
//
//   Matched fraction (reco side) = N(reco jets with jet_truth_match_idx
//   >= 0) / N(reco jets), binned in reco pT and centrality. Its
//   complement is the unmatched -- i.e. combinatorial / UE -- reco jet
//   fraction, the quantity that turns into fakes downstream.
//
// Both are reported twice, over two different denominators, because the
// two answer different questions and are easy to confuse:
//
//   "all"      -- every jet in the collection. Includes truth jets that
//                 no reco jet COULD have matched because they are outside
//                 the calorimeter acceptance, so this folds the eta
//                 acceptance into the number.
//   "accepted" -- only jets passing the acceptance the analysis actually
//                 uses (truth_jet_accept_eta / jet_accept_eta, plus
//                 E > 0 on the reco side). This is the matching
//                 efficiency proper, and the one to quote.
//
// The efficiency is also split by truth jet RANK (0 = leading, 1 =
// subleading, >=2 = the rest), since the dijet selection only ever asks
// about the two hardest jets and their efficiencies are not the same.
//
// CONSISTENCY CHECKS on the matcher itself, counted and printed at the
// end -- all four should be exactly zero:
//   - asymmetric pairs: truth_jet_reco_match_idx[it] == ir but
//     jet_truth_match_idx[ir] != it (or vice versa);
//   - double-claimed reco jets: one reco index appearing as the match of
//     more than one truth jet;
//   - out-of-range indices;
//   - matched pairs with dR >= 0.225, recomputed here from eta/phi.
// The dR of every matched pair is histogrammed as well, so the working
// point is visible rather than assumed.
//
// Input: the tree written by match_standalone.C (truth_jet_pT/eta/phi,
// truth_jet_accept_eta, truth_jet_reco_match_idx, jet_pT/eta/phi, jet_E,
// jet_accept_eta, jet_truth_match_idx; cent is optional). infile may be a
// single .root file or a .list of them.
//
// Output: one ROOT file holding the numerator/denominator TH1s, the
// TEfficiency objects built from them, and the 2D pT-vs-centrality
// versions, plus a summary table on stdout.
namespace MatchEff
{
    // Coarse centrality bins, plus an inclusive entry at index 0 that is
    // filled for every event regardless of centrality (and is the only
    // one filled when the input has no cent branch).
    const int   NCENT = 5;
    const float cent_lo[NCENT] = {  0,  0, 10, 30, 50 };
    const float cent_hi[NCENT] = { 90, 10, 30, 50, 90 };

    // std::string rather than the const char* Form() hands back: Form
    // uses a rotating static buffer, and these are used INSIDE the
    // argument list of another Form, which can clobber it.
    inline std::string cent_tag( const int i )
    {
        return ( i == 0 ) ? "incl"
                          : std::string( Form( "%d_%d", (int) cent_lo[i], (int) cent_hi[i] ) );
    }
    inline std::string cent_label( const int i )
    {
        return ( i == 0 ) ? "0-90%"
                          : std::string( Form( "%d-%d%%", (int) cent_lo[i], (int) cent_hi[i] ) );
    }

    // The dR working point of match_standalone.C, repeated here so the
    // matched pairs can be checked against it independently.
    const float max_dr = 0.75f * 0.3f;

    // A numerator/denominator pair plus the TEfficiency built from them.
    // Kept as two TH1s rather than only a TEfficiency so the raw counts
    // stay inspectable in the output file.
    struct Eff
    {
        TH1F * den = nullptr;
        TH1F * num = nullptr;

        void make( const std::string & name, const std::string & title,
                   const int nbins, const double * edges )
        {
            den = new TH1F( ( name + "_den" ).c_str(), title.c_str(), nbins, edges );
            num = new TH1F( ( name + "_num" ).c_str(), title.c_str(), nbins, edges );
        }
        void make( const std::string & name, const std::string & title,
                   const int nbins, const double lo, const double hi )
        {
            den = new TH1F( ( name + "_den" ).c_str(), title.c_str(), nbins, lo, hi );
            num = new TH1F( ( name + "_num" ).c_str(), title.c_str(), nbins, lo, hi );
        }

        void fill( const double x, const bool matched )
        {
            den -> Fill( x );
            if ( matched ) num -> Fill( x );
        }

        // Ratio of the integrals, with the binomial error. Used for the
        // stdout summary; the per-bin version is the TEfficiency.
        double integral_eff( double & err ) const
        {
            const double d = den -> Integral( 0, den -> GetNbinsX() + 1 );
            const double n = num -> Integral( 0, num -> GetNbinsX() + 1 );
            if ( d <= 0.0 ) { err = 0.0; return 0.0; }
            const double e = n / d;
            err = TMath::Sqrt( TMath::Max( e * ( 1.0 - e ), 0.0 ) / d );
            return e;
        }
        double n_den() const { return den -> Integral( 0, den -> GetNbinsX() + 1 ); }

        // TEfficiency needs num <= den bin by bin; that holds by
        // construction here, but guard anyway so a pathological input
        // cannot abort the job.
        void write( const std::string & name ) const
        {
            den -> Write();
            num -> Write();
            if ( TEfficiency::CheckConsistency( *num, *den ) )
            {
                auto * te = new TEfficiency( *num, *den );
                te -> SetName( ( name + "_eff" ).c_str() );
                te -> SetTitle( den -> GetTitle() );
                te -> Write();
            }
            else
            {
                std::cerr << "Warning: " << name
                          << " failed TEfficiency::CheckConsistency, "
                          << "writing the TH1s only." << std::endl;
            }
        }
    };
}

int matching_efficiency(
    const std::string & infile  = "output.root",
    const std::string & outfile = "matching_efficiency.root"
)
{
    using namespace MatchEff;

    auto * t = new TChain( "T" );
    // infile can be either a single merged .root file (as written by
    // match_standalone.C) or a .list file naming several of them --
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

    const int nentries = t -> GetEntries();
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

    const bool has_all = has_branch( "truth_jet_pT" )
                      && has_branch( "truth_jet_eta" )
                      && has_branch( "truth_jet_phi" )
                      && has_branch( "truth_jet_accept_eta" )
                      && has_branch( "truth_jet_reco_match_idx" )
                      && has_branch( "jet_pT" )
                      && has_branch( "jet_eta" )
                      && has_branch( "jet_phi" )
                      && has_branch( "jet_E" )
                      && has_branch( "jet_accept_eta" )
                      && has_branch( "jet_truth_match_idx" );
    if ( !has_all )
    {
        std::cerr << "Error: input tree is missing the jet or matching branches "
                  << "written by match_standalone.C (truth_jet_reco_match_idx / "
                  << "jet_truth_match_idx)." << std::endl;
        return -1;
    }

    int cent = -1;
    const bool has_cent = has_branch( "cent" );
    if ( has_cent ) enable( "cent", &cent );

    std::vector< float > * truth_jet_pT  = nullptr;
    std::vector< float > * truth_jet_eta = nullptr;
    std::vector< float > * truth_jet_phi = nullptr;
    std::vector< int >   * truth_jet_accept_eta = nullptr;
    std::vector< int >   * truth_jet_reco_match_idx = nullptr;

    std::vector< float > * jet_pT  = nullptr;
    std::vector< float > * jet_eta = nullptr;
    std::vector< float > * jet_phi = nullptr;
    std::vector< float > * jet_E   = nullptr;
    std::vector< int >   * jet_accept_eta = nullptr;
    std::vector< int >   * jet_truth_match_idx = nullptr;

    enable( "truth_jet_pT", &truth_jet_pT );
    enable( "truth_jet_eta", &truth_jet_eta );
    enable( "truth_jet_phi", &truth_jet_phi );
    enable( "truth_jet_accept_eta", &truth_jet_accept_eta );
    enable( "truth_jet_reco_match_idx", &truth_jet_reco_match_idx );
    enable( "jet_pT", &jet_pT );
    enable( "jet_eta", &jet_eta );
    enable( "jet_phi", &jet_phi );
    enable( "jet_E", &jet_E );
    enable( "jet_accept_eta", &jet_accept_eta );
    enable( "jet_truth_match_idx", &jet_truth_match_idx );

    auto * fout = new TFile( outfile.c_str(), "RECREATE" );

    //----------------------------------------------------------------
    // binning
    //----------------------------------------------------------------
    // Fine at low pT, where the efficiency actually turns on, and coarse
    // above it where the statistics run out.
    // Edges at 7, 8, 14 and 20 are the analysis thresholds (truth
    // subleading / reco subleading / truth leading / reco leading), so an
    // efficiency integrated above any of them is an exact bin sum rather
    // than a partial-bin estimate.
    const double pt_edges[] = {
        0, 2, 4, 6, 7, 8, 10, 12, 14, 16, 18, 20,
        25, 30, 35, 40, 50, 60, 80, 100
    };
    const int n_pt = sizeof( pt_edges ) / sizeof( double ) - 1;

    const std::string pt_ax_t   = ";truth jet p_{T} [GeV];matching efficiency";
    const std::string pt_ax_r   = ";reco jet p_{T} [GeV];matched fraction";
    const std::string eta_ax    = ";truth jet #eta;matching efficiency";
    const std::string cent_ax   = ";centrality [%];matching efficiency";

    // [cent][0 = all jets, 1 = jets in the analysis acceptance]
    Eff eff_truth_pt[NCENT][2];
    Eff eff_reco_pt [NCENT][2];
    // truth efficiency by rank, accepted denominator only
    Eff eff_truth_pt_rank[NCENT][3];

    for ( int ic = 0; ic < NCENT; ++ic )
    {
        for ( int ia = 0; ia < 2; ++ia )
        {
            const char * acc = ( ia == 0 ) ? "all" : "acc";
            const char * accl = ( ia == 0 ) ? "all truth jets" : "truth jets in acceptance";
            const char * accr = ( ia == 0 ) ? "all reco jets"  : "reco jets in acceptance";

            eff_truth_pt[ic][ia].make(
                Form( "h_eff_truth_pt_%s_cent%s", acc, cent_tag( ic ).c_str() ),
                Form( "truth->reco matching efficiency, %s, %s%s",
                      accl, cent_label( ic ).c_str(), pt_ax_t.c_str() ),
                n_pt, pt_edges );

            eff_reco_pt[ic][ia].make(
                Form( "h_frac_reco_pt_%s_cent%s", acc, cent_tag( ic ).c_str() ),
                Form( "reco jets with a truth match, %s, %s%s",
                      accr, cent_label( ic ).c_str(), pt_ax_r.c_str() ),
                n_pt, pt_edges );
        }

        for ( int ir = 0; ir < 3; ++ir )
        {
            const char * rl = ( ir == 0 ) ? "leading" : ( ir == 1 ) ? "subleading" : "rank #geq 3";
            eff_truth_pt_rank[ic][ir].make(
                Form( "h_eff_truth_pt_acc_rank%d_cent%s", ir, cent_tag( ic ).c_str() ),
                Form( "truth->reco matching efficiency, %s truth jet, %s%s",
                      rl, cent_label( ic ).c_str(), pt_ax_t.c_str() ),
                n_pt, pt_edges );
        }
    }

    // vs eta and vs centrality, for accepted truth jets above a pT floor
    // (below it the efficiency is dominated by the pT turn-on and the eta
    // and centrality dependence is unreadable).
    const float eta_pt_min = 10.0f;
    Eff eff_truth_eta;
    eff_truth_eta.make( "h_eff_truth_eta_acc",
        Form( "truth->reco matching efficiency, truth jets in acceptance, p_{T} > %.0f GeV%s",
              eta_pt_min, eta_ax.c_str() ),
        30, -1.5, 1.5 );

    Eff eff_truth_cent;
    eff_truth_cent.make( "h_eff_truth_cent_acc",
        Form( "truth->reco matching efficiency, truth jets in acceptance, p_{T} > %.0f GeV%s",
              eta_pt_min, cent_ax.c_str() ),
        18, 0.0, 90.0 );

    // 2D versions, so the pT and centrality dependence can be sliced
    // arbitrarily afterwards rather than only in the coarse bins above.
    auto * h2_truth_den = new TH2F( "h2_eff_truth_ptcent_acc_den",
        "truth jets in acceptance;truth jet p_{T} [GeV];centrality [%]",
        n_pt, pt_edges, 18, 0.0, 90.0 );
    auto * h2_truth_num = new TH2F( "h2_eff_truth_ptcent_acc_num",
        "matched truth jets in acceptance;truth jet p_{T} [GeV];centrality [%]",
        n_pt, pt_edges, 18, 0.0, 90.0 );
    auto * h2_reco_den = new TH2F( "h2_frac_reco_ptcent_acc_den",
        "reco jets in acceptance;reco jet p_{T} [GeV];centrality [%]",
        n_pt, pt_edges, 18, 0.0, 90.0 );
    auto * h2_reco_num = new TH2F( "h2_frac_reco_ptcent_acc_num",
        "matched reco jets in acceptance;reco jet p_{T} [GeV];centrality [%]",
        n_pt, pt_edges, 18, 0.0, 90.0 );
    for ( auto * h : { h2_truth_den, h2_truth_num, h2_reco_den, h2_reco_num } ) h -> Sumw2();

    // matcher diagnostics
    auto * h_dr = new TH1F( "h_match_dr",
        "dR of matched truth-reco pairs;#DeltaR(truth, reco);pairs", 100, 0.0, 0.5 );
    auto * h_dpt = new TH1F( "h_match_dpt",
        "p_{T} response of matched pairs;(reco p_{T} - truth p_{T}) / truth p_{T};pairs",
        120, -1.2, 1.2 );
    h_dr  -> Sumw2();
    h_dpt -> Sumw2();

    //----------------------------------------------------------------
    // event loop
    //----------------------------------------------------------------
    long n_events = 0;
    long n_truth_jets = 0, n_reco_jets = 0;
    long n_bad_asym = 0, n_bad_double = 0, n_bad_range = 0, n_bad_dr = 0;

    for ( int i = 0; i < nentries; ++i )
    {
        t -> GetEntry( i );
        ++n_events;

        const int n_truth = static_cast< int >( truth_jet_pT -> size() );
        const int n_reco  = static_cast< int >( jet_pT -> size() );

        // Which coarse centrality bins this event belongs to. Index 0
        // (inclusive) is always on; without a cent branch it is the only
        // one, so the macro still works on pp-like input.
        bool in_cent[NCENT] = { true, false, false, false, false };
        if ( has_cent )
        {
            for ( int ic = 1; ic < NCENT; ++ic )
            {
                in_cent[ic] = ( cent >= cent_lo[ic] && cent < cent_hi[ic] );
            }
            // keep the inclusive bin to the same 0-90% range the
            // sub-bins cover, so the numbers are comparable
            in_cent[0] = ( cent >= cent_lo[0] && cent < cent_hi[0] );
        }
        const double cent_val = has_cent ? cent : -1.0;

        //--- truth side -------------------------------------------------
        std::vector< int > reco_claimed( n_reco, 0 );

        for ( int it = 0; it < n_truth; ++it )
        {
            ++n_truth_jets;

            const int  mi      = truth_jet_reco_match_idx -> at( it );
            const bool matched = ( mi >= 0 );
            const bool in_acc  = ( truth_jet_accept_eta -> at( it ) != 0 );
            const double pt    = truth_jet_pT -> at( it );

            //--- matcher consistency ---
            if ( matched )
            {
                if ( mi >= n_reco )
                {
                    ++n_bad_range;
                }
                else
                {
                    ++reco_claimed[mi];
                    if ( reco_claimed[mi] > 1 ) ++n_bad_double;
                    if ( jet_truth_match_idx -> at( mi ) != it ) ++n_bad_asym;

                    const float dr = AnaUtils::calc_dr(
                        truth_jet_eta -> at( it ), truth_jet_phi -> at( it ),
                        jet_eta -> at( mi ), jet_phi -> at( mi ) );
                    h_dr -> Fill( dr );
                    if ( dr >= max_dr ) ++n_bad_dr;
                    if ( pt > 0.0 ) h_dpt -> Fill( ( jet_pT -> at( mi ) - pt ) / pt );
                }
            }

            //--- efficiency ---
            const int rank = ( it < 2 ) ? it : 2;
            for ( int ic = 0; ic < NCENT; ++ic )
            {
                if ( !in_cent[ic] ) continue;
                eff_truth_pt[ic][0].fill( pt, matched );
                if ( in_acc )
                {
                    eff_truth_pt[ic][1].fill( pt, matched );
                    eff_truth_pt_rank[ic][rank].fill( pt, matched );
                }
            }

            if ( in_acc )
            {
                h2_truth_den -> Fill( pt, cent_val );
                if ( matched ) h2_truth_num -> Fill( pt, cent_val );

                if ( pt > eta_pt_min )
                {
                    eff_truth_eta.fill( truth_jet_eta -> at( it ), matched );
                    if ( has_cent ) eff_truth_cent.fill( cent_val, matched );
                }
            }
        }

        //--- reco side --------------------------------------------------
        for ( int ir = 0; ir < n_reco; ++ir )
        {
            ++n_reco_jets;

            const int  ti      = jet_truth_match_idx -> at( ir );
            const bool matched = ( ti >= 0 );
            // the acceptance the analysis applies to reco jets
            const bool in_acc  = ( jet_accept_eta -> at( ir ) != 0 )
                              && ( jet_E -> at( ir ) > 0.0f );
            const double pt    = jet_pT -> at( ir );

            if ( matched )
            {
                if ( ti >= n_truth ) ++n_bad_range;
                else if ( truth_jet_reco_match_idx -> at( ti ) != ir ) ++n_bad_asym;
            }

            for ( int ic = 0; ic < NCENT; ++ic )
            {
                if ( !in_cent[ic] ) continue;
                eff_reco_pt[ic][0].fill( pt, matched );
                if ( in_acc ) eff_reco_pt[ic][1].fill( pt, matched );
            }

            if ( in_acc )
            {
                h2_reco_den -> Fill( pt, cent_val );
                if ( matched ) h2_reco_num -> Fill( pt, cent_val );
            }
        }
    }

    //----------------------------------------------------------------
    // summary
    //----------------------------------------------------------------
    std::cout << "\n---- match_standalone.C matching efficiency ----" << std::endl;
    std::cout << "Events: " << n_events
              << " | truth jets: " << n_truth_jets
              << " | reco jets: " << n_reco_jets << std::endl;
    if ( !has_cent )
    {
        std::cout << "NOTE: no cent branch in the input -- only the inclusive "
                  << "centrality bin is filled." << std::endl;
    }

    auto print_row = [ & ]( const char * what, const Eff & e )
    {
        double err = 0.0;
        const double v = e.integral_eff( err );
        printf( "  %-26s %8.4f +/- %.4f   (denominator %.0f)\n",
                what, v, err, e.n_den() );
    };

    std::cout << "\nIntegrated over all pT:" << std::endl;
    for ( int ic = 0; ic < NCENT; ++ic )
    {
        if ( !has_cent && ic > 0 ) break;
        std::cout << " cent " << cent_label( ic ) << ":" << std::endl;
        print_row( "truth eff (all)",       eff_truth_pt[ic][0] );
        print_row( "truth eff (in acc)",    eff_truth_pt[ic][1] );
        print_row( "reco matched (all)",    eff_reco_pt[ic][0] );
        print_row( "reco matched (in acc)", eff_reco_pt[ic][1] );
    }

    // The pT-differential table is what actually matters -- an efficiency
    // integrated over all pT is dominated by the soft jets that make up
    // the bulk of the collection.
    std::cout << "\nAccepted truth jets, inclusive centrality, vs truth pT:" << std::endl;
    std::cout << "   pT range [GeV]      efficiency        N(truth)" << std::endl;
    {
        const Eff & e = eff_truth_pt[0][1];
        for ( int b = 1; b <= e.den -> GetNbinsX(); ++b )
        {
            const double d = e.den -> GetBinContent( b );
            const double n = e.num -> GetBinContent( b );
            const double v = ( d > 0.0 ) ? n / d : 0.0;
            const double s = ( d > 0.0 ) ? TMath::Sqrt( TMath::Max( v * ( 1.0 - v ), 0.0 ) / d ) : 0.0;
            printf( "   %5.0f - %5.0f     %7.4f +/- %.4f     %9.0f\n",
                    e.den -> GetBinLowEdge( b ),
                    e.den -> GetBinLowEdge( b ) + e.den -> GetBinWidth( b ),
                    v, s, d );
        }
    }

    std::cout << "\nAccepted reco jets, inclusive centrality, vs reco pT:" << std::endl;
    std::cout << "   pT range [GeV]    matched frac        N(reco)" << std::endl;
    {
        const Eff & e = eff_reco_pt[0][1];
        for ( int b = 1; b <= e.den -> GetNbinsX(); ++b )
        {
            const double d = e.den -> GetBinContent( b );
            const double n = e.num -> GetBinContent( b );
            const double v = ( d > 0.0 ) ? n / d : 0.0;
            const double s = ( d > 0.0 ) ? TMath::Sqrt( TMath::Max( v * ( 1.0 - v ), 0.0 ) / d ) : 0.0;
            printf( "   %5.0f - %5.0f     %7.4f +/- %.4f     %9.0f\n",
                    e.den -> GetBinLowEdge( b ),
                    e.den -> GetBinLowEdge( b ) + e.den -> GetBinWidth( b ),
                    v, s, d );
        }
    }

    std::cout << "\nMatcher consistency (all four should be 0):" << std::endl;
    std::cout << "  asymmetric truth<->reco pairs : " << n_bad_asym   << std::endl;
    std::cout << "  double-claimed reco jets      : " << n_bad_double << std::endl;
    std::cout << "  out-of-range match indices    : " << n_bad_range  << std::endl;
    std::cout << "  matched pairs with dR >= "
              << Form( "%.3f", max_dr ) << "  : " << n_bad_dr << std::endl;
    const int last_dr_bin = h_dr -> FindLastBinAbove( 0.0 );
    std::cout << "  matched-pair dR: mean " << Form( "%.4f", h_dr -> GetMean() )
              << ", largest occupied bin up to "
              << ( ( last_dr_bin > 0 )
                   ? Form( "%.3f", h_dr -> GetXaxis() -> GetBinUpEdge( last_dr_bin ) )
                   : "(no matched pairs)" )
              << std::endl;

    //----------------------------------------------------------------
    // write
    //----------------------------------------------------------------
    fout -> cd();
    for ( int ic = 0; ic < NCENT; ++ic )
    {
        for ( int ia = 0; ia < 2; ++ia )
        {
            const char * acc = ( ia == 0 ) ? "all" : "acc";
            eff_truth_pt[ic][ia].write( Form( "h_eff_truth_pt_%s_cent%s", acc, cent_tag( ic ).c_str() ) );
            eff_reco_pt [ic][ia].write( Form( "h_frac_reco_pt_%s_cent%s", acc, cent_tag( ic ).c_str() ) );
        }
        for ( int ir = 0; ir < 3; ++ir )
        {
            eff_truth_pt_rank[ic][ir].write(
                Form( "h_eff_truth_pt_acc_rank%d_cent%s", ir, cent_tag( ic ).c_str() ) );
        }
    }
    eff_truth_eta.write( "h_eff_truth_eta_acc" );
    eff_truth_cent.write( "h_eff_truth_cent_acc" );

    h2_truth_den -> Write();
    h2_truth_num -> Write();
    h2_reco_den  -> Write();
    h2_reco_num  -> Write();
    h_dr  -> Write();
    h_dpt -> Write();

    fout -> Close();
    std::cout << "\nWrote histograms to " << outfile << std::endl;

    return 0;
}

#endif
