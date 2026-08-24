#ifndef PRIOR_REWEIGHT_QA_H
#define PRIOR_REWEIGHT_QA_H

// Shared construction, application and QA of the prior (truth) reweighting.
//
// createResponse_noempty_AA.cxx and testPriorReweight_AA.C both go through this
// header so the QA plots exercise the same index convention as the event loop.
// If the two ever diverge the QA is worthless, which is the failure mode this
// header exists to remove.

#include <iostream>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"

#include "dlUtility.h"
#include "histo_opps.h"

namespace prior_qa
{

  // ---------------------------------------------------------------------------
  // Index convention
  //
  // The flat pt1-pt2 histograms are TH1Ds with nbins*nbins unit-width bins
  // spanning [0, nbins*nbins) and are filled with Fill(k), k = i + nbins*j.
  // A value of k therefore lands in ROOT bin k+1.
  //
  // Every truth pair is filled in BOTH orderings, (i,j) and (j,i), so the flat
  // truth is symmetric and either ordering of the arguments returns the same
  // weight.  Reading GetBinContent(k) instead of GetBinContent(k+1) returns the
  // weight of the neighbouring truth bin (and, for j = 0, of a completely
  // unrelated (i-1, nbins-1) corner), which is exactly the off-by-one this
  // function is here to prevent.
  // ---------------------------------------------------------------------------
  inline int weight_bin(const int pt1_bin, const int pt2_bin, const int nbins)
  {
    return pt1_bin*nbins + pt2_bin + 1;
  }

  // ---------------------------------------------------------------------------
  // Space in which the prior fraction is applied.
  //
  // kBlendPerBin blends truth and unfolded data independently in each of the
  // ~100 populated (pt1, pt2) truth bins. It hits the target exactly, but it
  // also chases every per-bin statistical fluctuation of the unfolded data
  // (weights up to 23x in a 666-count corner bin at fraction 1).
  //
  // kBlendXJ blends the two xJ SHAPES instead and maps the resulting ratio back
  // onto the (pt1, pt2) bins. That is ~15 numbers rather than ~100, the pt1
  // spectrum is left alone, and the prior moves only in the observable the
  // systematic is quoted on.
  // ---------------------------------------------------------------------------
  enum blend_space { kBlendPerBin = 0, kBlendXJ = 1 };

  inline const char *blend_name(const int space)
  {
    return (space == kBlendXJ) ? "xj" : "perbin";
  }

  struct weights
  {
    TH1D *truth  = nullptr;  // normalized MC prior,   flat pt1pt2
    TH1D *unfold = nullptr;  // normalized unfold data, flat pt1pt2
    TH1D *target = nullptr;  // normalized truth + fraction*(unfold - truth)
    TH1D *weight = nullptr;  // target/truth, flat index k in ROOT bin k+1
    double fraction = 0.0;
    int    space = kBlendPerBin;
    int    n_zero_weight = 0;   // truth bins killed outright by the reweighting
    double zero_weight_truth_fraction = 0.0;
    double weight_min = 0.0;
    double weight_max = 0.0;
  };

  // Normalize a clone to unit integral; returns nullptr if the input is empty.
  inline TH1D *normalized_clone(const TH1D *h, const std::string &name)
  {
    if (!h) return nullptr;
    TH1D *c = (TH1D*) h->Clone(name.c_str());
    c->SetDirectory(nullptr);
    const double integral = c->Integral();
    if (integral != 0) { c->Scale(1.0/integral); }
    else { std::cerr << "prior_qa: warning, " << h->GetName() << " has zero integral." << std::endl; }
    return c;
  }

  // Range of the applied weights, and how much MC truth a zero weight deletes
  // from the response outright.
  inline void update_stats(weights &w, const int nbins)
  {
    w.n_zero_weight = 0;
    w.zero_weight_truth_fraction = 0.0;
    w.weight_min = w.weight_max = 0.0;
    bool first = true;
    for (int ibin = 1; ibin <= nbins*nbins; ibin++)
    {
      const double b = w.truth->GetBinContent(ibin);
      if (!(b > 0)) continue;
      const double r = w.weight->GetBinContent(ibin);
      if (!(r > 0))
      {
        w.n_zero_weight++;
        w.zero_weight_truth_fraction += b;
      }
      if (first) { w.weight_min = w.weight_max = r; first = false; }
      else
      {
        if (r < w.weight_min) w.weight_min = r;
        if (r > w.weight_max) w.weight_max = r;
      }
    }
  }

  // Build the prior target and the per-truth-bin weights.
  //   fraction = 0 -> target is the MC prior, every weight is 1
  //   fraction = 1 -> target is the unfolded data
  inline weights build(const TH1D *truth_in, const TH1D *unfold_in,
                       const double fraction, const int nbins,
                       const std::string &tag = "")
  {
    weights w;
    w.fraction = fraction;

    const int nflat = nbins*nbins;

    w.truth  = normalized_clone(truth_in,  "h_prior_truth_norm" + tag);
    w.unfold = normalized_clone(unfold_in, "h_prior_unfold_norm" + tag);
    if (!w.truth || !w.unfold) { return w; }

    w.target = (TH1D*) w.truth->Clone(("h_prior_target" + tag).c_str());
    w.target->SetDirectory(nullptr);
    w.target->Reset();
    w.weight = (TH1D*) w.truth->Clone(("h_prior_weight" + tag).c_str());
    w.weight->SetDirectory(nullptr);
    w.weight->Reset();

    for (int ibin = 1; ibin <= nflat; ibin++)
    {
      const double y1 = w.unfold->GetBinContent(ibin);
      const double y2 = w.truth ->GetBinContent(ibin);
      w.target->SetBinContent(ibin, y2 + fraction*(y1 - y2));
    }

    const double target_integral = w.target->Integral();
    if (target_integral != 0) { w.target->Scale(1.0/target_integral); }

    for (int ibin = 1; ibin <= nflat; ibin++)
    {
      const double v = w.target->GetBinContent(ibin);
      const double b = w.truth ->GetBinContent(ibin);
      // b <= 0 means the MC prior never populated this truth bin, so there is
      // nothing to reweight and the weight is irrelevant; 1 keeps it neutral.
      double rat = (b > 0) ? v/b : 1.0;
      if (rat < 0) { rat = 0.0; }
      w.weight->SetBinContent(ibin, rat);
    }

    update_stats(w, nbins);
    return w;
  }

  // Project a flat pt1pt2 distribution onto xJ over one leading/subleading window.
  // Shape only: normalized to unit integral, NOT divided by bin width, so a
  // ratio of two of these is directly a per-bin reweighting factor.
  inline TH1D *project_to_xj(const TH1D *flat, const int nbins,
                             const float *ipt_bins, const float *ixj_bins,
                             const int lead_lo, const int lead_hi,
                             const int sub_lo, const int sub_hi,
                             const std::string &name)
  {
    if (!flat) return nullptr;
    TH2D *h2 = new TH2D((name + "_pt1pt2").c_str(), "", nbins, ipt_bins, nbins, ipt_bins);
    h2->SetDirectory(nullptr);
    histo_opps::make_sym_pt1pt2((TH1D*) flat, h2, nbins);

    TH1D *hxj = new TH1D(name.c_str(), "", nbins, ixj_bins);
    hxj->SetDirectory(nullptr);
    histo_opps::project_xj(h2, hxj, nbins, lead_lo, lead_hi, sub_lo, sub_hi);
    delete h2;

    const double integral = hxj->Integral();
    if (integral != 0) { hxj->Scale(1.0/integral); }
    return hxj;
  }

  // Apply the weights to a flat truth histogram exactly as the event loop does:
  // a pair whose truth pt bins are (i,j) is scaled by weight_bin(i,j) and fills
  // both flat indices i + nbins*j and j + nbins*i.
  //
  // bin_offset is a diagnostic knob for testPriorReweight_AA.C: 0 reproduces the
  // production lookup, -1 reproduces the historical off-by-one.
  inline TH1D *apply(const TH1D *truth_in, const TH1D *weight, const int nbins,
                     const std::string &name, const int bin_offset = 0)
  {
    if (!truth_in || !weight) return nullptr;
    TH1D *out = (TH1D*) truth_in->Clone(name.c_str());
    out->SetDirectory(nullptr);
    out->Reset();

    const int nflat = nbins*nbins;
    for (int a = 0; a < nbins; a++)
    {
      for (int b = 0; b < nbins; b++)
      {
        const int flat = a + nbins*b;
        // The event loop indexes with (pt1, pt2) as read off the tree; the
        // weight is symmetric under the correct lookup, so ordering the pair
        // leading/subleading here is exact for bin_offset = 0 and is only a
        // convention for the broken lookup.
        const int i = std::max(a, b);
        const int j = std::min(a, b);
        const int wbin = weight_bin(i, j, nbins) + bin_offset;
        const double scale = (wbin >= 1 && wbin <= nflat) ? weight->GetBinContent(wbin) : 1.0;
        out->SetBinContent(flat + 1, truth_in->GetBinContent(flat + 1)*scale);
        out->SetBinError  (flat + 1, truth_in->GetBinError  (flat + 1)*scale);
      }
    }

    const double integral = out->Integral();
    if (integral != 0) { out->Scale(1.0/integral); }
    return out;
  }

  // One leading-pT window the xJ ratio is built in.
  struct lead_group { int lo; int hi; };

  // Blend the prior in xJ instead of per (pt1, pt2) bin, then map the xJ ratio
  // back onto the truth bins so it can still be applied as a per-pair weight.
  //
  // The target is unambiguous: projection onto xJ is linear, so projecting the
  // per-bin blend gives exactly the blend of the two xJ projections. What
  // changes is the WEIGHT -- it becomes a function of the leading-pT group and
  // the pT-bin separation d only, not of the individual bin's unfolded content.
  //
  // The ratio is built independently in each leading-pT group so the reweighted
  // prior reproduces the xJ target of every measurement range. A single group
  // spanning the whole measurement reproduces only the inclusive xJ shape, which
  // for the Run-24 binning pushes the 43-62 GeV range the wrong way by ~20%.
  //
  // Leading bins below the first group or above the last are in no projection
  // (the measurement does not cover them) and inherit the nearest group's ratio;
  // likewise pairs whose subleading bin is below sub_lo.
  inline weights build_xj(const TH1D *truth_in, const TH1D *unfold_in,
                          const double fraction, const int nbins,
                          const float *ipt_bins, const float *ixj_bins,
                          const std::vector<lead_group> &groups,
                          const int sub_lo, const int sub_hi,
                          const std::string &tag = "",
                          const int n_refine = 20)
  {
    weights w = build(truth_in, unfold_in, fraction, nbins, tag);
    if (!w.weight || groups.empty()) { return w; }
    w.space = kBlendXJ;

    const int ngroups = (int) groups.size();

    auto group_of = [&](const int lead_bin)
    {
      for (int g = 0; g < ngroups; g++)
      {
        if (lead_bin >= groups[g].lo && lead_bin < groups[g].hi) { return g; }
      }
      return (lead_bin < groups[0].lo) ? 0 : ngroups - 1;
    };

    std::vector<TH1D*> xj_truth(ngroups, nullptr);
    std::vector<TH1D*> xj_target(ngroups, nullptr);
    for (int g = 0; g < ngroups; g++)
    {
      xj_truth[g]  = project_to_xj(w.truth,  nbins, ipt_bins, ixj_bins,
                                   groups[g].lo, groups[g].hi, sub_lo, sub_hi,
                                   Form("h_prior_xjblend_truth%s_g%d", tag.c_str(), g));
      xj_target[g] = project_to_xj(w.target, nbins, ipt_bins, ixj_bins,
                                   groups[g].lo, groups[g].hi, sub_lo, sub_hi,
                                   Form("h_prior_xjblend_target%s_g%d", tag.c_str(), g));
      if (!xj_truth[g] || !xj_target[g])
      {
        std::cerr << "prior_qa: xJ projection failed, keeping the per-bin blend." << std::endl;
        for (auto *h : xj_truth)  { delete h; }
        for (auto *h : xj_target) { delete h; }
        return build(truth_in, unfold_in, fraction, nbins, tag);
      }
    }

    // histo_opps::project_xj splits a pair whose pT bins differ by d between xJ
    // bins nbins-d and nbins-d+1: the pT and xJ binnings share the same
    // geometric ratio alpha, so a pair's xJ = alpha^-d lands exactly on an xJ
    // bin edge. d = 0 (the balanced diagonal) goes entirely into bin nbins.
    // Averaging the two bins a given d feeds is the matching inverse.
    auto combine = [&](const int d, const std::function<double(int)> &value)
    {
      if (d == 0) { return value(nbins); }
      const double v_lo = value(nbins - d);
      const double v_hi = value(nbins - d + 1);
      if      (v_lo >= 0 && v_hi >= 0) { return 0.5*(v_lo + v_hi); }
      else if (v_lo >= 0)              { return v_lo; }
      else if (v_hi >= 0)              { return v_hi; }
      return -1.0;
    };

    std::vector<std::vector<double>> r(ngroups, std::vector<double>(nbins, 1.0));
    for (int g = 0; g < ngroups; g++)
    {
      // < 0 marks "no MC truth in this xJ bin", i.e. nothing constrains it.
      auto ratio_at = [&](const int b) -> double
      {
        if (b < 1 || b > nbins) return -1.0;
        const double t = xj_truth[g]->GetBinContent(b);
        if (!(t > 0)) return -1.0;
        const double v = xj_target[g]->GetBinContent(b);
        return (v > 0) ? v/t : 0.0;
      };
      for (int d = 0; d < nbins; d++)
      {
        const double val = combine(d, ratio_at);
        r[g][d] = (val >= 0) ? val : 1.0;
      }
    }

    // Back onto the flat truth bins. weight_bin() reads flat index
    // k = pt1_bin*nbins + pt2_bin, so bin k+1 belongs to (k/nbins, k%nbins).
    // max() and abs() are both symmetric under swapping the pair, so the weight
    // stays symmetric and the event loop may pass its pT bins in either order.
    auto write_weights = [&]()
    {
      for (int k = 0; k < nbins*nbins; k++)
      {
        const int a = k/nbins, b = k%nbins;
        w.weight->SetBinContent(k + 1, r[group_of(std::max(a, b))][std::abs(a - b)]);
      }
    };
    write_weights();

    // Averaging the two xJ bins a given d feeds is only an approximate inverse
    // of the half-half split, so refine by fixed-point iteration: reweight,
    // reproject, and fold the residual back into r. This stays inside the "one
    // weight per (leading-pT group, pT-bin separation)" family but converges to
    // the member of that family that best reproduces the xJ target.
    for (int iter = 0; iter < n_refine; iter++)
    {
      TH1D *rw = apply(w.truth, w.weight, nbins,
                       Form("h_prior_xjrefine%s_%d", tag.c_str(), iter));
      if (!rw) { break; }

      double max_shift = 0.0;
      for (int g = 0; g < ngroups; g++)
      {
        TH1D *xj_rw = project_to_xj(rw, nbins, ipt_bins, ixj_bins,
                                    groups[g].lo, groups[g].hi, sub_lo, sub_hi,
                                    Form("h_prior_xjrefine_proj%s_%d_g%d", tag.c_str(), iter, g));
        if (!xj_rw) { continue; }
        auto correction_at = [&](const int b) -> double
        {
          if (b < 1 || b > nbins) return -1.0;
          const double c = xj_rw->GetBinContent(b);
          const double t = xj_target[g]->GetBinContent(b);
          if (!(c > 0) || !(t > 0)) return -1.0;
          return t/c;
        };
        for (int d = 0; d < nbins; d++)
        {
          if (!(r[g][d] > 0)) continue;  // a bin the blend deliberately empties
          const double c = combine(d, correction_at);
          if (!(c > 0)) continue;
          r[g][d] *= c;
          max_shift = std::max(max_shift, std::fabs(c - 1.0));
        }
        delete xj_rw;
      }
      write_weights();
      delete rw;
      if (max_shift < 1e-6) { break; }
    }

    for (int g = 0; g < ngroups; g++)
    {
      std::cout << "prior_qa: xJ prior weights, leading bins [" << groups[g].lo
                << ", " << groups[g].hi << ") = " << ipt_bins[groups[g].lo] << " - "
                << ipt_bins[groups[g].hi] << " GeV, by pT-bin separation d:";
      for (int d = 0; d < nbins; d++) { std::cout << " " << d << ":" << r[g][d]; }
      std::cout << std::endl;
      delete xj_truth[g];
      delete xj_target[g];
    }

    update_stats(w, nbins);
    return w;
  }

  // Convenience overload: a single leading-pT window.
  inline weights build_xj(const TH1D *truth_in, const TH1D *unfold_in,
                          const double fraction, const int nbins,
                          const float *ipt_bins, const float *ixj_bins,
                          const int lead_lo, const int lead_hi,
                          const int sub_lo, const int sub_hi,
                          const std::string &tag = "",
                          const int n_refine = 20)
  {
    return build_xj(truth_in, unfold_in, fraction, nbins, ipt_bins, ixj_bins,
                    std::vector<lead_group>{{lead_lo, lead_hi}}, sub_lo, sub_hi,
                    tag, n_refine);
  }

  namespace detail
  {
    inline int &counter() { static int n = 0; return n; }

    inline TH1D *ratio(const TH1D *num, const TH1D *den, const std::string &name)
    {
      TH1D *r = (TH1D*) num->Clone(name.c_str());
      r->SetDirectory(nullptr);
      r->Reset();
      for (int ib = 1; ib <= num->GetNbinsX(); ib++)
      {
        const double d = den->GetBinContent(ib);
        if (d == 0) continue;
        r->SetBinContent(ib, num->GetBinContent(ib)/d);
        r->SetBinError  (ib, num->GetBinError(ib)/d);
      }
      return r;
    }

    // Two-panel "spectrum + ratio to target" plot shared by the flat and xJ QA.
    inline void draw_with_ratio(const TH1D *target, const TH1D *truth, const TH1D *reweighted,
                                const std::string &xtitle, const std::string &ytitle,
                                const std::vector<std::string> &captions,
                                const std::string &outfile,
                                const double ratio_min, const double ratio_max,
                                const bool logy, const int width, const int height,
                                const std::string &target_label = "target prior",
                                const bool ratio_errors = true)
    {
      const int id = detail::counter()++;
      TCanvas *c = new TCanvas(Form("c_prior_qa_%d", id), "", width, height);
      dlutility::ratioPanelCanvas(c, 0.35);

      TH1D *h_target = (TH1D*) target->Clone(Form("h_qa_target_%d", id));
      TH1D *h_truth  = (TH1D*) truth ->Clone(Form("h_qa_truth_%d",  id));
      TH1D *h_rw     = (TH1D*) reweighted->Clone(Form("h_qa_rw_%d", id));
      h_target->SetDirectory(nullptr);
      h_truth ->SetDirectory(nullptr);
      h_rw    ->SetDirectory(nullptr);

      c->cd(1);
      if (logy) gPad->SetLogy();
      dlutility::SetLineAtt(h_target, kBlack,    3, 1);
      dlutility::SetLineAtt(h_truth,  kBlue + 1, 2, 2);
      dlutility::SetLineAtt(h_rw,     kRed + 1,  2, 1);

      h_target->SetTitle(Form(";%s;%s", xtitle.c_str(), ytitle.c_str()));
      dlutility::SetFont(h_target, 42, 0.055, 0.055, 0.05, 0.05);
      h_target->GetYaxis()->SetTitleOffset(1.3);
      const double hi = std::max({h_target->GetMaximum(), h_truth->GetMaximum(), h_rw->GetMaximum()});
      if (logy)
      {
        double lo = 1e30;
        for (int ib = 1; ib <= h_target->GetNbinsX(); ib++)
        {
          const double v = h_target->GetBinContent(ib);
          if (v > 0 && v < lo) lo = v;
        }
        if (!(lo < 1e30) || lo <= 0) lo = hi*1e-8;
        // Leave headroom above the spectrum for the captions and the legend.
        h_target->GetYaxis()->SetRangeUser(0.5*lo, 300.0*hi);
      }
      else
      {
        h_target->GetYaxis()->SetRangeUser(0.0, 1.9*hi);
      }
      h_target->Draw("hist");
      h_truth ->Draw("hist same");
      h_rw    ->Draw("hist same");

      TLegend *leg = new TLegend(0.55, 0.68, 0.94, 0.90);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->AddEntry(h_target, target_label.c_str(), "l");
      leg->AddEntry(h_truth,  "MC truth (unweighted)", "l");
      leg->AddEntry(h_rw,     "reweighted truth", "l");
      leg->Draw();

      double y = 0.86;
      for (const auto &caption : captions)
      {
        dlutility::drawText(caption.c_str(), 0.23, y, 0, kBlack, 0.045);
        y -= 0.055;
      }

      c->cd(2);
      TH1D *r_truth = detail::ratio(h_truth, h_target, Form("h_qa_rtruth_%d", id));
      TH1D *r_rw    = detail::ratio(h_rw,    h_target, Form("h_qa_rrw_%d",    id));
      if (!ratio_errors)
      {
        // Across 225 sparse flat bins the propagated errors are metres tall and
        // hide the closure the panel exists to show.
        for (int ib = 1; ib <= r_truth->GetNbinsX(); ib++)
        {
          r_truth->SetBinError(ib, 0.0);
          r_rw   ->SetBinError(ib, 0.0);
        }
      }
      dlutility::SetLineAtt(r_truth, kBlue + 1, 2, 2);
      dlutility::SetLineAtt(r_rw,    kRed + 1,  2, 1);
      // Markers, not "hist": empty bins would otherwise be drawn as vertical
      // drops to zero and swamp the panel the QA is meant to be read from.
      dlutility::SetMarkerAtt(r_truth, kBlue + 1, 0.7, 24);
      dlutility::SetMarkerAtt(r_rw,    kRed + 1,  0.7, 20);
      r_truth->SetTitle(Form(";%s;X / target", xtitle.c_str()));
      dlutility::SetFont(r_truth, 42, 0.1, 0.1, 0.09, 0.09);
      r_truth->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
      r_truth->GetYaxis()->SetNdivisions(505);
      r_truth->GetYaxis()->SetTitleOffset(0.7);
      r_truth->GetXaxis()->SetTitleOffset(1.1);
      r_truth->Draw("p");
      r_rw   ->Draw("p same");

      TLine *unity = new TLine(r_truth->GetXaxis()->GetXmin(), 1.0,
                               r_truth->GetXaxis()->GetXmax(), 1.0);
      unity->SetLineStyle(2);
      unity->SetLineColor(kGray + 2);
      unity->Draw("same");

      // ROOT sizes the PDF page from gStyle, not from the canvas, so a wide
      // canvas would otherwise be printed into the corner of a portrait page.
      Float_t paper_x = 20, paper_y = 26;
      gStyle->GetPaperSize(paper_x, paper_y);
      gStyle->SetPaperSize(width*2.54/100.0, height*2.54/100.0);
      c->Print(outfile.c_str());
      gStyle->SetPaperSize(paper_x, paper_y);
      delete c;
    }
  }  // namespace detail

  // 1D flat pt1-pt2: target vs MC truth vs reweighted truth, plus ratios to target.
  inline void draw_flat(const weights &w, const TH1D *reweighted,
                        const std::vector<std::string> &captions,
                        const std::string &outfile)
  {
    if (!w.target || !w.truth || !reweighted) return;

    // Trim the empty tail: only the low-index corner of the 225-bin flat space
    // is ever populated, and the blank remainder squeezes the interesting part.
    int last = w.truth->GetNbinsX();
    while (last > 1 && w.truth->GetBinContent(last) <= 0 && w.target->GetBinContent(last) <= 0) { last--; }
    TH1D *target = (TH1D*) w.target->Clone("h_prior_flat_qa_target");
    TH1D *truth  = (TH1D*) w.truth ->Clone("h_prior_flat_qa_truth");
    TH1D *rw     = (TH1D*) reweighted->Clone("h_prior_flat_qa_rw");
    target->SetDirectory(nullptr); truth->SetDirectory(nullptr); rw->SetDirectory(nullptr);
    target->GetXaxis()->SetRange(1, last);
    truth ->GetXaxis()->SetRange(1, last);
    rw    ->GetXaxis()->SetRange(1, last);

    // With kBlendXJ the flat-space target is a REFERENCE, not something the
    // weights are meant to reproduce bin by bin -- say so on the legend so the
    // non-unity ratio is not read as a failure.
    const std::string target_label = (w.space == kBlendXJ)
      ? "per-bin blend (reference)" : "target prior";

    detail::draw_with_ratio(target, truth, rw,
                            "flat (#it{p}_{T,1}, #it{p}_{T,2}) truth bin",
                            "normalized counts", captions, outfile,
                            0.0, 2.0, true, 1100, 700, target_label, false);
    delete target; delete truth; delete rw;
  }

  // xJ projection of the same three distributions over one leading-pt window.
  inline void draw_xj(const weights &w, const TH1D *reweighted, const int nbins,
                      const float *ipt_bins, const float *ixj_bins,
                      const int lead_lo, const int lead_hi,
                      const int sub_lo, const int sub_hi,
                      const std::vector<std::string> &captions,
                      const std::string &outfile)
  {
    if (!w.target || !w.truth || !reweighted) return;

    const int id = detail::counter()++;
    auto project = [&](const TH1D *flat, const char *tag)
    {
      TH2D *h2 = new TH2D(Form("h2_prior_qa_%s_%d", tag, id), "", nbins, ipt_bins, nbins, ipt_bins);
      h2->SetDirectory(nullptr);
      histo_opps::make_sym_pt1pt2((TH1D*) flat, h2, nbins);
      TH1D *hxj = new TH1D(Form("hxj_prior_qa_%s_%d", tag, id), "", nbins, ixj_bins);
      hxj->SetDirectory(nullptr);
      histo_opps::project_xj(h2, hxj, nbins, lead_lo, lead_hi, sub_lo, sub_hi);
      histo_opps::normalize_histo(hxj, nbins);
      delete h2;
      return hxj;
    };

    TH1D *xj_target = project(w.target, "target");
    TH1D *xj_truth  = project(w.truth,  "truth");
    TH1D *xj_rw     = project(reweighted, "rw");

    detail::draw_with_ratio(xj_target, xj_truth, xj_rw,
                            "#it{x}_{J}", "(1/N) dN/d#it{x}_{J}",
                            captions, outfile, 0.0, 2.0, false, 700, 750);

    delete xj_target;
    delete xj_truth;
    delete xj_rw;
  }

  // One-line summary of how far the reweighted truth lands from the target.
  inline void print_summary(const weights &w, const TH1D *reweighted, const int nbins,
                            const std::string &label)
  {
    if (!w.target || !reweighted) return;
    double max_dev = 0.0, sum_dev = 0.0, sum_wt = 0.0;
    for (int ib = 1; ib <= nbins*nbins; ib++)
    {
      const double t = w.target->GetBinContent(ib);
      if (t <= 0) continue;
      const double dev = std::fabs(reweighted->GetBinContent(ib)/t - 1.0);
      if (dev > max_dev) max_dev = dev;
      sum_dev += dev*t;
      sum_wt  += t;
    }
    std::cout << "prior_qa[" << label << "]"
              << "  fraction=" << w.fraction
              << "  weights=[" << w.weight_min << ", " << w.weight_max << "]"
              << "  zero-weight truth bins=" << w.n_zero_weight
              << " (" << 100.0*w.zero_weight_truth_fraction << "% of truth)"
              << "  |reweighted/target - 1|: mean=" << (sum_wt > 0 ? sum_dev/sum_wt : 0.0)
              << " max=" << max_dev << std::endl;
  }

  // Closure in xJ: does the reweighted truth land on the xJ target? This is the
  // number that matters for kBlendXJ, where flat-space closure is not expected.
  inline void print_xj_summary(const weights &w, const TH1D *reweighted, const int nbins,
                               const float *ipt_bins, const float *ixj_bins,
                               const int lead_lo, const int lead_hi,
                               const int sub_lo, const int sub_hi,
                               const std::string &label)
  {
    if (!w.target || !w.truth || !reweighted) return;

    const std::string uid = Form("_%d", detail::counter()++);
    TH1D *xj_target = project_to_xj(w.target, nbins, ipt_bins, ixj_bins,
                                    lead_lo, lead_hi, sub_lo, sub_hi, "h_xjsum_target" + uid);
    TH1D *xj_truth  = project_to_xj(w.truth,  nbins, ipt_bins, ixj_bins,
                                    lead_lo, lead_hi, sub_lo, sub_hi, "h_xjsum_truth" + uid);
    TH1D *xj_rw     = project_to_xj(reweighted, nbins, ipt_bins, ixj_bins,
                                    lead_lo, lead_hi, sub_lo, sub_hi, "h_xjsum_rw" + uid);
    if (!xj_target || !xj_truth || !xj_rw) return;

    double max_dev = 0.0, sum_dev = 0.0, sum_wt = 0.0;
    for (int ib = 1; ib <= nbins; ib++)
    {
      const double t = xj_target->GetBinContent(ib);
      if (t <= 0) continue;
      const double dev = std::fabs(xj_rw->GetBinContent(ib)/t - 1.0);
      if (dev > max_dev) max_dev = dev;
      sum_dev += dev*t;
      sum_wt  += t;
    }
    std::cout << "prior_qa[" << label << "] xJ closure"
              << "  <xJ> truth=" << histo_opps::get_average_xj(xj_truth)
              << " target=" << histo_opps::get_average_xj(xj_target)
              << " reweighted=" << histo_opps::get_average_xj(xj_rw)
              << "  |reweighted/target - 1|: mean=" << (sum_wt > 0 ? sum_dev/sum_wt : 0.0)
              << " max=" << max_dev << std::endl;

    delete xj_target; delete xj_truth; delete xj_rw;
  }

}  // namespace prior_qa

#endif
