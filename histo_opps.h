#ifndef HISTO_OPPS_H
#define HISTO_OPPS_H

#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLatex.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TColor.h"
#include <cmath>

namespace histo_opps
{

  void get_xj_systematics(TGraphAsymmErrors *h1, TH1D *hsn,TH1D *hsp, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
    {
	    h1->SetPointError(i, hsn->GetBinWidth(i+1)/2.,hsn->GetBinWidth(i+1)/2., fabs(hsn->GetBinContent(i+1))*h1->GetPointY(i),fabs(hsp->GetBinContent(i+1))*h1->GetPointY(i));
    }
  }

  void get_xj_systematics_add(TGraphAsymmErrors *h1, TH1D *hsn,TH1D *hsp, TH1D *hsh, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
    {
      double o_neg_err = hsn->GetBinContent(i+1);
      double o_pos_err = hsp->GetBinContent(i+1);
      double sys_add = (hsh->GetBinContent(i+1));
      std::cout << "err: " << o_pos_err << " / " << o_neg_err <<  " / " << sys_add << std::endl;
      double neg_err = sqrt(o_neg_err*o_neg_err + sys_add * sys_add);
      double pos_err = sqrt(o_pos_err*o_pos_err + sys_add * sys_add);

	    if (!(neg_err > 0)) neg_err = 0;
	    if (!(pos_err > 0)) pos_err = 0;
	    std::cout << "err total : " <<  neg_err << " / " << pos_err << std::endl;	
	    h1->SetPointError(i, hsn->GetBinWidth(i+1)/2.,hsn->GetBinWidth(i+1)/2., neg_err*h1->GetPointY(i), pos_err*h1->GetPointY(i));
    }
  }

  TGraphAsymmErrors *get_xj_statistics(TH1D *h1, const int nbins)
  {
    TGraphAsymmErrors *gg = new TGraphAsymmErrors();
    for (int i = 0; i < nbins; i++)
    {
      gg->Set(i+1);
      gg->SetPoint(i, h1->GetBinCenter(i+1), h1->GetBinContent(i+1));
      gg->SetPointError(i, 0, 0,  h1->GetBinError(i+1), h1->GetBinError(i+1));
    }
    return gg;
  }

  void set_xj_errors(TH1D *h1, TProfile *hp, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
    {
      float old_error = h1->GetBinError(i+1);
      float stat_error = hp->GetBinError(i+1);
      float new_error = sqrt(TMath::Power(old_error, 2) + TMath::Power(stat_error, 2));
      h1->SetBinError(i+1, new_error);
    }
  }
  void set_xj_errors(TH1D *h1, TH1D *hp, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
    {
      float old_error = h1->GetBinError(i+1);
      float stat_error = hp->GetBinError(i+1);
      float new_error = sqrt(TMath::Power(old_error, 2) + TMath::Power(stat_error, 2));
      h1->SetBinError(i+1, new_error);
    }
  }


  void trim_tgraph(TGraphAsymmErrors *g, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
    { 
      if (g->GetPointX(nbins - 1 - i) < first_xj)
      {
        g->RemovePoint(nbins - 1 - i);
      }
    }
    for (int i = 0; i < g->GetN(); i++)
    {
      std::cout << "Trimmed : " << i << " / " << g->GetPointX(i) << " / " << g->GetPointY(i) << std::endl;	  
    }
  }
  void trim_tgraph(TGraph *g, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
    {
      if (g->GetPointX(nbins - 1 - i) < first_xj)
      {
        g->RemovePoint(nbins - 1 -i);
      }
    }
  }
  
  void finalize_xj(TH1D *h1, TH1D *h2, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
    {
      if (h1->GetBinLowEdge(i+1) >= first_xj)
      {
        h2->SetBinContent(i+1, h1->GetBinContent(i+1));
        h2->SetBinError(i+1, h1->GetBinError(i+1));
      }
    }
  }
  void finalize_xj(TProfile *h1, TH1D *h2, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
    {
      if (h1->GetBinLowEdge(i+1) >= first_xj)
      {
        h2->SetBinContent(i+1, h1->GetBinContent(i+1));
        h2->SetBinError(i+1, h1->GetBinError(i+1));
      }
    }
  }

  void normalize_histo( TH1D* h, const int nbins )
  {
    Double_t bin_contents[30] = {0};
    Double_t bin_errors[30] = {0};
    Double_t integral = 0;

    for (int i = 0; i < nbins; i++)
    {

      Double_t v = h->GetBinContent(i+1);
      Double_t w = h->GetBinWidth(i+1);
      Double_t err = h->GetBinError(i+1);

      Double_t vw = v/w;
      Double_t ew = err/w;

      integral += v;

      bin_contents[i] = vw;
      bin_errors[i] = ew;

    }

    if (!std::isfinite(integral) || integral == 0)
    {
      // Empty projections occur in proof samples and can occur in a very
      // sparse systematic bin.  Preserve the empty histogram instead of
      // manufacturing NaNs through a 0/0 normalization.
      return;
    }
    for (int i = 0; i < nbins; i++)
    {
      h->SetBinContent(i+1, bin_contents[i]/integral);
      h->SetBinError(i+1, bin_errors[i]/integral);
    }

    return;
  }

  void tprofile_to_histo(TProfile *hp, TH1D* h, const int nbins)
  {
    for (int ib = 0; ib < nbins; ib++)
    {
      h->SetBinContent(ib+1, hp->GetBinContent(ib+1));
      h->SetBinError(ib+1, hp->GetBinError(ib+1));
    }
  }

  void project_xj( 
    TH2D* h_pt1pt2, 
    TH1D* h_xj, 
    const int nbins, 
    const int start_leading_bin, 
    const int end_leading_bin, 
    const int start_subleading_bin, 
    const int end_subleading_bin
  )
  {

    TH1D *h_unc = (TH1D*) h_xj->Clone();
    h_unc->Reset();

    TH2D *h_asym_pt1pt2 = (TH2D*) h_pt1pt2->Clone();
    h_asym_pt1pt2->Reset();
    for (int ix = 0; ix < nbins; ix++)
    {
	    for (int iy = 0; iy < nbins; iy++)
	    {
	      int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

        if (ix == iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin));
          h_asym_pt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin));
        }
        else if (ix > iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin)*2.);
          h_asym_pt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin)*2);
        }
        else if (ix < iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, 0);
          h_asym_pt1pt2->SetBinError(bin, 0);
        }  
	    }
    }

    for (int ix = 0; ix < nbins; ix++)
    {
	    for (int iy = 0; iy <= ix; iy++)
      {
        int low =  iy - ix - 1;
        //int high = iy - ix + 1;
        int xjbin_low = nbins + low + 1;
        int xjbin_high = nbins + low + 2;

        
        int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

        if (ix < start_leading_bin) continue;
        if (iy < start_subleading_bin) continue;
        if (ix >= end_leading_bin) continue;
        if (iy >= end_subleading_bin) continue;

        if (ix == iy)
        {
          h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin));
          h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin), 2));
        }
        else
        {
          h_xj->Fill(h_xj->GetBinCenter(xjbin_high), h_asym_pt1pt2->GetBinContent(bin)/2.);
          h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin)/2.);
          h_unc->Fill(h_xj->GetBinCenter(xjbin_high),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
          h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
        }
      }
    }

    for (int ix = 0; ix < nbins; ix++)
    {
      h_xj->SetBinError(ix+1, sqrt(h_unc->GetBinContent(ix+1)));
    }

    return;
  }
  
  TH1D *project_pt1(TH2D* h_pt1pt2,  const int nbins)
  {

    TH2D *h_asym_pt1pt2 = (TH2D*) h_pt1pt2->Clone();
    h_asym_pt1pt2->Reset();
    for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
      {
        int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);
        if (ix == iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin));
          h_asym_pt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin));
        }
        else if (ix > iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin)*2.);
          h_asym_pt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin)*2);
        }
        else if (ix < iy)
        {
          h_asym_pt1pt2->SetBinContent(bin, 0);
          h_asym_pt1pt2->SetBinError(bin, 0);
        }      
      }
    }
    TH1D *h_pt1 = (TH1D*) h_asym_pt1pt2->ProjectionX();
    return h_pt1;
  }

  void make_sym_pt1pt2(TH1D *hflat, TH2D* hpt1pt2, const int nbins)
  {

    for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int xbin = ib/nbins;
      int ybin = ib%nbins;
          
      int b = hpt1pt2->GetBin(xbin+1, ybin+1);

      hpt1pt2->SetBinContent(b, hflat->GetBinContent(ib+1));

      hpt1pt2->SetBinError(b, hflat->GetBinError(ib+1));
    }
    return;
  }
  void make_asym_pt1pt2(TH2D *h_apt1pt2, TH2D* h_pt1pt2, const int nbins)
  {

    for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
      {
        int bin = h_apt1pt2->GetBin(ix+1, iy+1);

        if (ix == iy)
        {
          h_apt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin));
          h_apt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin));
        }
        else if (ix > iy)
        {
          h_apt1pt2->SetBinContent(bin, h_pt1pt2->GetBinContent(bin)*2.);
          h_apt1pt2->SetBinError(bin, h_pt1pt2->GetBinError(bin)*2);
        }
        else if (ix < iy)
        {
          h_apt1pt2->SetBinContent(bin, 0);
          h_apt1pt2->SetBinError(bin, 0);
        }
    
	    }
    }

    return;
  }
  void skim_down_histo(TH1D *h_skim, TH1D *h_full, TH1D *h_mapping)
  {

    int nbins = h_mapping->GetXaxis()->GetNbins();

    for (int ib = 0; ib < nbins*nbins; ib++)
    {
      if (!(h_mapping->GetBinContent(ib+1) == 0))
      {
        h_skim->SetBinContent(h_mapping->GetBinContent(ib+1), h_full->GetBinContent(ib+1));
        h_skim->SetBinError(h_mapping->GetBinContent(ib+1), h_full->GetBinError(ib+1));
      }
    }
  }

  void fill_up_histo(TH1D *h_skim, TH1D *h_full, TH1D *h_mapping)
  {
    
    int nbins = h_mapping->GetXaxis()->GetNbins();
 
    for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int bin = h_mapping->GetBinContent(ib+1);
      if (bin)
      {
        h_full->SetBinContent(ib+1, h_skim->GetBinContent(bin));
        h_full->SetBinError(ib+1, h_skim->GetBinError(bin));
      }
    }
  }

  double get_average_xj(TH1D *h_xj)
  {
    double mean = 0;
    double sum = 0;
    int point = 0;
    for (int b=0; b < h_xj->GetNbinsX(); b++)
    {
      double x = h_xj->GetBinCenter(b+1);
      double c = h_xj->GetBinContent(b+1);
      mean += c*x;
      sum += c;
    }
    mean /= sum;
    return mean;
  }
  double get_average_xj(TGraphAsymmErrors *h_xj)
  {
    double mean = 0;
    double sum = 0;
    int point = 0;
    for (int b=0; b < h_xj->GetN(); b++)
    {
      double x = h_xj->GetPointX(b);
      double c = h_xj->GetPointY(b);
      mean += c*x;
      sum += c;
    }
    mean /= sum;
    return mean;
  }

  double get_slope_dphi(TH1D *h_xj)
  {
    TF1 *fPower = new TF1("fPower", "[0]*pow(TMath::Pi()-x, [1]) - [2]", 0.0, TMath::Pi());
    fPower->SetParameters(1, -0.1, 1);
    h_xj->Fit(fPower, "NI");
    double rms = fPower->GetParameter(1);
    return rms;
  }

  std::pair<double, double> get_rms_dphi(TH1D *h_xj)
  {
    double weight = 0;
    double mean = TMath::Pi();
    double sum = 0;
    double var_sum = 0;
    int n = h_xj->GetNbinsX();
    for (int b=0; b < n; b++)
    {
      double x = (h_xj->GetBinCenter(b+1) - mean);
      double c = h_xj->GetBinContent(b+1);
      double dc = h_xj->GetBinError(b+1);
      weight += c*x*x;
      sum += c;
      var_sum += dc*dc*x*x*x*x;
    }

    weight /= sum;

    double dd = sqrt(var_sum) / sum;
    
    double rms = sqrt(weight);
    double drms = dd / (2.0 * rms);
    return std::make_pair(rms, dd);
  }
  TH1D *get_correction_smooth(TH1D *h1)
  {
    TH1D *h = (TH1D*) h1->Clone(Form("%s_smooth", h1->GetName()));

    TF1 *f = (TF1*) h1->GetFunction("p3");

    for (int ib = 0; ib < h1->GetNbinsX(); ib++)
    {
      double a = h1->GetBinLowEdge(ib+1);
      double w = h1->GetBinWidth(ib+1);;
      double b = a + w;
      double v = f->Integral(a, b)/w;
      h->SetBinContent(ib+1, v);
    }
    return h;
  }

  TH1D *average_results(TH1D *h_final, TH1D *h_add, const std::string new_name = "h_sys_add")
  {
    TH1D *h_orig = (TH1D*) h_final->Clone("h_orig");
    TH1D *h_sys = (TH1D*) h_final->Clone(new_name.c_str());

    for (int ibin = 0; ibin < h_final->GetNbinsX(); ibin++)
    {
      double og = h_orig->GetBinContent(ibin+1);
      double v = (h_orig->GetBinContent(ibin+1) + h_add->GetBinContent(ibin+1))/2.;
      double vd = (h_orig->GetBinContent(ibin+1) - h_add->GetBinContent(ibin+1))/(2.*v);
      double ve = (h_orig->GetBinError(ibin+1))*v/og;

      h_final->SetBinContent(ibin+1, v);
      h_final->SetBinError(ibin+1, ve);
      h_sys->SetBinContent(ibin+1, vd);
    }
    return h_sys;
  }
  
};

#endif
