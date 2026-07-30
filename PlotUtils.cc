#include "PlotUtils.h"

#include <TColor.h>
#include <TError.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <RVersion.h>

TStyle* PlotUtils::sphenix_style()
{
    auto* style = new TStyle("sPHENIX", "sPHENIX style");

    constexpr int color = 0;
    constexpr int font = 42;
    constexpr double text_size = 0.05;

    style->SetFrameBorderMode(color);
    style->SetFrameFillColor(color);
    style->SetCanvasBorderMode(color);
    style->SetCanvasColor(color);
    style->SetPadBorderMode(color);
    style->SetPadColor(color);
    style->SetStatColor(color);

    style->SetPaperSize(20, 26);
    style->SetPadTopMargin(0.05);
    style->SetPadRightMargin(0.05);
    style->SetPadBottomMargin(0.16);
    style->SetPadLeftMargin(0.16);
    style->SetTitleXOffset(1.4);
    style->SetTitleYOffset(1.4);

    style->SetTextFont(font);
    style->SetTextSize(text_size);
    for (const char* axis : {"x", "y", "z"}) {
        style->SetLabelFont(font, axis);
        style->SetTitleFont(font, axis);
        style->SetLabelSize(text_size, axis);
        style->SetTitleSize(text_size, axis);
    }

    style->SetMarkerStyle(20);
    style->SetMarkerSize(1.2);
    style->SetHistLineWidth(2.0);
    style->SetLineStyleString(2, "[12 12]");
    style->SetEndErrorSize(0.0);

    style->SetOptTitle(0);
    style->SetOptStat(0);
    style->SetOptFit(0);
    style->SetPadTickX(1);
    style->SetPadTickY(1);

    style->SetLegendBorderSize(0);
    style->SetLegendFillColor(0);
    style->SetLegendFont(font);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
    style->SetLegendTextSize(text_size);
    style->SetPalette(kBird);
#else
    constexpr int alpha = 0;
    double stops[9] = {0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    double red[9]   = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
    double green[9] = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
    double blue[9]  = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
    TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
#endif

    style->SetNumberContours(80);
    return style;
}

void PlotUtils::set_sphenix_style()
{
    static TStyle* style = sphenix_style();
    (void)style;
    gROOT->SetStyle("sPHENIX");
    gROOT->ForceStyle();
}

void PlotUtils::set_style()
{
    set_sphenix_style();
    gErrorIgnoreLevel = kWarning;
}

void PlotUtils::draw_text(double x, double y, int color, const char* text, float size)
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(size);
    latex.SetTextColor(color);
    latex.DrawLatex(x, y, text);
}
