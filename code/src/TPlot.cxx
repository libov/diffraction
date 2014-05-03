#include <TPlot.h>

TPlot::TPlot():
TCanvas()
{}

void TPlot::Initialise () {

    fHistogram = new TH1F ("","",1, fLowX,fUpX);
    fHistogram -> SetXTitle (fXAxisTitle);
    fHistogram -> SetYTitle (fYAxisTitle);
    if (fLogX) gPad -> SetLogx();
    if (fLogY) gPad -> SetLogy();
    fHistogram -> SetAxisRange(fLowY, fUpY, "Y");
    fHistogram -> SetStats(false);
    fHistogram -> Draw();

    fLegend = new TLegend (fLegendX1, fLegendY1, fLegendX2, fLegendY2);
    fLegend -> SetHeader(fLegendTitle);
    fLegend -> SetFillColor(0);
    fLegend -> SetBorderSize(0);
    fLegend -> SetTextSize(0.03);
    fLegend -> Draw();
}
