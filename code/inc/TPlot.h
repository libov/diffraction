#ifndef TPLOT_H
#define TPLOT_H

// ROOT headers
#include <TROOT.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

// my headers
#include <constants.h>

// system headers
#include <vector>
using namespace std;

class TPlot : public TCanvas {

    public:

        TPlot();
        ~TPlot(){};

        // setters
        void        set_x_axis_title (TString title)        { fXAxisTitle = title; }
        void        set_y_axis_title (TString title)        { fYAxisTitle = title; }
        void        set_log_x (bool log)                    { fLogX = log; }
        void        set_log_y (bool log)                    { fLogY = log; }
        void        set_x_axis_range (double low, double up){ fLowX = low; fUpX = up; }
        void        set_y_axis_range (double low, double up){ fLowY = low; fUpY = up; }
        void        set_legend_coordinates(double x1, double y1, double x2, double y2) { fLegendX1 = x1; fLegendY1 = y1; fLegendX2 = x2; fLegendY2 = y2; }
        void        set_legend_title (TString title)        { fLegendTitle = title; }

        // getters
        TLegend*    get_legend() {return fLegend;}

        void        Initialise();

        void        AddLegendEntry(TString legend);

    private:

        TH1F        *fHistogram;
        TLegend     *fLegend;
        TString     fXAxisTitle;
        TString     fYAxisTitle;
        bool        fLogX;
        bool        fLogY;
        Double_t    fLowX;
        Double_t    fUpX;
        Double_t    fLowY;
        Double_t    fUpY;
        TString     fLegendTitle;
        Double_t    fLegendX1;
        Double_t    fLegendX2;
        Double_t    fLegendY1;
        Double_t    fLegendY2;

};

extern TPlot *gTPlot;

#endif