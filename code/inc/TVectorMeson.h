#ifndef TVECTORMESON_H
#define TVECTORMESON_H

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

class TVectorMeson {

    public:

        TVectorMeson();
        ~TVectorMeson(){};

        // setters
        void        set_cms_energy (Double_t energy)        { fCMSEnergy = energy; }
        void        set_meson_type (MesonType type)         { fMesonType = type; }
        void        set_x_axis_title (TString title)        { fXAxisTitle = title; }
        void        set_y_axis_title (TString title)        { fYAxisTitle = title; }
        void        set_log_x (bool log)                    { fLogX = log; }
        void        set_log_y (bool log)                    { fLogY = log; }
        void        set_x_axis_range (double low, double up){ fLowX = low; fUpX = up; }
        void        set_y_axis_range (double low, double up){ fLowY = low; fUpY = up; }
        void        set_legend_coordinates(double x1, double y1, double x2, double y2) { fLegendX1 = x1; fLegendY1 = y1; fLegendX2 = x2; fLegendY2 = y2; }
        void        set_legend_title (TString title)        { fLegendTitle = title; }

        // getters
        Double_t    get_cms_energy ()                       { return fCMSEnergy;   }

        // initialiser
        void        Initialise();
        void        Print(TString filename);
        void        draw_theory (sigma_gamma_p_model sigma_model, process p, variable v, TString legend_entry, unsigned line_style, unsigned line_color);
        void        draw_data (TString filename, TString legend_entry);

        // formulas
        Double_t    sigma_gamma_p_power_law             ( Double_t W );
        Double_t    sigma_gamma_p_reggeometry           ( Double_t W );
        Double_t    photon_energy                       ( Double_t rapidity );
        Double_t    W                                   ( Double_t _photon_energy );
        Double_t    OMEGA                               ( Double_t _photon_energy);
        Double_t    photon_flux                         ( Double_t _photon_energy );
        
        Double_t    dSigma_dy                           ( Double_t *x, Double_t *par );
        Double_t    dSigma_dW                           ( Double_t *x, Double_t *par );
        
    private:
        
        Double_t    fCMSEnergy;
        Double_t    fLorentzGamma;
        MesonType   fMesonType;
        Double_t    fMesonMass;
        Double_t    fMesonBranchingRatio;
        Double_t    fYmin;
        Double_t    fYmax;
        Double_t    fWmin;
        Double_t    fWmax;

        TCanvas     *fCanvas;
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

extern TVectorMeson *gTVectorMeson;

#endif
