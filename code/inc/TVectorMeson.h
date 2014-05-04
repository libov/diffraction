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

        TVectorMeson(){};
        ~TVectorMeson(){};

        // setters
        void        set_cms_energy (Double_t energy)        { fCMSEnergy = energy; }
        void        set_meson_type (MesonType type)         { fMesonType = type; }
        void        set_sigma_gamma_p_model (sigma_gamma_p_model model) { fModel = model; }
        void        set_process (process _process) { fProcess = _process; }
        void        set_variable (variable _variable) { fVariable = _variable; }
        void        set_sigma_gamma_p_model_parameters (Double_t * par, unsigned n);

        // getters
        Double_t    get_cms_energy ()                       { return fCMSEnergy;   }

        // initialiser
        void        Initialise();

        void        draw_theory (sigma_gamma_p_model sigma_model, Double_t param[99], process p, variable v, TString legend_entry, unsigned line_style, unsigned line_color);
        void        draw_data (TString filename, TString legend_entry, unsigned marker_style);

        // formulas
        Double_t    sigma_gamma_p_power_law             ( Double_t W );
        Double_t    sigma_gamma_p_logarithmic           ( Double_t W );
        Double_t    sigma_gamma_p_reggeometry           ( Double_t W );
        Double_t    photon_energy                       ( Double_t rapidity );
        Double_t    W                                   ( Double_t _photon_energy );
        Double_t    OMEGA                               ( Double_t _photon_energy);
        Double_t    photon_flux                         ( Double_t _photon_energy );

        Double_t    dSigma_dy                           ( Double_t *x, Double_t *par );
        Double_t    dSigma_dW                           ( Double_t *x, Double_t *par );
        Double_t    sigma_gamma_p                       ( Double_t *x, Double_t *par );

        Double_t    theory_curve                        ( Double_t *x, Double_t *par );

    private:

        Double_t            fCMSEnergy;
        Double_t            fLorentzGamma;
        MesonType           fMesonType;
        Double_t            fMesonMass;
        Double_t            fMesonBranchingRatio;
        sigma_gamma_p_model fModel;
        process             fProcess;
        variable            fVariable;

        Double_t            fYmin;
        Double_t            fYmax;
        Double_t            fWmin;
        Double_t            fWmax;
        Double_t            fSigmaGammaPParameters[99];
};

#endif
