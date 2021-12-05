// WARNING!
// In order for this script to work, the psi2S-specific normalisation of the gamma-p
// crossection has to be deactivated! (set to 1)
// see TVectorMeson::sigma_gamma_p_reggeometry
//
#include <TXMLEngine.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include <TPlot.h>
#include <TVectorMeson.h>
#include <TXMLParser.h>
#include <TFitter.h>
#include <constants.h>

#include <iostream>
#include <fstream>
using namespace std;

TVectorMeson * jpsi;
TVectorMeson * psi2S;

Double_t R_regge(Double_t *x, Double_t *par) {

   return psi2S->theory_curve(x, par) / jpsi->theory_curve(x, par);
   // return psi2S->sigma_gamma_p(x, par) / jpsi->sigma_gamma_p(x, par);
}

int main (int argc, char **argv) {

    // create TPlot object - inherits from TCanvas
    TPlot plot;
    plot.cd();
    
    plot.set_x_axis_range(1, 1000);
    plot.set_x_axis_title("W [GeV]");

    plot.set_y_axis_range(0.45, 0.47);
    plot.set_y_axis_title("R_Regge = #sigma (#psi(2S))/#sigma (J/#psi)");
    
    plot.set_log_x(false);
    plot.set_log_y(false);

    // creates dummy histo, legend and applies the settings
    plot.Initialise();
    
    jpsi = new TVectorMeson;
    jpsi -> set_meson_type(kJPSI);

    psi2S = new TVectorMeson;
    psi2S -> set_meson_type(kPSI);

    jpsi  -> set_sigma_gamma_p_model (kREGGEOMETRY);
    psi2S -> set_sigma_gamma_p_model (kREGGEOMETRY);
    
    jpsi  -> set_process(kGammaP);
    psi2S  -> set_process(kGammaP);
    
    jpsi -> set_variable ( kW );
    psi2S -> set_variable ( kW );

    jpsi -> Initialise();
    psi2S -> Initialise();

    TF1 * f = new TF1("f", R_regge, 10, 1000, 0);

    // line style settings
    f -> SetNpx(10000);
    f -> SetLineWidth(1);
    f -> SetLineStyle(1);
    f -> SetLineColor(1);
    
    f -> Draw("same");

    plot.Print("../plots/R.eps");
    
    cout << "R_Regge at W=95 GeV: " << psi2S->sigma_gamma_p_reggeometry(95) / jpsi->sigma_gamma_p_reggeometry(95) << endl;
    cout << "R_Regge at W=53.2 GeV: " << psi2S->sigma_gamma_p_reggeometry(53.2) / jpsi->sigma_gamma_p_reggeometry(53.2) << endl;
    cout << "R_Regge at W=128.3 GeV: " << psi2S->sigma_gamma_p_reggeometry(128.3) / jpsi->sigma_gamma_p_reggeometry(128.3) << endl;
    
    cout << "f = R_exp(95)/R_Regge(95): " << 0.166 * (jpsi->sigma_gamma_p_reggeometry(95)/psi2S->sigma_gamma_p_reggeometry(95)) << endl;

    // done
    return 0;
}
