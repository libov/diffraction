// ROOT headers
#include <TF1.h>
#include <TLegend.h>
#include <TObjString.h>
#include <TGraphErrors.h>

// my headers
#include <TVectorMeson.h>
#include <constants.h>

// system headers
#include <iostream>
#include<fstream>
#include <vector>
#include <getopt.h>
using namespace std;

TVectorMeson::TVectorMeson() {
}

void    TVectorMeson::Initialise () {

    // calculate Lorentz gamma from the CMS energy
    fLorentzGamma =  fCMSEnergy / (2 * PROTON_MASS) ;
    
    // set vector meson mass
    if (fMesonType == kJPSI) {

        fMesonMass = JPSI_MASS;
        fMesonBranchingRatio = BR_JPSI_MUMU;

    } else if (fMesonType == kPSI) {

        fMesonMass = PSI_MASS;
        fMesonBranchingRatio = BR_PSI_MUMU;

    } else {

        cout << "ERROR: unknown meson type" << endl;
        abort();
    }

    // calculate max and min rapidity as well as the corresponding min max W
    Double_t y_abs_max = - TMath::Log(fMesonMass/fCMSEnergy);
    fYmin = -y_abs_max;
    fYmax =  y_abs_max;
    fWmin = sqrt (fMesonMass * fCMSEnergy * exp(fYmin));
    fWmax = sqrt (fMesonMass * fCMSEnergy * exp(fYmax));
    cout << fYmin << " " << fYmax << endl;
    cout << fWmin << " " << fWmax << endl;

    fCanvas = new TCanvas();

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

void TVectorMeson::draw_theory (sigma_gamma_p_model sigma_model, process p, variable v, TString legend_entry, unsigned line_style, unsigned line_color) {

    double par[1];
    
    if ( sigma_model == kPOWERLAW ) {
        par[0] = 1;
    } else if ( sigma_model == kREGGEOMETRY ) {
        par[0] = 2;
    } else {
        cout << "ERROR: unknown gamma-p cross section model" << endl;
        abort();
    }

    TF1 * f;
    if ( p == kPP ) {
        if ( v == kY ) f = new TF1("f", this, &TVectorMeson::dSigma_dy, fLowX, fUpX, 1, "","");
        if ( v == kW ) f = new TF1("f", this, &TVectorMeson::dSigma_dW, fLowX, fUpX, 1, "","");
    } else if ( p == kGammaP ) {
        f = new TF1("f", this, &TVectorMeson::sigma_gamma_p, fLowX, fUpX, 1, "","");
    }
    
    f -> SetParameters(par);
    f -> SetNpx(10000);
    f -> SetLineWidth(2);
    f -> SetLineStyle(line_style);
    f -> SetLineColor(line_color);
    f -> Draw("same");

    fLegend -> AddEntry(f, legend_entry, "l");
}

void TVectorMeson::draw_data (TString filename, TString legend_entry, unsigned marker_style) {

    // open the file
    ifstream f(filename);
    if (!f.is_open()) {
        cout << "ERROR: Unable to open file " << filename << endl;
        abort();
    }

    Double_t y[99];
    Double_t y_err[99];
    Double_t sigma[99];
    Double_t sigma_err[99];
    unsigned npoints = 0;
    
    string line;
    while ( f.good() ) {

        // read each line
        getline (f,line);

        // skip if an empty line
        if (line=="") continue;
        // tokenize
        TString line_str = line;
        TObjArray * tokens = line_str.Tokenize(" ");

        if (tokens -> GetEntries() == 0) continue;
        
        // check if this line is a comment
        TString first_word = ((TObjString*)tokens->At(0)) -> GetString();
        char first_char = first_word[0];
        if (first_char=='#') continue;

        y[npoints] = ((TObjString*)tokens->At(0)) -> GetString().Atof();
        y_err[npoints] = 0;
        sigma[npoints] = ((TObjString*)tokens->At(1)) -> GetString().Atof();
        sigma_err[npoints] = ((TObjString*)tokens->At(2)) -> GetString().Atof();
        
        npoints++;
    }
    cout << "INFO: npoints= " << npoints << endl;
    for (unsigned i=0;i<npoints; i++){
        cout << y[i] << " " << sigma[i] << " " << sigma_err[i] << endl;
    }
    
    TGraphErrors * g = new TGraphErrors(npoints, y, sigma, y_err, sigma_err);
    g -> SetMarkerStyle(marker_style);
    g -> Draw("p");
    
    fLegend -> AddEntry(g, legend_entry, "p");
}

void TVectorMeson::Print(TString filename) {
    fCanvas -> Print(filename);
}

Double_t TVectorMeson::sigma_gamma_p_power_law ( Double_t W ) {
    // sigma(gammap->Vp) = 1.5 nb *(W/W0)^0.8, W0 = 1 GeV
    return ( 1.5 * pow(W, 0.8) );
}

Double_t TVectorMeson::sigma_gamma_p_reggeometry ( Double_t W ) {

    double A0           = 29.8;
    double Q20          = 2.1;
    double n            = 1.37;
    double alpha0       = 1.2;
    double alphaprim    = 0.17;
    double a            = 1.01;
    double b            = 0.44;
    double W0           = 1.;

    double Q2tilde = pow(fMesonMass, 2);

    double term1 = A0*A0;
    double term2 = pow( W/W0, 2*2*(alpha0-1) ); // additional factor 2 because s=W^2
    double term3 = pow(1+Q2tilde/Q20, 2*n);
    double term4 = 2*2*alphaprim*log(W/W0) + 4*(a/Q2tilde+b/(2*pow(PROTON_MASS,2)) ) ; // additional factor 2 in the first term because s=W^2

    return ( (term1*term2) / (term3*term4) );
}

Double_t TVectorMeson::sigma_gamma_p (Double_t *x, Double_t *par) {
    if ( par[0] == 1 ) return sigma_gamma_p_power_law(x[0]);
    if ( par[0] == 2 ) return sigma_gamma_p_reggeometry(x[0]);
}

Double_t TVectorMeson::photon_energy ( Double_t rapidity ) {
    return ( 0.5 * fMesonMass * exp(rapidity) );
}

Double_t TVectorMeson::W ( Double_t _photon_energy ) {
    return sqrt( 2 * _photon_energy * fCMSEnergy );
}

Double_t TVectorMeson::OMEGA ( Double_t _photon_energy) {
    return ( 1 + 0.71 * pow(fLorentzGamma, 2.) / pow(_photon_energy, 2.) );
}

Double_t TVectorMeson::photon_flux ( Double_t _photon_energy ) {

    double _OMEGA = OMEGA ( _photon_energy );

    double prefactor = ALPHA / (2 * PI * _photon_energy );
    double term1 = 1 + TMath::Power( 1 - 2 * _photon_energy / fCMSEnergy, 2.);
    double term2 = log ( _OMEGA ) - 11./6. + 3. / _OMEGA - 3. / (2. * pow(_OMEGA, 2.) ) + 1./ ( 3. * pow(_OMEGA, 3.) );

    return ( prefactor * term1 * term2 );
}

Double_t TVectorMeson::dSigma_dy (Double_t *x, Double_t *par) {

    double rapidity = x[0];

    double photon_energy_plus  = photon_energy ( rapidity );
    double photon_energy_minus = photon_energy ( -rapidity );

    double W_plus  = W(photon_energy_plus);
    double W_minus = W(photon_energy_minus);

    double sigma_gamma_p_plus;
    double sigma_gamma_p_minus;

    if (par[0] == 1) {
        sigma_gamma_p_plus  = sigma_gamma_p_power_law(W_plus);
        sigma_gamma_p_minus = sigma_gamma_p_power_law(W_minus);
    } else if (par[0] == 2) {
        sigma_gamma_p_plus  = sigma_gamma_p_reggeometry(W_plus);
        sigma_gamma_p_minus = sigma_gamma_p_reggeometry(W_minus);
    } else {
        cout << "ERROR: par[0] not set for sigma_Y function. Please choose the model for sigma_gamma_p cross section" << endl;
        abort();
    }

    return S2GAP * (
             photon_energy_plus  * photon_flux(photon_energy_plus)  * sigma_gamma_p_plus
           + photon_energy_minus * photon_flux(photon_energy_minus) * sigma_gamma_p_minus
           );
}

Double_t TVectorMeson::dSigma_dW (Double_t *x, Double_t *par) {

    double W = x[0];

    double photon_energy = W*W/(2*fCMSEnergy);

    // my way (probably wrong)
//     double sigma_gamma_p;
//     if (par[0] == 1) {
//         sigma_gamma_p = sigma_gamma_p_power_law(W);
//     } else if (par[0] == 2) {
//         sigma_gamma_p = sigma_gamma_p_reggeometry(W);
//     } else {
//         cout << "ERROR: par[0] not set for sigma_W function. Please choose the model for sigma_gamma_p cross section" << endl;
//         abort();
//     }
//     return ( photon_flux(omega) * sigma_gamma_p * W/sqrt(S) );

    // Magno's way
    Double_t y[1] = {log(2*photon_energy/fMesonMass)};
    return ( dSigma_dy(y, par)/W );
}
