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
}

Double_t TVectorMeson::sigma_gamma_p_power_law ( Double_t W ) {
    // sigma(gammap->Vp) = a nb *(W/W0)^delta
    double a = fSigmaGammaPParameters[0];
    double delta = fSigmaGammaPParameters[1];
    double W0 = fSigmaGammaPParameters[2];
    return ( a * pow(W/W0, delta) );
}

Double_t TVectorMeson::sigma_gamma_p_logarithmic ( Double_t W ) {
    // sigma(gammap->Vp) = a nb *log(W/W0)
    double a = fSigmaGammaPParameters[0];
    double W0 = fSigmaGammaPParameters[1];
    return ( a * log(W/W0) );
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
    if ( fModel == kPOWERLAW ) return sigma_gamma_p_power_law(x[0]);
    if ( fModel == kLOGARITHMIC ) return sigma_gamma_p_logarithmic(x[0]);
    if ( fModel == kREGGEOMETRY ) return sigma_gamma_p_reggeometry(x[0]);

    // if nothing better found...
    return -1;
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

    Double_t dummy[99];
    Double_t _W[99];
    _W[0] = W_plus;
    double sigma_gamma_p_plus = sigma_gamma_p(_W, dummy);
    _W[0] = W_minus;
    double sigma_gamma_p_minus = sigma_gamma_p(_W, dummy);

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

void TVectorMeson::set_sigma_gamma_p_model_parameters (Double_t * par, unsigned n) {
    for (unsigned i=0; i<n; i++){
        fSigmaGammaPParameters[i] = par[i];
    }
}

Double_t TVectorMeson::theory_curve (Double_t *x, Double_t *par) {

    if ( fProcess == kPP ) {
        if ( fVariable == kY ) return dSigma_dy(x, par);
        if ( fVariable == kW ) return dSigma_dW(x, par);
    } else if ( fProcess == kGammaP ) {
        return sigma_gamma_p(x,par);
    }

    // if nothing better found...
    return -1;
}
