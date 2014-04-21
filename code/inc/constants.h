#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <TMath.h>

const double PROTON_MASS    = 0.938272;
const double JPSI_MASS      = 3.096916;
const double PSI_MASS       = 3.686109;

const double BR_JPSI_MUMU   = 0.0593;
const double BR_PSI_MUMU    = 0.0078;

const double PI             = TMath::Pi();
const double ALPHA          = 1./137;

const double S2GAP          = 0.8;
        
enum MesonType              { kJPSI, kPSI };
enum sigma_gamma_p_model    { kPOWERLAW, kREGGEOMETRY };
enum process                { kGammaP, kPP};
enum variable               { kY, kW }; 

#endif