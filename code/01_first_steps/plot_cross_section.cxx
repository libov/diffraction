#include<TF1.h>
#include<TMath.h>

const double S     = TMath::Power(14000., 2);
const double MP    = 0.94;
const double GAMMA = 0.5 * TMath::Sqrt(S) / MP;
const double MV    = 0.77;
const double DELTA = 0.8;
const double ALPHA = 1./137;
const double PI    = TMath::Pi();

double omega(double y);
double OMEGA(double y);
double W(double y);
double sigma_y(double y);
double photon_flux(double y);

int plot_cross_section () {

    cout << "sqrt(s) = " << sqrt(S) << " GeV " << endl;
    cout << "gamma = " << GAMMA << endl;
    cout << "vector meson mass MV = " << MV << " GeV " << endl;
    cout << "delta = " << DELTA << endl;

    TF1 * f = new TF1("f","sigma_y(x)", -10, 10);
    f -> SetNpx(1000);
    f->Draw();

    return 0;
}

double omega(double y) {

    // exponent argument should be negative
    //if (y>0) y *= (-1);

    return (0.5 * MV * exp(abs(y)) );

}

double OMEGA(double y) {

    return ( 1 + 0.71 * pow(GAMMA, 2) / pow(omega(y), 2.) );
}

double W(double y) {

    return sqrt( 2  * omega(y) * sqrt(S) );

}

double photon_flux(double y) {

    return ( (ALPHA/2/PI/omega(y)) * (1+TMath::Power((1-2*omega(y)/sqrt(S)), 2.)) * (TMath::Log(OMEGA(y))-11./6 + 3./OMEGA(y) - 3./2./TMath::Power(OMEGA(y),2.)+1./3./TMath::Power(OMEGA(y),3.)));
}

double sigma_y (double y) {

    return ( omega(y) * photon_flux(y) * pow( W(y), DELTA) );

}
