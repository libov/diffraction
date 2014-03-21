#include<TF1.h>
#include<TH1F.h>
#include<TMath.h>
#include<TCanvas.h>
#include<TLatex.h>
#include<TLegend.h>

// =============== constants  ================ //
const double S     = TMath::Power(7000., 2);
const double MP    = 0.938272;
const double GAMMA = 0.5 * TMath::Sqrt(S) / MP;
const double MV    = 3.096916;
const double DELTA = 0.8;
const double ALPHA = 1./137;
const double PI    = TMath::Pi();
const double S2GAP = 0.8;

// ========== helping functions ========== //
TString toStr(Double_t arg, Int_t decimals);
TString toStr(Int_t arg);

TString toStr(Double_t arg, Int_t decimals) {
    char tmp[256];
    TString format = "%."+toStr(decimals)+"f";
    sprintf(tmp, format, arg);
    TString result(tmp);
    return result;
}

TString toStr(Int_t arg) {
    TString result;
    result += arg;
    return result;
}

// ========== function definitions =========== //
double omega(double Y) { return (0.5 * MV * exp(Y) ); }
double omega_plus(double Y) { return omega(Y); }
double omega_minus(double Y) { return omega(-Y); }

double W(double om) { return sqrt( 2  * om * sqrt(S) ); }
double OMEGA(double om) { return ( 1 + 0.71 * pow(GAMMA, 2.) / pow(om, 2.) ); }

double photon_flux(double om) {

    double OM=OMEGA(om);

    double prefactor = ALPHA / (2 * PI * om );
    double term1 = 1 + TMath::Power( 1 - 2*om/sqrt(S), 2.);
    double term2 = log(OM) - 11./6. + 3./OM - 3. / (2. * pow(OM, 2.) ) + 1./ ( 3. * pow(OM, 3.) );

    return ( prefactor * term1 * term2 );
}

double sigma_gamma_p_simple(double W) {
    // sigma(gammap->Vp) = 1.5 nb *(W/W0)^0.8
    return ( 1.5 * pow(W, DELTA) );
}

double sigma_gamma_p_realistic(double W) {
    double A0=29.8;
    double Q20=2.1;
    double n=1.37;
    double alpha0=1.2;
    double alphaprim=0.17;
    double a=1.01;
    double b=0.44;
    double W0=1.;

    double Q2tilde=MV*MV;

    double term1 = A0*A0;
    double term2 = pow( W/W0, 2*2*(alpha0-1) ); // additional factor 2 because s=W^2
    double term3 = pow(1+Q2tilde/Q20, 2*n);
    double term4 = 2*2*alphaprim*log(W/W0) + 4*(a/Q2tilde+b/(2*MP*MP) ) ; // additional factor 2 in the first term because s=W^2

    return ( (term1*term2) / (term3*term4) );
}

double sigma_Y (Double_t *x, Double_t *par) {

    double y=x[0];

    double om_pl = omega_plus(y);
    double om_mi = omega_minus(y);

    double W_plus = W(om_pl);
    double W_minus = W(om_mi);

    double sigma_gamma_p_plus;
    double sigma_gamma_p_minus;

    if (par[0] == 1) {
        sigma_gamma_p_plus  = sigma_gamma_p_simple(W_plus);
        sigma_gamma_p_minus = sigma_gamma_p_simple(W_minus);
    } else if (par[0] == 2) {
        sigma_gamma_p_plus  = sigma_gamma_p_realistic(W_plus);
        sigma_gamma_p_minus = sigma_gamma_p_realistic(W_minus);
    } else {
        cout << "ERROR: par[0] not set for sigma_Y function. Please choose the model for sigma_gamma_p cross section" << endl;
        abort();
    }

    return S2GAP * (
             om_pl * photon_flux(om_pl) * sigma_gamma_p_plus
           + om_mi * photon_flux(om_mi) * sigma_gamma_p_minus
           );
}

double sigma_W (Double_t *x, Double_t *par) {

    double W = x[0];

    double omega = W*W/(2*sqrt(S));

    double sigma_gamma_p;

    if (par[0] == 1) {
        sigma_gamma_p = sigma_gamma_p_simple(W);
    } else if (par[0] == 2) {
        sigma_gamma_p = sigma_gamma_p_realistic(W);
    } else {
        cout << "ERROR: par[0] not set for sigma_W function. Please choose the model for sigma_gamma_p cross section" << endl;
        abort();
    }

    // return ( photon_flux(omega) * sigma_gamma_p * W/sqrt(S) );
    Double_t y[1] = {log(2*omega/MV)};
    return ( sigma_Y(y, par)/W );
}


// ============== plotting routine ============== //
void make_plot (double (*fcn)(Double_t *x, Double_t *par), Double_t low, Double_t up, Double_t * param, TString filename, TString x_axis_title, TString y_axis_title, TString legend_entry, bool logx, bool logy) {

    TCanvas *c = new TCanvas();

    TF1 * f = new TF1("f", fcn, low, up, 1);

    f -> SetParameters(param);

    f -> SetNpx(10000);
    f -> GetHistogram() -> SetXTitle (x_axis_title);
    f -> GetHistogram() -> SetYTitle (y_axis_title);
    f -> SetLineWidth(3);
    if (logx) gPad -> SetLogx();
    if (logy) gPad -> SetLogy();
    f -> Draw();
    ((TH1F*)f -> GetHistogram()) -> SetAxisRange(0, 1.25 * f -> GetMaximum(), "Y" );

    TLegend * leg = new TLegend (0.45, 0.8, 0.85, 0.9);
    leg -> AddEntry(f, legend_entry, "l");
    leg -> SetFillColor(0);
    leg -> SetBorderSize(0);
    leg -> SetTextSize(0.03);
    leg -> Draw();

    c -> Print(filename);
}

// ============== main function =============== //
int paper_plots () {

    // welcome info
    cout << "\nsqrt(s) = " << sqrt(S) << " GeV " << endl;
    cout << "gamma = " << GAMMA << endl;
    cout << "vector meson mass MV = " << MV << " GeV " << endl;
    cout << "delta = " << DELTA << endl << endl;

    TString legend_entry_prefix = "J/#psi at LHC (#sqrt{s} = " + toStr(sqrt(S)/1000, 0) + " TeV)";
    TString legend_entry_ending_simple = ", #sigma_{#gamma p#rightarrow V p} #propto W^{0.8}";
    TString legend_entry_ending_realistic = ", #sigma_{#gamma p#rightarrow V p} from eq.(2)";
    TString y_axis_title_prefix = "d#sigma(h_{1}+h_{2} #rightarrow h_{1}+V+h_{2})";

    Double_t par[99];

    Double_t y_abs_max = - TMath::Log(MV/sqrt(S));
    Double_t ymin = -y_abs_max;
    Double_t ymax =  y_abs_max;

    Double_t Wmin = sqrt (MV * sqrt(S) * exp(ymin));
    Double_t Wmax = sqrt (MV * sqrt(S) * exp(ymax));

    // simple case, sigma(gamma p->V p)=W^0.8
    par[0]=1;
    make_plot(sigma_Y, ymin, ymax, par, "sigma_y_simple.eps", "Rapidity Y",  y_axis_title_prefix+" / dY [nb]",     legend_entry_prefix+legend_entry_ending_simple, false, false);
    make_plot(sigma_W, Wmin, Wmax, par, "sigma_w_simple.eps",     "W [GeV]", y_axis_title_prefix+" / dW [nb/GeV]", legend_entry_prefix+legend_entry_ending_simple, true, false);

    // realistic case, Laszlo's model
    par[0]=2;
    make_plot(sigma_Y, ymin, ymax, par, "sigma_y_realistic.eps", "Rapidity Y", y_axis_title_prefix+" / dY [nb]",     legend_entry_prefix+legend_entry_ending_realistic, false, false);
    make_plot(sigma_W, Wmin, Wmax, par, "sigma_w_realistic.eps", "W [GeV]",    y_axis_title_prefix+" / dW [nb/GeV]", legend_entry_prefix+legend_entry_ending_realistic, true, false);

    return 0;
}
