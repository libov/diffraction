#include <TFitter.h>

TFitter *gTFitter;

TFitter::TFitter(void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t)):
TMinuit(),
fFitFunction(0)
{
    SetFCN(fcn);
    gTFitter = this;
}

void TFitter::set_data(unsigned npoints, Double_t *x, Double_t *y, Double_t *x_err, Double_t *y_err) {
    fNpoints = npoints;
    for (unsigned i=0; i<fNpoints; i++){
        fX[i] = x[i];
        fY[i] = y[i];
        fX_err[i] = x_err[i];
        fY_err[i] = y_err[i];
    }
}

void TFitter::set_fit_function(TF1 *f) {
    if (fFitFunction != 0) {
        cout << "ERROR: the fit function was already set" << endl;
        abort();
    }
    fFitFunction = f;
}

Double_t TFitter::chi2() {

    Double_t result = 0;

    for (unsigned i=0; i<fNpoints; i++) {
        double deviation = fFitFunction->Eval(fX[i]) - fY[i];
        result += pow( deviation/fY_err[i], 2 );
    }

    return result;
}
