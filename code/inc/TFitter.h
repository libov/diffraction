#ifndef TFITTER_H
#define TFITTER_H

// ROOT headers
#include <TROOT.h>
#include <TMinuit.h>
#include <TF1.h>

// my headers
#include <constants.h>

// system headers
using namespace std;

class TFitter : public TMinuit {

    public:

        TFitter(void (*fcn)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t));
        ~TFitter(){};

        Double_t    chi2 ();
        void        set_data(unsigned npoints, Double_t *x, Double_t *y, Double_t *x_err, Double_t *y_err);
        void        set_fit_function(TF1 * f);

    private:

        unsigned        fNpoints;
        Double_t        fX[99];
        Double_t        fX_err[99];
        Double_t        fY[99];
        Double_t        fY_err[99];
        TF1             *fFitFunction;
};

extern TFitter *gTFitter;

#endif