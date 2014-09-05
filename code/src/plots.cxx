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

TVectorMeson *gTVectorMeson;

void minimization_function(Int_t& npar, Double_t* grad, Double_t& f, Double_t par[], Int_t iflag) {
    gTVectorMeson -> set_sigma_gamma_p_model_parameters(par, gTVectorMeson -> get_n_sigma_gamma_p_model_parameters() );
    f = gTFitter -> chi2();
}

int main (int argc, char **argv) {

    // create TPlot object - inherits from TCanvas
    TPlot plot;
    plot.cd();

    // create TFitter object - needed for the fit (if requested in the xml file)
    TFitter fitter(minimization_function);

    // create parser instance, pass the xml file name
    TString  xmlfilename = argv[1];
    TXMLParser parser(xmlfilename);

    // top level properties
    double cms_energy = parser.getNodeContent("cms_energy").Atof();
    TString meson_type_str = parser.getNodeContent("meson_type");
    MesonType meson_type;
    if (meson_type_str == "jpsi") {
        meson_type = kJPSI;
    } else if (meson_type_str == "psi") {
        meson_type = kPSI;
    } else {
        cout << "ERROR: meson type " << meson_type_str << " not known " << endl;
        abort();
    }

    // legend
    parser.selectNode("legend");
    plot.set_legend_title(parser.getNodeContent("title"));
    double x1, x2, y1, y2;
    x1 = parser.getNodeContent("x1").Atof();
    x2 = parser.getNodeContent("x2").Atof();
    y1 = parser.getNodeContent("y1").Atof();
    y2 = parser.getNodeContent("y2").Atof();
    plot.set_legend_coordinates(x1, y1, x2, y2);

    // x-axis
    parser.selectMainNode();
    parser.selectNode("x_axis");
    plot.set_x_axis_title(parser.getNodeContent("title"));
    TString xlog = parser.getNodeContent("log");
    if (xlog == "false") plot.set_log_x(false);
    else if (xlog == "true") plot.set_log_x(true);
    double x_low_range = parser.getNodeContent("low_range").Atof();
    double x_up_range = parser.getNodeContent("up_range").Atof();
    plot.set_x_axis_range(x_low_range, x_up_range);

    // y-axis
    parser.selectMainNode();
    parser.selectNode("y_axis");
    plot.set_y_axis_title(parser.getNodeContent("title"));
    TString ylog = parser.getNodeContent("log");
    if (ylog == "false") plot.set_log_y(false);
    else if (ylog == "true") plot.set_log_y(true);
    double y_low_range = parser.getNodeContent("low_range").Atof();
    double y_up_range = parser.getNodeContent("up_range").Atof();
    plot.set_y_axis_range(y_low_range, y_up_range);

    parser.selectMainNode();

    // creates dummy histo, legend, applies settings
    plot.Initialise();

    // perform plotting - loop over plot tags, create TF1 object and draw it (canvas was selected above)
    parser.selectMainNode();
    parser.selectNode("plot");
    while (parser.getCurrentNode() != 0) {

        TVectorMeson * meson = new TVectorMeson;

        meson -> set_cms_energy(cms_energy);
        meson -> set_meson_type(meson_type);

        if ( parser.getNodeContent("sigma_gamma_p_model") == "powerlaw") {
            meson -> set_sigma_gamma_p_model (kPOWERLAW);
        } else if ( parser.getNodeContent("sigma_gamma_p_model") == "logarithmic") {
            meson -> set_sigma_gamma_p_model (kLOGARITHMIC);
        } else if ( parser.getNodeContent("sigma_gamma_p_model") == "reggeometry") {
            meson -> set_sigma_gamma_p_model (kREGGEOMETRY);
        }

        if (parser.getNodeContent("process") == "pp") {
            meson -> set_process ( kPP );
        } else if (parser.getNodeContent("process") == "gammap") {
            meson -> set_process ( kGammaP );
        }

        if (parser.getNodeContent("variable") == "y") {
            meson -> set_variable ( kY );
        } else if (parser.getNodeContent("variable") == "w") {
            meson -> set_variable ( kW );
        }

        TString par_str = parser.getNodeContent("parameters");
        TObjArray * tokens = par_str.Tokenize(" ");
        Double_t            par[99];
        for (int i=0; i<tokens -> GetEntries(); i++) {
            par[i] = ((TObjString*)tokens->At(i)) -> GetString().Atof();
        }
        meson -> set_sigma_gamma_p_model_parameters(par, tokens -> GetEntries());

        meson -> Initialise();

        TF1 * f = new TF1("f", meson, &TVectorMeson::theory_curve, x_low_range, x_up_range, 0, "","");

        // line style settings
        f -> SetNpx(10000);
        f -> SetLineWidth(1);
        f -> SetLineStyle(parser.getNodeContent("line_style").Atoi());
        f -> SetLineColor(parser.getNodeContent("line_color").Atoi());
        f -> Draw("same");

        plot.get_legend() -> AddEntry(f, parser.getNodeContent("legend_entry"), "l");

        // if this is a fit function
        if (parser.getNodeContent("use_in_the_fit") == "true") {
            fitter.set_fit_function(f);
            gTVectorMeson = meson;
            for (int i=0; i < tokens -> GetEntries(); i++) {
                fitter.DefineParameter(i, "", par[i], 0.1, 0, 0);
            }
            TString fixpar_str = parser.getNodeContent("fix_parameters");
            TObjArray * tokens = fixpar_str.Tokenize(" ");
            for (int i=0; i<tokens -> GetEntries(); i++) {
               if (((TObjString*)tokens->At(i)) -> GetString().Atof()) fitter.FixParameter(i);
            }
        }

        parser.selectNextNode("plot");
    }

    // plot data
    parser.selectMainNode();
    parser.selectNode("dataplot");
    while (parser.getCurrentNode() != 0) {

        // open the file
        ifstream f(parser.getNodeContent("filename"));
        if (!f.is_open()) {
            cout << "ERROR: Unable to open file " << parser.getNodeContent("filename") << endl;
            abort();
        }

        Double_t y[99];
        Double_t sigma[99];
        Double_t sigma_stat_err[99];
        Double_t sigma_tot_up[99];
        Double_t sigma_tot_down[99];
        Double_t zeros[99];
        unsigned npoints = 0;
        bool     systematics_given = false;
        for (unsigned i=0; i<99; i++) zeros[i] = 0;

        string line;
        while ( f.good() ) {

            // read each line
            getline (f,line);

            // skip if an empty line
            if (line=="") continue;
            // tokenize
            TString line_str = line;
            TObjArray * tokens = line_str.Tokenize(" ");

            unsigned nentries = tokens -> GetEntries();
            if ( nentries == 0) continue;

            // check if this line is a comment
            TString first_word = ((TObjString*)tokens->At(0)) -> GetString();
            char first_char = first_word[0];
            if (first_char=='#') continue;

            y[npoints] = ((TObjString*)tokens->At(0)) -> GetString().Atof();
            sigma[npoints] = ((TObjString*)tokens->At(1)) -> GetString().Atof();
            sigma_stat_err[npoints] = ((TObjString*)tokens->At(2)) -> GetString().Atof();
            sigma_tot_up[npoints] = sigma_stat_err[npoints];
            sigma_tot_down[npoints] = sigma_stat_err[npoints];

            //if systematic uncertainties are also given
            if ( nentries == 5) {
                systematics_given = true;
                Double_t unc1 = ((TObjString*)tokens->At(3)) -> GetString().Atof();
                Double_t unc2 = ((TObjString*)tokens->At(4)) -> GetString().Atof();
                Double_t unc_up, unc_down;
                if ( (unc1>=0) && (unc2<=0) ) {
                    unc_up = unc1;
                    unc_down = unc2;
                } else if ( (unc1<=0) && (unc2>=0) ) {
                    unc_up = unc2;
                    unc_down = unc1;
                } else {
                    cout << "ERROR: both systematic uncertainties are of same sign! " << endl;
                    abort();
                }
                sigma_tot_up[npoints]   = sqrt( pow(unc_up,2) + pow(sigma_stat_err[npoints],2) );
                sigma_tot_down[npoints] = sqrt( pow(unc_down,2) + pow(sigma_stat_err[npoints],2) );
            }

            npoints++;
        }
        cout << "INFO: npoints= " << npoints << endl;
        for (unsigned i=0;i<npoints; i++){
            cout << y[i] << " " << sigma[i] << " " << sigma_stat_err[i] << " +"<<sigma_tot_up[i] << " -" << sigma_tot_down[i] << endl;
        }

        if (systematics_given) {
            TGraphErrors * gstat = new TGraphErrors(npoints, y, sigma, zeros, sigma_stat_err);
            gstat -> SetLineColor(parser.getNodeContent("marker_color").Atoi());
            gstat -> Draw("||");
        }

        TGraphAsymmErrors * gtot = new TGraphAsymmErrors(npoints, y, sigma, zeros, zeros, sigma_tot_down, sigma_tot_up);
        gtot -> SetMarkerStyle(parser.getNodeContent("marker_style").Atoi());
        gtot -> SetMarkerSize(parser.getNodeContent("marker_size").Atof());
        gtot -> SetMarkerColor(parser.getNodeContent("marker_color").Atoi());
        gtot -> SetLineColor(parser.getNodeContent("marker_color").Atoi());
        gtot -> Draw("pz");

        plot.get_legend() -> AddEntry(gtot, parser.getNodeContent("legend_entry"), "p");

        // if want to use these data in the fitter
        if (parser.getNodeContent("use_in_the_fit") == "true") fitter.set_data(npoints, y, sigma, zeros, sigma_stat_err);

        parser.selectNextNode("dataplot");
    }

    // print to file
    parser.selectMainNode();
    plot.Print("../plots/" + parser.getNodeContent("filename"));

    // fit
    parser.selectMainNode();
    if ( parser.getNodeContent("perform_fit") == "true" ) {

        // run Minuit commands
        fitter.Command("MIGRAD");
        fitter.Command("HESSE");
        fitter.Command("MINOS");

        // print to eps
        plot.Print("../plots/fitted_" + parser.getNodeContent("filename"));
    }

    // done
    return 0;
}
