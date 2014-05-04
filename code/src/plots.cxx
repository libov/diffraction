#include <TXMLEngine.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include <TPlot.h>
#include <TVectorMeson.h>
#include <TXMLParser.h>
#include <constants.h>

#include <iostream>
#include <fstream>
using namespace std;

int main (int argc, char **argv) {

    // create TPlot object - inherits from TCanvas
    TPlot instance;
    instance.cd();

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
    instance.set_legend_title(parser.getNodeContent("title"));
    double x1, x2, y1, y2;
    x1 = parser.getNodeContent("x1").Atof();
    x2 = parser.getNodeContent("x2").Atof();
    y1 = parser.getNodeContent("y1").Atof();
    y2 = parser.getNodeContent("y2").Atof();
    instance.set_legend_coordinates(x1, y1, x2, y2);

    // x-axis
    parser.selectMainNode();
    parser.selectNode("x_axis");
    instance.set_x_axis_title(parser.getNodeContent("title"));
    TString xlog = parser.getNodeContent("log");
    if (xlog == "false") instance.set_log_x(false);
    else if (xlog == "true") instance.set_log_x(true);
    double x_low_range = parser.getNodeContent("low_range").Atof();
    double x_up_range = parser.getNodeContent("up_range").Atof();
    instance.set_x_axis_range(x_low_range, x_up_range);

    // y-axis
    parser.selectMainNode();
    parser.selectNode("y_axis");
    instance.set_y_axis_title(parser.getNodeContent("title"));
    TString ylog = parser.getNodeContent("log");
    if (ylog == "false") instance.set_log_y(false);
    else if (ylog == "true") instance.set_log_y(true);
    double y_low_range = parser.getNodeContent("low_range").Atof();
    double y_up_range = parser.getNodeContent("up_range").Atof();
    instance.set_y_axis_range(y_low_range, y_up_range);

    // creates dummy histo, legend, applies settings
    instance.Initialise();

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

        instance.get_legend() -> AddEntry(f, parser.getNodeContent("legend_entry"), "l");

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
        g -> SetMarkerStyle(parser.getNodeContent("marker_style").Atoi());
        g -> Draw("p");

        instance.get_legend() -> AddEntry(g, parser.getNodeContent("legend_entry"), "p");

        parser.selectNextNode("dataplot");
    }

    // print to file
    parser.selectMainNode();
    instance.Print("../plots/" + parser.getNodeContent("filename"));

    // done
    return 0;
}
