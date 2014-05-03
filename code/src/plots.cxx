#include <TXMLEngine.h>

#include <TVectorMeson.h>
#include <TXMLParser.h>
#include <constants.h>

#include <iostream>
using namespace std;

int main (int argc, char **argv) {

    TVectorMeson    instance;
    TString         filename;

    // create parser instance, pass the xml file name
    TString  xmlfilename = argv[1];
    TXMLParser parser(xmlfilename);

    // top level properties
    instance.set_cms_energy(parser.getNodeContent("cms_energy").Atof());
    TString meson_type = parser.getNodeContent("meson_type");
    if (meson_type=="jpsi") instance.set_meson_type(kJPSI);
    else if (meson_type=="psi")  instance.set_meson_type(kPSI);
    filename = parser.getNodeContent("filename");

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

    // calculates necessary parameters based on what was passed to the object since its creation
    instance.Initialise();

    // perform plotting - loop over plot tags and pass settings to the instance
    parser.selectMainNode();
    parser.selectNode("plot");
    while (parser.getCurrentNode() != 0) {

        // these are parameters we need for a theory plot
        sigma_gamma_p_model model;
        process             p;
        variable            v;
        unsigned            line_style;
        unsigned            line_color;
        TString             legend_entry;

        if      (parser.getNodeContent("sigma_gamma_p_model") == "powerlaw") model = kPOWERLAW;
        else if (parser.getNodeContent("sigma_gamma_p_model") == "reggeometry") model = kREGGEOMETRY;

        if      (parser.getNodeContent("process") == "pp") p = kPP;
        else if (parser.getNodeContent("process") == "gammap") p = kGammaP;

        if      (parser.getNodeContent("variable") == "y") v = kY;
        else if (parser.getNodeContent("variable") == "w") v = kW;

        line_style = parser.getNodeContent("line_style").Atoi();
        line_color = parser.getNodeContent("line_color").Atoi();
        legend_entry = parser.getNodeContent("legend_entry");

        instance.draw_theory(model, p, v, legend_entry, line_style, line_color);

        parser.selectNextNode("plot");
    }

    // plot data
    parser.selectMainNode();
    parser.selectNode("dataplot");
    while (parser.getCurrentNode() != 0) {
        TString filename;
        TString legend_entry;
        unsigned marker_style;

        filename = parser.getNodeContent("filename");
        marker_style = parser.getNodeContent("marker_style").Atoi();
        legend_entry = parser.getNodeContent("legend_entry");

        instance.draw_data(filename, legend_entry, marker_style);

        parser.selectNextNode("dataplot");
    }

    // print to file
    instance.Print("../plots/" + filename);

    // done
    return 0;
}
