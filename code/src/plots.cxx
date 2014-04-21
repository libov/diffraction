#include <TXMLEngine.h>

#include <TVectorMeson.h>
#include <constants.h>

#include <iostream>
using namespace std;

int main (int argc, char **argv) {

    TVectorMeson    instance;
    TString         filename;
    
    // open xml file
    TXMLEngine * xml = new TXMLEngine;
    TString  xmlfilename = argv[1];
    XMLDocPointer_t xmldoc = xml->ParseFile(xmlfilename);
    if (xmldoc==0) {
        delete xml;
        cout << "ERROR: could not open XML file" << endl;
        abort();
    }
    
    // various structs to hold compound objects
    struct axis {
        double  low_range;
        double  up_range;
        TString title;
        bool    log;
    } x, y;

    // take access to main node   
    XMLNodePointer_t mainnode = xml -> DocGetRootElement(xmldoc);
    
    // loop over children
    XMLNodePointer_t child = xml -> GetChild(mainnode);
    while (child!=0) {

        TString child_name    = xml -> GetNodeName(child);
        TString child_content = xml -> GetNodeContent(child);
        
        if ( child_name == "cms_energy" ) instance.set_cms_energy(child_content.Atof());
        if ( child_name == "meson_type" ) {
            if (child_content=="jpsi") instance.set_meson_type(kJPSI);
            if (child_content=="psi")  instance.set_meson_type(kPSI);
        }
        if ( child_name == "filename" ) filename = child_content;

        if ( child_name == "legend" ) {

            XMLNodePointer_t child2 = xml -> GetChild(child);
            double x1, x2, y1, y2;
            while (child2!=0) {
 
                TString child2_name = xml -> GetNodeName(child2);
                TString child2_content = xml -> GetNodeContent(child2);

                if (child2_name=="title") instance.set_legend_title(child2_content);
                if (child2_name=="x1") x1 = child2_content.Atof();
                if (child2_name=="y1") y1 = child2_content.Atof();
                if (child2_name=="x2") x2 = child2_content.Atof();
                if (child2_name=="y2") y2 = child2_content.Atof();
                child2 = xml->GetNext(child2);
            }
            
            instance.set_legend_coordinates(x1, y1, x2, y2);
        } 

        if ( child_name == "x_axis" ) {
            XMLNodePointer_t child2 = xml -> GetChild(child);
            while (child2!=0) {
 
                TString child2_name = xml -> GetNodeName(child2);
                TString child2_content = xml -> GetNodeContent(child2);

                if (child2_name=="title") x.title = child2_content;
                if (child2_name=="log") {
                    if (child2_content == "false") x.log = false;
                    if (child2_content == "true") x.log = true;
                }
                if (child2_name=="low_range") x.low_range = child2_content.Atof();
                if (child2_name=="up_range")  x.up_range  = child2_content.Atof();
                child2 = xml->GetNext(child2);
            }
        }
        
        if ( child_name == "y_axis" ) {
            XMLNodePointer_t child2 = xml -> GetChild(child);
            while (child2!=0) {
 
                TString child2_name = xml -> GetNodeName(child2);
                TString child2_content = xml -> GetNodeContent(child2);

                if (child2_name=="title") y.title = child2_content;
                if (child2_name=="log") {
                    if (child2_content == "false") y.log = false;
                    if (child2_content == "true") y.log = true;
                }
                if (child2_name=="low_range") y.low_range = child2_content.Atof();
                if (child2_name=="up_range")  y.up_range  = child2_content.Atof();
                child2 = xml->GetNext(child2);
            }
        }
        
        child = xml->GetNext(child);
    }

    instance.set_x_axis_title(x.title);
    instance.set_y_axis_title(y.title);
    instance.set_x_axis_range(x.low_range, x.up_range);
    instance.set_y_axis_range(y.low_range, y.up_range);
    
    instance.set_log_x(x.log);
    instance.set_log_y(y.log);
    
    // calculates necessary parameters based on what was passed to the object since its creation
    instance.Initialise();

    // loop over children
    child = xml -> GetChild(mainnode);
    while (child!=0) {

        TString child_name    = xml -> GetNodeName(child);
        TString child_content = xml -> GetNodeContent(child);
    
        if ( child_name == "plot" ) {
            XMLNodePointer_t child2 = xml -> GetChild(child);
            
            sigma_gamma_p_model model;
            process             p;
            variable            v;
            unsigned            line_style;
            unsigned            line_color;
            TString             legend_entry;
            
            while (child2!=0) {
 
                TString child2_name = xml -> GetNodeName(child2);
                TString child2_content = xml -> GetNodeContent(child2);

                if (child2_name == "sigma_gamma_p_model") {
                    if ( child2_content == "powerlaw" )     model = kPOWERLAW;
                    if ( child2_content == "reggeometry" )  model = kREGGEOMETRY;
                }
                if (child2_name == "process") {
                    if ( child2_content == "pp" )      p = kPP;
                    if ( child2_content == "gammap" )  p = kGammaP;
                }
                if (child2_name == "variable") {
                    if ( child2_content == "y" )  v = kY;
                    if ( child2_content == "w" )  v = kW;
                }
                
                if (child2_name == "line_style") line_style = child2_content.Atoi();
                if (child2_name == "line_color") line_color = child2_content.Atoi();
                if (child2_name == "legend_entry") legend_entry = child2_content;

                child2 = xml->GetNext(child2);
            }
            
            instance.draw_theory(model, p, v, legend_entry, line_style, line_color);
        }
        child = xml->GetNext(child);
    }

    // loop over children
    child = xml -> GetChild(mainnode);
    while (child!=0) {

        TString child_name    = xml -> GetNodeName(child);
        TString child_content = xml -> GetNodeContent(child);
    
        if ( child_name == "dataplot" ) {
            XMLNodePointer_t child2 = xml -> GetChild(child);
            
            TString filename;
            TString legend_entry;
            
            while (child2!=0) {
 
                TString child2_name = xml -> GetNodeName(child2);
                TString child2_content = xml -> GetNodeContent(child2);

                if (child2_name == "filename") filename = child2_content;
                if (child2_name == "legend_entry") legend_entry = child2_content;

                child2 = xml->GetNext(child2);
            }
            
            instance.draw_data(filename, legend_entry);
        }
        child = xml->GetNext(child);
    }
    
    instance.Print("../plots/" + filename);
    
    // Release memory before exit
    xml->FreeDoc(xmldoc);
    delete xml;
   
    return 0;
}