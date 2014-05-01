//ROOT headers
#include <TMath.h>

// own headers
#include <TXMLParser.h>

// system headers
#include <iostream>
using namespace std;

TXMLParser::TXMLParser(TString filename):
TXMLEngine()
{
    XMLDocPointer_t xmldoc = ParseFile(filename);
    if (xmldoc==0) {
        cout << "ERROR: could not open XML file" << endl;
        abort();
    }

    // take access to main node
    fMainNode = DocGetRootElement(xmldoc);
}

XMLDocPointer_t   TXMLParser::getNodeByParentAndName(XMLDocPointer_t parent, TString name) {

    // loop over children
    XMLNodePointer_t child = xml -> GetChild(parent);
    while (child!=0) {

        TString child_name    = xml -> GetNodeName(child);

	if (child_name)

        if ( child_name == "cms_energy" ) instance.set_cms_energy(child_content.Atof());
        if ( child_name == "meson_type" ) {
            if (child_content=="jpsi") instance.set_meson_type(kJPSI);
            if (child_content=="psi")  instance.set_meson_type(kPSI);
        }
        if ( child_name == "filename" ) filename = child_content;
}

