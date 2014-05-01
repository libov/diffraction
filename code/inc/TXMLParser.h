#ifndef TXMLPARSER_H
#define TXMLPARSER_H

// ROOT headers
#include <TROOT.h>
#include <TXMLEngine.h>

// system headers
#include <vector>
using namespace std;

class TXMLParser : public TXMLEngine {

    public:

        TXMLParser(TString filename);
        ~TXMLParser(){};

        XMLDocPointer_t   getNodeByParentAndName(XMLDocPointer_t parent, TString name);
        TString	          getNodeContentByName(TString name);
        TString	          getSubNodeContentByName(TString name);
        
    private:
        
        XMLDocPointer_t   fMainNode;
};

#endif
