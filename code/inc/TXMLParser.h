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

        void                selectMainNode();
        void                selectNode(TString node_name);
        void                selectNextNode(TString node_name);
        XMLDocPointer_t     getNode(TString node_name);
        TString             getNodeContent(TString node_name);

        XMLDocPointer_t     getCurrentNode() {return fCurrentNode;}

    private:

        XMLDocPointer_t     fMainNode;
        XMLDocPointer_t     fCurrentNode;
};

#endif
