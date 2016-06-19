#ifndef DataContainer_h 
#define DataContainer_h
#include "DataBody.h"
#include "Header.h"
#include <vector>
using namespace std;

//class Header;
//class DataBody;
//class vector;

class DataContainer
{
 public:
   DataContainer(){;}
   ~DataContainer();
  
   Header head;
   vector<DataBody*> pages;
   bool symetry, binning;
   void Get(const char* file, int block=1, bool sym=true, bool bins=true);
   void Print();
};
#endif 
