#ifndef DataBody_h
#define DataBody_h

//#include <ctype.h>
//#include <vector>
//#include <stdlib.h>
#include <stdio.h>
//#include "math.h"
//#include "string.h"
//#include <iostream>

class DataBody // Abstarct Data Type
{
  public:
   DataBody(){isconverted=0;N=0;P=-1111;X=Y=dYstat=dYtotal=NULL;};
  ~DataBody(){clean();};
 
   virtual void PrintBody(const char* header=NULL)=0; 
   int isconverted; 
   int N; // Number of raws.
   virtual bool ReadTable(FILE* file, char* header);
  
   float  P; 
   float *X,*Y,*dYstat,*dYtotal;

  protected:
   bool isSpace(const char* str);
   void clear(float*x){if(!x)delete[]x;x=NULL;}
   void clean();
};
#endif
