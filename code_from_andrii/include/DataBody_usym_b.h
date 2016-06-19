#ifndef DataBody_usym_b_h    
#define DataBody_usym_b_h           

//#include <vector>
//#include <stdlib.h>
//#include <stdio.h>
//#include "math.h"
//#include "string.h"
//#include <iostream>

#include "DataBody.h"

 class DataBody_usym_b:public DataBody
  {
   public:
    DataBody_usym_b(){X1=X2=dYsist_up=dYtotal_up=dYsist_low=dYtotal_low=dYtotal=NULL;};
    ~DataBody_usym_b(){clean();};
    bool ReadTable(FILE* file, char* header);
    void PrintBody(const char* header=NULL); 

    float *X1,*X2;
    float *dYsist_up,*dYtotal_up,*dYsist_low,*dYtotal_low;

   private:
    void clean() 
     {DataBody::clean();
      clear(X1);clear(X2);
      clear(dYsist_up); clear(dYtotal_up);
      clear(dYsist_low);clear(dYtotal_low);}
  };

#endif
