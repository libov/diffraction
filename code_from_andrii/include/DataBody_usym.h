#ifndef DataBody_usym_h    
#define DataBody_usym_h           
#include "DataBody.h"

 class DataBody_usym:public DataBody
  {
   public:
    DataBody_usym(){dYsist_up=dYtotal_up=dYsist_low=dYtotal_low=dYtotal=NULL;};
    ~DataBody_usym(){clean();};
    bool ReadTable(FILE* file, char* header);
    void PrintBody(const char* header=NULL); 

    float *dYsist_up,*dYtotal_up,*dYsist_low,*dYtotal_low;

   private:
    void clean() 
     {DataBody::clean(); 
      clear(dYsist_up); clear(dYtotal_up);
      clear(dYsist_low);clear(dYtotal_low);}
  };

#endif
