#ifndef DataBody_sym_b_h    
#define DataBody_sym_b_h           

#include "DataBody.h"
 class DataBody_sym_b:public DataBody
  {
   public:
    DataBody_sym_b(){X1=X2=dYsist=dYtotal=NULL;};
    ~DataBody_sym_b(){clean();};
    bool ReadTable(FILE* file, char* header);
    void PrintBody(const char* header=NULL); 

    float *X1,*X2;
    float *dYsist;

   private:
    void clean() 
     {DataBody::clean(); 
      clear(X1);clear(X2);
      clear(dYsist);}
  };
#endif
