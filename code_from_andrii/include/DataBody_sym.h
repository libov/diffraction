#ifndef DataBody_sym_h    
#define DataBody_sym_h           

#include "DataBody.h"
 class DataBody_sym:public DataBody
  {
   public:
    DataBody_sym(){dYsist=NULL;};
    ~DataBody_sym(){clean();};
    bool ReadTable(FILE* file, char* header);
    void PrintBody(const char* header=NULL); 

    float *dYsist;

   private:
    void clean() 
      {DataBody::clean(); 
       clear(dYsist);}
  };
#endif
