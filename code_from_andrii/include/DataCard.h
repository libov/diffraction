#ifndef DataCard_h
#define DataCard_h

#include "DataBody_sym.h"     
#include "DataContainer.h"

typedef double (*Ftype)(double x,  double* par);
typedef double (*FtypeI)(double s, double q2, double Mv2, double* par);
typedef double (*FtypeII)(double t,double s, double q2, double Mv2, double* par);
class DataCard
{public:
  DataCard(){DataType=WQ2pos=0;metric=1;}
  DataCard(DataCard& rhs)
   {Set(rhs.file, rhs.block, rhs.sym, rhs.bins, rhs.W, rhs.Q2, rhs.Mv2);metric=1;}
  DataCard(const char* f, int bk, bool sm, bool bn, double w, double q2, double mv2)
   {Set(f, bk, sm, bn, w, q2, mv2);metric=1;}
  ~DataCard(){;}

   void Set(const char* f, int bk, bool sm, bool bn, double w, double q2, double mv2);
   void SetMetric(double m){metric=m;};
   void SetFunc(Ftype   ff){func0 =ff;};
   void SetFunc(FtypeI  ff){funcI =ff;};
   void SetFunc(FtypeII ff){funcII=ff;};
   double Chi2(double* par);
   double Chi2(double* par,int iBlock);//Doesn't clean Xi2 and iii

   double Xi,Xi2,iii;
   DataContainer dc;
   double W,Q2,Mv2;
  private:
   const char* file;
   int block;
   bool sym, bins;
   double metric;

   Ftype   func0;
   FtypeI  funcI;
   FtypeII funcII;
  
  private:
   float *_X,*_Y,*_dY,_P;//Data
   int DataType, WQ2pos;
   int i,j;
   int Npages,Nlines;  
}; 
#endif
