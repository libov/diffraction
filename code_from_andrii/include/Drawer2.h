#ifndef Drawer2_h
#define Drawer2_h
#include "TGraphErrors.h"
#include "TCanvas.h"
//class TGraphErrors;
//class TCanvas;
//class TF1;
class DataContainer;

class Drawer2
{
 public:
  Drawer2() {gr=NULL; cv=NULL; X=Y=dX=dY=NULL;}
  ~Drawer2(){if(!gr)delete[]gr;if(!cv)delete[]cv;if(!dX)delete[]dX;}
  void SetDataContainer(DataContainer* _dc) {dc=_dc;Npages=dc->pages.size();}
  void SetcvNmAdd(const char* CVaddNm=""){cvaddNm=(char*)CVaddNm;};
  void Draw(const char* file);
  void Draw2(const char* file,int w,int h);
  void Draw3(const char* file);

 private:
  void GetPage(int i=0);
  DataContainer* dc;
  
  int Npages;
  TGraphErrors* gr;
  TCanvas* cv;
  char* cvaddNm;

 private:
  int Nlines_;
  float *X,*Y,*dX,*dY; 
  
};
#endif
