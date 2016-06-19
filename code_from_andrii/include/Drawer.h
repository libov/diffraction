#ifndef Drawer_h
#define Drawer_h
typedef double (*Ftype)(double* t,double* par);
//#include "math.h"
//#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
//#include "TFitter.h"
//#include "TFile.h"
#include "TCanvas.h"
//#include "DataBody_sym.h"
#include "DataContainer.h"
//class TGraphErrors;
//class TF1;
//class DataContainer;

class Drawer
{
 public:
  Drawer() {gr=NULL; f1=NULL; cv=NULL; par=vpar=NULL;X=Y=dX=dY=NULL; parNm=vparNm=NULL; Npar=vNpar=0; PNum=-1;cvaddNm=(char*)"";fx1=0;fx2=2;x1=x2=0;y1=y2=0; xLog=yLog=false;}
  ~Drawer(){if(!gr)delete[]gr;if(!cv)delete[]cv;if(!f1)delete[]f1;if(!dX)delete[]dX;}

  void SetDataContainer(DataContainer* _dc) {dc=_dc;Npages=dc->pages.size();}
  void SetFunc(Ftype func_){func=func_;}
  void SetParameters(double* Par, const char** ParNm, int NPar)      {par=Par; parNm=ParNm; Npar=NPar;}
  void SetVarParameters(double* vPar, const char** vParNm, int vNPar){vpar=vPar; vparNm=vParNm; vNpar=vNPar;}
  void SetPnum(int Pnum){PNum=Pnum;}

  void SetcvNmAdd(const char* CVaddNm=""){cvaddNm=(char*)CVaddNm;};
  void SetFuncRange(double _fx1, double _fx2){fx1=_fx1;fx2=_fx2;};
  void SetAxisLog(bool xlog, bool ylog){xLog=xlog; yLog=ylog;};
  void SetPlotXRange(double _x1, double _x2){x1=_x1;x2=_x2;};
  void SetPlotYRange(double _y1, double _y2){y1=_y1;y2=_y2;};
  void Draw(const char* file,int ii);
  void Draw1(const char* file, int block=-1);
  void Draw2(const char* file,int w,int h, int bk=-1);
  void Draw3(const char* file, int l);
  void DrawDisp(const char* file, int ii);

 protected:
  void SetGraphs(int l);
  void GetPage(int i=0);
  DataContainer* dc;
  
  int PNum;// P->vpar[Pnum] (!) variable parameter that are taken from table page.
  double* par, *vpar;
  const char**  parNm, **vparNm;
  int Npar, vNpar;

  Ftype func;
  double fx1,fx2;
  double x1,x2;
  double y1,y2;
  bool xLog, yLog;
  int Npages;
  TCanvas* cv;
  TGraphErrors* gr;
  TF1* f1;
  char* cvaddNm;

 protected:
  int Nlines_;
  float *X,*Y,*dX,*dY; 
  
};
#endif
