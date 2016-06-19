#include <string.h> 
#include <iostream> 
using namespace std;
#include <TFile.h>                  
#include <TROOT.h>                        
#include <TStyle.h>                              
#include <TSystem.h>
#include <math.h>                                  
#include <TMath.h>                 
#ifndef __CINT__              
//#include "Math/SpecFuncMathCore.h"
//#include "Math/SpecFuncMathMore.h"
  
#include "Math/SpecFunc.h"
#endif                            
//#include <gsl/gsl_sf_hyperg.h>
//#include "Math/SpecFunc.h"           
#include <TH1F.h>
#include <TDatime.h>       
#include <TF1.h>      
#include <TF2.h>        
#include <TH2D.h>
#include <TPad.h>                                                                
#include <TCanvas.h>
#include <TLegend.h>
#include <TFitter.h>
#include <TMinuit.h>
#include <TLegend.h>
#include <TGraph2D.h>
//#include <TGraphErrors.h>                                                      
#include "DataContainer.h"
#include "TGraphErrors.h"
double Chi(DataContainer& dc,const char* fname,double low, double up);


int main()
{
 //~~ DataContainers initialization ~~~
  int block; bool sym, bins;

  //## pp ##
  const char* f7000pp="../Data/mx/ppel.TOTEM.datx";
  DataContainer DC7000pp[2];
   DC7000pp[0].Get(f7000pp, block=1, sym=true,  bins=false);
   DC7000pp[1].Get(f7000pp, block=2, sym=false, bins=false);

  const char* pp_dat="../Data/mx/ppbar.sqrtsGt100.dc.Data";
  DataContainer DC1800pp[2];
   DC1800pp[0].Get(pp_dat, block=7, sym=true,  bins=false);
   DC1800pp[1].Get(pp_dat, block=8, sym=true,  bins=false);
  DataContainer DC546pp[4];
   DC546pp[0].Get(pp_dat, block=2, sym=true,  bins=false);
   DC546pp[1].Get(pp_dat, block=3, sym=true,  bins=false);
   DC546pp[2].Get(pp_dat, block=5, sym=true,  bins=false);
   DC546pp[3].Get(pp_dat, block=4, sym=true,  bins=false);
  DataContainer DC540pp;
   DC540pp.Get(pp_dat, block=1, sym=true,  bins=false);
//   DC540pp.Print();

  //## SD ##
  const char* fdcsdtSD="../Data/mx/dcsdt.UA4.Data";
  DataContainer dcsdtSD;
   dcsdtSD.Get(fdcsdtSD, block=1, sym=true,  bins=false);
//   dcsdtSD.Print();
   printf("SD dcsdt:  "); Chi(dcsdtSD,"dcsdt.SD.log",0,0.4);


  const char* fd2csSD="../Data/mx/d2cs.SD.Data";
  DataContainer d2csSD[3];
   d2csSD[0].Get(fd2csSD, block=1, sym=true,  bins=false);
   d2csSD[1].Get(fd2csSD, block=2, sym=true,  bins=false);
   d2csSD[2].Get(fd2csSD, block=3, sym=true,  bins=false);
   
  // d2csSD[0].Print();
   printf("SD d2cs|t=0.05 |546GeV :"); Chi(d2csSD[0],"d2cs.SD.546.log",1e3,1e5);

  // d2csSD[1].Print();
   printf("SD d2cs|t=0.05 |1800GeV:"); Chi(d2csSD[1],"d2cs.SD.1800.log",1e3,1e5);

  

 //## cs(s) ##
   //~~ sd ~~
   const char* fcsSD="../Data/mx/cs.SD.Data";
   DataContainer csSD_m2, csSD_xi;
    csSD_xi.Get(fcsSD, block=1, sym=true,  bins=false);
    csSD_m2.Get(fcsSD, block=2, sym=true,  bins=false);

   // csSD_m2.Print();
    printf("SD  <200: "); Chi(csSD_m2,"SD.200.log",100,10000);
   // csSD_xi.Print();
    printf("SD  <0.05:"); Chi(csSD_xi,"SD.0.05.log",100,10000);

   //~~ dd ~~
   const char* fcsDD="../Data/mx/cs.DD.Data";
   DataContainer csDD;
    csDD.Get(fcsDD, block=2, sym=true,  bins=false);
   // csDD.Print();
   printf("DD:       "); Chi(csDD,"DD.log",100,10000);
  
   printf("\n################\n");
   printf("SD  <200: "); Chi(csSD_m2,"SD.200.Kf.log",100,10000);
   printf("SD  <0.05:"); Chi(csSD_xi,"SD.0.05.Kf.log",100,10000);
   printf("DD :      "); Chi(csDD,"DD.Kf.log",100,10000);
   
   
}

  double sq(double a){return a*a;}
  //double sq(double a){return a*a};
  
  double Chi(DataContainer& dc,const char* fname,double low, double up)
  {
    TGraphErrors dd(fname, "%lg %lg");
    
    int   N=dc.pages[0]->N;
    float* X=dc.pages[0]->X;//(sqrt(s))
    float* Y=dc.pages[0]->Y;
    float* dY=dc.pages[0]->dYtotal;
  
    int nn=0;
    double chi=0;
    for(int i=0;i<N;i++)
     if((X[i]>low)&&(X[i]<up))
      { chi+=sq((Y[i]-dd.Eval(X[i],0,"S"))/dY[i]);nn++; } 
   //   { chi+=sq((Y[i]-dd.Eval(X[i]))/dY[i]);nn++; } 
  
    printf("Chi2/Ndf:%6.3f; Ndf:%3d \n", chi/nn,nn );
   return chi/nn;
  }
  
