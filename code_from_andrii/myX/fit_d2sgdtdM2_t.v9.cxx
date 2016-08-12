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

#include "plot_style_utils.h"

double Theta_function(double s1, double s2);
double Re_alpha(double s);
double Im_alpha(double s);
inline double sq(double a){return a*a;};
//**********************************
const double m_proton = 0.938272013;
const double mpi= 0.13498;
const double mp = m_proton;
const double m_proton2 = m_proton*m_proton;
const double mp2 = m_proton*m_proton;
const double MpMpi = m_proton + mpi;
const double sqMpMpi = MpMpi*MpMpi; 
const double Mp2Mpi2 = sq(mp)+sq(mpi); 

const double MRoper = 1.44;
const double MRoper2 = MRoper*MRoper;
const double WRoper = 0.325;

const double Pi = TMath::Pi();

const double alfa0 =-0.410839, delta = -0.460021;
const double c1_par[3] = {0.508309, 4.01083, 4558.44};
const double s_par[3]  = {sqMpMpi, 2.44, 11.7};
const double lambda[3] = {0.840083, 2.09707, 11.1778};

//const double _alpha =  0.25;

 bool _PNG_=false;//print also in png format
// bool _PNG_=true;//print also in png format

//(*Kf*)
//double alpha0=0.075, alpha_x=0.34,alpha2=0;//0.08;
//double dA=exp(-2.*alpha0*log(100));
//double db=-1/2.*0.44*log(100);
//double  k = 1.03377271988205033e+00;
// double _eta = 1; 
// double Ares = 15.5*(2*mp);//k*8*dA;
// double  Cbg = k*4*dA;
// double  nR0 = 0.45;//0.025//0.01
// double b_res= -0.41;//0.+db/2.; //f(x)
// double b_bg = -0.91;//0.+db;//3 //b goulianos
//
//double alpha0=0.075, alpha_x=0.34,alpha2=0;//0.08;
//double dA=exp(-2.*alpha0*log(100));
//double db=0;//-1/2.*0.44*log(100);
//double  k = 1.03377271988205033e+00;
// double _eta = 1; 
// double Ares = k*8*dA;
// double  Cbg = k*4*dA;
// double  nR0 = 0.45;//0.025//0.01
// double b_res= 0.+db/2.; //f(x)
// double b_bg = 0.+db;//3 //b goulianos
//
//const double kx = 1;//0.5;
//const double xi = 0.8;
//double Norm_el = 67 *dA;
//double Norm_DD = 1./4./Norm_el, s0=1.;//00;
//double bel  = 2.95+db;

double alpha0=0.075, alpha_x=0.34,alpha2=0;//0.08;
double dA=exp(-2.*alpha0*log(100));
double db=0;//-1/2.*0.44*log(100);
double  k = 1.03377271988205033e+00;
 double _eta = 1; 
 double Ares = k*8*dA;
 double  Cbg = k*4*dA;
 double  nR0 = 0.45;//0.025//0.01
 double b_res= -0.507+db/2.; //f(x)
 double b_bg = -1.013+db;//3 //b goulianos

const double kx = 1;//0.5;
const double xi = 0.8;
double Norm_el = 67 *dA;
double Norm_DD = 1./4./Norm_el, s0=1.;//00;
double bel  = 1.937+db;



const double M2_min =1, M2_max=5e8;
//const double M2_min =1, M2_max=12;
const double t_min = 0., t_max=1.;

int i;
int it,im2;
const double peak_M2_min = 2, peak_M2_max=3.1;


// August 2016 Jenkovszky, Libov
bool nonlinear_pomeron_trajectory = true;

void prestage_tlin(double*t, double* dtarr, int Nt, double tL, double tU)
 {double dt=(tU-tL)/Nt;
  for(int it=0;it<Nt;it++)
    {t[it]=-(tL + dt*(it+0.5)); dtarr[it]=dt;}
 }

void prestage_M2(double*M2, double* dm2arr, int Nm2, double M2L, double M2U, const char* opt="lin")
 {string a="";
  if((a+"lin")==opt)
   {
    double dM2=(M2U-M2L)/Nm2;
    for(int im2=0;im2<Nm2;im2++)
     {M2[im2]=M2L+dM2*(im2+0.5); dm2arr[im2]=dM2;}
   }
  if((a+"log")==opt)
   {
    const double x_min=log(M2L), x_max=log(M2U), dx=(x_max-x_min)/Nm2;
    for(int im2=0;im2<Nm2;im2++)
      {M2[im2]     = exp(x_min+dx*(im2+0.5));
       dm2arr[im2] = exp(x_min+dx*(im2+1))-exp(x_min+dx*im2);}
   }
 }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double alpha_t(double t) //alpha_t= alpha(t)-1  (!!!)
{
  // August 2016 Jenkovszky, Libov
  double alpha1=0;
  if (nonlinear_pomeron_trajectory) alpha1 = alpha_x / 10;
  return ( alpha0 + alpha_x*t - alpha1 * sqrt(4 * sq(mpi) - t) );
//return 0.08+0.25*t;
//return alpha0+alpha_x*t+alpha2*t*t;
//return 0.087+0.44*t+0.078*t*t;//** 
//return 0.09+0.44*t; //+
//return 0.09+0.44*t+0.078*t*t; //*
//return 0.05+0.4*t+0.078*t*t;
//return 0.08+0.25*t+0.078*t*t;
//return 0.08+0.25*t+0.25*t*t;
}

void stageI_t(double*Fp2_t,double*f4_res,double*f_bg, double* alpha_IP, double*t, int Nt)
{
//  printf("t[it],alpha_IP[it],f4_res[it],Fp2_t[it],f_bg[it]\n");
  double _t;
  double f2_tmp,fx; 

  for(it=0;it<Nt;it++)
   {
    _t=t[it];
    alpha_IP[it] = 2*alpha_t(_t);
    //fx =(4*mp2-2.8*_t)/(4*mp2-_t)*pow(1.-_t/0.71,-2.);fx*=fx;
//    fx =pow(1.-_t/0.71,-2.);fx*=fx;
//    fx*=exp(-2*2.1*_t);
//      fx=exp(1*_t*_t);//exp(0.5*_t+1*_t*_t);
    fx=1;
//    Fp2_t[it] = pow(1.-_t/0.71,-4.);
//     Fp2_t[it] =sq((4*mp2-2.8*_t)/(4*mp2-_t)*pow(1.-_t/0.71,-2.));
//    Fp2_t[it]   = exp(2*2.81*_t);
 //   Fp2_t[it] = exp(2*2.1*_t);
    Fp2_t[it] = exp(2*bel*_t);
    f_bg[it]  = exp(b_bg*t[it])*fx;

    f2_tmp = pow(exp(b_res*t[it])*fx, 2);
//    f4_res[it] = pow(f2_tmp,2);
//    f4_res[it+Nt] = pow(f2_tmp,3);
//    f4_res[it+2*Nt] = pow(f2_tmp,4);

    f4_res[it]      = pow(f2_tmp,2);
    f4_res[it+Nt]   = pow(f2_tmp,2);
    f4_res[it+2*Nt] = pow(f2_tmp,2);
  //  printf("%7.3f| %7.3f %7.3f %7.3f %7.3f \n",t[it],alpha_IP[it],f4_res[it],Fp2_t[it],f_bg[it]);
   }
}
 
void stageI_M2(double*ImA,double*ImBg,double*ImR, double*M2,int Nm2_, double eta=_eta)
{
  double im_tmp,re_tmp; int j;
  for(im2=0;im2<Nm2_;im2++)
   {im_tmp=Im_alpha(M2[im2]); re_tmp=Re_alpha(M2[im2]);
    for(j=1;j<4;j++)
      {ImA[(j-1)*Nm2_+im2] = im_tmp/( sq(2*j+0.5-re_tmp) + sq(im_tmp));}
       ImR[im2] = MRoper*(WRoper/2.)/( sq(M2[im2]-MRoper2) + sq(WRoper/2) );
    if(M2[im2]>sqMpMpi) 
        ImBg[im2] = 1./(1/kx/(pow(M2[im2]-sqMpMpi,xi))
                       +pow(M2[im2], eta));
     else ImBg[im2] = 0;
   }
}

inline double sigma_tot_pPN_m (int _it,int _im2, double* inF,double* imA, int Nt,int Nm2);
double kinematic_factor(double t, double M2);
void SetErrGraph(TGraphErrors& g, DataContainer* dc);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Plot_IIIel(const int Nt, double* t, double*dcs1, double*dcs2, double*dcs3,TString* leg, const char* fout);
void Plot_IIIx_el(const int Nt, double* t, double*dcs1,TString sleg, const char* fout);
void Plot_IIIxx_el(const int Nt, double* t, double*dcs7000,double*dcs1800, double*dcs546, double*dcs540,TString sleg, const char* fout);
void Plot_Vel(const int Nt, double* t,double* B1, double* B2, double* B3,TString* s, const char* fout);

//const double  norm_fr2 = 2.89623974996086773e+01;//28.3192954448929584;
//const double  norm_fr  = 2.89623974996086773e+01;//28.3192954448929584;
//const double  norm_r21 = 2.57839661588770674e+01;//25.2023042451354726;
//const double  norm_r28 = 3.00292458437867502e+01;//29.3636369693752570;

//const double Norm_el=norm_fr*1.63484936024432215;
//const double Norm_el=4.73491570233770034e+01;
void stageII_el(double*t, double*dtarr, int Nt )
{
 double dcsdt_el7000[Nt], dcsdt_el1800[Nt], dcsdt_el546[Nt], dcsdt_el540[Nt]; 
 double cs_el7000 = 0,  cs_el1800 = 0,  cs_el546 = 0; 
 double B_el7000[Nt-1], B_el1800[Nt-1], B_el546[Nt-1];

 double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
 stageI_t(Fp2_t, f4_res, f_bg,alpha_IP, t, Nt);
// TString str0="#frac{d#sigma_{el}}{dt} theoretical";
 TString str0="theoretical";

 for(it=0;it<Nt;it++)
  {
   dcsdt_el7000[it] = Norm_el*sq(Fp2_t[it])*pow(49e6     ,alpha_IP[it]);// /100
   dcsdt_el1800[it] = Norm_el*sq(Fp2_t[it])*pow(1800*1800,alpha_IP[it]);// /100
   dcsdt_el546[it]  = Norm_el*sq(Fp2_t[it])*pow(546*546  ,alpha_IP[it]);// /100 
   dcsdt_el540[it]  = Norm_el*sq(Fp2_t[it])*pow(540*540  ,alpha_IP[it]);// /100 
   cs_el7000 += dcsdt_el7000[it]*dtarr[it];  
   cs_el1800 += dcsdt_el1800[it]*dtarr[it];
   cs_el546  += dcsdt_el546 [it]*dtarr[it];
  }
 Plot_IIIx_el(Nt,t,dcsdt_el7000,str0, "dcsdt_xel.eps");
 printf("Hallo\n");
 Plot_IIIxx_el(Nt,t,dcsdt_el7000,dcsdt_el1800,dcsdt_el546,dcsdt_el540,str0, "dcsdt_xxel.eps");
// //### test el ###
// printf("%6s %6s (I)\n","t","dcsdt");
// for(it=0;it<Nt;it++) {printf("%6.3f %6.3f\n", t[it],dcsdt_el1800[it]);}
// //###############
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TString strx[4]={ "elastic pp",
   "#sqrt{s} = 7000GeV", 
   "#sqrt{s} = 1800GeV", 
   "#sqrt{s} =  546GeV"};
 Plot_IIIel (Nt,t,dcsdt_el7000,dcsdt_el1800,dcsdt_el546,strx, "dcsdt_x3el.eps");

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 for(it=0;it<Nt-1;it++)
  {B_el7000[it] =(log(dcsdt_el7000[it+1])-log(dcsdt_el7000[it]))/(t[it+1]-t[it]);
   B_el1800[it] =(log(dcsdt_el1800[it+1])-log(dcsdt_el1800[it]))/(t[it+1]-t[it]);
   B_el546 [it] =(log(dcsdt_el546 [it+1])-log(dcsdt_el546 [it]))/(t[it+1]-t[it]); }
 Plot_Vel(Nt,t,B_el7000, B_el1800, B_el546,strx, "Bel_t.ins.eps");

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 it=Nt/6;
 printf("### elastic pp ###\n");
 printf("cs_el7000:%6.3f cs(25.4+/-1.0mb); t:%6.4f; B(7000):%6.3f;B(19.9+/- 0.26)\n", cs_el7000,t[it],B_el7000[it]);
 printf("cs_el1800:%6.3f cs(16.6+/-1.6mb); t:%6.4f; B(1800):%6.3f;B(17,     17.9)\n", cs_el1800,t[it],B_el1800[it]);
 printf("cs_el546 :%6.3f cs(13.6+/-0.6mb); t:%6.4f; B(546 ):%6.3f;B(15,    15.35)\n", cs_el546 ,t[it],B_el546 [it]);
 printf("stageII_el  [passed]\n\n");
}

//#############################################
TString strLeg[4];
void Plot_IGoulFlat(const char* fout, const int Nm2, double* M2, double* d2A, double* d2B, double x1, double x2, double y1,double y2, bool xLog, bool yLog);
void Plot_IGoulian(const int Nm2, double* M2, double* d2A, double* d2B, double* d2C);
void Plot_IDDex7(const int Nm2, double* M2, double* d2A,double* d2B,double* d2C);
void Plot_IAdditionalData(const int Nm2, double* M2, double* d2A, double* d2B);
//#############################################
void stage_Goul(double*M2,const int Nm2, double* ImA,double* ImBg,double* ImR)
{
 const int Nt=2;
 double _t[2];
 double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
 _t[0]=-0.05; _t[1]=-0.55;
 stageI_t(Fp2_t, f4_res, f_bg,alpha_IP, _t, 2);
  
 double Kf;
 double Sg_b02, Bg_b02;
 double Propag_al20_s20,  Propag_al20_s1800,  Propag_al20_s546;
 double gd2cs_dtdM2_A[2*Nm2], gd2cs_dtdM2_B[2*Nm2], gd2cs_dtdM2_C[2*Nm2];   
 double tmp;
 strLeg[0]="SD, B-slope ";
 strLeg[1]="#sqrt{s} = 1.8TeV";
 strLeg[2]="#sqrt{s} = 546GeV";
 strLeg[3]="#sqrt{s} = 20 GeV";
  
 double _s1 =1800*1800;
 double _s2 = 546*546;
 double _s3 =  20*20;

 for(it=0;it<Nt;it++)
  for(im2=0;im2<Nm2;im2++)
   {
    Propag_al20_s1800 = pow(_s1/M2[im2], alpha_IP[it]);// /100
    Propag_al20_s546  = pow(_s2/M2[im2], alpha_IP[it]);// /100
    Propag_al20_s20   = pow(_s3/M2[im2], alpha_IP[it]);// /100
    if(Propag_al20_s1800<0)Propag_al20_s1800 = 0;
    if(Propag_al20_s546 <0)Propag_al20_s546  = 0;
    if(Propag_al20_s20  <0)Propag_al20_s20   = 0;

    Kf     = kinematic_factor(_t[it], M2[im2]);
    Sg_b02 = sigma_tot_pPN_m(it,im2, f4_res, ImA,Nt,Nm2)
                +nR0*sqrt(f4_res[it])*ImR[im2];
    Bg_b02 = f_bg[it]*ImBg[im2];

    tmp = Fp2_t[it]*(Ares*Sg_b02*Kf/(2*mp) + Cbg*Bg_b02);
//    if(it==0)
    tmp*=0.5;//single arm
    gd2cs_dtdM2_A[it*Nm2+im2]= tmp *Propag_al20_s1800;
    gd2cs_dtdM2_B[it*Nm2+im2]= tmp *Propag_al20_s546 ;
    gd2cs_dtdM2_C[it*Nm2+im2]= tmp *Propag_al20_s20  ;
   }
// //### test SD Goul ###
// printf("%6s %6s %6s (SD_Goul:I)\n","t","M2","dcsdt");
// it=0;
// for(im2=0;im2<50;im2++)
//  {printf("%6.3f %6.3f %6.3f\n", _t[it],M2[im2],gd2cs_dtdM2_A[it*Nm2+im2]);}
//######################################
 Plot_IGoulian(Nm2,M2, gd2cs_dtdM2_A,    gd2cs_dtdM2_B,     gd2cs_dtdM2_C);
 Plot_IDDex7  (Nm2,M2, gd2cs_dtdM2_A+Nm2, gd2cs_dtdM2_B+Nm2, gd2cs_dtdM2_C+Nm2);
 Plot_IAdditionalData(Nm2,M2, gd2cs_dtdM2_A,     gd2cs_dtdM2_B);
//######################################

 for(im2=0;im2<Nm2;im2++)
  {
    gd2cs_dtdM2_A[im2]*=pow(M2[im2],1.05);
    gd2cs_dtdM2_B[im2]*=pow(M2[im2],1.05);
    gd2cs_dtdM2_C[im2]*=pow(M2[im2],1.05);
  }
 //  for(int i=0;i<30;i++)printf("M2:%6.3f, cs:%6.3f\n",M2[i],gd2cs_dtdM2_A[i]);
   // Plot_IGoulFlat("FlatGoul.png",Nm2,M2, gd2cs_dtdM2_A, gd2cs_dtdM2_B, 1,1e5,  1e-1,3,true,false);
 Plot_IGoulFlat("FlatGoul.eps",Nm2,M2, gd2cs_dtdM2_A, gd2cs_dtdM2_B, 0,12,  0,14,false,false);
}

void stage_IIcernel(double* d2cs, double s, double*t,double*M2,
                    int Nt,int Nm2,
                    double*Fp2_t,double*f4_res,double*f_bg,double*alpha_IP,
                    double*ImA,double*ImBg,double*ImR,
                    double nA=Ares, double nCC=Cbg)
{
 // printf("%6s %6s %6s %6s %6s %6s\n", "t","M2", "Kf","Sg_b02","Bg_b02","d2cs");
 double Kf, propag, Sg_b02, Bg_b02;
 for(it=0;it<Nt;it++)
  for(im2=0;im2<Nm2;im2++)
   {
    propag = pow(s/M2[im2], alpha_IP[it]);// /100
    if(propag<0)printf("AgoV!");

      Kf     = kinematic_factor(t[it], M2[im2]);
      Sg_b02 = sigma_tot_pPN_m(it,im2, f4_res, ImA,Nt,Nm2)
                +nR0*sqrt(f4_res[it])*ImR[im2];
      Bg_b02 = f_bg[it]*ImBg[im2];
      {d2cs[it*Nm2+im2]=Fp2_t[it]*(nA*Sg_b02*Kf/(2*mp)+nCC*Bg_b02)*propag;}
   // else {d2cs[it*Nm2+im2]=0;}
//   printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", t[it],M2[im2], Kf,Sg_b02,Bg_b02, d2cs[it*Nm2+im2]/2.);
   }
}

//~~~~~ Technical ~~~~~~~~~~~~~~~~~~
double Integ2D(double* d2s,int Nt,int Nm2,double*M2, double *tarr,double *m2arr, double s_downlim=1, double s_uplim=49e6)
 {
   int jt,jm2;
//   printf("Here\n");
   double cs=0;
   int id=0;
   for(jt=0;jt<Nt;jt++) 
    for(jm2=0;(jm2<Nm2);jm2++) 
     {
      if((M2[jm2]>s_downlim)&&(M2[jm2]<s_uplim))
       {cs+=d2s[id]*tarr[jt]*m2arr[jm2];} 

      id++;
     }
   return cs;
 }
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Plot_VIcs_sII(int _Ns, double *_s_, double* cs1, double* cs2, string* leg, double _s_min, double _s_max, double ymin, double ymax, const char* label);
void Plot_VIcs_s(int _Ns, double *_s_, double* cs, double _s_min, double _s_max, double ymin, double ymax, const char* label);

void stage_IIs_dependense()
{
 const int Ns=20;
 double sqrts_min=10, sqrts_max=14000;
 double sarr[Ns]={0},sqarr[Ns]={0};

 const int Nt=200, Nm2=1000;//6000;
// const int Nt=200, Nm2=3000;//6000;
 double M2[Nm2], t[Nt];
 double dtarr[Nt]={0}, dm2arr[Nm2]={0};

 double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
 double  ImA[3*Nm2], ImBg[Nm2], ImR[Nm2];
//~~~~~~~~~~~~~~~~~
 prestage_tlin(t, dtarr, Nt, 0,1);
 prestage_M2(M2, dm2arr, Nm2, 1,1e8,"log");
 stageI_M2(ImA,ImBg,ImR,M2,Nm2);
 stageI_t(Fp2_t, f4_res, f_bg,alpha_IP, t, Nt);

 //~~~ s-array ~~~~
 const double x_min=log(sqrts_min), x_max=log(sqrts_max), dx=(x_max-x_min)/(Ns-1);
 for(int i=0;i<Ns;i++) {sqarr[i]=exp((x_min+dx*i));sarr[i]=sqarr[i]*sqarr[i];}
 //~~~ calculate d2cs matix ~~~~
 double _s=49e6;
 double cd2cs_dtdM2[Nt*Nm2], tmp_d2cs_dtdM2[Nt*Nm2];
 stage_IIcernel(cd2cs_dtdM2,_s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP,ImA,ImBg,ImR, Ares,Cbg);

 //~~~ propagate s-dependense + integrate ~~~
 int ind, im2;
 double propag;
 double  cs_s[Ns]={0}, cs_s200[Ns]={0}, cs_s005[Ns]={0};
 for(int is=0;is<Ns;is++)
  {
   ind=0;
   for(it=0;it<Nt;it++)
    {propag=pow(sarr[is]/_s, alpha_IP[it]);
     for(im2=0;im2<Nm2;im2++)
      {tmp_d2cs_dtdM2[ind]=propag*cd2cs_dtdM2[ind]; ind++;}
    }
// stage_IIcernel(cd2cs_dtdM2,sarr[is],t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);

   cs_s005[is]=Integ2D(tmp_d2cs_dtdM2,Nt,Nm2,M2,dtarr,dm2arr,1,0.05*sarr[is]);
 //  cs_s005[is]=Integ2D(tmp_d2cs_dtdM2,Nt,Nm2,M2,dtarr,dm2arr,1,(0.05*sarr[is]>200*200)?(0.05*sarr[is]):200*200);
   cs_s200[is]=Integ2D(tmp_d2cs_dtdM2,Nt,Nm2,M2,dtarr,dm2arr,1,200*200);
  }
  Plot_VIcs_s(Ns, sqarr, cs_s,   sqrts_min, sqrts_max,0,25,"SD");
  string leg[2]={"#xi < 0.05","M < 200GeV"};
  Plot_VIcs_sII(Ns, sqarr, cs_s005,cs_s200,leg, sqrts_min, sqrts_max,0,25,"SD");
//  Plot_VIcs_s(Ns, sqarr, cs_s005, sqrts_min, sqrts_max,0,25,"SD0.05");
//  Plot_VIcs_s(Ns, sqarr, cs_s200, sqrts_min, sqrts_max,0,25,"SD200");
//  Plot_VIcs_sII(int _Ns, double *_s_, double* cs1, double* cs2, string* leg, double _s_min, double _s_max, double ymin, double ymax, const char* label);
}

//********************************
void stage_IInormalize(double*t,double* dtarr,double*M2,double* dm2arr,const int Nt, const int Nm2,double* Fp2_t, double* f4_res, double* f_bg, double* alpha_IP,double* ImA, double* ImBg,double*ImR)
{
//~~~ calculate d2cs matix ~~~~
//  double _s=sq(13000);
////  double Wl=4;
//  double Wl=3.6;

 double _s=sq(8000);
 double Wl=3.6;

// double _s=sq(7000);
// double Wl=3.4;
// double Wl=4;

 double sl=Wl*Wl;
 double sh=_s*0.05; //GeV2

 double n_d2cs_dtdM2[Nt*Nm2];
 double csNorm;
 //~~ integrate ~~
 stage_IIcernel(n_d2cs_dtdM2, _s,t,M2, Nt,Nm2, Fp2_t,f4_res,  f_bg, alpha_IP,ImA,ImBg,ImR, Ares,Cbg);
 printf("***\n");
 printf("Integral dcs/(dtdM2) #sqrt(s)=%6.1f\n", sqrt(_s));
 csNorm=Integ2D(n_d2cs_dtdM2,Nt,Nm2,M2, dtarr, dm2arr,1,sh);
 printf("I full:  %6.4f  (ksi<0.05s)\n", csNorm);
 csNorm=Integ2D(n_d2cs_dtdM2, Nt,Nm2,M2, dtarr, dm2arr,1,sl);
 printf("I<%4.2fGeV:%6.4f\n",Wl, csNorm);
 csNorm=Integ2D(n_d2cs_dtdM2,Nt,Nm2,M2, dtarr, dm2arr,sl,sh);
 printf("I>%4.2fGeV:%6.4f   (ksi<0.05s)\n",Wl, csNorm);

 stage_IIcernel(n_d2cs_dtdM2,_s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares, 0);
 csNorm=Integ2D(n_d2cs_dtdM2, Nt,Nm2,M2, dtarr, dm2arr,1,sh);
 printf("I(Sg): %6.4f\n", csNorm);
 stage_IIcernel(n_d2cs_dtdM2, _s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, 0,Cbg);
 csNorm=Integ2D(n_d2cs_dtdM2, Nt,Nm2,M2, dtarr, dm2arr,1,sh);
 printf("I(Bg): %6.4f\n", csNorm);
 printf("***\n\n");
}

void stage_II_SDpredict()
{
 double w_max=15000;
const int Nt=300, Nm2=500;//canonical
// const int Nt=200, Nm2=6000;
//  const int Nt=500, Nm2=3000;

 double M2[Nm2], t[Nt];
 double dt[Nt]={0}, dm2[Nm2]={0};
 double f2[3*Nt], Fp2[Nt], f_bg[Nt], alpha_IP[Nt];
 double ImA[3*Nm2], ImBg[Nm2], ImR[Nm2];

 prestage_tlin(t,dt,Nt, 0,1);
 prestage_M2(M2,dm2,Nm2, 1,pow(w_max,2),"log"); //(!)1e8(!)
 stageI_t(Fp2, f2, f_bg, alpha_IP, t, Nt);
 stageI_M2(ImA,ImBg,ImR,M2,Nm2);

 double s[3]={sq(7000),sq(8000),sq(13000)};
 double n_d2cs_dtdM2[Nt*Nm2];
 double cs005, cs200;

 //~~ integrate ~~~~~~~~~~~~~
 printf("***\n");
 printf("                        cs200     0.05s\n",i,sqrt(s[i]),cs200 ,cs005);
 for(int i=0;i<3;i++)
  {stage_IIcernel(n_d2cs_dtdM2, s[i],t,M2, Nt,Nm2, Fp2,f2,f_bg, alpha_IP,ImA,ImBg,ImR, Ares,Cbg);
   cs005=Integ2D(n_d2cs_dtdM2, Nt,Nm2,M2, dt, dm2,1,0.05*s[i]);
   cs200=Integ2D(n_d2cs_dtdM2, Nt,Nm2,M2, dt, dm2,1,sq(200));
   printf("[%2d] sqrt(s):%6.3f  %10.6f %10.6f\n",i,sqrt(s[i]),cs200 ,cs005);}

  printf("#########\n\n");
}


void integr(double* d2cs,int Nt, int Nm2, double* M2, double *dt, double *dm2, double* dcsdt, double* dcsdM2, double* peakdcsdt, double* integdcsdM2, double sl=1, double su=12)
{
 for(it=0;it<Nt;it++)
  for(im2=0;im2<Nm2;im2++)
   if((M2[im2]<su)&&(M2[im2]>sl))
    {dcsdt[it]  +=d2cs[it*Nm2+im2]*dm2[im2];
     dcsdM2[im2]+=d2cs[it*Nm2+im2]*dt[it];}

 integdcsdM2[0]=0;
 for(im2=1;im2<Nm2;im2++)
  if((M2[im2]<su)&&(M2[im2]>sl))
   {integdcsdM2[im2] = integdcsdM2[im2-1]+dcsdM2[im2]*dm2[im2];}

 for(it=0;it<Nt;it++)
  for(im2=0;im2<Nm2;im2++)
   if((M2[im2]>peak_M2_min)&&(M2[im2]<peak_M2_max))
    {peakdcsdt[it]  +=d2cs[it*Nm2+im2]*dm2[im2];}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Plot_Iint_t (const int Nt, double* t,double* dcs1,double* dcs2,double* dcs3, TString* stLeg,double _m1,double _m2, const char* fout);
void Plot_Iint_m2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, TString* stLeg, const char* fout);
void Plot_Iintint_m2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, TString* stLeg, const char* fout);
void Plot_III_Mx(const int Nt, double* t,double **acs, double* fcs ,const double* Mx);
void Plot_Im2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, double* ta, const char* fout);
void Plot_It(const int Nt, double* t,double* dcs1,double* dcs2,double* dcs3, double* m2a, const char* fout);

void Plot_Vint(const int Nt, double* t, TString fout, double* B1_t, double* B2_t, double* B3_t, TString *legs, double hmin=0, double hmax=20);
void Plot_Vdif(const int Nt, double* t, TString fout, double* B1_t, double* B2_t, double* B3_t, TString *legs, double hmin=0, double hmax=25);
void Bslope(double* cs_t, double* B,int Nt, double*t)
 {for(it=0; it<Nt-1;it++){B[it]=(log(cs_t[it+1])-log(cs_t[it]))/(t[it+1]-t[it]);} } 
void mBslope(double* cs_t1,double* cs_t2, double* B,int Nm2, double dt)
 {for(int im=0; im<Nm2;im++){B[im]=(log(cs_t2[im])-log(cs_t1[im]))/dt;} } 
void zero(double*arr, int N){for(int i=0;i<N;i++)arr[i]=0;}

void test_SD(double* t,double* M2);
void stageII(double*t,double* dtarr,double*M2,double* dm2arr,const int Nt,const int Nm2, double*Fp2_t,double*f4_res,double*f_bg,double*alpha_IP, double*ImA,double*ImBg, double* ImR)
{
 double integ_dcsdt_A[Nt];        //,integ_dcsdt_B[Nt],        integ_dcsdt_C[Nt];
 double integ_dcsdM2_A[Nm2];      //,integ_dcsdM2_B[Nm2],      integ_dcsdM2_C[Nm2];
 double integ_dcsdt_A_peak[Nt];   //,integ_dcsdt_B_peak[Nt],   integ_dcsdt_C_peak[Nt];
 double integinteg_dcsdM2_A[Nm2]; //,integinteg_dcsdM2_B[Nm2] ,integinteg_dcsdM2_C[Nm2];
 zero(integ_dcsdt_A, Nt);       // zero(integ_dcsdt_B, Nt);        zero(integ_dcsdt_C, Nt);
 zero(integ_dcsdM2_A, Nm2);     // zero(integ_dcsdM2_B, Nm2);      zero(integ_dcsdM2_C, Nm2);
 zero(integ_dcsdt_A_peak, Nt);  // zero(integ_dcsdt_B_peak, Nt);   zero(integ_dcsdt_C_peak, Nt);
 zero(integinteg_dcsdM2_A, Nm2);// zero(integinteg_dcsdM2_B, Nm2); zero(integinteg_dcsdM2_C, Nm2);

 double integ_dcsdt_A_1[Nt]; zero(integ_dcsdt_A_1, Nt);  
 double integ_dcsdt_A_2[Nt]; zero(integ_dcsdt_A_2, Nt); 
 double integ_dcsdt_A_3[Nt]; zero(integ_dcsdt_A_3, Nt); 
 double integ_dcsdt_A_4[Nt]; zero(integ_dcsdt_A_4, Nt); 
 double integ_dcsdt_A_5[Nt]; zero(integ_dcsdt_A_5, Nt); 
 double d2cs_dtdM2[Nt*Nm2];

// strLeg[1]="#alpha'=0.25, b_{in}=0.2";
// for(int i=0;i<4;i++){printf("[%4i]%s; ",i ,(char*)strLeg[i].Data());} printf("\n");
 double _s=sq(7000);
 stage_IIcernel(d2cs_dtdM2, _s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);
 integr(d2cs_dtdM2, Nt,Nm2,M2, dtarr, dm2arr,
        integ_dcsdt_A, integ_dcsdM2_A, integ_dcsdt_A_peak, integinteg_dcsdM2_A, 1,0.05*_s);

 double tmp;
 const double Mx[5]={2.5, 3.5, 6, 12.0};                               
 for(im2=0;im2<Nm2;im2++)
  for(it=0;it<Nt;it++)
   {tmp = d2cs_dtdM2[it*Nm2+im2]*dm2arr[im2];
    if(M2[im2]<=Mx[0]) integ_dcsdt_A_1[it]+=tmp;
    if(M2[im2]<=Mx[1]) integ_dcsdt_A_2[it]+=tmp;
    if(M2[im2]<=Mx[2]) integ_dcsdt_A_3[it]+=tmp;
    if(M2[im2]<=Mx[3]) integ_dcsdt_A_4[it]+=tmp;}
 double* acs[4]={ integ_dcsdt_A_1, integ_dcsdt_A_2,
 integ_dcsdt_A_3, integ_dcsdt_A_4};

// //### test SD ###
// printf("%6s %6s %6s (SD:I)\n","t","M2","dcsdt");
// for(it=0;it<50;it++)
//  for(im2=0;im2<50;im2++)
//   {printf("%6.3f %6.3f %6.3f\n", t[it],M2[im2],0.5*d2cs_dtdM2[it*Nm2+im2]);}
//
// test_SD(t,M2);
//######################################

//~~~ B-slopes ~~~~
 double B_A_t[Nt-1];
 double B_A_peak_t[Nt-1];
 Bslope(integ_dcsdt_A, B_A_t, Nt,t);
 Bslope(integ_dcsdt_A_peak, B_A_peak_t, Nt,t);
// Bslope(integ_dcsdt_B, B_B_t, Nt,t);
// Bslope(integ_dcsdt_B_peak, B_B_peak_t, Nt,t);
{ TString legends[4];
  legends[0]="SD, #sqrt{s}=7TeV";
  legends[1].Form("M^{2}#in[%4.2f: 0.05s]",M2_min);
  legends[2].Form("M^{2}#in[%4.2f: %4.2f]",peak_M2_min,peak_M2_max);

  Plot_Vint(Nt,t,"B.SD.integ.eps", B_A_t, B_A_peak_t, NULL, legends);
  if(_PNG_) Plot_Vint(Nt,t,"B.SD.integ.png", B_A_t, B_A_peak_t, NULL, legends);
}

//~~~~~~~~~~~~~~~~~~~~~~~~
 //for(int i=0; i<Nt; i++)printf("t: %5.3f| %6.4f\n",t[i], integ_dcsdt_A[i]);
 //for(int i=0; i<Nm2; i++)printf("M2: %5.3f| %6.4f %6.4f\n",M2[i], integ_dcsdM2_B[i], integ_dcsdM2_C[i]);
 strLeg[0]="SD, #sigma(M_{x}^{2})";
 Plot_Iintint_m2(Nm2,M2,integinteg_dcsdM2_A,NULL ,NULL,strLeg,"intdcsdM2.eps");
 strLeg[0]="SD, integrated";
 Plot_Iint_m2(Nm2,M2,integ_dcsdM2_A, integ_dcsdM2_A, integ_dcsdM2_A, strLeg,"dcsdM2.eps");

 double I1, I2, I01,Ix; Ix=I1=I2=I01=0;
 double Ax,dIx;
 for(int i=0;i<Nm2,M2[i]<10;i++)
  { Ax=integ_dcsdM2_A[i];
   // printf("%.4f %8.4f\n",M2[i],Ax);
    dIx=Ax*dm2arr[i];
    if(M2[i]<2.44){I01+=dIx;} 
    if(M2[i]>1.67&&M2[i]<2.44){I1+=dIx;} 
    if(M2[i]>2.44&&M2[i]<3.5){I2+=dIx;} 
    Ix+=dIx;
  }
  printf("I01:%.3f\n",I01);
  printf("I1:%.3f\n",I1);
  printf("I2:%.3f\n",I2);
  printf("Ix:%.3f\n",Ix);


// Plot_Iint_t(Nt,t,integ_dcsdt_A, NULL,NULL, strLeg,M2_min,12,"dcsdt.full.eps");
// Plot_Iint_t(Nt,t,integ_dcsdt_A_peak, NULL,NULL, strLeg,peak_M2_min,peak_M2_max,"dcsdt.peak.eps");

if(_PNG_)
 {Plot_Iintint_m2(Nm2,M2,integinteg_dcsdM2_A, NULL,NULL,strLeg,"intdcsdM2.png");
 Plot_Iint_m2(Nm2,M2,integ_dcsdM2_A, integ_dcsdM2_A, integ_dcsdM2_A, strLeg,"dcsdM2.png");
// Plot_Iint_t(Nt,t,integ_dcsdt_A, NULL, NULL, strLeg,M2_min,12,"dcsdt.full.png");
// Plot_Iint_t(Nt,t,integ_dcsdt_A_peak, NULL, NULL, strLeg,peak_M2_min,peak_M2_max,"dcsdt.peak.png");
 }

 Plot_III_Mx(Nt,t,acs, integ_dcsdt_A, Mx);
}

void Plot_IIAdditionalDataII(const int Nt, double* t, double*dcsdt);
void stageIIAddonData()
{
  const int Nt=400, Nm2=2000;//6000;
  const double s=546*546;
  const double M2min=1, M2max=0.05*s;

  double M2[Nm2], t[Nt];
  double dtarr[Nt]={0}, dm2arr[Nm2]={0};
  double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
  double  ImA[3*Nm2], ImBg[Nm2], ImR[Nm2];
 //~~~~~~~~~~~~~~~~~
  prestage_tlin(t, dtarr, Nt, 0,1.5);
  prestage_M2(M2, dm2arr, Nm2, M2min, M2max,"log");
  stageI_M2(ImA,ImBg,ImR,M2,Nm2);
  stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, t, Nt);

  double d2cs_dtdM2[Nt*Nm2];
  double dcsdt[Nt]; zero(dcsdt, Nt);  
  // double _alf;
 stage_IIcernel(d2cs_dtdM2, s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);

 double cs546=0,tmp;
 for(int it=0;it<Nt;it++)
  for(int im2=0;im2<Nm2;im2++)
  if((M2[im2]<0.05*s)&&(M2[im2]>1))
    {dcsdt[it] +=tmp=d2cs_dtdM2[it*Nm2+im2]*dm2arr[im2];
     cs546+=tmp*dtarr[it];} 

 double _cs546=0;
 for(int it=0;it<Nt;it++)
  {_cs546+=dcsdt[it]*dtarr[it];} 

// for(int it=0;it<Nt;it+=40){printf("[%2d] t:%5.3f %5.3e\n",it,t[it],dcsdt[it]);}
 printf("### cs(546|ksi<0.05): %5.3f | %5.3f\n",cs546,_cs546);
 Plot_IIAdditionalDataII(Nt,t,dcsdt);
}


void Plot_VmSlope(const int Nm2,double* M2, TString fout, double*B1,double*B2,double*B3, TString *legs, double hmin=0,double hmax=20);
void stageIIother(double*t,double* dtarr,double*M2,double* dm2arr,const int Nt,const int Nm2)
{
 double d2cs1_dtdM2[3*Nm2];
//#########################################
 double dcs_dtdM2_A_M2_1[Nm2];
 double dcs_dtdM2_A_M2_2[Nm2];
 double dcs_dtdM2_A_M2_3[Nm2];
 double _s=49e6;//GeV2

 double  f4_res[3*Nt],Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
 double  ImA[3*Nm2], ImBg[Nm2],ImR[Nm2];
 double _t[3];
 _t[0]=-0.01; _t[1]=-0.05; _t[2]=-0.1; 
 stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, _t, 3);
 stageI_M2(ImA,ImBg,ImR,M2,Nm2);
 stage_IIcernel(d2cs1_dtdM2, _s,_t,M2, 3,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);

 for(im2=0;im2<Nm2;im2++)
  {dcs_dtdM2_A_M2_1[im2] = d2cs1_dtdM2[0*Nm2+im2];
   dcs_dtdM2_A_M2_2[im2] = d2cs1_dtdM2[1*Nm2+im2];
   dcs_dtdM2_A_M2_3[im2] = d2cs1_dtdM2[2*Nm2+im2];}

// for(int i=0;i<Nm2;i++) printf("[%6.3f], %6.3f  %6.3f\n",M2[i],dcs_dtdM2_A_M2_1[i],dcs_dtdM2_A_M2_3[i]);
 Plot_Im2(Nm2,M2,dcs_dtdM2_A_M2_1, dcs_dtdM2_A_M2_2, dcs_dtdM2_A_M2_3, _t, "d2cs_m2.eps");
 if(_PNG_) Plot_Im2(Nm2,M2,dcs_dtdM2_A_M2_1, dcs_dtdM2_A_M2_2, dcs_dtdM2_A_M2_3, _t, "d2cs_m2.png");

// /**/ Plot_Iint_m2(Nm2,M2,dcs_dtdM2_A_M2_1, dcs_dtdM2_A_M2_2, dcs_dtdM2_A_M2_3, strLeg,"xxxx.eps");
//#########################################
 double d2cs_dtdM2[Nt*3];
 double dcs_dtdM2_A_t_1[Nt];
 double dcs_dtdM2_A_t_2[Nt];
 double dcs_dtdM2_A_t_3[Nt];
 double _M2[3];

 prestage_tlin(t, dtarr, Nt, t_min,t_max);
 _M2[0]=2.0; _M2[1]=2.8; _M2[2]=8; 
 stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, t, Nt);
 stageI_M2(ImA,ImBg,ImR,_M2,3);
 stage_IIcernel(d2cs_dtdM2, _s,t,_M2, Nt,3, Fp2_t,f4_res,f_bg, alpha_IP,  ImA,ImBg,ImR, Ares,Cbg);

 for(it=0;it<Nt;it++)
  {dcs_dtdM2_A_t_1[it] = d2cs_dtdM2[it*3+0];
   dcs_dtdM2_A_t_2[it] = d2cs_dtdM2[it*3+1];
   dcs_dtdM2_A_t_3[it] = d2cs_dtdM2[it*3+2];
//   printf("%4.2f| %6.2f, %6.2f, %6.2f\n",t[it],dcs_dtdM2_A_t_1[it],dcs_dtdM2_A_t_2[it],dcs_dtdM2_A_t_3[it]);
  }

 if(_PNG_) Plot_It(Nt,t, dcs_dtdM2_A_t_1, dcs_dtdM2_A_t_2, dcs_dtdM2_A_t_3, _M2,"d2cs_t.png");
 Plot_It(Nt,t, dcs_dtdM2_A_t_1, dcs_dtdM2_A_t_2, dcs_dtdM2_A_t_3, _M2,"d2cs_t.eps");
// ~~~ Slopes ~~~
 double B_d2cs_A_t_1[Nt-1];
 double B_d2cs_A_t_2[Nt-1];
 double B_d2cs_A_t_3[Nt-1];

 Bslope(dcs_dtdM2_A_t_1, B_d2cs_A_t_1, Nt,t);
 Bslope(dcs_dtdM2_A_t_2, B_d2cs_A_t_2, Nt,t);
 Bslope(dcs_dtdM2_A_t_3, B_d2cs_A_t_3, Nt,t);
//~~~~~~~~~~~~~~~~~~~~~~~~
{ TString legends[4];
  for(int i=1;i<=3;i++)legends[i].Form("M^{2} = %5.2f",_M2[i-1]);
  legends[0]="SD, B-slope";//+strLeg[1];
 if(_PNG_) Plot_Vdif(Nt,t,"B.1.dMdt.png", B_d2cs_A_t_1,B_d2cs_A_t_2,B_d2cs_A_t_3, legends); 
  Plot_Vdif(Nt,t,"B.1.dMdt.eps", B_d2cs_A_t_1,B_d2cs_A_t_2,B_d2cs_A_t_3, legends); }

//#########################################

 double d2cs_M2_t1_7000[Nm2];
 double d2cs_M2_t2_7000[Nm2];
 double d2cs_M2_t1_1800[Nm2];
 double d2cs_M2_t2_1800[Nm2];
 double d2cs_M2_t1_546[Nm2];
 double d2cs_M2_t2_546[Nm2];
 double B7000[Nm2];
 double B1800[Nm2];
 double B546[Nm2];
 double _s7000=sq(7000);//GeV2
 double _s1800=sq(1800);//GeV2
 double _s546 =sq(546); //GeV2

// double f4_res[3*Nt],Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
// double ImA[3*Nm2], ImBg[Nm2],ImR[Nm2];
// double _t[u3];
 _t[0]=-0.1; _t[1]=-0.11; 
 stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, _t, 2);
 stageI_M2(ImA,ImBg,ImR,M2,Nm2);

 //## 7000
 stage_IIcernel(d2cs1_dtdM2, _s7000,_t,M2, 2,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);
 for(im2=0;im2<Nm2;im2++)
  {d2cs_M2_t1_7000[im2] = d2cs1_dtdM2[0*Nm2+im2];
   d2cs_M2_t2_7000[im2] = d2cs1_dtdM2[1*Nm2+im2];}
  mBslope(d2cs_M2_t1_7000,d2cs_M2_t2_7000, B7000,Nm2,-0.01);
 //## 1800
 stage_IIcernel(d2cs1_dtdM2, _s1800,_t,M2, 2,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);
 for(im2=0;im2<Nm2;im2++)
  {d2cs_M2_t1_1800[im2] = d2cs1_dtdM2[0*Nm2+im2];
   d2cs_M2_t2_1800[im2] = d2cs1_dtdM2[1*Nm2+im2];}
  mBslope(d2cs_M2_t1_1800,d2cs_M2_t2_1800, B1800,Nm2,-0.01);
 //## 546
 stage_IIcernel(d2cs1_dtdM2, _s546, _t,M2, 2,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);
 for(im2=0;im2<Nm2;im2++)
  {d2cs_M2_t1_546[im2] = d2cs1_dtdM2[0*Nm2+im2];
   d2cs_M2_t2_546[im2] = d2cs1_dtdM2[1*Nm2+im2];}
  mBslope(d2cs_M2_t1_546,d2cs_M2_t2_546, B546,Nm2,-0.01);

TString Legs[4]=
 {"",
  "#sqrt(s) =  7  TeV",
  "#sqrt(s) = 1.8 TeV",
  "#sqrt(s) = 546 GeV" };
 Plot_VmSlope(Nm2, M2, "B.SD.vs.M2.eps", B7000,B1800,B546, Legs, 0, 20);
 if(_PNG_) 
   Plot_VmSlope(Nm2, M2, "B.SD.vs.M2.png", B7000,B1800,B546, Legs, 0, 20);
}   

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//#############   Double Difraction  #############
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//const double Norm_DD = 1./4./Norm_el, s0=1./100;//00;
//int im1;

inline double DD_d3sg(int it, int im1, int im2, double s,
                    double* alpha,
                    double*t,double*M2, int Nt,int Nm2,
                    double*f4_res,double*f_bg,
                    double*ImA,double*ImBg,double*ImR,
                    double nA=Ares, double nCC=Cbg)
{
 double s_M1M2=s*s0/(M2[im1]*M2[im2]); 
// if(s_M1M2>1)
 return Norm_DD
  *(nA*(sigma_tot_pPN_m(it,im1, f4_res, ImA,Nt,Nm2)
        +nR0*sqrt(f4_res[it])*ImR[im1])
       *kinematic_factor(t[it], M2[im1])/(2*mp)
    +nCC*f_bg[it]*ImBg[im1])
  *(nA*(sigma_tot_pPN_m(it,im2, f4_res, ImA,Nt,Nm2)
        +nR0*sqrt(f4_res[it])*ImR[im2])
       *kinematic_factor(t[it], M2[im2])/(2*mp)
    +nCC*f_bg[it]*ImBg[im2])
  *pow(s_M1M2, alpha[it]);
 return 0;
}
/*
void stage_IIIcernel(double* d3cs_dtdM2dM2, double s, double* alpha,
                    double*t,double*M2, int Nt,int Nm2,
                    double*f4_res,double*f_bg,
                    double*ImA,double*ImBg,double*ImR,
                    double nA=Ares, double nCC=Cbg)
{
 int im1, ind=0;
 double propag;
 double Kf1, Sg1_b02, Bg1_b02, sg1_tot;
 double Kf2, Sg2_b02, Bg2_b02, sg2_tot;
 double M1M2, s_M1M2;
 
 for(it=0;it<Nt;it++)
  for(im1=0;im1<Nm2;im1++)
  for(im2=0;im2<Nm2;im2++)
   {
    ind=(it*Nm2+im1)*Nm2+im2;

    s_M1M2 = s*s0/(M2[im1]*M2[im2]); 
    propag = pow(s_M1M2, alpha[it]);
   // if(s_M1M2>1)
     {Kf1     = kinematic_factor(t[it], M2[im1]);
      Kf2     = kinematic_factor(t[it], M2[im2]);
      Sg1_b02 = sigma_tot_pPN_m(it,im1, f4_res, ImA,Nt,Nm2)
                +nR0*sqrt(f4_res[it])*ImR[im1];
      Sg2_b02 = sigma_tot_pPN_m(it,im2, f4_res, ImA,Nt,Nm2)
                +nR0*sqrt(f4_res[it])*ImR[im2];

      sg1_tot = (nA*Sg1_b02*Kf1/(2*mp) + nCC*f_bg[it]*ImBg[im1]);
      sg2_tot = (nA*Sg2_b02*Kf2/(2*mp) + nCC*f_bg[it]*ImBg[im2]);
  
      d3cs_dtdM2dM2[ind] = Norm_DD*sg1_tot*sg2_tot*propag;}
   // else {d3cs_dtdM2dM2[ind] = 0;}
//    ind++;
   }
}
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void stage_IIIs_dependese()
{
 const int Ns=10;//20;
 double w_min=10, w_max=1e4;
 double cs_s[Ns], sarr[Ns], sqarr[Ns];

 const int NtDD=100, Nm2DD=300;
// const int NtDD=200, Nm2DD=1000;//canonical
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,w_max*w_max,"log"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("cs(s) prestage.\n");
 //(!!!)
 const double x_min=log(w_min), x_max=log(w_max), dx=(x_max-x_min)/(Ns-1);
 for(int i=0;i<Ns;i++) {sqarr[i]=exp(x_min+dx*i); sarr[i]=sqarr[i]*sqarr[i];}

 int im1;
 double cs,cs_ut,cs_d;
 for(int is=0;is<Ns;is++)
  { cs=0;
    for(it=0;it<NtDD;it++)
     {cs_ut=cs_d=0;
      for(im1=0;im1<Nm2DD;im1++)
       for(im2=im1+1;im2<Nm2DD;im2++)
        cs_ut += DD_d3sg(it,im1,im2,sarr[is],alpha_IP,tDD,M2DD,NtDD,Nm2DD,f2DD,f_bgDD,ImADD,ImBgDD,ImRDD,Ares,Cbg)*dm2DD[im1]*dm2DD[im2];

      for(im1=0;im1<Nm2DD;im1++)
        cs_d  += DD_d3sg(it,im1,im1,sarr[is],alpha_IP,tDD,M2DD,NtDD,Nm2DD,f2DD,f_bgDD,ImADD,ImBgDD,ImRDD,Ares,Cbg)*dm2DD[im1]*dm2DD[im1];
      cs+=(2*cs_ut+cs_d)*dtDD[it];
      //printf("%6.3f| %6.3f %6.3f\n",tDD[it],cs_ut,cs_d);
     }
   cs_s[is]=cs;
   printf("sqrt(s):%6.3f  %10.6f\n",sqarr[is],cs);
  }
  Plot_VIcs_s(Ns, sqarr, cs_s, w_min, w_max,0,15,"DD");
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void stage_IIIs_dependese2()
{
 const int Ns=10;//3;//6;//10;
 double w_min=10, w_max=15000;
// w_max=7e3;
 double cs_s[Ns]={0}, sarr[Ns], sqarr[Ns];
 double cs_sm005[Ns]={0};
 double cs_s005[Ns]={0};
 double cs_s200[Ns]={0};
 const int NtDD=100, Nm2DD=300;
// const int NtDD=200, Nm2DD=500;//canonical
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,0.05*pow(w_max,2),"log"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("cs(s) prestage I.\n");

 const double x_min=log(w_min), x_max=log(w_max), dx=(x_max-x_min)/(Ns-1);
 for(int i=0;i<Ns;i++) {sqarr[i]=exp(x_min+dx*i);sarr[i]=sqarr[i]*sqarr[i];}

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
    { Vtx[it][im]= (Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                         +nR0*sqrt(f2DD[it])*ImRDD[im])
                       *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                  +Cbg*f_bgDD[it]*ImBgDD[im]);
//     printf("%6.3f %6.3f| %6.3f %6.3f %6.3f\n ",tDD[it],M2DD[im],Vtx[it][im],f_bgDD[it],ImBgDD[im]);
    }

 printf("cs(s) prestage II.\n");
//~~~~~ Integral ~~~~~~~~~~~~~~~~~~~~
 int im1,im2;
 double csm005,cs_utm005,cs_dm005,dcs;
 double cs005,cs_ut005,cs_d005;
 double cs200,cs_ut200,cs_d200,dcs200;
 double M1M2, s_M1M2;
 for(int is=0;is<Ns;is++)
  { csm005=cs005=cs200=0;
//sarr[is]=49e6; 
    for(it=0;it<NtDD;it++)
     {cs_utm005=cs_dm005=0;
      cs_ut005 =cs_d005 =0;
      cs_ut200 =cs_d200 =0;
      for(im1=0;im1<Nm2DD;im1++)
       for(im2=im1+1;im2<Nm2DD;im2++)
        {s_M1M2=sarr[is]*s0/(M2DD[im1]*M2DD[im2]);
         //s_M1M2=sarr[is]/sqrt(M2DD[im1]*M2DD[im2]);
          dcs =Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it])
                   *dm2DD[im1]*dm2DD[im2];
         M1M2=M2DD[im1]*M2DD[im2];
         if(log(s_M1M2)>3)      cs_utm005 += dcs;
         if(M1M2<0.05*sarr[is]) cs_ut005 += dcs;
         if(M2DD[im1]<4e4&&M2DD[im2]<4e4) cs_ut200 += dcs;
//          printf("%6.5f| %6.5f %6.5e\n",tDD[it],cs_ut,dcs);
        }
      for(im1=0;im1<Nm2DD;im1++)
        {s_M1M2=sarr[is]*s0/(M2DD[im1]*M2DD[im1]);
        // s_M1M2=sarr[is]/sqrt(M2DD[im1]*M2DD[im2]);
          dcs = sq(Vtx[it][im1]*dm2DD[im1])*pow(s_M1M2,alpha_IP[it]);
         M1M2=sq(M2DD[im1]);
         if(log(s_M1M2)>3)      cs_dm005 += dcs;
         if(M1M2<0.05*sarr[is]) cs_d005 += dcs;
         if(M2DD[im1]<4e4&&M2DD[im2]<4e4) cs_d200 += dcs;
        }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      csm005+=2*Norm_DD*(2*cs_utm005+cs_dm005)*dtDD[it];
      cs005 +=2*Norm_DD*(2*cs_ut005 +cs_d005 )*dtDD[it];
      cs200 +=2*Norm_DD*(2*cs_ut200 +cs_d200 )*dtDD[it];
//      printf("%6.3f| %6.3f %6.3f\n",tDD[it],Norm_DD*cs_ut,Norm_DD*cs_d);
     }
   cs_sm005[is]=csm005;
   cs_s005[is] =cs005 ;
   cs_s200[is] =cs200 ;
   printf("[%2d] sqrt(s):%6.3f  %10.6f\n",is,sqarr[is],csm005);
  }
//  Plot_VIcs_s(Ns, sqarr, cs_sm005, w_min, w_max,0,15,"DD");
  string leg[2]={"#Delta#eta>3",""};
  Plot_VIcs_sII(Ns, sqarr, cs_sm005,NULL,leg, w_min, w_max,0,15,"DD");
 // string leg2[2]={"#xi < 0.05","ln(s/(M_{1}^{2}M_{2}^{2}))>3"};
 // Plot_VIcs_sII(Ns, sqarr, cs_s005,cs_sm005,leg2, w_min, w_max,0,15,"DD");
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void stage_III_DDprediction()
{
 const int Ns=3;//6;//10;
 double w_max=14000;
 double sarr[Ns], sqarr[Ns]={7000, 8000,13000};
 for(int i=0;i<3;i++)sarr[i]=sq(sqarr[i]);

 double cs_s[Ns]={0};
 //const int NtDD=200, Nm2DD=1500;//canonical
 const int NtDD=200, Nm2DD=3000;
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,pow(w_max,2),"log"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("cs(s) prestage I.\n");

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
    { Vtx[it][im]= (Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                         +nR0*sqrt(f2DD[it])*ImRDD[im])
                       *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                  +Cbg*f_bgDD[it]*ImBgDD[im]);
//     printf("%6.3f %6.3f| %6.3f %6.3f %6.3f\n ",tDD[it],M2DD[im],Vtx[it][im],f_bgDD[it],ImBgDD[im]);
    }

 printf("cs(s) prestage II.\n");
//~~~~~ Integral ~~~~~~~~~~~~~~~~~~~~
 int im1,im2;
 double cs, dcs, cs_u,cs_d;

 double M1M2, s_M1M2;
 printf("\nDD predictions: \n");
 for(int is=0;is<Ns;is++)
  { cs=0;
    for(it=0;it<NtDD;it++)
     {cs_u =cs_d =0;
      for(im1=0;im1<Nm2DD;im1++)
       for(im2=im1+1;im2<Nm2DD;im2++)
        {s_M1M2=sarr[is]*s0/(M2DD[im1]*M2DD[im2]);
          dcs =Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it])
                   *dm2DD[im1]*dm2DD[im2];
         M1M2=M2DD[im1]*M2DD[im2];
         if(log(s_M1M2)>3) cs_u+=dcs;
        }
      for(im1=0;im1<Nm2DD;im1++)
        {s_M1M2=sarr[is]*s0/(M2DD[im1]*M2DD[im1]);
          dcs = sq(Vtx[it][im1]*dm2DD[im1])*pow(s_M1M2,alpha_IP[it]);
         M1M2=sq(M2DD[im1]);
         if(log(s_M1M2)>3) cs_d+=dcs;
        }
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      cs+=2*Norm_DD*(2*cs_u+cs_d)*dtDD[it];
     }
   cs_s[is]=cs;
   printf("[%2d] sqrt(s):%6.3f  %10.6f\n",is,sqarr[is],cs);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PlotDD_Iint_t (const int Nt, double* t, double* dcs1,double* dcs2, TString* stLeg, const char* fout);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void stage_IIIs_normalize2()
{
// double s=sq(1800);//49e6;
 double s=sq(7000);//49e6;
// const double Ares=0;
// const double Cbg=0;
// const int NtDD=200, Nm2DD=500;
// const int NtDD=200, Nm2DD=1000;//canonical
 const int NtDD=200, Nm2DD=3000;
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,s,"log"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("Normalize prestage I.\n");

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
    {Vtx[it][im]=(Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                       +nR0*sqrt(f2DD[it])*ImRDD[im])
                     *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                +Cbg*f_bgDD[it]*ImBgDD[im]);
//     printf("%6.3f %6.3f| %6.3f %6.3f %6.3f\n ",tDD[it],M2DD[im],Vtx[it][im],f_bgDD[it],ImBgDD[im]);
//            printf("[t:%6.3f m:%6.3f] | Kf:%6.3f; Sg:%6.3f; Bg:%6.3f;\n",
//    tDD[it],M2DD[im],
//kinematic_factor(tDD[it], M2DD[im]),
//sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)+nR0*sqrt(f2DD[it])*ImRDD[im],
//f_bgDD[it]*ImBgDD[im]);
}

 printf("Normalize prestage II.\n");
//~~~~~ Integral ~~~~~~~~~~~~~~~~~~~~
// printf("\n(!)>>>s0: %f\n\n",s0);
 int im1,im2;
 double dcsdt_sh[NtDD], dcsdtL[NtDD];
 double s_M1M2;
 double propag;

 double dcs, cs,cs_u,cs_d,  cs_xl,cs_xh, cs_xlu,cs_xld, cs_xhu,cs_xhd;
 double M2x=16;//GeV2
 cs=cs_xl=cs_xh=0;
 for(it=0;it<NtDD;it++)
  {cs_u=cs_d=cs_xlu=cs_xld=cs_xhu=cs_xhd=0;
   for(im1=0;im1<Nm2DD;im1++)
    for(im2=im1+1;im2<Nm2DD;im2++)
      {s_M1M2=s*s0/(M2DD[im1]*M2DD[im2]);
     //  propag=pow(s/(sqrt(M2DD[im1])+sqrt(M2DD[im2])),alpha_IP[it]); 
     //  propag=pow(s/sqrt(M2DD[im1]*M2DD[im2]),alpha_IP[it]); 
       propag=pow(s_M1M2,alpha_IP[it]); 
     //  if(s_M1M2>1)
       dcs=Vtx[it][im1]*Vtx[it][im2]*propag*dm2DD[im1]*dm2DD[im2];
//       printf("%6.5f| %6.5f %6.5e\n",tDD[it],cs_u,dcs);
       if(log(s_M1M2)>3)
        {if((M2DD[im1]*M2DD[im2])<M2x) cs_xlu+=dcs; else cs_xhu+=dcs;
         cs_u+=dcs;}
  //       printf(">[t:%6.3f m1:%6.3f m2:%6.3f] Ares:%6.3f,Cbg:%6.3f|d2cs:%6.3f, Vtx1:%6.3f, propag:%6.3f\n", tDD[it],M2DD[im1],M2DD[im2],Ares/(2*mp),Cbg,Vtx[it][im1]*Vtx[it][im2]*propag,Vtx[it][im1],propag);
//         printf("[t:%6.3f m1:%6.3f m2:%6.3f] propag:%6.3f| Vtx1:%6.3f1;\n",
 //        tDD[it],M2DD[im1],M2DD[im2],propag,Vtx[it][im1]);
      }
   for(im1=0;im1<Nm2DD;im1++)
      {s_M1M2=s*s0/(M2DD[im1]*M2DD[im1]);
      // propag=pow(s/(sqrt(M2DD[im1])+sqrt(M2DD[im2])),alpha_IP[it]); 
      // propag=pow(s/sqrt(M2DD[im1]*M2DD[im2]),alpha_IP[it]); 
       propag=pow(s_M1M2,alpha_IP[it]); 
      // if(s_M1M2>1)
        dcs=sq(Vtx[it][im1]*dm2DD[im1])*propag;

       if(log(s_M1M2)>3)
        {if(sq(M2DD[im1])<M2x) cs_xld+=dcs; else cs_xhd+=dcs;
         cs_d+=dcs;}
      }

//   printf("%6.3f| %6.3f %6.3f\n",tDD[it],Norm_DD*cs_u,Norm_DD*cs_d);
   dcsdtL[it] = Norm_DD*(2*cs_u+cs_d); //Long (full cs, deta>3 )
   dcsdt_sh[it]= Norm_DD*(2*cs_xlu+cs_xld);
   cs    += dcsdtL[it] *dtDD[it];
   cs_xl += dcsdt_sh[it]*dtDD[it];
   cs_xh += Norm_DD*(2*cs_xhu+cs_xhd)*dtDD[it];
  }

 cs    *=2;
 cs_xl *=2;
 cs_xh *=2;
 printf("Integral cs(DD) #sqrt(s)=%6.1f\n", sqrt(s));
 printf("I full:  %6.4f\n",  cs);
 printf("I<%4.1fGeV:%6.4f\n",  sqrt(M2x),cs_xl);
 printf("I>%4.1fGeV:%6.4f\n\n",sqrt(M2x),cs_xh);

//~~~ dcs/dt (t) ~~~~ 
 TString str[3];
 str[0].Form("DD, #sqrt{s} = %2.0fTeV",sqrt(s)/1000);
 str[1].Form("M_{1}^{2}#upointM_{2}^{2} < %2.0fGeV^{2}",M2x);;
 str[2]="#Delta#eta > 3";
 if(_PNG_) PlotDD_Iint_t (NtDD, tDD, dcsdt_sh, dcsdtL, str, "dcsdt.DD.png");
 PlotDD_Iint_t (NtDD, tDD, dcsdt_sh, dcsdtL, str, "dcsdt.DD.eps");

//~~~ B-slope ~~~~ 
 TString strB[3];
 strB[0]="B-slope DD, #sqrt{s} = 7 TeV";
 strB[1]=Form("M_{1}^{2}#upointM_{2}^{2}<%2.0fGeV^{2}",M2x);
 strB[2]="#Delta#eta>3";
 double Bs_t_DD[NtDD-1], Bl_t_DD[NtDD-1];
 Bslope(dcsdt_sh,  Bs_t_DD, NtDD, tDD);
 Bslope(dcsdtL, Bl_t_DD, NtDD, tDD);
 Plot_Vint(NtDD,tDD,"intB.DD.eps", Bs_t_DD, Bl_t_DD, NULL, strB,0,10);
 }
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void stage_IIIs_Risto_Invis()
{

 double Mx=3.4;//GeV
// double Mx=2*3.4;//GeV
 double s=sq(1800);//49e6;
// double s=sq(7000);//49e6;
// double Mx=2*3.4;//GeV
// double Mx=3.4;//GeV
// double s=sq(8000);//49e6;
// double Mx=2*3.6;//GeV
// double Mx=3.6;//GeV
//   double s=sq(13000);//49e6;
//  // double Mx=2*4;//GeV
//   double Mx=4;//GeV
// const double Ares=0;
// const double Cbg=0;
// const int NtDD=200, Nm2DD=500;
// const int NtDD=200, Nm2DD=1000;//canonical
 const int NtDD=200, Nm2DD=3000;
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,s,"log"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("Normalize prestage I.\n");

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
    {Vtx[it][im]=(Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                       +nR0*sqrt(f2DD[it])*ImRDD[im])
                     *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                +Cbg*f_bgDD[it]*ImBgDD[im]);
//     printf("%6.3f %6.3f| %6.3f %6.3f %6.3f\n ",tDD[it],M2DD[im],Vtx[it][im],f_bgDD[it],ImBgDD[im]);
//            printf("[t:%6.3f m:%6.3f] | Kf:%6.3f; Sg:%6.3f; Bg:%6.3f;\n",
//    tDD[it],M2DD[im],
//kinematic_factor(tDD[it], M2DD[im]),
//sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)+nR0*sqrt(f2DD[it])*ImRDD[im],
//f_bgDD[it]*ImBgDD[im]);
}

 printf("Normalize prestage II.\n");
//~~~~~ Integral ~~~~~~~~~~~~~~~~~~~~
// printf("\n(!)>>>s0: %f\n\n",s0);
 int im1,im2;
 double dcsdt_sh[NtDD], dcsdtL[NtDD];
 double s_M1M2;
 double propag;

 double dcs, cs,cs_u,cs_d,  cs_xl,cs_xh, cs_xlu,cs_xld, cs_xhu,cs_xhd;
 cs=cs_xl=cs_xh=0;
 for(it=0;it<NtDD;it++)
  {cs_u=cs_d=cs_xlu=cs_xld=cs_xhu=cs_xhd=0;
   for(im1=0;im1<Nm2DD;im1++)
    for(im2=im1+1;im2<Nm2DD;im2++)
      {s_M1M2=s*s0/(M2DD[im1]*M2DD[im2]);
     //  propag=pow(s/(sqrt(M2DD[im1])+sqrt(M2DD[im2])),alpha_IP[it]); 
     //  propag=pow(s/sqrt(M2DD[im1]*M2DD[im2]),alpha_IP[it]); 
       propag=pow(s_M1M2,alpha_IP[it]); 
     //  if(s_M1M2>1)
       dcs=Vtx[it][im1]*Vtx[it][im2]*propag*dm2DD[im1]*dm2DD[im2];
//       printf("%6.5f| %6.5f %6.5e\n",tDD[it],cs_u,dcs);
       if(log(s_M1M2)>3)
        {if(sqrt(M2DD[im1])+sqrt(M2DD[im2])<Mx) cs_xlu+=dcs; else cs_xhu+=dcs;
         cs_u+=dcs;}
  //       printf(">[t:%6.3f m1:%6.3f m2:%6.3f] Ares:%6.3f,Cbg:%6.3f|d2cs:%6.3f, Vtx1:%6.3f, propag:%6.3f\n", tDD[it],M2DD[im1],M2DD[im2],Ares/(2*mp),Cbg,Vtx[it][im1]*Vtx[it][im2]*propag,Vtx[it][im1],propag);
//         printf("[t:%6.3f m1:%6.3f m2:%6.3f] propag:%6.3f| Vtx1:%6.3f1;\n",
 //        tDD[it],M2DD[im1],M2DD[im2],propag,Vtx[it][im1]);
      }
   for(im1=0;im1<Nm2DD;im1++)
      {s_M1M2=s*s0/(M2DD[im1]*M2DD[im1]);
      // propag=pow(s/(sqrt(M2DD[im1])+sqrt(M2DD[im2])),alpha_IP[it]); 
      // propag=pow(s/sqrt(M2DD[im1]*M2DD[im2]),alpha_IP[it]); 
       propag=pow(s_M1M2,alpha_IP[it]); 
      // if(s_M1M2>1)
        dcs=sq(Vtx[it][im1]*dm2DD[im1])*propag;

       if(log(s_M1M2)>3)
        {if(2*sqrt(M2DD[im1])<Mx) cs_xld+=dcs; else cs_xhd+=dcs;
         cs_d+=dcs;}
      }

//   printf("%6.3f| %6.3f %6.3f\n",tDD[it],Norm_DD*cs_u,Norm_DD*cs_d);
   dcsdtL[it] = Norm_DD*(2*cs_u+cs_d); //Long (full cs, deta>3 )
   dcsdt_sh[it]= Norm_DD*(2*cs_xlu+cs_xld);
   cs    += dcsdtL[it] *dtDD[it];
   cs_xl += dcsdt_sh[it]*dtDD[it];
   cs_xh += Norm_DD*(2*cs_xhu+cs_xhd)*dtDD[it];
  }

 cs    *=2;
 cs_xl *=2;
 cs_xh *=2;
 printf("Integral cs(DD) #sqrt(s)=%6.1f\n", sqrt(s));
 printf("I full:  %6.4f\n",  cs);
 printf("M1+M2<%4.1fGeV:%6.4f\n",  (Mx),cs_xl);
 printf("M1+M2>%4.1fGeV:%6.4f\n\n",(Mx),cs_xh);
 }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Plot3D_DD_int_t (const int Nm2, double* M2, double* d2cs, const char* fout, const char* Ztitle);

void stage_IIIint_other() //2D
{
 double s=sq(7000);
// const double Ares=0;
// const double Cbg=0;
 const double M2x=10;//GeV2
// const int NtDD=200, Nm2DD=50;//1000
 const int NtDD=200, Nm2DD=200;//1000
// const int NtDD=200, Nm2DD=6000;
 double M2DD[Nm2DD], tDD[NtDD];
 double dtDD[NtDD]={0}, dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD], alpha_IP[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];

 prestage_tlin(tDD,dtDD,NtDD, 0,1);
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,M2x,"lin"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD,alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("Normalize prestage I.\n");

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
     Vtx[it][im]=(Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                       +nR0*sqrt(f2DD[it])*ImRDD[im])
                     *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                +Cbg*f_bgDD[it]*ImBgDD[im]);

 printf("Normalize prestage II.\n");
//~~~~~ Integral (int|t)~~~~~~~~~~~~~~~~~~~~
 double d2cs[Nm2DD][Nm2DD];
 double s_M1M2, dm2dm2;
 double d2s;
 int im1,im2;

 for(im1=0;im1<Nm2DD;im1++)
  for(im2=im1;im2<Nm2DD;im2++)
    {s_M1M2  = s*s0/(M2DD[im1]*M2DD[im2]);
     d2s=0;
     for(it=0;it<NtDD;it++)
      d2s+=Norm_DD*Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it])*dtDD[it];
     d2cs[im1][im2]=d2cs[im2][im1]=d2s;
//     printf("%6.4f %6.4f %12.8f\n",M2DD[im1], M2DD[im2], d2cs[im1][im2]);
//     printf("%6.4f %6.4f %12.8f\n",M2DD[im2], M2DD[im1], d2cs[im1][im2]);
    }

Plot3D_DD_int_t(Nm2DD, M2DD, &d2cs[0][0], "int_t_DD_2D.png", "d^{2}#sigma_{DD}/dM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-4}}");
//Plot3D_DD_int_t(Nm2DD, M2DD, &d2cs[0][0], "int_t_DD_2D.eps", "d^{2}#sigma_{DD}/dM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-4}}");
//Plot3D_DD_int_t(Nm2DD, M2DD, &d2cs[0][0], "int_t.DD.2D.root", "d^{2}#sigma_{DD}/dM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-4}}");
// PlotDD_Iint_t (NtDD, tDD, dcsdt, dcsdtL, str, "dcsdt.DD.eps");
//##################################################
}

void stage_IIIdiff_other() // 2D
{
 double s=49e6;
// const double Ares=0;
// const double Cbg=0;
 const double M2x=10;//GeV2
// const int NtDD=200, Nm2DD=50;//1000
 const int NtDD=3, Nm2DD=200;//1000
 double M2DD[Nm2DD], tDD[NtDD]={-0.1, -0.2, -0.3};
 double dm2DD[Nm2DD]={0};
 double f2DD[3*NtDD], Fp2DD[NtDD], f_bgDD[NtDD];
 double ImADD[3*Nm2DD], ImBgDD[Nm2DD], ImRDD[Nm2DD];
 double alpha_IP[NtDD];
 
 prestage_M2(M2DD,dm2DD,Nm2DD, 1,M2x,"lin"); //(!)1e8(!)
 stageI_t(Fp2DD, f2DD, f_bgDD, alpha_IP, tDD, NtDD);
 stageI_M2(ImADD,ImBgDD,ImRDD,M2DD,Nm2DD);
 printf("prestage I.\n");

//~~~~~ PreStage ~~~~~~~~~~~~~~~~~~~~
 double Vtx[NtDD][Nm2DD];
// double Norm_DD = 1./4./norm_el;
 int is,it,im;
 for(it=0;it<NtDD;it++)
   for(im=0;im<Nm2DD;im++)
     Vtx[it][im]=(Ares*(sigma_tot_pPN_m(it,im,f2DD,ImADD,NtDD,Nm2DD)
                       +nR0*sqrt(f2DD[it])*ImRDD[im])
                     *kinematic_factor(tDD[it], M2DD[im])/(2*mp)
                +Cbg*f_bgDD[it]*ImBgDD[im]);
 printf("Calculate.\n");
//~~~~~ Integral (int|t)~~~~~~~~~~~~~~~~~~~~
 double d3cs0[Nm2DD][Nm2DD];
 double d3cs1[Nm2DD][Nm2DD];
 double d3cs2[Nm2DD][Nm2DD];
 double s_M1M2, dm2dm2;
 double d2s;
 int im1,im2;

 for(im1=0;im1<Nm2DD;im1++)
  for(im2=im1;im2<Nm2DD;im2++)
    {s_M1M2  = s*s0/(M2DD[im1]*M2DD[im2]);
     {  
      it=0; d2s=Norm_DD*Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it]);
      d3cs0[im1][im2]=d3cs0[im2][im1]=d2s;
      it=1; d2s=Norm_DD*Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it]);
      d3cs1[im1][im2]=d3cs1[im2][im1]=d2s;
      it=2; d2s=Norm_DD*Vtx[it][im1]*Vtx[it][im2]*pow(s_M1M2, alpha_IP[it]);
      d3cs2[im1][im2]=d3cs2[im2][im1]=d2s;
     }
//     printf("%6.4f %6.4f %12.8f\n",M2DD[im1], M2DD[im2], d2cs[im1][im2]);
//     printf("%6.4f %6.4f %12.8f\n",M2DD[im2], M2DD[im1], d2cs[im1][im2]);
    }

Plot3D_DD_int_t(Nm2DD, M2DD, &d3cs0[0][0], "DD_2D_t01.eps", "d^{3}#sigma_{DD}/dtdM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-6}}|_{t=-0.1}");
Plot3D_DD_int_t(Nm2DD, M2DD, &d3cs1[0][0], "DD_2D_t02.eps", "d^{3}#sigma_{DD}/dtdM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-6}}|_{t=-0.2}");
Plot3D_DD_int_t(Nm2DD, M2DD, &d3cs2[0][0], "DD_2D_t03.eps", "d^{3}#sigma_{DD}/dtdM^{2}_{1}dM^{2}_{2} #(){mb#upointGeV^{-6}}|_{t=-0.3}");
// PlotDD_Iint_t (NtDD, tDD, dcsdt, dcsdtL, str, "dcsdt.DD.eps");
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

double re_al;
double th1, th2;
double GGG,f21;
double piscot, GGGsf21;
double A;
//***********************************
double Theta_function(double s1, double s2)
{ if(s1 > s2) return 1.; else if(s1 == s2) return 0.5; return 0;}

const double piscot_re= Pi/TMath::Tan(TMath::Pi()*(1.-delta));
const double GGG_re[3]=
{ROOT::Math::tgamma(1.-delta) * ROOT::Math::tgamma(lambda[0]+1.)/ROOT::Math::tgamma(lambda[0]-delta+2.) / pow(s_par[0], 1. - delta),
 ROOT::Math::tgamma(1.-delta) * ROOT::Math::tgamma(lambda[1]+1.)/ROOT::Math::tgamma(lambda[1]-delta+2.) / pow(s_par[1], 1. - delta),
 ROOT::Math::tgamma(1.-delta) * ROOT::Math::tgamma(lambda[2]+1.)/ROOT::Math::tgamma(lambda[2]-delta+2.) / pow(s_par[2], 1. - delta)};
const double GGGsf21_re[3] =
{ROOT::Math::tgamma(-1.*delta)*ROOT::Math::tgamma(lambda[0]+1.)*pow(s_par[0],delta)/ROOT::Math::tgamma(lambda[0]-delta+1.),
 ROOT::Math::tgamma(-1.*delta)*ROOT::Math::tgamma(lambda[1]+1.)*pow(s_par[1],delta)/ROOT::Math::tgamma(lambda[1]-delta+1.),
 ROOT::Math::tgamma(-1.*delta)*ROOT::Math::tgamma(lambda[2]+1.)*pow(s_par[2],delta)/ROOT::Math::tgamma(lambda[2]-delta+1.)};

int ire;
const int n_sum = 3;
double Re_alpha(double s)
{ re_al=alfa0;
  for(int ire=0; ire<n_sum; ire++)
    {      
      f21=0;GGG=0;
      th1 = Theta_function(s_par[ire], s);
      if(th1 > 0) {
        GGG=GGG_re[ire];
	if(s/s_par[ire] != 1.) {//See Abramovitz Stegun p. 556 
	  if(lambda[ire]-delta+2.<10)
	    f21 = ROOT::Math::hyperg(1., 1.-delta, lambda[ire]-delta+2., s/s_par[ire]);
	  else
	    f21 = ROOT::Math::hyperg(1., 1.-delta, 9.999, s/s_par[ire]);
	}
	else 
         {f21 = ROOT::Math::tgamma(lambda[ire]-delta+2.) *  ROOT::Math::tgamma(lambda[ire]) ;
	  f21 /= ROOT::Math::tgamma(lambda[ire]-delta+1.) *  ROOT::Math::tgamma(lambda[ire] + 1.); }
       } 

      piscot=0; GGGsf21 = 0.;
      th2 = Theta_function(s, s_par[ire]);
      if(th2 > 0) {
	piscot  = piscot_re*pow(s, delta - 1.);
	piscot *= pow((s-s_par[ire])/s, lambda[ire]);

	GGGsf21 =GGGsf21_re[ire] /s;

	if(s_par[ire]/s != 1.) {//See Abramovitz Stegun p. 556
	  if(delta-lambda[ire]>-10) GGGsf21 *= ROOT::Math::hyperg(delta-lambda[ire], 1., delta+1., s_par[ire]/s);
	  else GGGsf21 *= ROOT::Math::hyperg(-9.9999, 1., delta+1., s_par[ire]/s);
	}
	else 
         {GGGsf21 *= ROOT::Math::tgamma(delta+1.) * ROOT::Math::tgamma(lambda[ire]);
	  GGGsf21 /= ROOT::Math::tgamma(lambda[ire] + 1.) * ROOT::Math::tgamma(delta);}	
      }

      A = GGG * f21 * th1 + (piscot - GGGsf21) * th2;
      re_al += s/Pi*c1_par[ire]*A ;
    }
  return re_al;
}

//~~~
int ial;
double im_al;
double Im_alpha(double s)
{ im_al = 0.;
  for(ial=0; ial<3; ial++)
   // if(s>s_par[ial])im_al += c1_par[ial]*pow((s-s_par[ial])/s, Re_alpha(s_par[ial])); 
    if(s>s_par[ial])im_al += c1_par[ial]*pow((s-s_par[ial])/s, lambda[ial]); 
    if(s==0)return 0;
  return pow(s,delta)*im_al;
}
//~~~
int n;
inline double sigma_tot_pPN_m (int _it,int _im2, double* inF,double* imA, int Nt,int Nm2)
 {return  inF[_it]*imA[_im2] +inF[Nt+_it]*imA[Nm2+_im2] +inF[2*Nt+_it]*imA[2*Nm2+_im2];}
//**********************************
double factor, x_bj;
double kinematic_factor(double t, double M2)
{ 
if(M2<sqMpMpi) return 0; 
if(M2>sqMpMpi) return 1/sq(M2);else return 0; 

//x_bj=-t/(M2 - m_proton2 - t);
//  factor=pow((1-4*m_proton2*x_bj*x_bj/t),1.5)*(M2 - m_proton2);
//  if(factor!=0) return  x_bj*(1-x_bj)*(1-x_bj)/factor;
  return 0; 
}

//######################################
//############# Drawing ################
void SetGraphStyle(TGraph& gr, int mkSyle,int mkSize, int color=kBlack,int mkcolor=-1);
/*
void Plot_Vall(const int Nt, double* t)
{
//### B(dsg/dt) (SD) [all] ####  
  printf("\nB(dsg/dt) (SD) [all]\n");
 
  double _t[Nt-1]; for(it=0;it<Nt-1;it++) _t[it]=-(t[it+1]+t[it])/2;

//  TGraph gCross_1(Nt-1, _t, B_el);     gCross_1.SetLineWidth(2);
//  TGraph gCross_2(Nt-1, _t, B_el_r21); gCross_2.SetLineWidth(2);
//  TGraph gCross_3(Nt-1, _t, B_el_r28); gCross_3.SetLineWidth(2);
  TGraph gCross_4(Nt-1, _t, B_A_t);    gCross_4.SetLineWidth(1);
  TGraph gCross_5(Nt-1, _t, B_B_t);    gCross_5.SetLineWidth(1);
  TGraph gCross_6(Nt-1, _t, B_C_t);    gCross_6.SetLineWidth(1);
//  gCross_1.SetLineColor(kBlack);       gCross_1.SetLineStyle(9); 
//  gCross_2.SetLineColor(kGreen);// gCross_2.SetLineStyle(9);  
//  gCross_3.SetLineColor(kCyan); // gCross_3.SetLineStyle(9);
  gCross_4.SetLineColor(kBlack); 
  gCross_5.SetLineColor(kBlue); 
  gCross_6.SetLineColor(kRed);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString title; 
  title.Form("#sqrt{s} = 7 TeV, M2:[%2.1f;%2.1f] ",M2_min,M2_max);

  TLegend leg( 0.45, 0.25, 0.95, 0.45);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_4, title, "");
  leg.AddEntry(&gCross_4, "b=2.1, Kf", "l");
  leg.AddEntry(&gCross_5, "b=2.1", "l");
  leg.AddEntry(&gCross_6, "b=2.81", "l");

//  TLegend legUp( 0.45, 0.80, 0.95, 0.95);
//  legUp.SetBorderSize(0);
//  legUp.SetFillColor(0);
//  legUp.AddEntry(&gCross_1, "elastic pp", "");
//  legUp.AddEntry(&gCross_1, "(1-t/0.71)^{-8}","l");
//  legUp.AddEntry(&gCross_2, "exp(4*2.1*t)", "l");
//  legUp.AddEntry(&gCross_3, "exp(4*2.81*t)", "l");
//
  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 1, 1000,-0,30);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetTicky();
  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","B-slope #(){GeV^{-2}}", "");
  h_windowCross1_xx.Draw();
//  gCross_1.Draw("C SAME");
//  gCross_2.Draw("C SAME");
//  gCross_3.Draw("C SAME");
  gCross_4.Draw("C SAME");
  gCross_5.Draw("C SAME");
  gCross_6.Draw("C SAME");
  leg.Draw("SAME");
//  legUp.Draw("SAME");

  cCross_xx.Print("Ball_t.eps");
}
*/

//### B(dsg/dt) (SD) [orig] ####  
void Plot_Vel(const int Nt, double* t, double* B1, double* B2, double* B3, TString* sleg, const char* fout)
{
  double _t[Nt-1]; for(it=0;it<Nt-1;it++) _t[it]=-(t[it+1]+t[it])/2;
  TGraph gCross_1(Nt-1, _t, B1); gCross_1.SetLineWidth(2);
  TGraph gCross_2(Nt-1, _t, B2); gCross_2.SetLineWidth(2);
  TGraph gCross_3(Nt-1, _t, B3); gCross_3.SetLineWidth(2);
  gCross_2.SetLineColor(kRed);  
  gCross_3.SetLineColor(kBlue);  
  ////////////////////////
  //        DRAW!
  ////////////////////////

  TLegend leg( 0.55, 0.3, 0.93, 0.45);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  if(sleg[0]!="")leg.AddEntry(&gCross_1, sleg[0],"");//"#sqrt{s} = 7 TeV, elastic pp", "");
  leg.AddEntry(&gCross_1, sleg[1],"l");//"(1-t/0.71)^{-8}","l");
  leg.AddEntry(&gCross_2, sleg[2],"l");//"exp(4*2.1*t)", "l");
  leg.AddEntry(&gCross_3, sleg[3],"l");//"exp(4*2.81*t)", "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 1, 1000,-0,30);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetTicky();
  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","B #(){GeV^{-2}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//### (int|M2) dsg/dt (t) (SD) [orig, mod] ####  
void Plot_III_Mx(const int Nt, double* t, double **acs, double* fcs ,const double* Mx)
{
  double _t[Nt]; for(it=0;it<Nt;it++) _t[it]=-t[it];
  TGraph gCross_1 (Nt, _t, acs[0]); gCross_1.SetLineWidth(1);
  TGraph gCross_2 (Nt, _t, acs[1]); gCross_2.SetLineWidth(1);
  TGraph gCross_3 (Nt, _t, acs[2]); gCross_3.SetLineWidth(1);
  TGraph gCross_4 (Nt, _t, acs[3]); gCross_4.SetLineWidth(1);
//  TGraph gCross_5 (Nt, _t, acs[4]); gCross_5.SetLineWidth(1);
  TGraph gCross_6 (Nt, _t, fcs);    gCross_6.SetLineWidth(1);
  gCross_6.SetLineColor(kBlack);  
//  gCross_5.SetLineColor(kRed);  
  gCross_4.SetLineColor(kMagenta);  
  gCross_3.SetLineColor(kBlue);  
  gCross_2.SetLineColor(kCyan);  
  gCross_1.SetLineColor(kGreen);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString title; title.Form("#sqrt{s} = 7 TeV");
  TString sleg[6];
  for(int i=0;i<5;i++) sleg[i].Form("M^{2}#in [%3.1f; %3.1f]",M2_min,Mx[i]);
  sleg[5].Form("M^{2}#in [%3.1f; %3.1e]",M2_min,M2_max);

  TLegend leg( 0.65, 0.7, 0.95, 0.95);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, title, "");
  leg.AddEntry(&gCross_1, sleg[0], "l");
  leg.AddEntry(&gCross_2, sleg[1], "l");
  leg.AddEntry(&gCross_3, sleg[2], "l");
  leg.AddEntry(&gCross_4, sleg[3], "l");
//  leg.AddEntry(&gCross_5, sleg[4], "l");
  leg.AddEntry(&gCross_6, "M^{2}#in [1.0; 0.05s]", "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 1, 1000,1e-4,1e3);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{SD}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  gCross_4.Draw("C SAME");
//  gCross_5.Draw("C SAME");
  gCross_6.Draw("C SAME");
  leg.Draw("SAME");

  if(_PNG_) cCross_xx.Print("int.fr21_Mx_t.png");
  cCross_xx.Print("int.fr21_Mx_t.eps");
}

//### (int|M2) dsg/dt (t) (el) online  ####  
void Plot_IIIx_el(const int Nt, double* t, double*dcs1,TString sleg, const char* fout)
{
//"dcsdt_el.eps"
  int block;
  bool sym, bins;
  float dx[100]={0};

//### TOTEM #######
  const char* fdat="../Data/mx/ppel.TOTEM.datx";
  DataContainer dc[3];
  dc[0].Get(fdat, block=1, sym=true,  bins=false);
  dc[1].Get(fdat, block=2, sym=false, bins=false);

  TGraphErrors gTOTEM1(dc[0].pages[0]->N); SetErrGraph(gTOTEM1,dc);
  gTOTEM1.SetLineColor(kGreen+3);
  gTOTEM1.SetMarkerColor(kGreen+3);
  gTOTEM1.SetMarkerSize(0.9); gTOTEM1.SetMarkerStyle(24);

  TGraphErrors gTOTEM2(dc[1].pages[0]->N); SetErrGraph(gTOTEM2,dc+1);
//  gTOTEM2.SetLineColor(kGreen+3);
//  gTOTEM2.SetMarkerColor(kGreen+3);
  gTOTEM2.SetMarkerSize(0.9); gTOTEM2.SetMarkerStyle(24);

  const char* fdat2="../Data/mx/ppel.TOTEM.full.datx";
  dc[2].Get(fdat2, block=2, sym=false,  bins=true);
  int N=dc[2].pages[0]->N;
  for(int i=0;i<N;i++) // mkb -> mb
   {dc[2].pages[0]->Y[i]/=1.e3;
    dc[2].pages[0]->dYtotal[i]/=1.e3;}

  TGraphErrors gTOTEM3(dc[2].pages[0]->N); SetErrGraph(gTOTEM3,dc+2);
  gTOTEM3.SetLineColor(kRed+2);
  gTOTEM3.SetMarkerColor(kRed+2);
  gTOTEM3.SetMarkerSize(.9); gTOTEM3.SetMarkerStyle(8);
//### TOTEM #######

  double _t[Nt]; for(it=0;it<Nt;it++) _t[it]=-t[it];
  TGraph gCross_1 (Nt, _t, dcs1); gCross_1.SetLineWidth(1);
  gCross_1.SetLineStyle(7);
  gCross_1.SetLineColor(kBlue);  
  gCross_1.SetLineWidth(2);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TLegend leg( 0.38, 0.65, 0.95, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "#sqrt{s} = 7 TeV, elastic pp", "");
  leg.AddEntry(&gCross_1, " "," ");
  leg.AddEntry(&gCross_1, sleg,"l");
  leg.AddEntry(&gTOTEM1,  " ", "");
  leg.AddEntry(&gTOTEM1,  "TOTEM 2012-002, CERN-PH-EP-2012-239", "P");
  leg.AddEntry(&gTOTEM3,  "TOTEM 2011-001, CERN-PH-EP-2011-101", "P");

// TH2D h_windowCross1_xx("h_windowCross1", "",  600, 0,0.5, 1000,1e-2,600);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,  1, 1000,1e-5,1e3);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,0.5, 1000,1e-2,500);//*0:1
//  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,2.6, 1000,1e-5,500);//0:2.6
  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,0.5, 1000,5e-4,500);//0:2.6
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{el}^{pp}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();

  gTOTEM1.Draw("P Same");  
//  gTOTEM2.Draw("P Same");  
  gTOTEM3.Draw("P Same");  
  gCross_1.Draw("C SAME");
  leg.Draw("SAME");
  cCross_xx.Print(fout);
}

//### (int|M2) dsg/dt (t) (el) online  ####  
void Plot_IIIxx_el(const int Nt, double* t,
  double*dcs7000,double*dcs1800, double*dcs546, double*dcs540,TString sleg, const char* fout)
{
//"dcsdt_el.eps"
  int block;
  bool sym, bins;
  float dx[100]={0};

// printf("1.Hallo\n");
//### TOTEM #######
  const char* fdat="../Data/mx/ppel.TOTEM.datx";
  DataContainer dc[3];
  dc[0].Get(fdat, block=1, sym=true,  bins=false);
  dc[1].Get(fdat, block=2, sym=false, bins=false);
//  dc[0].Print(); //  dc[1].Print();

  TGraphErrors gTOTEM1(dc[0].pages[0]->N); SetErrGraph(gTOTEM1,dc);
  gTOTEM1.SetLineColor(kGreen+3);
  gTOTEM1.SetMarkerColor(kGreen+3);
  //gTOTEM1.SetLineColor(kMagenta+1);
  gTOTEM1.SetMarkerSize(.5); gTOTEM1.SetMarkerStyle(24);

  TGraphErrors gTOTEM2(dc[1].pages[0]->N); SetErrGraph(gTOTEM2,dc+1);
  gTOTEM2.SetLineColor(kGreen+3);
  gTOTEM2.SetMarkerColor(kGreen+3);
  gTOTEM2.SetMarkerSize(.5); gTOTEM2.SetMarkerStyle(24);

//  printf("1.5 Hallo\n");
  const char* fdat2="../Data/mx/ppel.TOTEM.full.datx";
  dc[2].Get(fdat2, block=2, sym=false,  bins=true);
  int N=dc[2].pages[0]->N;
  for(int i=0;i<N;i++) // mkb -> mb
   {dc[2].pages[0]->Y[i]/=1.e3;
    dc[2].pages[0]->dYtotal[i]/=1.e3;}

  TGraphErrors gTOTEM3(dc[2].pages[0]->N); SetErrGraph(gTOTEM3,dc+2);
  gTOTEM3.SetLineColor(kRed+2);
  gTOTEM3.SetMarkerColor(kRed+2);
  gTOTEM3.SetMarkerSize(.5); gTOTEM3.SetMarkerStyle(8);

// printf("2.Hallo\n");

//### pp(540|546, 630, 1800 ) #######
  int Pallete[8]= {kMagenta+2, kYellow+3, kYellow,kOrange-3, kSpring,  kRed, kBlue, kCyan+1};
  DataContainer dc2[8]; //  dc2[0].Print(); 
  const char* ppdat="../Data/mx/ppbar.sqrtsGt100.dc.Data";

  float dn=1e-2; float dn3=dn*dn*dn, dn2=dn*dn;
  float norm[8] = { dn3, dn2,dn2,dn2,dn2, dn2,  dn,dn};
  TGraphErrors gpp[8];
  for(int i=0;i<8;i++){
     dc2[i].Get(ppdat, block=i+1, sym=true,  bins=false);
     N=dc2[i].pages[0]->N;
     for(int ii=0;ii<N;ii++)
        {dc2[i].pages[0]->Y[ii]*=norm[i];
         dc2[i].pages[0]->dYtotal[ii]*=norm[i];}
//     printf("3.Hallo\n");
     gpp[i]=TGraphErrors(dc2[i].pages[0]->N,
                         dc2[i].pages[0]->X,
                         dc2[i].pages[0]->Y,
                      dx,dc2[i].pages[0]->dYtotal);
     gpp[i].SetLineColor(Pallete[i]);
     gpp[i].SetMarkerColor(Pallete[i]);
     gpp[i].SetMarkerSize(.5); gTOTEM1.SetMarkerStyle(24);
   }

//dc2[7].Print();

//###################################

  int N7000, N1800, N546, N540;
  N7000=N1800=N546=N540=0;

  double _t[Nt], _dcs1800[Nt],_dcs546[Nt],_dcs540[Nt];
  for(it=0;it<Nt;it++) 
    {_t[it]=-t[it];
    _dcs1800[it] =dcs1800[it]*dn;
    _dcs546[it]  =dcs546[it] *dn2;
    _dcs540[it]  =dcs540[it] *dn3;
    if(_t[it]<0.6){N7000++;N540++;}
    if(_t[it]<0.7){N1800++;}
    if(_t[it]<0.86){N546++;}
    
       }
  TGraph g7000 (N7000, _t, dcs7000); g7000.SetLineWidth(2);
  TGraph g1800 (N1800, _t,_dcs1800); g1800.SetLineWidth(2);
  TGraph g546  (N546 , _t,_dcs546 );  g546.SetLineWidth(2);
  TGraph g540  (N540 , _t,_dcs540 );  g540.SetLineWidth(2);
  //gCross_1.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TLegend leg( 0.6, 0.75, 0.95, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&g7000, "#sqrt{s} = 7 TeV, elastic pp", "");
  leg.AddEntry(&g7000, " "," ");
  leg.AddEntry(&g7000, sleg,"l");
  leg.AddEntry(&gTOTEM1,  " ", "");
  leg.AddEntry(&gTOTEM1,  "TOTEM 2012-002, CERN-PH-EP-2012-239", "PLE1");
  leg.AddEntry(&gTOTEM3,  "TOTEM 2011-001, CERN-PH-EP-2011-101", "PLE1");

// TH2D h_windowCross1_xx("h_windowCross1", "",  600, 0,0.2, 1000,1,600);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,  1, 1000,1e-5,1e3);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,0.5, 1000,1e-2,500);//*0:1

//*  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,2.6, 1000,1e-5,500);//0:2.6
  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,2.6, 1000,1e-8,500);//0:2.6
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{el}^{pp}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();

   gTOTEM1.Draw("P Same");  
 //  gTOTEM2.Draw("P Same");  
   gTOTEM3.Draw("P Same");  
 g7000.Draw("C SAME");
 g1800.Draw("C SAME");
 g546 .Draw("C SAME");
 g540 .Draw("C SAME");
 
 gpp[0].Draw("P Same");

 gpp[1].Draw("P Same");
 gpp[2].Draw("P Same");
 gpp[3].Draw("P Same");
 gpp[4].Draw("P Same");

// //--  gpp[5].Draw("P Same");

  gpp[6].Draw("P Same");
  gpp[7].Draw("P Same");

  leg.Draw("SAME");
 // printf("5.Hallo\n");
  cCross_xx.Print(fout);
}

//### (int|M2) dsg/dt (t) (el)  ####  
void Plot_IIIel(const int Nt, double* t, double*dcs1, double*dcs2, double*dcs3,TString* sleg, const char* fout)
{
//"dcsdt_el.eps"
  int block;
  bool sym, bins;
  float dx[100]={0};

  const char* fdat="../Data/mx/ppel.TOTEM.datx";
  DataContainer dc[3];
  dc[0].Get(fdat, block=1, sym=true,  bins=false);
  dc[1].Get(fdat, block=2, sym=false, bins=false);
//  dc[0].Print(); //  dc[1].Print();
  TGraphErrors gTOTEM1 (dc[0].pages[0]->N); SetErrGraph(gTOTEM1,dc);
  gTOTEM1.SetLineColor(kGreen+3);
  gTOTEM1.SetMarkerColor(kGreen+3);
  gTOTEM1.SetMarkerSize(.5); gTOTEM1.SetMarkerStyle(24);

  TGraphErrors gTOTEM2 (dc[1].pages[0]->N); SetErrGraph(gTOTEM2,dc+1);
  gTOTEM2.SetLineColor(kGreen+3);
  gTOTEM2.SetMarkerColor(kGreen+3);
  gTOTEM2.SetMarkerSize(.5); gTOTEM2.SetMarkerStyle(24);

//  printf("1.Here!\n");
  const char* fdat2="../Data/mx/ppel.TOTEM.full.datx";
  dc[2].Get(fdat2, block=2, sym=false,  bins=true);
  int N=dc[2].pages[0]->N;
  for(int i=0;i<N;i++) // mkb -> mb
   {dc[2].pages[0]->Y[i]/=1.e3;
    dc[2].pages[0]->dYtotal[i]/=1.e3;}
//  printf("2.Here!\n");
  TGraphErrors gTOTEM3(dc[2].pages[0]->N); SetErrGraph(gTOTEM3,dc+2);
  gTOTEM3.SetLineColor(kRed+2);
  gTOTEM3.SetMarkerColor(kRed+2);
  gTOTEM3.SetMarkerSize(.5); gTOTEM3.SetMarkerStyle(8);

  double _t[Nt]; for(it=0;it<Nt;it++) _t[it]=-t[it];
  TGraph gCross_1 (Nt, _t, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nt, _t, dcs2); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nt, _t, dcs3); gCross_3.SetLineWidth(2);
//  gCross_2.SetLineColor(kRed);  
  gCross_2.SetLineColor(kRed);  
  gCross_3.SetLineColor(kBlue);  
//  gCross_3.SetLineWidth(1);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TLegend leg( 0.6, 0.70, 0.95, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "#sqrt{s} = 7 TeV, elastic pp", "");
  leg.AddEntry(&gCross_1, " "," ");
  leg.AddEntry(&gCross_1, sleg[1],"l");
  leg.AddEntry(&gCross_1, " "," ");
  leg.AddEntry(&gCross_2, sleg[2],"l");
  leg.AddEntry(&gCross_3, sleg[3],"l");
  leg.AddEntry(&gTOTEM1,  " ", "");
  leg.AddEntry(&gTOTEM1,  "TOTEM 2012-002, CERN-PH-EP-2012-239", "PLE1");
  leg.AddEntry(&gTOTEM3,  "TOTEM 2011-001, CERN-PH-EP-2011-101", "PLE1");

// TH2D h_windowCross1_xx("h_windowCross1", "",  600, 0,0.5, 1000,1e-2,600);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,  1, 1000,1e-5,1e3);
// TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,0.5, 1000,1e-2,500); //**0:1
  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0,2.6, 1000,1e-5,500);//**0:2.6
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{el}^{pp}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");

  gTOTEM1.Draw("P Same");  
//  gTOTEM2.Draw("P Same");  
  gTOTEM3.Draw("P Same");  
  leg.Draw("SAME");
  cCross_xx.Print(fout);
}

//### (int|t) dsg/dM2 (M2) (SD) ####  
void Plot_Iint_m2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, TString* stLeg, const char* fout)
{
  TGraph gCross_1 (Nm2-1, M2, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nm2-1, M2, dcs2); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nm2-1, M2, dcs3); gCross_3.SetLineWidth(2);

  gCross_1.SetMarkerSize(1);  
  gCross_1.SetMarkerStyle(4);//20  
  gCross_1.SetMarkerColor(kBlack);  
  gCross_1.SetLineColor(kBlack);  
  gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString title; 
  title.Form("t:[%2.1f;%2.1f] ",t_min,t_max);
  TLegend leg( 0.6, 0.72, 0.95, 0.87);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, strLeg[0], "");
  leg.AddEntry(&gCross_1, title, "");
  leg.AddEntry(&gCross_1, strLeg[1], "");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 10, 1000,1e-4,2.7);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})","#frac{d#sigma}{dM^{2}}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("L SAME");
//  gCross_2.Draw("LP SAME");
//  gCross_3.Draw("CP SAME");
  leg.Draw("SAME");

// for(int i=0; i<Nm2; i++)printf("M2: %5.3f| %6.4f %6.4f\n",M2[i], dcs2[i], dcs3[i]);
  cCross_xx.Print(fout);
}

//### (int|t) dsg/dM2 (M2) (SD) ####  
void Plot_Iintint_m2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, TString* stLeg, const char* fout)
{
  TGraph gCross_1 (Nm2-1, M2, dcs1); gCross_1.SetLineWidth(2);
//  TGraph gCross_2 (Nm2-1, M2, dcs2); gCross_2.SetLineWidth(2);
//  TGraph gCross_3 (Nm2-1, M2, dcs3); gCross_3.SetLineWidth(2);
//  gCross_2.SetLineColor(kRed);  
//  gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString title; 
  title.Form("t:[%2.1f;%2.1f] ",t_min,t_max);
  TLegend leg( 0.6, 0.34, 0.95, 0.50);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, strLeg[0], "");
  leg.AddEntry(&gCross_1, title, "");
  leg.AddEntry(&gCross_1, strLeg[1], "l");
//  leg.AddEntry(&gCross_2, strLeg[2], "l");
//  leg.AddEntry(&gCross_3, strLeg[3], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 10, 1000,1e-4,5);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2}_{x} (GeV^{2})", "#sigma(M^{2}_{x}) (mb)", "");
//"#sigma(M^{2}_{x})=#int^{M^{2}_{x}}_{1}#frac{d#sigma_{SD}}{dM^{2}}dM^{2} (mb)", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
 // gCross_2.Draw("C SAME");
 // gCross_3.Draw("C SAME");
  leg.Draw("SAME");

// for(int i=0; i<Nm2; i++)printf("M2: %5.3f| %6.4f %6.4f\n",M2[i], dcs2[i], dcs3[i]);
  cCross_xx.Print(fout);
}

//### (int|M2) dsg/dt (t) (SD) [orig, mod] ####  
void Plot_Iint_t (const int Nt, double* t, double* dcs1,double* dcs2,double* dcs3, TString* stLeg,double _m1,double _m2, const char* fout)
{ 
  double _t[Nt]; for(it=0;it<Nt;it++) _t[it]=-t[it];
  TGraph gCross_1 (Nt, _t, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nt, _t, dcs2); gCross_2.SetLineWidth(2);
//  TGraph gCross_3 (Nt, _t, dcs3); gCross_3.SetLineWidth(2);
  gCross_2.SetLineColor(kRed);  
//  gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString title; 
  title.Form("M2:[%2.1f;%2.1f] ",_m1,_m2);

  TLegend leg( 0.6, 0.75, 0.95, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, stLeg[0], "");
  leg.AddEntry(&gCross_1, title, "");
  leg.AddEntry(&gCross_1, stLeg[1], "l");
 // leg.AddEntry(&gCross_2, stLeg[2], "l");
///  leg.AddEntry(&gCross_3, stLeg[3], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 1, 100,1e-4,1e3);
//  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 10, 100,1e-24,1e3);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{SD}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
//  gCross_2.Draw("C SAME");
//  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//### d2sg/dtdM2 (M2) (SD) [orig m2={1,2,3}] ####  
void Plot_It(const int Nt, double* t, double* dcs1,double* dcs2,double* dcs3, double* m2a, const char* fout)
{ double _t[Nt];
  for(it=0;it<Nt;it++) _t[it]=-t[it];

  TGraph gCross_1 (Nt, _t, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nt, _t, dcs2); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nt, _t, dcs3); gCross_3.SetLineWidth(2);
  gCross_1.SetLineStyle(1);  
  gCross_2.SetLineStyle(7);  
  gCross_3.SetLineStyle(3);  
 // gCross_2.SetLineColor(kRed);  
 // gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString s1; s1.Form("M^{2} = %4.1f GeV^{2}",m2a[0]);
  TString s2; s2.Form("M^{2} = %4.1f GeV^{2}",m2a[1]);
  TString s3; s3.Form("M^{2} = %4.1f GeV^{2}",m2a[2]);

  TLegend leg( 0.56, 0.67, 0.93, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "SD, #sqrt{s} = 7 TeV", "");
  leg.AddEntry(&gCross_1, s1, "l");
  leg.AddEntry(&gCross_2, s2, "l");
  leg.AddEntry(&gCross_3, s3, "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 1, 1000,2e-5,100);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d^{2}#sigma_{SD}}{dtdM^{2}}  #(){#frac{mb}{GeV^{4}}}", "");
  h_windowCross1_xx.Draw();
  h_windowCross1_xx.GetYaxis()->SetTitleOffset(1.9);
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}


//### d2sg/dtdM2 (M2) (SD) [orig t={1,2,3}] ####  
void Plot_Im2(const int Nm2, double* M2, double* dcs1,double* dcs2,double* dcs3, double* ta, const char* fout)
{
  TGraph gCross_1 (Nm2, M2, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nm2, M2, dcs2); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nm2, M2, dcs3); gCross_3.SetLineWidth(2);
  gCross_1.SetLineStyle(1);  
  gCross_2.SetLineStyle(7);  
  gCross_3.SetLineStyle(3);  
//  gCross_2.SetLineColor(kRed);  
//  gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString s1; s1.Form("t = %6.2f GeV^{2}", ta[0]);
  TString s2; s2.Form("t = %6.2f GeV^{2}", ta[1]);
  TString s3; s3.Form("t = %6.2f GeV^{2}", ta[2]);

  TLegend leg( 0.56, 0.67, 0.93, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "SD, #sqrt{s} = 7 TeV", "");
  leg.AddEntry(&gCross_1, s1, "l");
  leg.AddEntry(&gCross_2, s2, "l");
  leg.AddEntry(&gCross_3, s3, "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.001, 10, 1000,0,35);
//  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.001, 10, 1000,0,18);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})","#frac{d^{2}#sigma_{SD}}{dtdM^{2}}  #(){#frac{mb}{GeV^{4}}}", "");
  h_windowCross1_xx.GetYaxis()->SetTitleOffset(1.6);
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//############ DD part ############################################
/*
//*************************
void PlotDD_I()
{
//### d3sg/dtdM2dM2 (M2j) (DD) [t={.., .., ..}, M2k ={set_M2[1]}] ####  
  printf("\nd3sg/dtdM2dM2 (M2j) (DD) [t={%6.3f, %6.3f, %6.3f}, M2k=%6.3f]\n",
          t[set_t[0]],t[set_t[1]],t[set_t[2]], M2[set_M2[1]]);

  TGraph gCross_1 (Nm2, M2, d3cs_dtdM2dM2_M2_1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nm2, M2, d3cs_dtdM2dM2_M2_2); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nm2, M2, d3cs_dtdM2dM2_M2_3); gCross_3.SetLineWidth(2);
  gCross_2.SetLineColor(kRed);  
  gCross_3.SetLineColor(kBlue);  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TString s0; s0.Form("M_{2}^{2} = %6.3f",M2[set_M2[1]]);
  TString s1; s1.Form("t = %6.3f",t[set_t[0]]);
  TString s2; s2.Form("t = %6.3f",t[set_t[1]]);
  TString s3; s3.Form("t = %6.3f",t[set_t[2]]);

  TLegend leg( 0.6, 0.7, 0.93, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "DD, #sqrt{s} = 7 TeV", "");
  leg.AddEntry(&gCross_1, s0, "");
  leg.AddEntry(&gCross_1, s1, "l");
  leg.AddEntry(&gCross_2, s2, "l");
  leg.AddEntry(&gCross_3, s3, "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.001, 10., 10000,-0.001,40);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  cCross_xx.GetPad(1)->SetGrid();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "(M_{1})^{2} (GeV^{2})","#frac{d^{3}#sigma}{dtd(M_{1})^{2}d(M_{2})^{2}}  #(){#frac{mb}{GeV^{6}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print("DD_M2.eps");
//  cCross_xx.Print("DD_M2.root");
  ////////////////////////
}
*/
// ########### Slopes #######################
void Plot_VmSlope(const int Nm2,double* M2, TString fout, double*B1,double*B2,double*B3, TString *legs, double hmin,double hmax)
{
//for(int i=0;i<Nm2;i++) printf("%6.2f|%6.2f,%6.2f,%6.2f\n",M2[i],B1[i],B2[i],B3[i]);
  TGraph gcs_1(Nm2, M2, B1); gcs_1.SetLineWidth(2);
  TGraph gcs_2(Nm2, M2, B2); gcs_2.SetLineWidth(2);
  TGraph gcs_3(Nm2, M2, B3); gcs_3.SetLineWidth(2);
  gcs_2.SetLineColor(kRed);  
  gcs_3.SetLineColor(kGreen);  

  TLegend leg( 0.6, 0.23, 0.9, 0.4);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  if(legs[0]!="")leg.AddEntry(&gcs_1, legs[0], "");
  leg.AddEntry(&gcs_1, legs[1], "l");
  leg.AddEntry(&gcs_2, legs[2], "l");
  leg.AddEntry(&gcs_3, legs[3], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 10, 1000,hmin,hmax);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetTicky();

  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})","B #(){GeV^{-2}}", "");
  h_windowCross1_xx.Draw();
  gcs_1.Draw("C SAME");
  gcs_2.Draw("C SAME");
  gcs_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//##############################################
void Plot_Vdif(const int Nt, double* t, TString fout, double* B1_t, double* B2_t, double* B3_t, TString *legs, double hmin, double hmax)
{
//  printf((char*)(fout+" | "+title));
  double _t[Nt-1]; for(it=0;it<Nt-1;it++) _t[it]=-(t[it+1]+t[it])/2;
// for(it=0;it<Nt-1;it++) printf("[%4d] %8.6f | %7.3f\n",it,_t[it], B_A_t[it]);

  TGraph gCross_1(Nt-1, _t, B1_t); gCross_1.SetLineWidth(2);
  TGraph gCross_2(Nt-1, _t, B2_t); gCross_2.SetLineWidth(2);
  TGraph gCross_3(Nt-1, _t, B3_t); gCross_3.SetLineWidth(2);
  gCross_2.SetLineColor(kRed);  
  gCross_3.SetLineColor(kGreen);  
//  TString title; 
//  title.Form("#sqrt{s} = 7 TeV, M2:[%2.1f;%2.1f] ",M2_min,M2_max);

  TLegend leg( 0.6, 0.23, 0.9, 0.4);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, legs[0], "");
  leg.AddEntry(&gCross_1, legs[1], "l");
  leg.AddEntry(&gCross_2, legs[2], "l");
  leg.AddEntry(&gCross_3, legs[3], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 1, 1000,-0,20);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetTicky();

  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","B #(){GeV^{-2}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

void Plot_Vint(const int Nt, double* t, TString fout, double* B1_t, double* B2_t, double* B3_t, TString *legs, double hmin, double hmax)
{
//  printf((char*)(fout+" | "+title));
  double _t[Nt-1]; for(it=0;it<Nt-1;it++) _t[it]=-(t[it+1]+t[it])/2;
// for(it=0;it<Nt-1;it++) printf("[%4d] %8.6f | %7.3f\n",it,_t[it], B_A_t[it]);

  TGraph gCross_1(Nt-1, _t, B1_t); gCross_1.SetLineWidth(2);
  TGraph gCross_2(Nt-1, _t, B2_t); gCross_2.SetLineWidth(2);
  TGraph gCross_3(Nt-1, _t, B3_t); gCross_3.SetLineWidth(2);
  gCross_1.SetLineColor(kBlue);  
  if(B2_t!=NULL)  gCross_2.SetLineColor(kBlack);  
  if(B3_t!=NULL)  gCross_3.SetLineColor(kGreen);  
//  TString title="#sqrt{s} = 7 TeV"; 

  double dx=0;
   if(hmax<12)dx=0.26;
  TLegend leg( 0.58, 0.75-dx, 0.94, 0.95-dx);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
//  leg.AddEntry(&gCross_1, title, "");
  leg.AddEntry(&gCross_1, legs[0], "");
  leg.AddEntry(&gCross_1, legs[1], "l");
  if(B2_t!=NULL)  leg.AddEntry(&gCross_2, legs[2], "l");
  if(B3_t!=NULL)  leg.AddEntry(&gCross_3, legs[3], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.02, 1, 1000, hmin,hmax);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetTicky();

  //cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","B #(){GeV^{-2}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  if(B2_t!=NULL)gCross_2.Draw("C SAME");
  if(B3_t!=NULL)gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//##### Appendix(Background) ##############################################
void Plot_IDDex7(const int Nm2, double* M2, double* d2A, double* d2B, double* d2C)
{
  TGraph gCross_1(Nm2, M2, d2A); gCross_1.SetLineWidth(2);
  TGraph gCross_2(Nm2, M2, d2B); gCross_2.SetLineWidth(2);
  TGraph gCross_3(Nm2, M2, d2C); gCross_3.SetLineWidth(2);
//  gCross_1.SetLineColor(kRed);  
//  gCross_2.SetLineColor(kBlue);  
////  gCross_3.SetLineColor(kBlack);  

//  gCross_1.SetMarkerStyle(22);      gCross_2.SetMarkerStyle(24);
  gCross_1.SetMarkerSize(1.2);      gCross_2.SetMarkerSize(1.2);
  gCross_1.SetMarkerColor(kBlack);  gCross_2.SetMarkerColor(kBlack);
  gCross_1.SetLineColor  (kBlack);  gCross_2.SetLineColor  (kBlack);
  gCross_1.SetLineStyle(7); 

  TF1 gl1("g1","0.19/x",1,M2_max);
  gl1.SetLineColor(kGreen);
  
  TGraph g540("../Data/DDex7.540UA4.dat", "%lg %lg");
  g540.SetMarkerStyle(25);
  g540.SetMarkerSize(1.2);
  g540.SetMarkerColor(kBlack);
  g540.SetLineColor  (kBlack);
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TLegend leg( 0.55, 0.7, 0.93, 0.9);
  leg.SetBorderSize(0); leg.SetFillColor(0);
//  leg.AddEntry(&gCross_1, strLeg[0], "");
  leg.AddEntry(&gCross_1, "(single arm) " ,"");
  leg.AddEntry(&gCross_1, strLeg[1], "l");
  leg.AddEntry(&gCross_2, strLeg[2], "l");
 // leg.AddEntry(&gCross_3, strLeg[3], "l");
 // leg.AddEntry(&gl2, "Fit D/(M^{2})^{#varepsilon}" ,"l");
  leg.AddEntry(&gl1, " " ,"");
//  leg.AddEntry(&gl1, "0.19/M^{2}","l");
  leg.AddEntry(&g540,"540 GeV UA4","p");

  TH2D h_windowCross1_xx("h_windowCross1", "|t|=0.5", 10*Nm2, 2, 1e7, 1000,1e-8,1e-1);

  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})",
        "d^{2}#sigma/dtdM^{2}|_{t=-0.5} (mb GeV^{-4})", "");

  h_windowCross1_xx.Draw();

 // gCross_3.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_1.Draw("C SAME");
//  gl1.Draw("SAME");
  g540.Draw("P SAME");

  leg.Draw("SAME");

 // TF1 f2("ddexp7_540GeV", ddexp7_M2,0,1e7,2); f2.SetParameters(0.5, 540*540);  f2.SetLineColor(kBlue);
 // TF1 f3("ddexp7_45GeV",  ddexp7_M2,0,1e7,2); f3.SetParameters(0.5, 45*45);    f3.SetLineColor(kCyan);
 // TF1 f4("ddexp7_23.5GeV",ddexp7_M2,0,1e7,2); f4.SetParameters(0.5, 23.5*23.5);f4.SetLineColor(kGreen);
 // f2.Draw("same");
 // f3.Draw("same");
 // f4.Draw("same");

  if(_PNG_)cCross_xx.Print("cl_t.png");
  cCross_xx.Print("cl_t.eps");
}

//*************************
void Plot_IIAdditionalDataII(const int Nt, double* t, double*dcsdt)
//(const int Nm2, double* M2, double* d2A, double* d2B, double* d2C)
{
    int Nm2=100;//!!!/ 

    TGraphErrors gt_s546("../Data/New/Needed/dcsdt.ksi<0.05.ppbar_Xpbar.UA4.1987", "%lg %lg %lg");
    gt_s546.SetMarkerStyle(22);   
    gt_s546.SetMarkerSize(2);   
//    gt_s546.SetMarkerColor(kBlue);
//    gt_s546.SetLineColor  (kBlue);
    
    double _t[Nt];for(int i=0;i<Nt;i++)_t[i]=-t[i];  
    double _dcsdt[Nt];for(int i=0;i<Nt;i++)_dcsdt[i]=dcsdt[i]/2.;  
    TGraph gdcsdt(Nt, _t, _dcsdt); gdcsdt.SetLineWidth(2); 
//    gdcsdt.Print(); 

    TF1 gl1("gl1","33.5*exp(8.0*(-x)+2.3*x*x)",0.0,1.5);
    gl1.SetLineWidth(2);
    gl1.SetLineColor(kGray+3);
//    gl1.SetLineColor(kMagenta);
    gl1.SetLineStyle(7); 
    /////////////////////////
    TLegend leg( 0.45, 0.67, 0.93, 0.9);
    leg.SetBorderSize(0); leg.SetFillColor(0);
    leg.AddEntry(&gt_s546,"(single arm)"," ");
    leg.AddEntry(&gt_s546,"546GeV, UA4 SPS (1987)","P");
    leg.AddEntry(&gl1,    "33.5*exp(8.0*t + 2.3*t^{2})", "l");
    leg.AddEntry(&gdcsdt, "d#sigma/dt, #xi<0.05", "l");

  
    TH2D h_windowCross1_xx("h_windowCross1", "", 500, 0, 1.5, 500,2.5e-2,40);
//    TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, -0.01, 0.12, 1000,10,1e5);
    TCanvas cCross_xx("cCross", "cCross", 800, 800);
    cCross_xx.Divide(1, 1);
    cCross_xx.Draw();
  
    make_clean_pads(&cCross_xx, 1, 0, 0);
    cCross_xx.GetPad(1)->cd();
   // cCross_xx.GetPad(1)->SetGrid();
    cCross_xx.GetPad(1)->SetLogy();
//    cCross_xx.GetPad(1)->SetLogx();
    sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})",
          "#frac{d#sigma_{SD}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
    h_windowCross1_xx.GetYaxis()->SetTitleOffset(1.3);
  
    h_windowCross1_xx.Draw();
    gl1.Draw("C SAME");
    gdcsdt.Draw("C SAME");
    gt_s546.Draw("P SAME");
    leg.Draw("SAME");

    cCross_xx.Print("AdditionalDataII.eps");
}


//*************************
void Plot_xGoul()
{
  int block; bool sym, bins;
  DataContainer dc[2];
  dc[0].Get("Data/ggSD1.Data", block=1, sym=false,  bins=false);
  dc[1].Get("Data/ggSD2.Data", block=1, sym=false,  bins=false);
//  dc[0].Print();
  TGraphErrors gr1(dc[0].pages[0]->N); SetErrGraph(gr1,dc);
  gr1.SetLineColor(kBlack); gr1.SetMarkerColor(kBlack);
  gr1.SetMarkerSize(1); gr1.SetMarkerStyle(8);
  TGraphErrors gr2(dc[1].pages[0]->N); SetErrGraph(gr2,dc+1);
  gr2.SetLineColor(kBlack); gr2.SetMarkerColor(kBlack);
  gr2.SetMarkerSize(1); gr2.SetMarkerStyle(24);
  
  /////////////////////////
  TLegend leg( 0.58, 0.75, 0.93, 0.9);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gr2, "|t|=0.035", "P");
  leg.AddEntry(&gr1, "|t|=0.05", "P");

  TH2D h_windowCross1_xx("h_windowCross1", "", 100,1,30, 1000,0,6);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1,1);
  cCross_xx.Draw();
  
  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
//  cCross_xx.GetPad(1)->SetGrid();
//  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})",
        "#frac{d^{2}#sigma}{dtdM^{2}}  #(){#frac{mb}{GeV^{4}}}", "");
  
  h_windowCross1_xx.Draw();
  gr1.Draw("P SAME");
  gr2.Draw("P SAME");
  leg.Draw("SAME");
  
  cCross_xx.Print("xGoul.eps");
}

//*************************
void Plot_IAdditionalData(const int Nm2, double* M2, double* d2A, double* d2B)
{
//  int Nm2=100;//!!!/ 
    TGraph gCross_1(Nm2, M2, d2A); gCross_1.SetLineWidth(2);
    TGraph gCross_2(Nm2, M2, d2B); gCross_2.SetLineWidth(2);
  //  TGraph gCross_3(Nm2, M2, d2C); gCross_3.SetLineWidth(2);
    gCross_2.SetLineColor(kBlue);  
  //  gCross_3.SetLineColor(kRed);  
//  
//    TF1 gl1("g1","[0]/pow(x,[1])",1,M2_max);
//    TF1 gl2("g2","[0]/pow(x,[1])",1,M2_max);
//    TF1 gl3("g3","[0]/pow(x,[1])",1,M2_max);
//    double norm=5.0;
//    gl1.SetParameters(norm, 1.05);gl1.SetParNames("D", "epsilon");gl1.SetLineColor(kBlue);
//    gl2.SetParameters(norm, 1.1 );gl2.SetParNames("D", "epsilon");gl2.SetLineColor(kGreen);
//    gl3.SetParameters(norm, 1.15);gl3.SetParNames("D", "epsilon");gl3.SetLineColor(kRed);  
//    gl1.SetLineStyle(7); 
//    gl2.SetLineStyle(7); 
//    gl3.SetLineStyle(7); 

  TGraphErrors g1800D("../Data/goul.1800GeV.dat", "%lg %lg %lg %lg");
  TGraphErrors g546D("../Data/goul.546GeV.dat", "%lg %lg %lg %lg");
  g546D.SetMarkerStyle(22);     g1800D.SetMarkerStyle(24);
  g546D.SetMarkerSize(1.2);     g1800D.SetMarkerSize(1.2);
  g546D.SetMarkerColor(kBlue);  g1800D.SetMarkerColor(kBlue);
  g546D.SetLineColor  (kBlue);  g1800D.SetLineColor  (kBlue);

  
    TGraphErrors g1800dat("../Data/New/Needed/dcsdtdksi.FromCDFfit.Goul.1.8TeV.t=0.05.ksi<0.05.dat", "%lg %lg %lg");
    TGraphErrors g546dat ("../Data/New/Needed/dcsdtdksi.FromCDFfit.Goul.546GeV.t=0.05.ksi<0.05.dat", "%lg %lg %lg");
    g546dat.SetMarkerStyle(22);      g1800dat.SetMarkerStyle(24);
    g546dat.SetMarkerSize(1.2);      g1800dat.SetMarkerSize(1.2);
    g546dat.SetMarkerColor(kBlack);  g1800dat.SetMarkerColor(kBlack);
    g546dat.SetLineColor  (kBlack);  g1800dat.SetLineColor  (kBlack);
//~~~ Convert d2cs(ksi)/dtdksi -> d2cs(M2)/dtdM2 ~~~
    Double_t* x =g546dat.GetX();
    Double_t* ex=g546dat.GetEX();
    Double_t* y =g546dat.GetY();
    Double_t* ey=g546dat.GetEY();
    Int_t     N =g546dat.GetN();
    Double_t  s =pow(546,2);
    for(int i=0;i<N;i++) {y[i]/=s; x[i]*=s; ey[i]/=s; ex[i]*=s;}

    x =g1800dat.GetX();
    ex=g1800dat.GetEX();
    y =g1800dat.GetY();
    ey=g1800dat.GetEY();
    N =g1800dat.GetN();
    s =pow(1800,2);
    for(int i=0;i<N;i++) {y[i]/=s; x[i]*=s; ey[i]/=s; ex[i]*=s;}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    /////////////////////////
    TLegend leg( 0.58, 0.75, 0.93, 0.9);
    leg.SetBorderSize(0); leg.SetFillColor(0);
  //  leg.AddEntry(&gCross_1, strLeg[0], "");
    leg.AddEntry(&gCross_1, strLeg[1], "l");
    leg.AddEntry(&gCross_2, strLeg[2], "l");
//  //  leg.AddEntry(&gCross_3, strLeg[3], "l");
//   // leg.AddEntry(&gl2, "Fit D/(M^{2})^{#varepsilon}" ,"l");
//    leg.AddEntry(&gl2, " " ,"");
//    leg.AddEntry(&gl3, "5/(M^{2})^{1.15}","l");
//    leg.AddEntry(&gl2, "5/(M^{2})^{1.10}","l");
//    leg.AddEntry(&gl1, "5/(M^{2})^{1.05}","l");
//    leg.AddEntry(&gl2, " " ,"");
    leg.AddEntry(&g546dat, "546 GeV,  |t|=0.05, CDF","P");
    leg.AddEntry(&g1800dat,"1800 GeV, |t|=0.05, CDF","P");

  
//    TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, 1, 1e6, 1000,1e-5,30);
    TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, 1e2, 1e6, 1000,1e-5,1e-1);
    TCanvas cCross_xx("cCross", "cCross", 800, 800);
    cCross_xx.Divide(1, 1);
    cCross_xx.Draw();
  
    make_clean_pads(&cCross_xx, 1, 0, 0);
    cCross_xx.GetPad(1)->cd();
//    cCross_xx.GetPad(1)->SetGrid();
    cCross_xx.GetPad(1)->SetLogy();
    cCross_xx.GetPad(1)->SetLogx();
    sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})",
          "#frac{d^{2}#sigma}{dtdM^{2}}  #(){#frac{mb}{GeV^{4}}}", "");
  
    h_windowCross1_xx.Draw();
  //  gCross_3.Draw("C SAME");
    gCross_2.Draw("C SAME");
    gCross_1.Draw("C SAME");
//    gl1.Draw("SAME");
//    gl2.DrawClone("SAME");
//    gl3.Draw("SAME");
  

    g546dat.Draw("P SAME");
    g1800dat.Draw("P SAME");

    g546D.Draw("P SAME");
    g1800D.Draw("P SAME");
//  /*  gCross_1.Fit(&gl2,"R","SAME",3.1,7.7);
//    s1.Form(           "D = %5.2f +/- %5.2f",gl2.GetParameter(0),gl2.GetParError(0));
//    leg.AddEntry(&gl2,s1,"");
//    s1.Form("#varepsilon  = %5.2f +/- %5.2f",gl2.GetParameter(1),gl2.GetParError(1));
//    leg.AddEntry(&gl2,s1,"");
//  */
    leg.Draw("SAME");
  
  // if(_PNG_)
     cCross_xx.Print("AdditionalData.eps");
    // cCross_xx.Print("AdditionalData.png");
  //  cCross_xx.Print("gl_t.eps");
  //  cCross_xx.Print("gl_t.root");
}

//*************************
void Plot_IGoulian(const int Nm2, double* M2, double* d2A, double* d2B, double* d2C)
{
  TGraph gCross_1(Nm2, M2, d2A); gCross_1.SetLineWidth(2);
  TGraph gCross_2(Nm2, M2, d2B); gCross_2.SetLineWidth(2);
//gCross_1.Print();
//  printf("A\n");
//  for(int i=0;i<Nm2;i++) if(M2[i]>1000) {printf("%f\t%f\n",M2[i], d2A[i]);}
//  printf("B\n");
//  for(int i=0;i<Nm2;i++) if(M2[i]>1000) {printf("%f\t%f\n",M2[i], d2B[i]);}

//  TGraph gCross_3(Nm2, M2, d2C); gCross_3.SetLineWidth(2);
  gCross_1.SetLineStyle(1);  
  gCross_2.SetLineStyle(4);  
  gCross_2.SetLineColor(kBlack);  
//  gCross_2.SetLineColor(kBlue);  
//  gCross_3.SetLineColor(kRed);  

  TF1 gl1("g1","[0]/pow(x,[1])",1,M2_max);
  TF1 gl2("g2","[0]/pow(x,[1])",1,M2_max);
  TF1 gl3("g3","[0]/pow(x,[1])",1,M2_max);
  double norm=5.0;
  gl1.SetParameters(norm, 1.05);gl1.SetParNames("D", "epsilon");gl1.SetLineColor(kBlue);
  gl2.SetParameters(norm, 1.1 );gl2.SetParNames("D", "epsilon");gl2.SetLineColor(kGreen);
  gl3.SetParameters(norm, 1.15);gl3.SetParNames("D", "epsilon");gl3.SetLineColor(kRed);  
  gl1.SetLineStyle(7); 
  gl2.SetLineStyle(7); 
  gl3.SetLineStyle(7); 

  TGraphErrors g1800D("../Data/goul.1800GeV.dat", "%lg %lg %lg %lg");
  TGraphErrors g546D("../Data/goul.546GeV.dat", "%lg %lg %lg %lg");
  g546D.SetMarkerStyle(22);      g1800D.SetMarkerStyle(24);
  g546D.SetMarkerSize(1.2);      g1800D.SetMarkerSize(1.2);
  g546D.SetMarkerColor(kBlack);  g1800D.SetMarkerColor(kBlack);
  g546D.SetLineColor  (kBlack);  g1800D.SetLineColor  (kBlack);

  
  /////////////////////////
  //        DRAW!
  ////////////////////////
  TLegend leg( 0.5, 0.67, 0.93, 0.9);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "(single arm)", "");
  leg.AddEntry(&gCross_1, strLeg[1], "l");
  leg.AddEntry(&gCross_2, strLeg[2], "l");
//  leg.AddEntry(&gCross_3, strLeg[3], "l");
 // leg.AddEntry(&gl2, "Fit D/(M^{2})^{#varepsilon}" ,"l");
  leg.AddEntry(&gl2, " " ,"");
//  leg.AddEntry(&gl3, "5/(M^{2})^{1.15}","l");
//  leg.AddEntry(&gl2, "5/(M^{2})^{1.10}","l");
//  leg.AddEntry(&gl1, "5/(M^{2})^{1.05}","l");
//  leg.AddEntry(&gl2, " " ,"");
  leg.AddEntry(&g546D, "546 GeV,  CDF","P");
  leg.AddEntry(&g1800D,"1800 GeV, CDF","P");

  TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, 1, 1e6, 1000,1e-5,30);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
//  cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
//  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})",
//        "#frac{d^{2}#sigma}{dtdM^{2}}|_{t=0.05}  #(){#frac{mb}{GeV^{4}}}", "");
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M_{x}^{2} (GeV^{2})",
        "#frac{d^{2}#sigma}{dtdM_{x}^{2}}|_{t=0.05}  #(){#frac{mb}{GeV^{4}}}", "");

  h_windowCross1_xx.Draw();
//  gCross_3.Draw("C SAME");
  gCross_2.Draw("C");
  gCross_1.Draw("C SAME");
//  gl1.Draw("SAME");
//  gl2.DrawClone("SAME");
//  gl3.Draw("SAME");

  g546D.Draw("P SAME");
  g1800D.Draw("P SAME");
/*  gCross_1.Fit(&gl2,"R","SAME",3.1,7.7);
  s1.Form(           "D = %5.2f +/- %5.2f",gl2.GetParameter(0),gl2.GetParError(0));
  leg.AddEntry(&gl2,s1,"");
  s1.Form("#varepsilon  = %5.2f +/- %5.2f",gl2.GetParameter(1),gl2.GetParError(1));
  leg.AddEntry(&gl2,s1,"");
*/
  leg.Draw("SAME");

 // TF1 df1("ddexp7_1800GeV",ddexp7_M2,1,2e5,2);df1.SetParameters(0.05,1800*1800);df1.SetLineColor(kMagenta);
 // TF1 df2("ddexp7_546GeV", ddexp7_M2,1,2e5,2);df2.SetParameters(0.05,  546*546);df2.SetLineColor(kBlue);
 // TF1 df3("ddexp7_20GeV",  ddexp7_M2,1,2e5,2);df3.SetParameters(0.05,    20*20);df3.SetLineColor(kCyan);
 // TF1 df4("ddexp7_14GeV",  ddexp7_M2,1,2e5,2);df4.SetParameters(0.05,    14*14);df4.SetLineColor(kGreen);
 // df1.SetLineStyle(1);
 // df2.SetLineStyle(1);
 // df3.SetLineStyle(1);
 // df4.SetLineStyle(1);

 // TF1 dP1("ddexp7P_1800GeV",ddexp7P_M2,1,2e5,2);dP1.SetParameters(0.05,1800*1800);dP1.SetLineColor(kMagenta);
 // TF1 dP2("ddexp7P_546GeV", ddexp7P_M2,1,2e5,2);dP2.SetParameters(0.05,  546*546);dP2.SetLineColor(kBlue);
 //  dP1.SetLineStyle(10);
 //  dP2.SetLineStyle(10);
//  df1.Draw("same");
//  df2.Draw("same");
//  df3.Draw("same");
//  df4.Draw("same");
//  dP1.Draw("same");
//  dP2.Draw("same");

  if(_PNG_) cCross_xx.Print("gl_t.png");
  cCross_xx.Print("gl_t.eps");
//  cCross_xx.Print("gl_t.root");
}

void Plot_IGoulFlat(const char* fout, const int Nm2, double* M2, double* d2A, double* d2B, double x1=1, double x2=1e6, double y1=1e-5,double y2=30,bool xLog=true, bool yLog=true)
{
  int N;
  for(N=0;M2[N]<x2;N++);
  TGraph gCross_1(N, M2, d2A); gCross_1.SetLineWidth(2);
  TGraph gCross_2(N, M2, d2B); gCross_2.SetLineWidth(2);
  gCross_2.SetLineColor(kBlue);  

//  gCross_1.Print();
  TLegend leg( 0.65, 0.75, 0.93, 0.9);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, strLeg[1], "l");
  leg.AddEntry(&gCross_2, strLeg[2], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, x1, x2, 1000,y1,y2);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  if(xLog)cCross_xx.GetPad(1)->SetLogx();
  if(yLog)cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})",
        "(M^{2})^{1.05}#frac{d^{2}#sigma}{dtdM^{2}} ", "");

  h_windowCross1_xx.Draw();
  gCross_2.Draw("C SAME");
  gCross_1.Draw("C SAME");

  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//*************************
void Plot_XXX_s()
{
  const int _Ns=20;
  double _s[_Ns]= { 100.000,   127.427,   162.378,   206.914,   263.665,   335.982,   428.133,   545.559,   695.193,   885.867,   1128.838,  1438.450,  1832.981,  2335.721,  2976.351,  3792.690,  4832.930,  6158.482,  7847.600,  10000.000 };

  double cs[_Ns]= {1.177179,  1.312843, 1.460185, 1.620322, 1.794501, 1.984116, 2.190740, 2.416140, 2.662330, 2.931591, 3.226532, 3.550147, 3.905870, 4.297658, 4.730114, 5.208538, 5.739134, 6.329109, 6.986921, 7.722473}; 
 
  TGraph gCross_1(_Ns, _s, cs); gCross_1.SetLineWidth(2);
  gCross_1.SetLineColor(kBlue);  

  TH2D h_windowCross1_xx("h_windowCross1", "", 10*_Ns, 1e2, 1e4, 1000,1,100);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
  char aa[20];
  const char* label="DD";
  sprintf(aa,"#sigma_{%s} (mb)",label);
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx,"#sqrt{s} (GeV)",aa,"");

  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");

//  const int Nx=13;
//  double ss2[Nx]={6e2,7e2,1.5e3,4e3,3e5,3e5,3.5e6,8e5,1e6,3.5e6,3.5e6,9e6, 49e6};
//  double ss[Nx]; for(int i=0;i<Nx;i++)ss[i]=sqrt(ss2[i]);
//  double css[Nx] ={ 5,6.5,  7.5,7.5,5.5,9.5,11.75,7.9,10, 8    ,11.75,13.5,16.7};

//  TGraph gx(Nx, ss, css);
//  gx.SetMarkerStyle(20);  
//  gx.SetMarkerSize(1);  
//  gx.SetMarkerColor(kBlack);  
//  gx.SetLineColor(kBlack);  
//  gx.Draw("CP SAME");
 
//  cCross_xx.Print("cs_s.root");
  if(_PNG_)cCross_xx.Print("cs_s.png");
  cCross_xx.Print("cs_s.eps");
}

//*************************
void Plot_VIcs_sII(int _Ns, double *_s_, double* cs1, double* cs2, string* leg, double _s_min, double _s_max, double ymin, double ymax, const char* label)//,const char*file)
{
//*************************************
 int aType=-1; //1-SD, 2-DD;
 string a="";
 if((a+label[0]+label[1])=="SD") {aType=1;}
 else if((a+label[0]+label[1])=="DD") {aType=2;}
 //printf("%s,%d\n",label,aType);

//   printf("We here\n");
  int pass=0;
  if(aType==1)for(pass=0;_s_[pass]<130;pass++); // _s_ -- [GeV]! not [GeV2] 
  TGraph gCross_1(_Ns, _s_, cs1); gCross_1.SetLineWidth(2);
//  printf(">>pass:%3d, type:%2d ,%p\n",pass,aType,cs2+pass);
  TGraph gCross_2(_Ns-pass, _s_+pass, cs2+pass); gCross_2.SetLineWidth(2);
//   gCross_1.Print();
//   gCross_2.Print();
  gCross_1.SetLineColor(1);//kBlue);  
  gCross_2.SetLineColor(1);//kBlue);  
  gCross_2.SetLineStyle(7);  

 // printf("Ok\n");

  TH2D h_windowCross1_xx("h_windowCross1", "", 10*_Ns, _s_min,_s_max, 1000,ymin,ymax);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
//  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
  char aa[100];
  if(aType==1) sprintf(aa,"2#sigma_{%s} (mb)",label);
  if(aType==2) sprintf(aa,"#sigma_{%s} (mb)",label);
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx,"#sqrt{s} (GeV)",aa,"");

  h_windowCross1_xx.GetYaxis()->SetTitleOffset(1.3);
  h_windowCross1_xx.Draw();

//!!!//
  gCross_1.Draw("C SAME");
  if(cs2!=NULL) gCross_2.Draw("C SAME");

//  const int Nx=13;
//  double sx[Nx]={6e2,7e2,1.5e3,4e3,3e5,3e5,3.5e6,8e5,1e6,3.5e6,3.5e6,9e6,49e6};
//  double Wx[Nx]; for(int i=0;i<Nx;i++)Wx[i]=sqrt(sx[i]);
//  double csx[Nx]={ 5,6.5,  7.5,7.5,5.5,9.5,11.75,7.9,10, 8   ,11.75,13.5,16.7};
//  TGraph gx(Nx, Wx, csx);   
  TGraphErrors gSD("../Data/SD.Data.dat", "%lg %lg %lg %lg");
  TGraphErrors gDD("../Data/DD.Data.dat", "%lg %lg %lg %lg");
  //SetGraphStyle(gSD,7, 1, kBlack);
  //SetGraphStyle(gDD,20, 1, kBlack);
  //SetGraphStyle(gx, 20, 1, kGreen);
  SetGraphStyle(gSD, 7, 1,1);// kBlue);
  SetGraphStyle(gDD,20, 1,1);// kBlue);

  TGraphErrors gAbeCDF("../Data/SD.sigma/Abe.CDF.1993.0.15s.dat", "%lg %lg %lg %lg");
  TGraphErrors gAlbrowISR("../Data/SD.sigma/Albrow.CHLM_ISR_CERN.1976.dat", "%lg %lg %lg %lg");
  TGraphErrors gALICE005("../Data/SD.sigma/ALICE_prelim.005.dat", "%lg %lg %lg %lg");
  TGraphErrors gALICE200("../Data/SD.sigma/ALICE_prelim.M200.dat", "%lg %lg %lg %lg");
  TGraphErrors gAmosE710("../Data/SD.sigma/Amos.E710.1992.dat", "%lg %lg %lg %lg");
  TGraphErrors gAnsorgeUA5("../Data/SD.sigma/Ansorge.UA5.1986.dat", "%lg %lg %lg %lg");
  TGraphErrors gArmitageISR("../Data/SD.sigma/Armitage.CHLM_ISR_CERN.1976.dat","%lg %lg %lg %lg");
  TGraphErrors gBernardUA4("../Data/SD.sigma/Bernard.UA4.1987.dat", "%lg %lg %lg %lg");
  TGraphErrors gCoolFNAL("../Data/SD.sigma/Cool.FERMILAB.1981.dat", "%lg %lg %lg %lg");
  TGraphErrors gSchambergerFNAL("../Data/SD.sigma/Schamberger.FERMILAB.1978.dat","%lg %lg %lg %lg");

  double *x=gALICE005.GetX();
  int nn=gALICE005.GetN();
  for(int i=0;i<nn;i++){x[i]*=0.98;}

  SetGraphStyle(gCoolFNAL       ,2,  1, kBlack);//kBlue-7);   //);//
  SetGraphStyle(gSchambergerFNAL,3,  1, kBlack);//kBlue-7);   //);//
  SetGraphStyle(gAlbrowISR      ,5,  1, kBlack);//kBlue);     //);//
  SetGraphStyle(gArmitageISR    ,28, 1, kBlack);//kCyan+2);   //);//

  SetGraphStyle(gAnsorgeUA5     ,26, 2, kBlack);//kGreen+2);  //);//
  SetGraphStyle(gBernardUA4     ,22, 2, kBlack);//kGreen+2);  //);//
  SetGraphStyle(gAmosE710       ,3,  2, kBlack);//kGreen-6);  //);//
  SetGraphStyle(gAbeCDF         ,24, 2, kBlack);//kOrange+3); //);//
  SetGraphStyle(gALICE200       ,21, 2, kBlack);//kRed);      //);//
  SetGraphStyle(gALICE005       ,25, 2, kBlack);//kRed-7);//, kMagenta);

  TGraph gGoul_st ("Data/csgSD1st.txt", "%lg %lg");
  TGraph gGoul_ren("Data/csgSD2ren.txt","%lg %lg");
  TGraph gGoul_GLM("Data/csgSD3GLM.txt","%lg %lg");

  gGoul_st.SetLineColor(kBlack);  gGoul_st.SetLineStyle(2);
  gGoul_ren.SetLineColor(kBlack); gGoul_ren.SetLineStyle(1);
  gGoul_GLM.SetLineColor(kBlack); gGoul_GLM.SetLineStyle(1);

  TLegend gLeg( 0.63, 0.21, 0.98, 0.395);
  gLeg.SetBorderSize(0); gLeg.SetFillColor(0);
  gLeg.AddEntry(&gCoolFNAL       ,"*Cool, FERMILAB (1981)",       "EP");
  gLeg.AddEntry(&gSchambergerFNAL,"*Schamberger,FERMILAB(1978)",  "EP");
  gLeg.AddEntry(&gAlbrowISR      ,"Albrow, CHLM(ISR)CERN(1976)",  "EP");
  gLeg.AddEntry(&gArmitageISR    ,"Armitage, CHLM(ISR) (1976)",   "EP");

  TLegend jLeg(0.22,0.55,0.65,0.94);
  jLeg.SetBorderSize(0); jLeg.SetFillColor(0);
//  TLegend Leg( 0.7, 0.85, 0.95, 0.95);
//  Leg.SetBorderSize(0); Leg.SetFillColor(0);
  jLeg.AddEntry(&gCross_1, leg[0].c_str(), "l");
  jLeg.AddEntry(&gCross_2, leg[1].c_str(), "l");
  jLeg.AddEntry(&gSD, " "," ");
  jLeg.AddEntry(&gAnsorgeUA5     ,"Ansorge, UA5(1986)|#xi<0.05",      "P");
  jLeg.AddEntry(&gBernardUA4     ,"Bernard, UA4 (1987)|#xi<0.05",     "P");
  jLeg.AddEntry(&gAmosE710       ,"Amos,E710(1992)|2<M2<0.05s",       "P");
  jLeg.AddEntry(&gAbeCDF         ,"Abe, CDF(1993) M2[1.4: 0.15s]",    "P"); 
  jLeg.AddEntry(&gALICE005       ,"ALICE(prelim) (extrapol #xi<0.05)","P");
  jLeg.AddEntry(&gALICE200       ,"ALICE(prelim) M<200",              "P");

  TGraphErrors gCDF_DD("../Data/DD.sigma/CDF.DD.dat", "%lg %lg %lg %lg");
  TGraphErrors gUA5_DD("../Data/DD.sigma/UA5.DD.dat", "%lg %lg %lg %lg");
  TGraphErrors gLow_DD("../Data/DD.sigma/LowEnergy.DD.dat", "%lg %lg %lg %lg");
  TGraphErrors gALICE_DD("../Data/DD.sigma/ALICE.DD.dat", "%lg %lg %lg %lg");
  double* xdd=gALICE_DD.GetX(); xdd[0]*=0.97;
  SetGraphStyle(gLow_DD,   3,  2,kBlack);// kBlue);     //);//
  SetGraphStyle(gUA5_DD,   23, 2,kBlack);// kBlue-7);   //);//
  SetGraphStyle(gCDF_DD,   24, 2,kBlack);// kOrange+10);//);//
  SetGraphStyle(gALICE_DD ,21, 2,kBlack);// kRed);      //);//
  TLegend dLeg( 0.22, 0.65, 0.65, 0.94);
  dLeg.SetBorderSize(0); dLeg.SetFillColor(0);
  dLeg.AddEntry(&gCross_1, leg[0].c_str(), "l");
  dLeg.AddEntry(&gCross_1, " ", " ");
  dLeg.AddEntry(&gLow_DD,   "Low energy data", "P");
  dLeg.AddEntry(&gUA5_DD,   "UA5",             "P");
  dLeg.AddEntry(&gCDF_DD,   "CDF",             "P");
  dLeg.AddEntry(&gALICE_DD ,"ALICE",           "P");


//  double x1,y1,x2,y2;
//  if(cs2==NULL) {x1=0.30; y1=0.75; x2=0.50; y2=0.85;} 
//           else {x1=0.22; y1=0.8; x2=0.58; y2=0.95;} 
//  TLegend xLeg(x1,y1,x2,y2);
//  xLeg.SetBorderSize(0); xLeg.SetFillColor(0);
////  TLegend Leg( 0.7, 0.85, 0.95, 0.95);
////  Leg.SetBorderSize(0); Leg.SetFillColor(0);
//  xLeg.AddEntry(&gCross_1, leg[0].c_str(), "l");
//  if(cs2!=NULL) xLeg.AddEntry(&gCross_2, leg[1].c_str(), "l");

  string suf;
  if(aType==1)//SD 
   {
    gCoolFNAL.Draw("P SAME");
    gSchambergerFNAL.Draw("P SAME");
    gArmitageISR.Draw("P SAME");
    gAlbrowISR.Draw("P SAME");
    gALICE005.Draw("P SAME");
    gALICE200.Draw("P SAME");
    gAbeCDF.Draw("P SAME");
    gAmosE710.Draw("P SAME");
    gAnsorgeUA5.Draw("P SAME");
    gBernardUA4.Draw("P SAME");
//    gSD.Draw("P SAME");

    gLeg.Draw("SAME");
    jLeg.Draw("SAME");

//    gGoul_st.Draw("L,SAME");   
//    gGoul_ren.Draw("L,SAME"); 
//    gGoul_GLM.Draw("L,SAME"); 
   }
  if(aType==2)//"DD"
   {
//    gDD.Draw("P SAME");
    gCDF_DD.Draw("P SAME");
    gUA5_DD.Draw("P SAME");
    gLow_DD.Draw("P SAME");
    gALICE_DD.Draw("P SAME");
 //   xLeg.Draw("SAME");
    dLeg.Draw("SAME");
    
   }

  a="cs_s.";
//  cCross_xx.Print((a+label+".root").c_str());
  cCross_xx.Print((a+label+".eps").c_str());
  if(_PNG_) cCross_xx.Print((a+label+".png").c_str());
}

void SetGraphStyle(TGraph& gr, int mkSyle,int mkSize, int color,int mkcolor)
{ gr.SetMarkerStyle(mkSyle);    
  gr.SetMarkerSize(mkSize);      
  gr.SetMarkerColor(color);
  if(mkcolor<0) gr.SetLineColor(color);
           else gr.SetLineColor(mkcolor);
}

//*************************
void Plot_VIcs_s(int _Ns, double *_s_, double* cs, double _s_min, double _s_max, double ymin, double ymax, const char* label)//,const char*file)
{
//   printf("We here\n");
  TGraph gCross_1(_Ns, _s_, cs); gCross_1.SetLineWidth(2);
  gCross_1.SetLineColor(kBlue);  

  TH2D h_windowCross1_xx("h_windowCross1", "", 10*_Ns, _s_min,_s_max, 1000,ymin,ymax);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
//  cCross_xx.GetPad(1)->SetLogy();
  cCross_xx.GetPad(1)->SetLogx();
  char aa[100]; sprintf(aa,"#sigma_{%s} (mb)",label);
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx,"#sqrt{s} (GeV)",aa,"");

  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");

//  const int Nx=13;
//  double sx[Nx]={6e2,7e2,1.5e3,4e3,3e5,3e5,3.5e6,8e5,1e6,3.5e6,3.5e6,9e6,49e6};
//  double Wx[Nx]; for(int i=0;i<Nx;i++)Wx[i]=sqrt(sx[i]);
//  double csx[Nx]={ 5,6.5,  7.5,7.5,5.5,9.5,11.75,7.9,10, 8   ,11.75,13.5,16.7};
//  TGraph gx(Nx, Wx, csx);   
  TGraphErrors gSD("../Data/SD.Data.dat", "%lg %lg %lg %lg");
  TGraphErrors gDD("../Data/DD.Data.dat", "%lg %lg %lg %lg");

  gSD.SetMarkerStyle(20);     gDD.SetMarkerStyle(20);   //gx.SetMarkerStyle(20);     
  gSD.SetMarkerSize(1);       gDD.SetMarkerSize(1);     //gx.SetMarkerSize(1);       
  gSD.SetMarkerColor(kBlack); gDD.SetMarkerColor(kBlack);//gx.SetMarkerColor(kGreen); 
  gSD.SetLineColor  (kBlack); gDD.SetLineColor  (kBlack);//gx.SetLineColor  (kGreen); 

  string a="",suf;
  if((a+label[0]+label[1])=="SD") gSD.Draw("P SAME");
  if((a+label[0]+label[1])=="DD") gDD.Draw("P SAME");

  a="cs_s";
//  cCross_xx.Print((a+label+".root").c_str());
  cCross_xx.Print((a+label+".eps").c_str());
  if(_PNG_) cCross_xx.Print((a+label+".png").c_str());
}

//*************************
/*
void Plot_imal(double _t_)
{
  TGraph gCross_1 (Nm2, M2, Im); gCross_1.SetLineWidth(2);
  TString s1; s1.Form("original t = %6.3f",_t_);

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.001, 10., 10000, 0,1);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  cCross_xx.GetPad(1)->SetGrid();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})","Im #alpha(M^{2})", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");

  cCross_xx.Print("Im.eps");
}
*/
//*************************
/*
void Plot_imA()
{
  TGraph gCross_1 (Nm2, M2, ImA[0]); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nm2, M2, ImA[1]); gCross_2.SetLineWidth(2);
  TGraph gCross_3 (Nm2, M2, ImA[2]); gCross_3.SetLineWidth(2);
  gCross_1.SetLineColor(kRed);  
  gCross_2.SetLineColor(kBlue);  
  gCross_3.SetLineColor(kGreen);  

  TLegend leg( 0.7, 0.8, 0.93, 0.9);
  leg.SetBorderSize(0); leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, "Bw[0]", "l");
  leg.AddEntry(&gCross_2, "Bw[1]", "l");
  leg.AddEntry(&gCross_3, "Bw[2]", "l");

 // TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0.001, 10., 10000, 0,5);
  TH2D h_windowCross1_xx("h_windowCross1", "", 10*Nm2, 1, 15, 1000,5e-4,5);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  cCross_xx.GetPad(1)->SetGrid();
//  cCross_xx.GetPad(1)->SetLogy();
//  cCross_xx.GetPad(1)->SetLogx();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "M^{2} (GeV^{2})","Bw(i)", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
  gCross_3.Draw("C SAME");
  leg.Draw("SAME");
  cCross_xx.Print("Bw.eps");
}
*/
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//      DD   
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//### (int|M2) dsg/dt (t) (DD) ####  
void PlotDD_Iint_t (const int Nt, double* t, double* dcs1,double* dcs2, TString* stLeg, const char* fout)
{ 
  double _t[Nt]; for(it=0;it<Nt;it++) _t[it]=-t[it];
  TGraph gCross_1 (Nt, _t, dcs1); gCross_1.SetLineWidth(2);
  TGraph gCross_2 (Nt, _t, dcs2); gCross_2.SetLineWidth(2);
  gCross_1.SetLineStyle(1);  
  gCross_2.SetLineStyle(7);  
//  gCross_1.SetLineColor(kBlue);  
//  gCross_2.SetLineColor(kBlack);  
  /////////////////////////
  //        DRAW!
  ////////////////////////

  TLegend leg( 0.56, 0.50, 0.94, 0.65);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.AddEntry(&gCross_1, stLeg[0], "");
  leg.AddEntry(&gCross_1, stLeg[1], "l");
  leg.AddEntry(&gCross_2, stLeg[2], "l");

  TH2D h_windowCross1_xx("h_windowCross1", "", 1000, 0, 1, 100,1e-3,5e1);
  TCanvas cCross_xx("cCross", "cCross", 800, 800);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
  cCross_xx.GetPad(1)->cd();
  //cCross_xx.GetPad(1)->SetGrid();
  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1), &h_windowCross1_xx, "t (GeV^{2})","#frac{d#sigma_{DD}}{dt}  #(){#frac{mb}{GeV^{2}}}", "");
  h_windowCross1_xx.Draw();
  gCross_1.Draw("C SAME");
  gCross_2.Draw("C SAME");
//  gCross_3.Draw("C SAME");
  leg.Draw("SAME");

  cCross_xx.Print(fout);
}

//### (int|t) d2sg/dM1dM2  (DD) 3D-graph ####  
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

void Plot3D_DD_int_t (const int Nm2, double* M2, double* d2cs, const char* fout, const char* Ztitle)
{ 
//  printf("Start plotter\n");
  double M2arr[Nm2+1];
  double delta=(M2[Nm2-1]-M2[0])/(Nm2-1);
  double xl=M2[0]-delta, xu=M2[Nm2-1]+delta;
  TString str();
  TH2D h2D("h2D","Plot;M1;M2;#frac{d#sigma_{DD}}{dt}  #(){#frac{mb}{GeV^{2}}}",Nm2,xl,xu,Nm2,xl,xu); 
  h2D.SetZTitle(Ztitle); 
//  TGraph2D gGr2D(Nm2*Nm2);
  int id=-1;
  int im1,im2;
  for(im1=0;im1<Nm2;im1++)
   for(im2=0;im2<Nm2;im2++)
    {id++; h2D.SetBinContent(im1,im2,d2cs[id]);}
  
  TCanvas cCross_xx("cCross", "cCross", 1400, 1000);
  cCross_xx.Divide(1, 1);
  cCross_xx.Draw();

  make_clean_pads(&cCross_xx, 1, 0, 0);
//  cCross_xx.GetPad(1)->cd();
//  cCross_xx.GetPad(1)->SetGrid();
//  cCross_xx.GetPad(1)->SetLogy();
  sign_window(cCross_xx.GetPad(1),&h2D,"M_{1}^{2} (GeV^{2})","M_{2}^{2} (GeV^{2})","");
  h2D.GetXaxis()->SetTitleOffset(2.);
  h2D.GetYaxis()->SetTitleOffset(2.);
  h2D.GetZaxis()->SetTitleOffset(1.4);
//  cCross_xx.GetPad(1)->SetLeftMargin(50);
//  cCross_xx.GetPad(1)->SetRightMargin(50);

//  h_windowCross1_xx.Draw();
 // gGr2D.Draw("A surf1, Z");
// const UInt_t Number = 3;
// Double_t Red[Number]    = { 1.00, 0.00, 0.00};
// Double_t Green[Number]  = { 0.00, 1.00, 0.00};
// Double_t Blue[Number]   = { 0.00, 0.00, 1.00};
// Double_t Length[Number] = { 0.10, 0.10, 0.10 };
// Int_t nb=10;
// TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
//TColor::SetGrayscale( kTRUE);
// gStyle->CreateColorsGray();
// gStyle->SetPalette(53,0);  

 int N=10;
 TColor* Col[N];
 int colors[N];
 double d=1./(N+2), d0=1-d;
// for(int i=0;i<10;i++) Col[i]=new TColor(255+i,d0-d*i,d0-d*i,d0-d*i);//colors[i]= 
// gStyle->SetPalette(N,colors);  

 double Zmin=0.001;
 double Zmax=h2D.GetMaximum();
// printf("Min:%f, Max:%f\n",Zmin,Zmax);
 double x1=log(Zmin);
 double x2=log(Zmax);
 double dx=(x2-x1)/(N-1);
 double cont[N];
 printf("Counurs: ");
 for(int i=0;i<N;i++){cont[i]=exp(x1+i*dx); printf("%5.3f ",cont[i]);}printf("\n");
 cont[0]=0;
 h2D.SetContour(N-1, cont);
 h2D.SetLineWidth(0.1);
 h2D.SetLineColor(kBlack);

//########################################
//#### Set Gray Pallette #####
 const int Ncol=9;                                                               
 int stColor=101, curColor;                                                                 
 double grScale;
 int palette[Ncol];
 TColor *color;
                                                                                  
 for(int i=0; i<Ncol; i++)
  {
   curColor= stColor+i;
   grScale = i*1./Ncol;
   if(i==6) {grScale=grScale+0.5/Ncol;} 
   if(i==7) {grScale=grScale+1.2/Ncol;} 
   if(i+1==Ncol) {grScale=grScale/2.;} 
  

   if(!(color=gROOT->GetColor(curColor)))                                        
    {color = new TColor(curColor, grScale,grScale,grScale, "");}
   else
    {color->SetRGB(grScale,grScale,grScale);}
   palette[Ncol-i]=curColor;                                                          
  }
 
 gStyle->SetPalette(Ncol, palette); 

// for(int i=0;i<3;i++)
// for(int j=0;j<1;j++)
//  h2D.SetBinContent(i,Nm2-j,Zmax);
  
 // h2D.Draw("Surf1, Z");

  h2D.Draw("Surf3");
  h2D.Draw("surf1 Z,Same");

  //gGr2D.Draw("Surf1, Z");
  //leg.Draw("SAME");
//  printf("Hallo1 (printig)\n");
//  cCross_xx.Print("a.root");
  //cCross_xx.Print("a.gif");
  cCross_xx.Print(fout);
//  printf("Hallo2\n");
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@ Fitting parameters @@@
 double   Ael_,bel_,al0_,al_,Ares_,Cbg_,R_,b_res_,b_bg_,s0_,eta_;
 double  dAel,dbel,dal0,dal,dAres,dCbg,dR,db_res,db_bg,ds0,deta;

//~~ elastic pp ~~~~~~~~~~~~~~
double dcsdt_el(double t, double s) //4-param +(s0)
 {return  Ael_*exp(4*bel_*t) *pow(s/s0_ ,2*(al0_+al_*t) ); } 

//~~ SD ~~~~~~~~~~~~~~
double Sg_res(double t, double M2);
double Sg_bg(double t, double M2);
double d2cs_SD(double t, double M2, double s)//single arm 
{ return 0.5* exp(2*bel_*t)
          *(Ares_*kinematic_factor(t,M2)*Sg_res(t,M2) + Cbg_*Sg_bg(t,M2))
          *pow(s/M2, 2*(al0_+al_*t) );
}

double im_,re_,sg; int jj;
double Sg_res(double t, double M2)
{sg=0;
 im_=Im_alpha(M2);
 re_=Re_alpha(M2);

 for(jj=1;jj<4;jj++) {sg+=im_/(sq(2*jj+0.5-re_) + sq(im_));}
 sg*=exp(4*b_res_*t); 
 sg+=R_*exp(2*b_res*t)*MRoper*WRoper/2./(sq(M2-MRoper2) + sq(WRoper/2));
 return sg;
}

double Sg_bg(double t, double M2)
{ if(M2>sqMpMpi) return exp(b_bg_*t)/(1/pow(M2-sqMpMpi,xi) +pow(M2, eta_));
  return 0; }
//~~ SD integ ~~~~~~~~~~~~~~~~~~~~~~~~~
double s_max=0.05*sq(7000);
const int xNt=200, xNm2=500;
//const int xNt=100, xNm2=300;
//const int xNt=100, xNm2=500;//canonical
//double xd2cs[xNt*xNm2]; 
double xM2[xNm2], xdm2[xNm2]={0};
double xIRes[3*xNm2], xIBg[xNm2], xIRoper[xNm2];

double xt[20];
double st[xNt], sdt[xNt]={0};
double x1f2[xNt],x2f2[xNt],x3f2[xNt],xFp2[xNt],xf_bg[xNt], xalpha_IP[xNt];


//~~~~~~~~~~~~~~~~
//
//double dcsdt_SD(double M2, double s) 
//{ 
//  double cs,t,dt_; int j;
//  int N=200;
//  dt_=1./N; 
//  t=dt*0.5;
//  for(jjj=0;jjj<N;jjj+++)
//   {cs;t+=dt;}
//  return 0.5* exp(2*bel_*t)
//          *(Ares_*kinematic_factor(t,M2)*Sg_res(t,M2) + Cbg_*Sg_bg(t,M2))
//          *pow(s/M2, 2*(al0_+al_*t) );
//}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@  Staff  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void SetErrGraph(TGraphErrors& g, DataContainer* dc)
{
 int i, N=dc[0].pages[0]->N;
 float *X=dc[0].pages[0]->X;
 float *Y=dc[0].pages[0]->Y;
 float*dY=dc[0].pages[0]->dYtotal;
 for(i=0;i<N;i++)
  {g.SetPoint(i,X[i],Y[i]);
   g.SetPointError(i,0,dY[i]);} 
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@  Fitting  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#include "TFitter.h"
  DataContainer DC7000pp[2];
  DataContainer DC1800pp[2];
  DataContainer DC546pp[4];
  DataContainer DC540pp;

  DataContainer d2csSD[3];
  DataContainer dcsdtSD;
  DataContainer csSD_m2, csSD_xi;
  DataContainer csDD;

//~~~ Local ~~~~
int ii,N;
float *X,*Y,*dY,P;
double iChi2;
//~~~~~~~~~~~~~~
double Chi2ppDC(DataContainer&dc, int& Nn)
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;
  P=dc.pages[0]->P;//(!)s
//  printf(">>%f\n",P); 
//  printf("Ael:%7.3f, bel:%6.3f, al0:%4.2f, al:%6.3f\n",Ael_, bel_,al0_,al_);
 iChi2=0;
 for(ii=0;ii<N;ii++)
  {iChi2+=sq( (dcsdt_el(-X[ii],sq(P))-Y[ii])/dY[ii] );
//  printf("Theor:%7.3f, Exp: %7.3f+/-%7.3f\n ",dcsdt_el(X[ii],sq(540)),Y[ii],dY[ii]);
}
 Nn+=N;
 return iChi2; 
}
//~~~~~~~~~~~~~~
double Chi2_d2SD_DC(DataContainer&dc,double t, int& Nn)// (t>0)
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;
  P=dc.pages[0]->P;//(!)s
 iChi2=0;
 for(ii=0;ii<N;ii++)
  {iChi2+=sq((d2cs_SD(-t, X[ii],sq(P))-Y[ii])/dY[ii] );}
 Nn+=N;
 return iChi2; 
}
//~~~~~~~~~~~~~~

double Chi2_dcsdtSD_DC(DataContainer&dc, int& Nn)// (t>0)
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;
  P=dc.pages[0]->P;//(!)s
 double s_lim=sq(546)*0.05;
 
 double dcsdt; 
 int Nx=0;

 iChi2=0;
// printf("\n#################\n");
 for(ii=0;ii<N;ii++)
  if(X[ii]<0.6) // X==t
    {dcsdt=0;Nx++;
     for(int i=0;xM2[i]<s_lim;i++)
        dcsdt+=d2cs_SD(-X[ii],xM2[i],sq(P))*xdm2[i];
//        printf("%5.2f  %7.3f\n",X[ii],dcsdt);
     iChi2+=sq((dcsdt-Y[ii])/dY[ii]);}
 Nn+=Nx;
 return iChi2; 
}

double Chi2_dcsdtSD_shDC(DataContainer&dc, int& Nn)
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;
  P=dc.pages[0]->P;//(!)s
 double s=sq(546);
 double s_lim=s*0.05;
// t -dependence
 for(ii=0;ii<N;ii++)
  {
    xalpha_IP[ii]=2*(al0_-al_*X[ii]);
    x1f2[ii]  = exp(-4*b_res_*X[ii]); 
    xFp2[ii]  = exp(-2*bel_*X[ii]);
    xf_bg[ii] = exp(-b_bg_*X[ii]);
  }
// printf("Hello!\n");
 
 double dcsdt,d2cs; 
 double Kf, propag, Sg_b02, Bg_b02;
 int Nx=0;

 iChi2=0;
// printf("\n#################\n");
 for(ii=0;ii<N;ii++)
  if(X[ii]<0.6)
    {dcsdt=0;Nx++;
     for(int i=0;xM2[i]<s_lim;i++)
       {
         propag = pow(s/xM2[i], 2*(al0_-al_*X[ii]));           
         Kf     = kinematic_factor(-X[ii], xM2[i]);
         Sg_b02 = x1f2[ii]*xIRes[i]
                 +x1f2[ii]*xIRes[i+xNm2]
                 +x1f2[ii]*xIRes[i+2*xNm2] 
                  +R_*sqrt(x1f2[ii])*xIRoper[i];   
         if(xM2[i]>sqMpMpi)
           Bg_b02 = xf_bg[ii]/(1/pow(xM2[i]-sqMpMpi,xi) +pow(xM2[i], eta_));
         else Bg_b02=0;
         d2cs = xFp2[ii]*(Ares_*Sg_b02*Kf+Cbg_*Bg_b02)*propag;
 
         printf("%5.2f  %7.3f\n",xM2[ii],d2cs);
         dcsdt+=d2cs *xdm2[i];
       }
      dcsdt*=0.5;// single arm;
//      printf("%5.2f  %7.3f\n",X[ii],dcsdt);
     iChi2+=sq((dcsdt-Y[ii])/dY[ii]);}
 Nn+=Nx;
 return iChi2; 
}

//### cs(s) SD ###
double Chi2_csSD_DC(DataContainer&dc, int& Nn,int type=1)// 1-- xi<0.05; 2 -- M<200
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;//(sqrt(s))
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;

 double s,s_lim;

// t -dependence
 for(ii=0;ii<xNt;ii++)
  {
    xalpha_IP[ii]=2*(al0_+al_*st[ii]);
    x1f2[ii]  = exp(4*b_res_*st[ii]); 
    xFp2[ii]  = exp(2*bel_*st[ii]);
    xf_bg[ii] = exp(b_bg_*st[ii]);
  }
// printf("Hello!\n");
 
 double dcsdt,d2cs,cs; 
 double Kf, propag, Sg_b02, Bg_b02;
 int i,Nx=0;
 iChi2=0;
// printf("\n#################\n");
 for(int is=0; is<N; is++)
  {
    s=sq(X[is]);
    if(type==1) s_lim=0.05*s;
    if(type==2) s_lim=sq(200);
     
    cs=0;
    for(ii=0;ii<xNt;ii++)// sum_t
       {dcsdt=0;
        for(i=0;xM2[i]<s_lim;i++)//sum_M2
          {
            propag = pow(s/xM2[i], 2*(al0_+al_*st[ii]));           
            Kf     = kinematic_factor(st[ii], xM2[i]);
            Sg_b02 = x1f2[ii]*xIRes[i]
                    +x1f2[ii]*xIRes[i+xNm2]
                    +x1f2[ii]*xIRes[i+2*xNm2] 
                     +R_*sqrt(x1f2[ii])*xIRoper[i];   
            if(xM2[i]>sqMpMpi)
              Bg_b02 = xf_bg[ii]/(1/pow(xM2[i]-sqMpMpi,xi) +pow(xM2[i], eta_));
            else Bg_b02=0;
            d2cs = xFp2[ii]*(Ares_*Sg_b02*Kf+Cbg_*Bg_b02)*propag;
    
            dcsdt+=d2cs *xdm2[i];
          }
//         dcsdt*=0.5;// single arm;
         cs+=dcsdt*sdt[ii];
       }
    iChi2+=sq((cs-Y[is])/dY[is]); Nx++;
  }
 Nn+=Nx;
 return iChi2; 
}

//### cs(s) DD ###
double Chi2_csDD_DC(DataContainer&dc, int& Nn)
{
  N=dc.pages[0]->N;
  X=dc.pages[0]->X;//(sqrt(s))
  Y=dc.pages[0]->Y;
 dY=dc.pages[0]->dYtotal;
// printf(">>%2d| %6.3f. %6.3f %6.3f\n",N,X[0],Y[0],dY[0]); 
// t -dependence
//   printf("xalpha_IP[ii],x1f2[ii],xFp2[ii],xf_bg[ii]\n");
 for(ii=0;ii<xNt;ii++)
  {
    xalpha_IP[ii]=2*(al0_+al_*st[ii]);
    x1f2[ii]  = exp(4*b_res_*st[ii]); 
    xFp2[ii]  = exp(2*bel_*st[ii]);
    xf_bg[ii] = exp(b_bg_*st[ii]);
//   printf("%7.3f|%7.3f %7.3f %7.3f %7.3f \n",st[ii],xalpha_IP[ii],x1f2[ii],xFp2[ii],xf_bg[ii]);
  }
// printf("Hello!\n");
 
 double ss0,ss0_lim;
 double dcsdt,d2cs,cs; 
 double Kf1,Kf2, propag, Sg1_b02, Sg2_b02, Bg1_b02, Bg2_b02;
 int i1,i2,Nx=0;
 iChi2=0;
// printf("\n#################\n");
 for(int is=0; is<N; is++)
  if(X[is]>200)
  {
//    printf(">>xNt:%3d, xNm2:%3d\n",xNt,xNm2);
    ss0=s0_*sq(X[is]);
    ss0_lim=ss0/exp(3.);
    cs=0;
    for(ii=0;ii<xNt;ii++)// sum_t
     {
        dcsdt=0;
        for(i1=0;i1<xNm2;i1++)//sum_M1
        for(i2=i1+1;i2<xNm2;i2++)//sum_M2
         if(xM2[i1]*xM2[i2]<ss0_lim)
          {
            propag = pow(ss0/(xM2[i1]*xM2[i2]), 2*(al0_+al_*st[ii]));           
            Kf1    = kinematic_factor(st[ii], xM2[i1]);
            Kf2    = kinematic_factor(st[ii], xM2[i2]);
            Sg1_b02 = x1f2[ii]*xIRes[i1] +x1f2[ii]*xIRes[i1+xNm2]
                     +x1f2[ii]*xIRes[i1+2*xNm2] +R_*sqrt(x1f2[ii])*xIRoper[i1];   
            Sg2_b02 = x1f2[ii]*xIRes[i2] +x1f2[ii]*xIRes[i2+xNm2]
                     +x1f2[ii]*xIRes[i2+2*xNm2] +R_*sqrt(x1f2[ii])*xIRoper[i2];   

            if(xM2[i1]>sqMpMpi) 
              Bg1_b02 = xf_bg[ii]/(1/pow(xM2[i1]-sqMpMpi,xi) +pow(xM2[i1], eta_));
            else Bg1_b02=0;
            if(xM2[i2]>sqMpMpi)
              Bg1_b02 = xf_bg[ii]/(1/pow(xM2[i2]-sqMpMpi,xi) +pow(xM2[i2], eta_));
            else Bg1_b02=0;

            d2cs=(Ares_*Sg1_b02*Kf1+Cbg_*Bg1_b02)*(Ares_*Sg2_b02*Kf2+Cbg_*Bg2_b02)*propag;
            dcsdt+=d2cs *xdm2[i1]*xdm2[i2];
 
  //          printf("#[t:%6.3f m1:%6.3f m2:%6.3f] Ares_:%6.3f,Cbg_:%6.3f| d2cs:%6.3f, Vtx1:%6.3f, propag:%6.3f\n",st[ii],xM2[i1],xM2[i2],Ares_,Cbg_ ,d2cs,(Ares_*Sg1_b02*Kf1+Cbg_*Bg1_b02),propag);
         //   printf("[t:%6.3f m1:%6.3f m2:%6.3f] propag:%6.3f| Kf1:%6.3f,Kf2:%6.3f; Sg1:%6.3f,Sg2:%6.3f;  Bg1:%6.3f,Bg2:%6.3f;\n",st[ii],xM2[i1],xM2[i2] ,propag,Kf1,Kf2,Sg1_b02,Sg2_b02, Bg1_b02,Bg1_b02);
          }
         cs+=2*dcsdt*sdt[ii];

        dcsdt=0;
        for(i1=0;i1<xNm2;i1++)//sum_M1
         if(sq(xM2[i1])<ss0_lim)
          {
            propag = pow(ss0/(xM2[i1]*xM2[i2]), 2*(al0_+al_*st[ii]));           
            Kf1    = kinematic_factor(st[ii], xM2[i1]);
            Sg1_b02 = x1f2[ii]*xIRes[i1] +x1f2[ii]*xIRes[i1+xNm2]
                     +x1f2[ii]*xIRes[i1+2*xNm2] +R_*sqrt(x1f2[ii])*xIRoper[i1];   
            if(xM2[i1]>sqMpMpi) 
              Bg1_b02 = xf_bg[ii]/(1/pow(xM2[i1]-sqMpMpi,xi) +pow(xM2[i1], eta_));
            else Bg1_b02=0;

            d2cs =sq(Ares_*Sg1_b02*Kf1+Cbg_*Bg1_b02)*propag;
            dcsdt+=d2cs *sq(xdm2[i1]);
          }
         cs+=dcsdt*sdt[ii];
     }
      cs*=1./(4*Ael_);
//    cs*=0.5;// single arm;
    cs*=2;  // two  arms;
//    printf(">>%5.1f| t:%7.3f; exp:%7.3f %7.3f\n",X[is],cs,Y[is],dY[is]);
    iChi2+=sq((cs-Y[is])/dY[is]); Nx++;
  }
 Nn+=Nx;
 return iChi2; 
}

//  printf("%5.2f  %7.3f\n",X[ii],dcsdt);

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
double ippChi2;
int Npp;
double Chi2pp()
{ Npp=0;   
  ippChi2=0;
// printf("##Chi2:%6.3f\n",iChi2);
  ippChi2 =
        Chi2ppDC(DC7000pp[0],Npp)
        +Chi2ppDC(DC7000pp[1],Npp);
     //  +Chi2ppDC(DC1800pp[0],Npp)
     //  +Chi2ppDC(DC1800pp[1],Npp)
//   //    +Chi2ppDC(DC546pp[0],Npp)
     //  +Chi2ppDC(DC546pp[1],Npp)
     //  +Chi2ppDC(DC546pp[2],Npp)
     //  +Chi2ppDC(DC540pp,Npp);
 return  ippChi2/Npp;
}

double Chi2_d2SD()
{ Npp=0;   
  ippChi2 =
        Chi2_d2SD_DC(d2csSD[0],0.05,Npp)
       +Chi2_d2SD_DC(d2csSD[1],0.05,Npp)
       +0.25*Chi2_d2SD_DC(d2csSD[2],0.55,Npp);
 return  ippChi2/Npp;
}

double Chi2_dcsdtSD()
{ Npp=0;   
//  ippChi2 = Chi2_dcsdtSD_DC(dcsdtSD,Npp); //41.s
  ippChi2 = Chi2_dcsdtSD_shDC(dcsdtSD,Npp); //9.7s
 return  ippChi2/Npp;
}

double Chi2_csSD()
{ Npp=0;   
  ippChi2 = Chi2_csSD_DC(csSD_xi,Npp,1); // 1-- xi<0.05; 2 -- M<200
  ippChi2 = Chi2_csSD_DC(csSD_m2,Npp,2); // 1-- xi<0.05; 2 -- M<200
 return  ippChi2/Npp;
}

double Chi2_csDD()
{ Npp=0;   
  ippChi2 = Chi2_csDD_DC(csDD,Npp);
//  printf(">>Chi2DD%7.3f; N:%2d \n",ippChi2,Npp);
 return  ippChi2/Npp;
}

double Chi2full()
 {return Chi2pp()+Chi2_d2SD()+Chi2_dcsdtSD()
        +5*(Chi2_csSD()+Chi2_csDD());}

void minuitFunction1(int& nDim, double* gout, double& result, double* par, int flg)
 { Ael_=par[0];bel_=par[1];al0_=par[2];al_=par[3];s0_=par[9];eta_=par[10];  
  Ares_=par[4];Cbg_=par[5]; R_ =par[6];b_res_=par[7];b_bg_=par[8];
  result=Chi2full();}

 void GetFitPar(int i, double& x, double& dx, TFitter* m)
  {x=m->GetParameter(i); dx=m->GetParError(i);}

//??? Tests  ????????????????????????????????
void test_el (double*t, int Nt )
{
 printf("%6s %6s (II)\n","t","dcsdt"); 
 for(int it=0;it<Nt;it++) {printf("%6.3f %6.3f\n", t[it],dcsdt_el(t[it], 1800*1800));}
}

void test_SD_Goul (double* M2) {
 printf("%6s %6s %6s (SD_Goul:II)\n","t","M2","dcsdt");
  for(int im2=0;im2<50;im2++)
   {printf("%6.3f %6.3f %6.3f\n", -0.05,M2[im2],d2cs_SD(-0.05,M2[im2],1800*1800));}
}

void test_SD(double* t,double* M2) {
  printf("%6s %6s %6s (SD:II)\n","t","M2","dcsdt");
  for(int it=0;it<50;it++)
  for(int im2=0;im2<50;im2++)
   {printf("%6.3f %6.3f %6.3f\n", t[it],M2[im2],d2cs_SD(t[it],M2[im2],7000*7000));}
}

void test_SD_I() {
  double s=1800*1800;
  const int Nt=10, Nm2=10;
  double d2cs[Nt*Nm2]; 
  double  t[Nt], M2[Nm2];
  double  dt[Nt], dM2[Nm2];
  double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
  double  ImA[3*Nm2], ImBg[Nm2], ImR[Nm2];

  prestage_tlin(t, dt, Nt, 0,1);
  prestage_M2(M2, dM2, Nm2, 1.45, 10,"log");
  stageI_M2(ImA,ImBg,ImR,M2,Nm2);
  stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, t, Nt);

  printf("(SD:I)\n");
  stage_IIcernel(d2cs, s,t,M2, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP, ImA,ImBg,ImR, Ares,Cbg);

  double Kf,Sg_b02,Bg_b02,d2cs_;
  printf("(SD:II)\n");
 printf("%6s %6s %6s %6s %6s %6s\n", "t","M2", "Kf","Sg_b02","Bg_b02","d2cs");
  for(int it=0;it<Nt;it++)
  for(int im2=0;im2<Nm2;im2++)
   { Kf  = kinematic_factor(t[it],M2[im2]);
     Sg_b02= Sg_res(t[it],M2[im2]);
     Bg_b02= Sg_bg(t[it],M2[im2]);
     d2cs_=d2cs_SD(t[it],M2[im2],1800*1800); 
     
     printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
      t[it],M2[im2], Kf,Sg_b02,Bg_b02, d2cs_);}
}

//~~~~~~~~~~~~~~~~~
 void Minimize()
  {
 //   Ael_=33.579; bel_=1.937; al0_=0.075; al_=0.34; Ares_=3.89; Cbg_=2.07; R_=0.45;
    Ael_=33.579; bel_=1.937; al0_=0.075; al_=0.34; Ares_=1.38; Cbg_=2.07; R_=0.45;
   dAel =0.01;   dbel=0.001; dal0=0.001; dal=0.001;dAres=0.001;dCbg=0.001;dR=0.001;

    b_res_=-0.507; b_bg_=-1.013; s0_=1;   eta_=1;
   db_res =0.001;  db_bg = 0.001; ds0 =0.01;deta=0.001;
   TFitter* minimizer =new TFitter(11);
   {double p1=-1;minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);}
 
   minimizer->SetFCN(minuitFunction1);
   minimizer->SetParameter(0,"Ael",  Ael_,  dAel,   0,0); //minimizer->FixParameter(0);
   minimizer->SetParameter(1,"bel",  bel_,  dbel,   0,0); //minimizer->FixParameter(1);
   minimizer->SetParameter(2,"al0",  al0_,  dal0,   0,0); //minimizer->FixParameter(2);
   minimizer->SetParameter(3,"al",   al_,   dal,    0,99);//minimizer->FixParameter(3);

   minimizer->SetParameter(4,"Ares", Ares_, dAres,  0,0); minimizer->FixParameter(4);
   minimizer->SetParameter(5,"Cbg",  Cbg_,  dCbg,   0,99);minimizer->FixParameter(5);
   minimizer->SetParameter(6,"R"   , R_,    dR,     0,0); minimizer->FixParameter(6);
   minimizer->SetParameter(7,"b_res",b_res_,db_res, 0,99);minimizer->FixParameter(7);
   minimizer->SetParameter(8,"b_bg", b_bg_, db_bg,  0,0); minimizer->FixParameter(8);

   minimizer->SetParameter(9,"s0",   s0_,   ds0,    0,5); minimizer->FixParameter(9);
   minimizer->SetParameter(10,"eta", eta_,  deta,   0,40);minimizer->FixParameter(10);
//   minimizer->ExecuteCommand("SIMPLEX",0,0);
 //  minimizer->ExecuteCommand("MIGRAD", 0,0);
 
   GetFitPar(0,  Ael_,  dAel,  minimizer);  
   GetFitPar(1,  bel_,  dbel,  minimizer);  
   GetFitPar(2,  al0_,  dal0,  minimizer);  
   GetFitPar(3,  al_,   dal,   minimizer);  
   GetFitPar(4,  Ares_, dAres, minimizer);  
   GetFitPar(5,  Cbg_,  dCbg,  minimizer);  
   GetFitPar(6,  R_,    dR,    minimizer);  
   GetFitPar(7,  b_res_,db_res,minimizer);  
   GetFitPar(8,  b_bg_, db_bg, minimizer);  
   GetFitPar(9,  s0_,   ds0,   minimizer);  
   GetFitPar(10, eta_,  deta,  minimizer);
 
//   double par[]={ Ael_,bel_,al0_,al_,Ares_,Cbg_,R_,b_res_,b_bg_,s0_,eta_};

   double chi2=Chi2full();
   printf("Ael:%7.3f %7.3f, bel:%6.3f %6.3f, al0:%4.2f %4.2f, al:%6.3f %6.3f, Ares:%6.3f %6.3f, Cbg:%6.3f %6.3f, R:%6.3f %6.3f, bres:%6.3f %6.3f, bbg:%6.3f %6.3f\n",Ael_,dAel, bel_,dbel, al0_,dal0, al_,dal, Ares_,dAres, Cbg_,dCbg, R_,dR, b_res_,db_res, b_bg_,db_bg);
  
//   printf("##Chi2:%6.3f\n",Chi2full());
   printf(" pp_el   chi2:%6.3f\n",Chi2pp());
   printf(" d2SD 1  chi2:%6.3f\n",Chi2_d2SD_DC(d2csSD[0],0.05,Npp));
   printf(" d2SD 2  chi2:%6.3f\n",Chi2_d2SD_DC(d2csSD[1],0.05,Npp));
   printf(" dsdt SD chi2:%6.3f\n",Chi2_dcsdtSD());
   printf(" cs SD   chi2:%6.3f\n",Chi2_csSD());
   printf(" cs SD <0.05  chi2:%6.3f\n",Chi2_csSD_DC(csSD_xi,Npp,1));
   printf(" cs SD <200   chi2:%6.3f\n",Chi2_csSD_DC(csSD_m2,Npp,2));

   printf(" cs DD   chi2:%6.3f\n",Chi2_csDD());

//  for(int i=0;i<xNt;i++)
//   for(int j=0;j<xNm2;j++)
//    {printf("%10.8f  %7.4f\n",xM2[j],d2cs_SD(0.05, xM2[j], sq(546)));
//    printf("%10.8f  %7.4f\n",xM2[j],d2cs_SD(0.05, xM2[j], sq(1800))); }
  }

//~~~~~~~~~~~~~~~~~~~~~~
void Fitting()
{
 //~~ DataContainers initialization ~~~
  int block; bool sym, bins;
  
  //## pp ##
  const char* f7000pp="../Data/mx/ppel.TOTEM.datx";
  //DataContainer DC7000pp[2];
   DC7000pp[0].Get(f7000pp, block=1, sym=true,  bins=false);
   DC7000pp[1].Get(f7000pp, block=2, sym=false, bins=false);

  const char* pp_dat="../Data/mx/ppbar.sqrtsGt100.dc.Data";
  //DataContainer DC1800pp[2];
   DC1800pp[0].Get(pp_dat, block=7, sym=true,  bins=false);
   DC1800pp[1].Get(pp_dat, block=8, sym=true,  bins=false);
  //DataContainer DC546pp[4];
   DC546pp[0].Get(pp_dat, block=2, sym=true,  bins=false);
   DC546pp[1].Get(pp_dat, block=3, sym=true,  bins=false);
   DC546pp[2].Get(pp_dat, block=5, sym=true,  bins=false);
   DC546pp[3].Get(pp_dat, block=4, sym=true,  bins=false);
  //DataContainer DC540pp;
   DC540pp.Get(pp_dat, block=1, sym=true,  bins=false);
  // DC540pp.Print();

  //## SD ##
  const char* fdcsdtSD="../Data/mx/dcsdt.UA4.Data";
  //DataContainer dcsdtSD;
   dcsdtSD.Get(fdcsdtSD, block=1, sym=true,  bins=false);
    
   dcsdtSD.Print();
 
  const char* fd2csSD="../Data/mx/d2cs.SD.Data";
  //DataContainer d2csSD[3];
   d2csSD[0].Get(fd2csSD, block=1, sym=true,  bins=false);
   d2csSD[1].Get(fd2csSD, block=2, sym=true,  bins=false);
   d2csSD[2].Get(fd2csSD, block=3, sym=true,  bins=false);
  
 //## cs(s) ##
   //~~ sd ~~
   const char* fcsSD="../Data/mx/cs.SD.Data";
   //DataContainer csSD_m2, csSD_xi;
    csSD_xi.Get(fcsSD, block=1, sym=true,  bins=false);
    csSD_m2.Get(fcsSD, block=2, sym=true,  bins=false);
    csSD_m2.Print();
   csSD_xi.Print();

   //~~ dd ~~
   const char* fcsDD="../Data/mx/cs.DD.Data";
   //DataContainer csDD;
    csDD.Get(fcsDD, block=1, sym=true,  bins=false);
    csDD.Print();

   s_max=1800*1800;
   prestage_M2(xM2, xdm2, xNm2, 1,s_max,"log");
   stageI_M2(xIRes,xIBg,xIRoper,xM2,xNm2);
   prestage_tlin(st, sdt, xNt, 0, 1);
//   stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, t, Nt);
 //~~ Minimization ~~~
  Minimize(); 
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //const int Nt=400, Nm2=1000;
   //const int Nt=100, Nm2=500;//6000;
   //const int Nt=100, Nm2=2000;//6000;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
bool isFitting=true;
int main()
{
 printf("Program is starting\n");
 if(isFitting){
   printf("Fitting...\n");
   al0_=alpha0; al_=alpha_x;
   Ael_=Norm_el; bel_=bel;
   s0_=s0;
   Norm_DD = 1./4./Norm_el;
   
   Ares_=Ares/(2*mp);
   Cbg_=Cbg;
   R_=nR0;
   b_res_=b_res;
   b_bg_=b_bg; 
   eta_=_eta;

 //  Fitting();

   alpha0=al0_; alpha_x=al_;
   Norm_el = Ael_; bel=bel_;
   s0=s0_;
   Norm_DD = 1./4./Norm_el;
   
   Ares=Ares_*(2*mp);
   Cbg=Cbg_;
   nR0=R_;
   b_res=b_res_;
   b_bg=b_bg_; 
   _eta=eta_;
  }
//printf()

//#### Print Parameters ############################
printf("Ael:%7.3f, bel:%6.3f, al0:%4.3f, al:%6.3f,Ares:%6.3f, Cbg:%6.3f, R:%6.3f, bres:%6.3f, bbg:%6.3f\n",Norm_el,bel,alpha0,alpha_x,Ares/(2*mp), Cbg,nR0,b_res,b_bg);
printf("al0:%4.3f, al:%6.3f \t ", alpha0,alpha_x);
printf("eta:%6.2f, xi:%6.2f\n",_eta,xi);

//   Plot_xGoul();

//#### Elastic part ############################
 const int Nt=300, Nm2=2000;//6000;
 double t[Nt], dtarr[Nt]={0};
 //~ Elastic  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 prestage_tlin(t, dtarr, Nt, 0,2.6);
 stageII_el(t,dtarr,Nt);
 //_ test_el(t,Nt);

///########## SD  ############################
  printf("####### SD #######\n");
  stage_IIs_dependense();//(!)1e8(!!!)
 //~~~~~~~~~~~~~~~~~
  double  M2[Nm2], dm2arr[Nm2]={0};
  double  f4_res[3*Nt], Fp2_t[Nt], f_bg[Nt], alpha_IP[Nt];
  double  ImA[3*Nm2], ImBg[Nm2], ImR[Nm2];
  prestage_tlin(t, dtarr, Nt, t_min,t_max);
  prestage_M2(M2, dm2arr, Nm2, M2_min,M2_max,"log");
  stageI_M2(ImA,ImBg,ImR,M2,Nm2);
  stageI_t(Fp2_t, f4_res, f_bg, alpha_IP, t, Nt);

  stage_II_SDpredict();
  stage_IInormalize  (t,dtarr,M2,dm2arr, Nt,Nm2, Fp2_t,f4_res,f_bg, alpha_IP,ImA,ImBg,ImR);
  stage_Goul(M2,Nm2, ImA,ImBg,ImR);
  //_ test_SD_Goul (M2);
  stageIIAddonData();
  stageII(t,dtarr,M2,dm2arr, Nt,Nm2, Fp2_t,f4_res,f_bg,alpha_IP,ImA,ImBg,ImR); //{integr}

 //_   test_SD_I();
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   int oNm2=1000;
   double oM2L=0.5, oM2U=12; 
   prestage_M2(M2,dm2arr, oNm2,oM2L,oM2U,"lin");
   stageIIother(t,dtarr, M2,dm2arr, Nt,oNm2);//{t1,t2,t3};{m1,m2,m3}; //{diff}
//    printf("stageII     [passed]\n");
//    //~~~~~~~~~~~~~~~~~
//    //  Plot_Vall(Nt,t); // stageII & stageII_el needed (!)
//

//########## DD ############################
//   printf("####### DD #######\n");
//   //stage_IIIs_normalize2();
//     stage_IIIs_Risto_Invis();// for Risto limits
////   //.... 2D-plots ....
//     stage_IIIint_other();  //3D (integrated)
//     stage_IIIdiff_other(); //3D (diff)
//
// //  stage_IIIs_dependese();// ~7min
//     stage_III_DDprediction();
   stage_IIIs_dependese2(); // ~12 sec
  
  
 //  PlotDD_I();
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  printf("Done\n");
  return 0;
}


//### Coments #######
//9Jan2013  dcs_dt_el is the same for fitting and for graphs (+)
//          stage_Goul(..) & d2cs_SD give the same results   (+)
//11Jan2013 stage_IIcernel & double d2cs_SD(t,M2,s) give the same results   (+) 

//### Remarks #######                     
// stage_IIcernel -- double arms
// d2cs_SD(t,M2,s) -- single arm (for fitting) 




