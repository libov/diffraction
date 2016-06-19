//#include "string.h"
//#include <iostream>

#include "math.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
//#include "TFitter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "DataBody_sym.h"
#include "DataContainer.h"
#include "Drawer.h"
//double FdCSdt(double* t,double* par)// par=={vpar1,...,vparn, par1, ... ,parm} 
//If PNum <0 then P isn't taken here
//vNpar could be 0;
void Drawer::GetPage(int i)
{
  Nlines_ = dc->pages[i]->N;
  if(!dX)delete[]dX; dX=new float[Nlines_];
  for(int j=0;j<Nlines_;j++)dX[j]=0;
  Y  = dc->pages[i]->Y;
  dY = dc->pages[i]->dYtotal;
  X  = dc->pages[i]->X;
//  if(PNum==0)vpar[PNum] = pow(dc->pages[i]->P,2);//s=W*W
  if(PNum>=0)vpar[PNum] = dc->pages[i]->P;//Q2,W
/*  printf("     ");for(int i=0;i<vNpar;  i++)printf("%6s ",vparNm[i]);printf("\n");
  printf("vPar:");for(int i=0;i<vNpar;  i++)printf("%6.3f ",vpar[i]);printf("\n");
  printf(" x: "); for(int i=0;i<Nlines_;i++)printf("%6.3f ", X[i]);printf("\n");
  printf(" y: "); for(int i=0;i<Nlines_;i++)printf("%6.3f ", Y[i]);printf("\n");
  printf("dy: "); for(int i=0;i<Nlines_;i++)printf("%6.3f ",dY[i]);printf("\n");
*/
}

void Drawer::SetGraphs(int l=-1) //l=-1 -all; l==0,...,(Npages-1) -only one plot
{
 if(f1)delete []f1;  
 if(gr)delete []gr;  
 int start,finish,j,jj;
 if(l==-1)
  {gr = new TGraphErrors[Npages];
   f1 = new TF1[Npages];
   start=0;finish=Npages;}
 else 
  {gr = new TGraphErrors[1];
   f1 = new TF1[1];
   start=l;finish=l+1;}

 char cstr[200];
 for(int j=start;j<finish;j++)
  {
//   printf("~!~[%d]\n",j);
   GetPage(j);
//   printf(">>%5.2f<<\n",dc->pages[j]->P);
//   printf(">>%s || %s || %s || %s <<\n",dc->head.Pnm.data(),dc->head.Punits.data(),(dc->head.Xnm+" "+dc->head.Xunits).data(),(dc->head.Ynm+" "+dc->head.Yunits).data());

   sprintf(cstr,"\"%s\"|%s %5.2f %s;%s;%s",
      dc->head.title,dc->head.Pnm.data(),
      dc->pages[j]->P,dc->head.Punits.data(),
     (dc->head.Xnm+" "+dc->head.Xunits).data(),
     (dc->head.Ynm+" "+dc->head.Yunits).data());
//   printf("### %s\n",cstr);
//   printf("%d. t_:%f, Y_:%f, dt_:%f, dY_:%f\n",j,*t_,*Y_,*dt_,*dY_);
   jj=j;
   if(l!=-1)jj=0;
    
   gr[jj] = TGraphErrors(Nlines_,X,Y,dX,dY);
   gr[jj].SetTitle(cstr);
   gr[jj].SetMarkerColor(4);
   gr[jj].SetMarkerStyle(21);
   if(y1<y2){gr[jj].SetMaximum(y2);gr[jj].SetMinimum(y1);}

//   s_=70*70;
//   s_=13.3*13.3;
   f1[jj] = TF1("f1",func,fx1,fx2,Npar+vNpar);//6 -- number of par.
//   f1[jj].SetNpx(200);
   f1[jj].SetNpx(20000);
   int i=0;
   for(i=0;i<vNpar;i++)
    {f1[jj].SetParameter(i,vpar[i]);
     f1[jj].SetParName(i,vparNm[i]);}
   for(i=0;i<Npar;i++)
    {f1[jj].SetParameter(i+vNpar,par[i]);
     f1[jj].SetParName(i+vNpar,parNm[i]);}
//   printf("~!~<[%d]\n",j);
  }
//   printf("~!~#[%d]\n",jj);
}
//###########################################################################################
void Drawer::DrawDisp(const char* file="x.root",int ii=0) //one plot 
{
 char cstr[50];
 SetGraphs(ii);//-1
// printf("#~!~#[%d]\n",ii);

 if(cv)delete []cv;  
 cv = new TCanvas[1];

 sprintf(cstr,"%s %f %s",dc->head.Pnm.data(),dc->pages[ii]->P,dc->head.Punits.data() );
 cv[0].SetTitle(cstr);
 sprintf(cstr,"cv%d %s %f %s",ii,dc->head.Pnm.data(),dc->pages[ii]->P, cvaddNm);
 cv[0].SetName(cstr); cv[0].SetWindowSize(800,600);
 cv[0].Divide(2,2); 
 TF1 f1_l[Npar],f1_u[Npar];
 printf("# %s\n",file); TFile f(file,"RECREATE");
 TLegend legend[Npar];

  for(int ipar=0;ipar<Npar;ipar++)
   {
    f1_l[ipar]=f1[0]; f1_l[ipar].SetLineColor(kBlue);f1_l[ipar].SetParameter(ipar+vNpar,par[ipar]*1.1);
    f1_u[ipar]=f1[0]; f1_u[ipar].SetLineColor(kRed );f1_u[ipar].SetParameter(ipar+vNpar,par[ipar]*0.1);
    f1[0].SetLineColor(kGreen);
    char ch[50];
    legend[ipar].SetX1NDC(0.63);
    legend[ipar].SetY1NDC(0.15);
    legend[ipar].SetX2NDC(0.90);
    legend[ipar].SetY2NDC(0.38);
    legend[ipar].SetHeader(parNm[ipar]);
      sprintf(ch,"up  %f",par[ipar]*1.2);legend[ipar].AddEntry(f1_u+ipar,ch,"l");
      sprintf(ch,"f1  %f",par[ipar]);    legend[ipar].AddEntry(f1,ch,"l");
      sprintf(ch,"low %f",par[ipar]*0.9);legend[ipar].AddEntry(f1_l+ipar,ch,"l");
    TPad* pad;
   //    for(int k=0;k<60;k+=5) {printf("%4.1f %8.2f\n",(double)k, f1[j].Eval(k));}
    pad=(TPad*)cv[0].cd(ipar+1);
    if(xLog)pad->SetLogx();
    if(yLog)pad->SetLogy();
    pad -> SetGridx();
    pad -> SetGridy();
    gr[0].Draw("AP");
    legend[ipar].Draw();
    f1[0].Draw("same");
    f1_l[ipar].Draw("same");
    f1_u[ipar].Draw("same");
   }
 cv[0].Write();
 f.Close();

}
//###########################################################################################
void Drawer::Draw(const char* file="x.root",int ii=0) //one plot 
{
 char cstr[50];
 SetGraphs(ii);//-1
// printf("#~!~#[%d]\n",ii);

 if(cv)delete []cv;  
 cv = new TCanvas[1];

  sprintf(cstr,"%s %f %s",dc->head.Pnm.data(),dc->pages[ii]->P,dc->head.Punits.data() );
  cv[0].SetTitle(cstr);
  sprintf(cstr,"cv%d %s %f %s",ii,dc->head.Pnm.data(),dc->pages[ii]->P, cvaddNm);
  cv[0].SetName(cstr); cv[0].SetWindowSize(800,600);
  

 printf("# %s\n",file); TFile f(file,"RECREATE");
 TPad* pad;
//    for(int k=0;k<60;k+=5) {printf("%4.1f %8.2f\n",(double)k, f1[j].Eval(k));}
 pad=(TPad*)cv[0].cd();
 if(xLog)pad->SetLogx();
 if(yLog)pad->SetLogy();
 pad -> SetGridx();
 pad -> SetGridy();
 gr[0].Draw("AP");
 f1[0].Draw("same");
 cv[0].Write();
 f.Close();

}

//######################################################################################
void Drawer::Draw1(const char* file="x.root",int bk) // all Npage-graphs one plot per canvas;
{
 char cstr[250];
 SetGraphs();//-1
 for(int j=0;j<Npages;j++) 
    {gr[j].SetMarkerSize(0.9);
     f1[j].SetLineWidth(1);}

 if(cv)delete []cv;  
 cv = new TCanvas[Npages];
//dc->head.title,
 for(int j=0;j<Npages;j++)
  {
   sprintf(cstr,"%s %s %f %s",dc->head.title,dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
//  printf(">>%s\n",cstr);
//  sprintf(cstr,"%s %f %s",dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
   cv[j].SetTitle(cstr);
   if(bk>=0) sprintf(cstr,"bk%dcv%d %s %f %s",bk,j,dc->head.Pnm.data(),dc->pages[j]->P, cvaddNm);
        else sprintf(cstr,"cv%d %s %f %s",j,dc->head.Pnm.data(),dc->pages[j]->P, cvaddNm);
   cv[j].SetName(cstr); cv[j].SetWindowSize(800,600);
  }

 printf("# %s\n",file); TFile f(file,"RECREATE");
 TPad* pad;
 for(int j=0;j<Npages;j++)
  {
//    for(int k=0;k<100;k+=5) {printf("%4.1f %8.2f\n",(double)k, f1[j].Eval(k));}
    pad=(TPad*)cv[j].cd();
    if(xLog)pad->SetLogx();
    if(yLog)pad->SetLogy();
    pad -> SetGridx();
    pad -> SetGridy();
    gr[j].Draw("AP");
    f1[j].Draw("same");
    cv[j].Write();
  }
 f.Close();
}

//######################################################################################
void Drawer::Draw2(const char* file,int w, int h, int bk)
{
 char cstr[50];
 SetGraphs();//-1
 for(int j=0;j<Npages;j++) 
    {gr[j].SetMarkerSize(0.8);
     f1[j].SetLineWidth(1);}

 int Nhw=h*w;
 int Ncv=ceil(1.0*Npages/Nhw),i;// n of Canvas;
 int j; //n of graph in Canvas;
 int fg; //n of graph

 if(cv)delete []cv;  
 cv = new TCanvas[Ncv];
 printf("# %s\n",file);
 TFile f(file,"RECREATE");
 TPad* pad;
 for(fg=i=0;fg<Npages;i++)
  { 
    sprintf(cstr,"cv%d",i);  cv[i].SetTitle(cstr);
    if(bk>=0) sprintf(cstr,"bk%dcv%d %s %5.2f %s",bk,i,dc->head.Pnm.data(),dc->pages[i]->P, cvaddNm);
         else sprintf(cstr,"cv%d %s %5.2f %s",i,dc->head.Pnm.data(),dc->pages[i]->P, cvaddNm);
  //  sprintf(cstr,"cv%d%s",i,cvaddNm); 
    cv[i].SetName(cstr);
    cv[i].SetWindowSize(800,600);
    cv[i].Divide(w,h);

    for(j=1;j<=Nhw&&fg<Npages;j++)
     {
      pad=(TPad*)cv[i].cd(j);
      if(xLog)pad->SetLogx();
      if(yLog)pad->SetLogy();
      pad -> SetGridx();
      pad -> SetGridy();
      gr[fg].Draw("AP");
//      if(x1<x2){TAxis* ax=((TAxis*)gr[fg].GetXaxis()); ax->SetRangeUser(x1,x2);}
      f1[fg].Draw("same");
      fg++;
     }
    cv[i].Write();
  }

 f.Close();
}

void Drawer::Draw3(const char* file, int l) //only one graph
{
 if(l<0) return;
 char cstr[50];
 SetGraphs(l);
 gr = new TGraphErrors[Npages];
 gr[0].SetMarkerSize(0.3);

// int fg; //n of graph
 printf("# %s\n",file);
 TFile f(file,"RECREATE");
 TCanvas canv("one",gr[0].GetTitle(),800,600,0,0);
 if(xLog)canv.SetLogx();
 if(yLog)canv.SetLogy();
 canv.SetLogy();
 canv.SetGridx();
 canv.SetGridy();
 gr[0].Draw("AP");
 f1[0].Draw("same");
 canv.Write();
 f.Close();
}

