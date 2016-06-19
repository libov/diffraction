#include "math.h"
#include "TF1.h"
#include "TFile.h"
#include "DataBody_sym.h"
#include "DataContainer.h"
#include "Drawer2.h"

class TGraphErrors;
class TCanvas;

void Drawer2::GetPage(int i)
{
  Nlines_ = dc->pages[i]->N;
  if(!dX)delete[]dX; dX=new float[Nlines_];
  for(int j=0;j<Nlines_;j++)dX[j]=0;
  Y  = dc->pages[i]->Y;
  dY = dc->pages[i]->dYtotal;
  X  = dc->pages[i]->X;
}

void Drawer2::Draw(const char* file="x.root")
{
 gr = new TGraphErrors[Npages];
 cv = new TCanvas[Npages];

 string str=""; 
 char cstr[50];
 for(int j=0;j<Npages;j++)
  {
   GetPage(j);
   sprintf(cstr,"%s %5.4f %s ;-t;d#sigma/dt",dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
   gr[j] = TGraphErrors(Nlines_,X,Y,dX,dY);
   gr[j].SetTitle(cstr);
   gr[j].SetMarkerColor(4);
   gr[j].SetMarkerStyle(21);
   //###########################################################
   sprintf(cstr,"%s %f %s",dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
   cv[j].SetTitle(cstr);
   sprintf(cstr,"cv%d %s %f %s",j,dc->head.Pnm.data(),dc->pages[j]->P, cvaddNm);
   cv[j].SetName(cstr);
   cv[j].SetWindowSize(800,600);
  }
 printf("# %s\n",file);
 TFile f(file,"RECREATE");
 TPad* pad;

 for(int j=0;j<Npages;j++)
  {
    pad=(TPad*)cv[j].cd();
    pad -> SetLogy();
    pad -> SetGridx();
    pad -> SetGridy();
    gr[j].Draw("AP");
    cv[j].Write();
  }
 f.Close();
}

void Drawer2::Draw2(const char* file,int w, int h)
{
 char cstr[50];
 string str=""; 
 gr = new TGraphErrors[Npages];
 for(int j=0;j<Npages;j++)
  {
   GetPage(j);
   sprintf(cstr,"%s %5.4f %s ;-t;d#sigma/dt",dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
   gr[j] = TGraphErrors(Nlines_,X,Y,dX,dY);
   gr[j].SetTitle(cstr);
   gr[j].SetMarkerColor(4);
   gr[j].SetMarkerStyle(21);
   gr[j].SetMarkerSize(0.5);
  }

 int Nhw=h*w;
 int Ncv=ceil(1.0*Npages/Nhw),i;// n of Canvas;
 int j; //n of graph in Canvas;
 int fg; //n of graph

 cv = new TCanvas[Ncv];
 printf("# %s\n",file);
 TFile f(file,"RECREATE");
 TPad* pad;
// printf("Npages:%d; Ncv:%d; Nhw:%d\n",Npages, Ncv, Nhw); 
 for(fg=i=0;fg<Npages;i++)
  { 
// printf("##%d\n",i);
    sprintf(cstr,"cv%d",i);            cv[i].SetTitle(cstr);
    sprintf(cstr,"cv%d%s",i,cvaddNm);  cv[i].SetName(cstr);
    cv[i].SetWindowSize(800,600);
    cv[i].Divide(w,h);

//    printf(">>");
    for(j=1;j<=Nhw&&fg<Npages;j++)
     {
      pad=(TPad*)cv[i].cd(j);
      pad -> SetLogy();
      pad -> SetGridx();
      pad -> SetGridy();
      gr[fg].Draw("AP");
      fg++;
     }
//    printf("cv[%d]",i);
    cv[i].Write();
//    printf(" Ok\n",i);
  }
// printf("##%d\n",i);

 f.Close();
}

void Drawer2::Draw3(const char* file)
{
 char cstr[50];
 string str=""; 
 gr = new TGraphErrors[Npages];
 for(int j=0;j<Npages;j++)
  {
   GetPage(j);
   sprintf(cstr,"%s %5.4f %s ;-t;d#sigma/dt",dc->head.Pnm.data(),dc->pages[j]->P,dc->head.Punits.data() );
   gr[j] = TGraphErrors(Nlines_,X,Y,dX,dY);
   gr[j].SetTitle(cstr);
   gr[j].SetMarkerColor(800+j);//4
   gr[j].SetLineColor(800+j);//4
   gr[j].SetMarkerStyle(21);
   gr[j].SetMarkerSize(0.3);
  }

 int fg; //n of graph
 printf("# %s\n",file);
 TFile f(file,"RECREATE");
 TCanvas canv("all","All",800,600,0,0);
 canv.SetLogy();
 canv.SetGridx();
 canv.SetGridy();
 gr[1].Draw("AP");

 for(fg=0;fg<Npages;fg++) gr[fg].Draw("P, same");
 canv.Write();
 f.Close();
}
