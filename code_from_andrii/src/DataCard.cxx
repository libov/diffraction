#include "DataCard.h"
#include "math.h"

void DataCard::Set(const char* f, int bk, bool sm, bool bn, double w, double q2, double mv2)
 {file=f; block=bk; sym=sm; bins=bn; W=w; Q2=q2; Mv2=mv2;
  if(w==-2||q2==-2)
    {DataType=1;//B-slope or cs
     if(w==-2&&q2==-1){WQ2pos=1;}
     else if(w==-1&&q2==-2){WQ2pos=-1;} else printf("\nError!! WQ2pos...\n");
    }
  else
    {DataType=2;//dcs/dt 
     if(w==-1&&q2>=0){WQ2pos=1;}
     else if(w>=0&&q2==-1){WQ2pos=-1;} else printf("\nError!! WQ2pos...\n");
    }
  dc.Get(file,block,sym,bins);
//  printf(">>>!<<<\n");
//  dc.Print();
}
  
double DataCard::Chi2(double* par)
{
 iii=Xi2=0;
 Npages=dc.pages.size(); 
 for(i=0;i<Npages;i++)
   Chi2(par,i);
 return Xi2;
}  

double DataCard::Chi2(double* par,int iBlock)// Does'n clean Xi and iii
{
 Nlines = dc.pages[iBlock]->N;
 _Y     = dc.pages[iBlock]->Y;
 _dY    = dc.pages[iBlock]->dYtotal;
 _X     = dc.pages[iBlock]->X;
 _P     = dc.pages[iBlock]->P;
 if(DataType==1)//B-slope, cs. No t dependance.
  {
   if(WQ2pos==1)//{W,Q2,Mv2,par}
    {for(j=0;j<Nlines;j++)
      {Xi=((double)_Y[j]-funcI(pow(_X[j],2), _P, Mv2, par)/metric)/_dY[j];Xi2+=Xi*Xi;iii++;}}
   else if(WQ2pos==-1)
    {for(j=0;j<Nlines;j++)
      {Xi=((double)_Y[j]-funcI(    _P*_P, _X[j], Mv2, par)/metric)/_dY[j];Xi2+=Xi*Xi;iii++;}}
  }
 else if(DataType==2)//dcsdt. |t|-dependant.
  {if(WQ2pos==1)//{t,W,Q2,Mv2,par}
    {for(j=0;j<Nlines;j++)
      {Xi=((double)_Y[j]-funcII(_X[j], _P*_P, Q2, Mv2, par)/metric)/_dY[j];Xi2+=Xi*Xi;iii++;}}
   else if(WQ2pos==-1)
    {for(j=0;j<Nlines;j++)
      {Xi=((double)_Y[j]-funcII(_X[j],   W*W, _P, Mv2, par)/metric)/_dY[j];Xi2+=Xi*Xi;iii++;}}
  }
 return Xi2;
}

