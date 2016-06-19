#include "math.h"
#include <vector>
#include "DataBody_sym_b.h"

using namespace std;                  

bool DataBody_sym_b::ReadTable(FILE* file, char* header)
{
  if(!DataBody::ReadTable(file,header))return false;
  fpos_t pos;
 
  float  tX1,tX2,tX,tY,tdYsist,tdYstat;
  vector<float> vX1,vX2,vX,vY,vdYstat;
  vector<float> vdYsist;
  fgetpos(file,&pos);
  while(fscanf(file,"%f %f %f %f %f %f",&tX1,&tX2,&tX,&tY,&tdYstat,&tdYsist)==6)
   {
    vX1.push_back(tX1);
    vX2.push_back(tX2);
    vX.push_back(tX);
    vY.push_back(tY);
    vdYstat.push_back(tdYstat);
    vdYsist.push_back(tdYsist);
    N++;
    fgetpos(file,&pos);
   }
  fsetpos(file,&pos); 
 
  clean();
  X1 = new float[N];
  X2 = new float[N];
  X  = new float[N];
  Y  = new float[N];
  dYstat  = new float[N];
  dYsist  = new float[N];
  dYtotal = new float[N];
 
  for(int i = 0;i<N; i++)
   {
    X1[i]= vX1[i];
    X2[i]= vX2[i];
    X[i] = vX[i];
    Y[i] = vY[i];
    dYstat[i] = vdYstat[i];
    dYsist[i] = vdYsist[i];
    dYtotal[i]= sqrt(vdYsist[i]*vdYsist[i]+vdYstat[i]*vdYstat[i]);
   }
  return true;
}

void DataBody_sym_b::PrintBody(const char* header) 
{
 if(P!=-1111)printf("## %8.2f\n",P);
 if(header!=NULL)printf("   %s  %s\n",header,"  dtotal");
 for(int i=0;i<N;i++) 
   printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
           X1[i],X2[i],X[i],Y[i],dYstat[i],dYsist[i],dYtotal[i]);
 printf("\n");
}
