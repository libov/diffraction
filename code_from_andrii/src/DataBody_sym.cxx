#include "math.h"
#include <vector>
#include "DataBody_sym.h"

using namespace std;                  

bool DataBody_sym::ReadTable(FILE* file, char* header)
{
  if(!DataBody::ReadTable(file,header))return false;
  fpos_t pos;
 
  float tX,tY,tdYsist,tdYstat;
  vector<float> vX,vY,vdYstat;
  vector<float> vdYsist;
  fgetpos(file,&pos);
  while(fscanf(file,"%f %f %f %f",&tX,&tY,&tdYstat,&tdYsist)==4)
   {
    vX.push_back(tX);
    vY.push_back(tY);
    vdYstat.push_back(tdYstat);
    vdYsist.push_back(tdYsist);
    N++;
    fgetpos(file,&pos);
   }
  fsetpos(file,&pos); 
 
  clean();
  X  = new float[N];
  Y  = new float[N];
  dYstat  = new float[N];
  dYsist  = new float[N];
  dYtotal = new float[N];
 
  for(int i = 0;i<N; i++)
   {
    X[i] = vX[i];
    Y[i] = vY[i];
    dYstat[i] = vdYstat[i];
    dYsist[i] = vdYsist[i];
    dYtotal[i]= sqrt(vdYsist[i]*vdYsist[i]+vdYstat[i]*vdYstat[i]);
   }
  return true;
}

void DataBody_sym::PrintBody(const char* header) 
{
 if(P!=-1111)printf("## %8.2f\n",P);
 if(header!=NULL)printf("   %s  %s\n",header,"  dtotal");
 for(int i=0;i<N;i++) 
   printf("%7.2f %7.2f %7.2f %7.2f %7.2f\n",
           X[i],Y[i],dYstat[i],dYsist[i],dYtotal[i]);
 printf("\n");
}
