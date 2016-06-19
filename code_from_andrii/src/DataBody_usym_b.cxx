#include <vector>
#include "math.h"
#include "DataBody_usym_b.h"
using namespace std;   

bool DataBody_usym_b::ReadTable(FILE* file, char* header)
{
 if(!DataBody::ReadTable(file,header))return false;
 fpos_t pos;

 float  tX1,tX2,tX,tY,tdYstat,tdYsist_up,tdYsist_low;
 vector<float> vX1,vX2,vX,vY,vdYstat;
 vector<float> vdYsist_up,vdYsist_low;
 fgetpos(file,&pos);
 while(fscanf(file,"%f %f %f %f %f %f %f",&tX1,&tX2,&tX,&tY,&tdYstat,&tdYsist_low,&tdYsist_up)==7)
  {
   vX1.push_back(tX1);
   vX2.push_back(tX2);
   vX.push_back(tX);
   vY.push_back(tY);
   vdYstat.push_back(tdYstat);
   vdYsist_up.push_back(tdYsist_up);
   vdYsist_low.push_back(tdYsist_low);
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
 dYsist_up   = new float[N];
 dYsist_low  = new float[N];
 dYtotal     = new float[N];
 dYtotal_up  = new float[N];
 dYtotal_low = new float[N];
 for(int i = 0;i<N; i++)
  {
   X1[i]= vX1[i];
   X2[i]= vX2[i];
   X[i] = vX[i];
   Y[i] = vY[i];
   dYstat[i] = vdYstat[i];
   dYsist_up[i]  = vdYsist_up[i];
   dYsist_low[i] = vdYsist_low[i];
   dYtotal_up[i] = sqrt(vdYsist_up[i]*vdYsist_up[i]  +vdYstat[i]*vdYstat[i]);
   dYtotal_low[i]= sqrt(vdYsist_low[i]*vdYsist_low[i]+vdYstat[i]*vdYstat[i]);
   isconverted=1;
   dYtotal[i]= sqrt(pow((vdYsist_up[i]+vdYsist_low[i])/2.,2.)+vdYstat[i]*vdYstat[i]);
  }
}

void DataBody_usym_b::PrintBody(const char* header) 
{
 if(P!=-1111)printf("## %8.2f\n",P);
 if(header!=NULL) printf("    %s %s\n",header,"(-)dtotal (+)dtotal dtotal");
 for(int i=0;i<N;i++)
   printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
	   X1[i],X2[i],X[i],Y[i],dYstat[i],dYsist_low[i],dYsist_up[i],dYtotal_low[i],
	   dYtotal_up[i],dYtotal[i]);
 printf("\n");
}
