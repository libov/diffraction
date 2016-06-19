#ifndef Header_h
#define Header_h

#include <string> // Щоб використовувати string
using namespace std;

class Header
{
 public:
  Header(){Pnm=Punits="";header[0]='\0';}
  ~Header(){;}
  int GetHeader(FILE* file, int block=1);
  void PrintHeader();

  char header[100];// header of table X[..] Y[...] dY[...] 
  char title[100];
  char info[100];

  string Xnm,Xunits,Ynm,Yunits;
  string Pnm,Punits;
};

#endif

