#include <stdlib.h>
#include "string.h"
#include <string>
#include "DataBody.h"
using namespace std;

bool DataBody::isSpace(const char* str)
{
  int i=0;
  while(isspace(str[i])&&str[i]!='\0') {i++;}
  if(str[i]=='\0') return true;
  return false;
}

bool DataBody::ReadTable(FILE* file, char* header)
{
 int pass,pass2=0;
 char buffer [100];
 fpos_t pos;
 string s;

// printf(">>###\n"); 
 do{
//    printf(">>     inside\n"); 
    pass=0;
    fgetpos(file,&pos);
    fgets(buffer,100,file);
    if(isSpace(buffer)){pass++;}

    if(memcmp("## ",buffer,3)==0)
      {pass++;pass2++; P=atof(strtok(buffer+3," "));}

    if(buffer[0]=='%')
     if(buffer[1]=='%'){pass++;}
     else{pass++;pass2++; if(header!=NULL)strcpy(header,strtok(buffer+1,"\n"));}
 //   printf("<>%d|%s",isSpace(buffer),buffer); 
    if(feof(file))pass=0;
  }while(pass);
 fsetpos(file,&pos); 
// printf(">>pass2:%d\n",pass2); 
 if(pass2==0)return false;
 return true; 
}

 void DataBody::clean()
  {clear(X);clear(Y);clear(dYstat);clear(dYtotal);}
