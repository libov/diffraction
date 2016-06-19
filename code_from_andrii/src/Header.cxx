#include <stdio.h>
#include <string.h> // Щоб використ. strlen, memcmp, strtok
#include "Header.h"

int Header::GetHeader(FILE* file, int block)
{
  char buffer [100];

  char sprobe[5]; sprintf(sprobe,"%d)",block);
  int  slen=strlen(sprobe);
  fpos_t pos;
  while (!feof(file))
   {
    if(fgets(buffer,100,file)!=NULL)
      {
        if(memcmp(sprobe,buffer,slen)==0)
          {
            strcpy(title,buffer+slen);
            title[strlen(title)-1]='\0';

            fgets(buffer,100,file);
            if(buffer[0]!='#'){printf("Wrong format!");return 1;}
            strcpy(info,buffer+1);
            info[strlen(info)-1]='\0';
            int pass;
            char* ch;
             do{
                pass=0;
                fgetpos(file,&pos);
                fgets(buffer,100,file);
                if(memcmp("#%",buffer,2)==0) {pass++;}

                if(memcmp("#X:",buffer,3)==0)
                  {pass++;
                   Xnm=strtok(buffer+3," ");
                   Xunits=strtok(NULL,"\n ");}

                if(memcmp("#Y:",buffer,3)==0)
                  {pass++;
                   Ynm=strtok(buffer+3," ");
                   Yunits=strtok(NULL,"\n ");}

                if(memcmp("#P:",buffer,3)==0)
                  {pass++;
                   Pnm    = strtok(buffer+3," ");
                   Punits = strtok(NULL,"\n ");}
              }while(pass);
            fsetpos(file,&pos);
	    return 0;
          }
      }
    }
return 0;
}

void Header::PrintHeader()
 {
  printf("T: %s\n",title);
  printf("I: %s\n",info);
  printf("x: %s; %s\n",Xnm.data(),Xunits.data());
  if(Pnm!="")
  printf("P: %s; %s\n",Pnm.data(),Punits.data());
  printf("y: %s; %s\n",Ynm.data(),Yunits.data());
 }
