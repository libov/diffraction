#include "Header.h"
#include "DataBody_sym.h"
#include "DataBody_sym_b.h"
#include "DataBody_usym.h"
#include "DataBody_usym_b.h"
#include "DataContainer.h"

void DataContainer::Get(const char* file, int block, bool sym, bool bins)
{
  bool bRead=false;
  symetry=sym;
  binning=bins;
  FILE* pFile=fopen(file,"r");
  if (pFile == NULL) {perror ("Error opening file");return;}
 
  head.GetHeader(pFile, block);
  DataBody* tmpPage;
  do{
     if(sym) 
       {if(binning) tmpPage = new DataBody_sym_b;
          else      tmpPage = new DataBody_sym; } 
     else
       {if(binning) tmpPage = new DataBody_usym_b;
          else      tmpPage = new DataBody_usym; } 
 
//     printf("Hallo\n"); 
     bRead = tmpPage->ReadTable(pFile,head.header);
     if(bRead) pages.push_back(tmpPage);
    }while(bRead);

  if(sym) 
    {if(binning) delete ((DataBody_sym_b*)tmpPage);
       else      delete ((DataBody_sym*)  tmpPage);} 
  else
    {if(binning) delete ((DataBody_usym_b*)tmpPage);
       else      delete ((DataBody_usym*)  tmpPage);} 
  
  fclose(pFile); 
  return;
}

 void DataContainer::Print()
 {
   head.PrintHeader(); 
   int L=pages.size(); 
   for(int i=0;i<L;i++)
     pages[i]->PrintBody(head.header);
 }

 DataContainer::~DataContainer()
 {int L=pages.size(); 
  if(symetry) 
    {if(binning) for(int i=0;i<L;i++) delete ((DataBody_sym_b*)pages[i]);
       else      for(int i=0;i<L;i++) delete ((DataBody_sym*)  pages[i]);} 
  else
    {if(binning) for(int i=0;i<L;i++) delete ((DataBody_usym_b*)pages[i]);
       else      for(int i=0;i<L;i++) delete ((DataBody_usym*)  pages[i]);} 
 }
