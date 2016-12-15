#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <vector>
using namespace std;
void manual()
{
     cout<< "\n\nSXFilterRepeatGap\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXFilterRepeatGap <File_Name> \n\n";
}

bool FilterGapsnps(int i, size_t offset)
{
    string filename = "../../syhwah/Manish/gapdata/chr";
    std::stringstream out;
    out << i;
    string s = filename;
    s  += out.str();
    s  += ".txt";
    size_t start,stop;
    ifstream in(s.c_str());
    char* pch;
    char line[1000];
    pch = strtok(line,"\t");
    while(in)
    {
      in.getline(line,1000);
      if(in)
      {
         char* pch;
         vector <const char*> vec;
         pch = strtok(line,"\t");
         while(pch != NULL)
         {
            vec.push_back(pch);
            pch = strtok(NULL,"\t");
         }
         size_t start = (size_t) atol(vec[0]);
         size_t stop = (size_t) atol(vec[1]);

         if(offset >= start && offset <= stop)
         {
           in.close();
           return 1;
         }
      }
    }
    in.close();
    return 0;
}

unsigned int GetRepeatregion(int i, unsigned int offset)
{
    string filename;
    filename = "Reptmask/Chr"; 
    std::stringstream out;
    out << i;       
    string s = filename;
    s  += out.str();
    s += "_mask_den";
    FILE* pf = fopen(s.c_str(),"rb");
    if (!pf) return 0;
    unsigned int g_unReads;   
    unsigned char g_acBuf[2];   
    fseek(pf,2*offset,SEEK_SET);   
    fread(g_acBuf,2,1,pf);
    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|=g_acBuf[0];
    fclose(pf);
    return g_unReads;
}
 
int main(int argc, char* argv[])
{
   if(argc != 2)
   {
     manual();
     exit(1);
   }
   ifstream in(argv[1]);
   char outfile[1000];
   strcpy(outfile,argv[1]);
   strcat(outfile,".filtrpt.filtgap");
   ofstream out(outfile);
   char line[5000], recline[5000];
   while(in)
   {
      in.getline(line,5000);
      strcpy(recline,line);
      if(in)
      {
         char* pch;
         pch = strtok(line,"\t");
         vector<const char*> vec;
         while(pch != NULL)
         {   
             vec.push_back(pch);
             pch = strtok(NULL,"\t"); 
         }
           int chrnum;
           if(vec[0]==  "X")
           chrnum = 23;
           else if (vec[0] ==  "Y")
           chrnum = 24;
           else
           chrnum = atoi(vec[0]);

         int rpc = GetRepeatregion(chrnum, atoi(vec[1]));
         int gap = FilterGapsnps(chrnum, atoi(vec[1])); 
         if(rpc == 0 )
           out  << recline << endl;
      }
   }
   out.close();
  return 0;
}
