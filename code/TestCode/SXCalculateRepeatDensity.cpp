#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdio.h>
using namespace std;

int main(int argc, char* argv[])
{
   string filename;
   
   filename = "/archive/project/oneclick/syhwah/GS000117/37_1_CG1/raw_read/";
   string s = filename;
   s += argv[1];
   s += "_rept_den";
   cout << s <<endl;
   FILE *pf = fopen(s.c_str(),"rb");
   if (!pf) exit(1);
   unsigned long sum  = 0;
   double stdev = 0;
   size_t chrlength = (size_t)atol(argv[2]);
   size_t nrpt_count =0;
   for (unsigned int i =1; i <= chrlength; i++)
   {
       unsigned int g_unReads;
       unsigned char g_acBuf[2];
       fseek(pf,(i-1)*2,SEEK_SET);
       fread(g_acBuf,2,1,pf);
       g_unReads = g_acBuf[1];g_unReads<<=8;
       g_unReads|= g_acBuf[0];
       int rpt = g_unReads;
       if(rpt != 0)
         nrpt_count++;    
       sum += rpt;
  }
    double avg = double(sum)/(double)nrpt_count;
    for (unsigned int j =1; j <= chrlength; j++)
    {
       unsigned int g_unReads2;
       unsigned char g_acBuf2[2];
       fseek(pf,(j-1)*2,SEEK_SET);
       fread(g_acBuf2,2,1,pf);
       g_unReads2 = g_acBuf2[1];g_unReads2<<=8;
       g_unReads2|= g_acBuf2[0];
       int rpt2 = g_unReads2;
       if(rpt2 != 0)
       {
         double diff = rpt2 - avg ;
         double square  = diff*diff;
         stdev += square;
       }
   }
  fclose(pf);
  cout << argv[1]<< "\t" << avg <<"\t"<<sqrt(stdev/nrpt_count)<< endl;
  return 0;
}

