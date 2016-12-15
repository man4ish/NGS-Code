#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdio.h>
using namespace std;
unsigned int GetReadDensity(int i, unsigned int uioffset)
{
    string filename;
    filename = "raw_read/Chr";
    std::stringstream out;
    out << i;
    string s = filename;
    s  += out.str();
    s += "_read_den"; 
    FILE* pf = fopen(s.c_str(),"rb");
    if (!pf) return 0;
    unsigned int g_unReads;
    unsigned char g_acBuf[2];
    fseek(pf,uioffset*2,SEEK_SET);
    fread(g_acBuf,2,1,pf);
    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|= g_acBuf[0];
    fclose(pf);
    return g_unReads;
}

int main(int argc, char* argv[])
{
   unsigned long sum  = 0;
   double stdev = 0;
   int chrnum = atoi(argv[1]);
   size_t chrlength = (size_t)atol(argv[2]);
  
  for (unsigned int i =1; i<=chrlength; i++)
  {
      int rd = GetReadDensity(chrnum,i);
      sum += rd;
  }
  double avg = double(sum)/(double)chrlength;
  for (unsigned int i =1; i<=chrlength; i++)
  {
      int rd1 = GetReadDensity(chrnum,i);
      double diff = rd1 - avg ;
      double square  = diff*diff;
      stdev += square;
  }
  double cutoff = avg + 3*sqrt(stdev/chrlength);
  cout << argv[1]<< "\t" << avg <<"\t"<<sqrt(stdev/chrlength)<<"\t"<< cutoff << endl;
   return 0;
}

