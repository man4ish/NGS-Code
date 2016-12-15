#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
using namespace std;
unsigned int GetReadDensity(int i, unsigned int uioffset)
{
    string filename;
    filename = "CG5ReadDensity/Chr";
    std::stringstream out;
    out << i;
    string s = filename;
    s  += out.str();
    s += "_read_den";
    FILE* pf = fopen(s.c_str(),"rb");
    if (!pf) return 0;
    unsigned int g_unReads;
    unsigned char g_acBuf[2];
    fseek(pf,(uioffset-1)*2,SEEK_SET);
    fread(g_acBuf,2,1,pf);
    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|= g_acBuf[0];
    fclose(pf);
    return g_unReads;
}

map <int,float> m;
int main(int argc, char* argv[])
{
     ifstream in(argv[1]);
     char line1[5000],line2[5000]; 
     while(in)
     {
       in.getline(line1,5000);
       strcpy(line2,line1);
       if(in)
       {  
          char* pch;
          pch = strtok(line1,"\t");
          vector <string> records;
          while(pch != NULL)
          {
              records.push_back(pch);
              pch = strtok(NULL,"\t");
          }
         
          int tot_rd = 0;
          int chrnum;
          if(records[0]== "X")
          chrnum = 23;
          else if(records[0]== "Y")
          chrnum = 24;
          else
          chrnum = atoi(records[0].c_str());
          size_t offset = (size_t)atol(records[1].c_str());
          for (int index = 0; index <= 200;index++)
          tot_rd += GetReadDensity(chrnum,offset-100+index);
          float avg_rd_unique2 = ((float)tot_rd/201);
          float cutoff = (float) avg_rd_unique2/(float)(40);
          cout << line2 <<"\t"<< GetReadDensity(chrnum,offset)+ atoi(records[3].c_str())<< endl;
       }
     }
     return 0;
}
