#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


void manual()
{
     cout<< "\n\nSXFilterSNPRSRS\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXFilterSNPRSRS <File_Name> <Supporting_Reads> <Read_Strength> <Read_Density> \n\n";
}

void readsnpfile1(char* line1,int nSR, int nRS, int nRD)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chr"))
           return;
        else
        {
           char recline1[10000];
           strcpy(recline1,line1);
           char* pch;
           vector<string> sdf;
           string s;
           pch = strtok (recline1,"\t");
           while(pch != NULL) 
           {
              sdf.push_back(pch);
              pch = strtok (NULL,"\t");
           }
           if(atoi(sdf[8].c_str()) > 0)
           {
             int rs =  (int)(atoi(sdf[5].c_str())*100/atoi(sdf[8].c_str())+0.5);
             if( atoi(sdf[5].c_str()) >= nSR && rs >= nRS && atoi(sdf[8].c_str()) >= nRD)
               cout << line1<< endl;
          }
         }
     }
}

int main(int argc, char* argv[])
{
  if (argc != 5)
  {
     manual();
     exit(1);
  }

  int nSR = atoi(argv[2]);
  int nRS = atoi(argv[3]);
  int nRD = atoi(argv[4]);
     
  ifstream in(argv[1]);
  char line[5000],recline[5000];
  while(in)
  {
     in.getline(line,5000);
     if(in)
     {
        readsnpfile1(line,nSR,nRS,nRD);
     }
  }
  in.close();
}
