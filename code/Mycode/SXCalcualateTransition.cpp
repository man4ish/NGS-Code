#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
#include <string>
using namespace std;

int transitioncount = 0;
int transversioncount = 0;
int reccount =0;
vector <const char*> split(const char* line)
{
   char record[5000];
   vector <const char*> vec;
   strcpy(record,line);
   char* pch;
   pch = strtok(record,">");
   while(pch != NULL)
   {
      vec.push_back(pch);
      pch = strtok(NULL,">");
   }
   return vec;
}

void readsnpfile1(char* line1, int suppreads)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chromosome"))
           return;
        else
        {
           reccount++;
           char recline1[10000];
           strcpy(recline1,line1);
           char* pch;
           vector<const char*> sdf;
           string s;
           pch = strtok (recline1,"\t");
           while(pch != NULL) 
           {
              sdf.push_back(pch);
              pch = strtok (NULL,"\t");
           }
           vector <const char*> basevec = split(sdf[2]);
           if((strcmp(basevec[0],"A") == 0 && strcmp(basevec[1],"G") == 0)||(strcmp(basevec[0],"G") == 0 && strcmp(basevec[1],"A") == 0)||(strcmp(basevec[0],"C") == 0 && strcmp(basevec[1],"T") == 0)||(strcmp(basevec[0],"T") == 0 && strcmp(basevec[1],"C") == 0) )
           {
               transitioncount++;
           }
           else
           {
              transversioncount++;
           }
         }
     }
}

int main(int argc ,char* argv[]) 
{
   ifstream snpfile1 (argv[1]);
   if(!snpfile1.is_open())
   {
      cout<<"Unable to open the snpfile1\n";
      exit(1);
   }
   while(snpfile1) 
   {
          char line1[10000];
          snpfile1.getline(line1,10000);
          if(snpfile1)
            readsnpfile1(line1,atoi(argv[1]));
   }   
   snpfile1.close();
   char outfile[1000];
   strcpy(outfile,argv[1]);
   strcat(outfile,"TransTranv_summary.txt");
   ofstream out (outfile);
   out << "Transition = "<< (double) transitioncount*100 /(double) reccount << endl;
   out << "Transversion = "<< (double) transversioncount*100 /(double) reccount << endl;
   out.close();
   return 0;
}
