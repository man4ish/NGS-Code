#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chromosome"))
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
           cout << (float) atoi(sdf[0].c_str())/float(1000) <<"\t" << (float)atoi(sdf[1].c_str())/float(1000) <<"\t"<< sdf[2]<< endl;
         }
     }
}

int main(int argc, char* argv[])
{
  ifstream in(argv[1]);
  char line[5000],recline[5000];
  while(in)
  {
     in.getline(line,5000);
     if(in)
     {
        readsnpfile1(line);
     }
  }
  in.close();
}
