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
           if(atoi(sdf[3].c_str()) > 1)
          // int rs = (int)((double) (atoi(sdf[5].c_str())*100)/(double)(atoi(sdf[8].c_str())));  
           cout << line1<< endl;
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
