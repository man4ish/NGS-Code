#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
using namespace std;
map <string,int> m;

void ClusterGenes(const char* genenames)
{
         char gene[10000];
         strcpy(gene, genenames);
          char* pch;
          pch = strtok(gene,", ");
          vector <string> records;
          while(pch != NULL)
          {
             if(m.count(pch) != 0)
              m[pch] += 1;
             else 
              m[pch]= 1;
              pch = strtok(NULL,", ");
          }
}
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
          if(records[9] != "-")
            ClusterGenes(records[9].c_str());
         }
     }
     for(map<string,int>:: iterator it = m.begin(); it != m.end(); it++)
      cout << it->first<<"\t"<<it->second<<endl;
     return 0;
}
