#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
map<int,int> m;
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
          // cout << m[3] << endl;
           if( atoi(sdf[8].c_str()) <= m[atoi(sdf[0].c_str())] )
           // cout << atoi(sdf[5].c_str())<<"\t"<< m[atoi(sdf[0].c_str())]<< "\t"<< atoi(sdf[0].c_str())<<endl;
            cout << line1<< endl;
         }
     }
}

int main(int argc, char* argv[])
{
  int chr =1;
  ifstream in1(argv[1]);
  char line1[5000],recline2[5000];
  while(in1)
  {
     in1.getline(line1,5000);
     if(in1)
     {
   /*        char recline1[10000];
           strcpy(recline1,line1);
           char* pch1;
           vector<const char*> vec1;
           pch1 = strtok (recline1,"\t");
           while(pch1 != NULL)
           {
              vec1.push_back(pch1);
              pch1 = strtok (NULL,"\t");
           } */
          m[chr] = atoi(line1);
          
     }
      chr++;
  }
  in1.close();


  ifstream in(argv[2]);
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
