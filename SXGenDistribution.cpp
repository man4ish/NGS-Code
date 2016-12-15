#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <map>
using namespace std;
void manual()
{
    cout<< "\n\nSXGenDistribution\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXGenDistribution <query file>\n\n";
}

struct ltstr
{
    bool operator()(string  _lhs, string _rhs) const
    {
        return strcmp(_lhs.c_str(), _rhs.c_str()) < 0;
    }
};

typedef map<string ,int,ltstr > Dist_Map;
Dist_Map dmap;

void SetDistribution(char* Seq)
{
   int len = strlen(Seq);
   if(len >= 18 and len <= 28)
   {
      std::stringstream ss;
      ss << len;
      if(dmap.count(ss.str())==0)
      { 
         dmap[ss.str()] =1;
      }
      else
      {
        dmap[ss.str()] += 1;
      }
   }
  else if(len > 29)
  {
    if(dmap.count("29-36")== 0)
    {
       dmap["29-36"] = 1;
    }
    else 
    {
       dmap["29-36"] +=1;
    }
  }
}

void GetDistribution(const char* flname)
{
    ofstream out(flname);
    for(Dist_Map :: iterator it = dmap.begin(); it != dmap.end(); ++it)
        out<< (*it).first<<" : "<<(*it).second<<endl;
}
int main(int argc ,char* argv[])
{

         if(argc != 2)
         {
           manual();
           exit(1);
         }
         char qline[1024] = "\0";
         char qHeader[1024] = "\0";
         char qSequence[65536] = "\0";
         FILE * queryFile = fopen(argv[1], "r");
         char flname[200]; 
         strcpy(flname,argv[1]);
         strcat(flname,".dbn");
         while (!feof(queryFile))
         {
                fgets(qline, 1024, queryFile);
                if((qline[0] == '>') || (feof(queryFile)))
                {
                    if (qSequence[0] == '\0')
                    {
                        strcpy(qHeader, qline);
                    }
                    else
                    {      
                           if(qSequence[strlen(qSequence)-1]=='\n')
                           qSequence[strlen(qSequence)-1]='\0'; 
                           SetDistribution(qSequence); 
                           qSequence[0] = '\0';
                           strcpy(qHeader, qline);
                     }
                  }
                 else
                 {
                       strcat(qSequence, qline);
                 }
          }
          GetDistribution(flname);  
          fclose(queryFile);
          return 0;
}
