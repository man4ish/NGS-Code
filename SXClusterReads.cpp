#include <iostream>
#include <map>
#include <fstream>
using namespace std;
using namespace std;

struct ltstr
{
    bool operator()(string _lhs, string _rhs) const
    {
        return strcmp(_lhs.c_str(), _rhs.c_str()) < 0;
    }
};

typedef map<string, int,ltstr> readstats;
readstats rd;
typedef multimap<string,string,ltstr> Cluster_Map;
Cluster_Map  C_Map; 
int rep_count = 0;
int uni_count = 0;


void manual()
{
    cout<< "\n\nSXClusterReads\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXClusterReads <query file>\n\n";
}


const char* Getsequence(const char* seq)
{
      int i =0;
      int seqend;
      for(;i<strlen(seq);i++)
      {
          if(*(seq+i)=='\n' && *(seq+i+1)== '+')
          {
             seqend = i;
          }
      }
      char* sequence = new char[seqend+1];
      memcpy(sequence, seq,seqend);
      sequence[seqend]='\0'; 
      return sequence;     
}

void Loadreads(string reads)
{
    if(rd[reads]==0)
       rd[reads]=1;
    else 
       rd[reads]++;
}
 
void GenerateStatistics()
{
     for(readstats::iterator it = rd.begin(); it !=rd.end();++it)
     {
         if(it->second == 1)
            uni_count++;
          else
           rep_count++; 
         cout<<it->second<<"\t"<<it->first<<"\t"; 
         pair<multimap<string, string,ltstr>::iterator, multimap<string, string,ltstr>::iterator> ppp;
         ppp = C_Map.equal_range(it->first);
         for (multimap<string,string,ltstr>::iterator it2 = ppp.first; it2 != ppp.second;++it2)
         {
             cout << (*it2).second <<" ";
         } 
         cout<<endl; 
     }
}
 
int  GetRepeatcount()
{
     return rep_count;
}

int  GetUniquecount()
{
     return uni_count;
}

int main(int argc, char* argv[])
{
    if(argc!= 2)
    {
       manual();
       exit(1);
    } 

               char qline[1024] = "\0";
               char qHeader[1024] = "\0";
               char qSequence[65536] = "\0";
               FILE * queryFile = fopen(argv[1], "r");
               char flname[200];
               strcpy(flname, argv[1]);
               strcat(flname,".log");
               ofstream out(flname);
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
 
                          
                               Loadreads(qSequence);
                               if(qHeader[strlen(qHeader)-1]=='\n')
                                  qHeader[strlen(qHeader)-1]='\0';
                              
                              C_Map.insert(pair<string,string>(qSequence,qHeader));     
                              qSequence[0] = '\0';
                              strcpy(qHeader, qline);
                          }
                     }
                     else
                     {
                         strcat(qSequence, qline);
                     }
              }
       fclose(queryFile);
       GenerateStatistics();
out<<"File Name              : "<<argv[1]<<endl;
out<<"No of repeat sequences : "<<GetRepeatcount()<<endl;
out<<"No of Unique sequences : "<<GetUniquecount()<<endl;
C_Map.clear();
rd.clear();
return 0; 
}

