#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
class RangeSearch
{
    char flname[500];
    char  binfilename[500];
    char indexfilename[5000];
    size_t mapquery;
    int reccount;
    public:
    void Setreccount()
    {
      reccount = 0;
    }
    int Getcount()
    {
      return  reccount ;
    }
    void SetMapquery(size_t query)
    {
         mapquery= query;
    }
    size_t GetMapQuery()
    {
         return mapquery;
    }

    vector<string> Get_Record(const char* filename, size_t start, size_t stop, int Recordcount )
    {
        int recsize ;
        ifstream snpfile(filename, ios::in | ios::binary);
        if(!snpfile.is_open())
        {
            cout<<"could not open feature file\n";
            exit(1);
        }
             snpfile.seekg (0, ios::end);
             int filesize;
             bool recflag = 1;
             filesize = snpfile.tellg();
             vector<string> Records;
             int inde;
             int interval;            
             interval =  Getcount()/(Recordcount + 50);
             if(interval == 0)
             interval = 1 ;
             inde = GetMapQuery();
             char line[1000],rec[1000];
             while(1)
             { 
                snpfile.seekg(inde, ios::beg);
                snpfile.getline(line, 1000);
                if(recflag == 1)
                {
                   recsize = strlen(line);
                   recflag = 0;
                }
                if(strlen(line) == 0)
                 {
                    if(Records.size()!= 0)
                    {
                       filesize -= recsize;
                       char lline[1000];
                       ifstream sfile(filename, ios::in | ios::binary);
                       sfile.seekg(filesize-1, ios::beg);
                       sfile.getline(lline, 1000);
                       Records.push_back(lline);
                       sfile.close(); 
                    }
                    return Records;
                 }
                strcpy(rec,line);
                inde += interval*(recsize + 1);  
                char* pch;
                pch =  strtok(line,"\t");
                vector <const char* > vec;
                while(pch!= NULL)
                {
                    vec.push_back(pch);
                    pch = strtok(NULL,"\t");
                }
                if((size_t)atol(vec[0]) >= stop )
                {
                    break;
                }
                else if((size_t)atol(vec[1]) >= start )   
                  Records.push_back(rec);
             }
             if(Records.size() != 0)
             {
                  int nextindex = inde - (recsize + 1) *interval;
                  int  lastindex = inde - 2*(recsize + 1) *interval;
                  if(nextindex > filesize - recsize)
                   nextindex = filesize -recsize;
                   int begin;
                   begin = nextindex;
                   while(begin > lastindex)
                    {
                        snpfile.seekg(begin, ios::beg);
                        char lastline[1000], lastrec[1000];
                        snpfile.getline(lastline, 1000);
                        strcpy(lastrec, lastline);
                        char* pchh;
                        pchh =  strtok(lastline,"\t");
                        vector <const char* > veclast;
                        while(pchh!= NULL)
                        {
                             veclast.push_back(pchh);
                             pchh = strtok(NULL,"\t");
                        }
                        if(atoi(veclast[1]) <= stop)
                        {
                          Records.push_back(lastrec);
                          return Records;
                        } 
                        begin -= (recsize +1);
                    } 
              }                          
        snpfile.close();        
        return Records;
   }

  void load_binfile(const char* binfilename, size_t start, size_t stop)
  {
       ifstream in(binfilename);
       if(!in.is_open())
       {
          cout<<"could not open bin file\n";
          exit(1);
       }
       size_t llimit,ulimit;
       size_t gid;
       int count;
       string mapquery;
       while(in>>llimit>>ulimit>>gid >> count )
       {
           size_t query = start;
           if(query >= llimit && query <= ulimit )
           {
                 SetMapquery(gid);
           }        
           if(start >= llimit || stop >= ulimit )
           {
                reccount += count;
           }
       }
       in.close();
   }

  void Setrecfile (const char* indexrecfilename)
  {
       strcpy(flname,indexrecfilename);
  }

  char* GetRecFile()
  {
       return flname;
   }
 
 void SetParameter(char* species, char* sampletype,char* feature, char* chr,size_t startpos, size_t stoppos)
 {
     char flname[500],binfilename[500];
     strcpy(flname,"../species/");
     strcat(flname,species);
     strcat(flname,"/");
     strcat(flname,sampletype);
     strcat(flname,"/");
     strcat(flname, feature);
     strcat(flname, "/");
     strcat(flname,species);
     strcat(flname,".");
     strcat(flname,sampletype);
     strcat(flname,".");
     strcat(flname, feature);
     strcat(flname, ".");
     strcat(flname,chr);
     strcat(flname,".txt");
     strcpy(binfilename,flname);
     strcat(binfilename,".bin");
     Setreccount();
     load_binfile(binfilename, startpos, stoppos);
     Setrecfile(flname);
  }
};

int main(int argc, char* argv[] )
{
     RangeSearch r; 
     r.SetParameter(argv[1],argv[2],argv[3],argv[4],(size_t) atol(argv[5]), (size_t)atol(argv[6]));
     vector<string> feature =  r.Get_Record(r.GetRecFile(), (size_t)atol(argv[5]), (size_t) atol(argv[6]), atoi(argv[7]));
     for(int i =0 ; i< feature.size(); i++)
     cout<<feature[i]<<endl;
     feature.clear();
     return 0;
}
