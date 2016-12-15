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
    public:
    void SetMapquery(size_t query)
    {
         mapquery= query;
    }
    size_t GetMapQuery()
    {
         return mapquery;
    }
    vector<string> Get_Record(const char* filename, size_t start, size_t stop)
    {
        ifstream snpfile(filename, ios::in | ios::binary);
        if(!snpfile.is_open())
        {
            cout<<"could not open feature file\n";
            exit(1);
        }
        vector<string> Records;
             int inde;
             inde = GetMapQuery();
             char line[1000],rec[1000];
             while(1)
             { 
                snpfile.seekg(inde, ios::beg);
                snpfile.getline(line, 1000);
                if(strlen(line) == 0)
                return Records;
                strcpy(rec,line);
                inde +=  strlen(rec);
                inde += 1; 
                char* pch;
                pch =  strtok(line,"\t");
                vector <const char* > vec;
                while(pch!= NULL)
               {
                    vec.push_back(pch);
                    pch = strtok(NULL,"\t");
               }
               if((size_t)atol(vec[1]) >= stop )
                  return Records;
               else if((size_t)atol(vec[2]) >= start )   // for others it is 1.
                  Records.push_back(rec);
             }      
        snpfile.close();
   }
  void load_binfile(const char* binfilename, size_t start)
  {
       ifstream in(binfilename);
       if(!in.is_open())
       {
          cout<<"could not open bin file\n";
          exit(1);
        }
       size_t llimit,ulimit;
       size_t gid;
       string mapquery;
       while(in>>llimit>>ulimit>>gid)
       {
           size_t query = start;
           if(query >= llimit && query <= ulimit )
           {
                 SetMapquery(gid);
           }
       }
       in.close();
   }

  void SetRecFile (char* indexrecfilename)
  {
       strcpy(flname,indexrecfilename);
  }

  char* GetRecFile()
  {
       return flname;
  }
 void SetParameter(char* species, char* sampletype,char* feature, char* chr,size_t startpos)
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
     load_binfile(binfilename, startpos);
     SetRecFile(flname);
  }
};

int main(int argc, char* argv[] )
{
     RangeSearch r; 
     r.SetParameter(argv[1],argv[2],argv[3],argv[4],(size_t) atol(argv[5]));
   vector<string> feature =  r.Get_Record(r.GetRecFile(), (size_t)atol(argv[5]), (size_t) atol(argv[6]));
   for(int i =0 ; i< feature.size(); i++)
   cout<<feature[i]<<endl;
     return 0;
}
