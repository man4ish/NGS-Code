#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <string.h>
using namespace std;
class RangeSearch
{
    char featurefname[500];
    char binfilename[500];
    char indexfilename[5000];
    size_t initindex;
    int reccount;
    public:
    void Setreccount()
    {
      reccount = 0;
    }

    int Getreccount()
    {
      return reccount ;
    }

    void SetInitialindex(size_t sindex)
    {
         initindex= sindex;
    }

    size_t GetInitialindex()
    {
         return initindex;
    }

    void Setrecfile (const char* indexrecfilename)
    {
       strcpy(featurefname,indexrecfilename);
    }

    char* Getrecfile()
    {
       return featurefname;
    }
     
    vector <string> scaling(vector <string> vec, int scalesize)
    {
        vector <string> scaled_vec;
        int diff = vec.size() - scalesize ;
        if(diff > 0)
        {
           int mid = (vec.size()-1)/2;
           int ii = 0;
           
           for(; ii < (mid-1); ii ++)
           {
              scaled_vec.push_back(vec[ii]);
           }
           ii += diff;
           for(; ii < vec.size(); ii++)
              scaled_vec.push_back(vec[ii]);
            return scaled_vec;
         }   
         else
           return  vec;
    }
 

   
    vector<string> Getrecord(const char* featurefilename, size_t start, size_t stop, int recordcount )
    {
        int recsize ;
        ifstream featurefile(featurefilename, ios::in | ios::binary);
        if(!featurefile.is_open())
        {
            cout<<"could not open feature file\n";
            exit(1);
        }
             featurefile.seekg (0, ios::end);
             int filesize;
             bool recflag = 1;
             filesize = featurefile.tellg();
             vector<string> records;
             int recindex;
             int interval;
             
             interval =  Getreccount()/(recordcount + 50  );
             if(interval == 0)
             interval = 1 ;
             recindex = GetInitialindex();
             char line[1000],rec[1000];
             while(1)
             { 
                featurefile.seekg(recindex, ios::beg);
                featurefile.getline(line, 1000);
                if(recflag == 1)
                {
                   recsize = strlen(line);
                   recflag = 0;
                }
                if(strlen(line) == 0)
                {
                    if(records.size()!= 0)
                    {
                       
                       filesize -= recsize;
                       char lline[1000];
                       ifstream sfile(featurefilename, ios::in | ios::binary);
                       sfile.seekg(filesize-1, ios::beg);
                       sfile.getline(lline, 1000);
                       records.push_back(lline);
                       sfile.close(); 
                    }
                    return records;
                }
                strcpy(rec,line);
                recindex += interval*(recsize + 1);  
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
                {  
                   records.push_back(rec);
                }
             }
             if(records.size() != 0)
             {
                   int nextindex = recindex - (recsize + 1) *interval;
                   int lastindex = recindex - 2*(recsize + 1) *interval;
                   if(nextindex > filesize - recsize)
                   nextindex = filesize -recsize;
                   int begin;
                   begin = nextindex;
                   while(begin > lastindex)
                    {
                        featurefile.seekg(begin, ios::beg);
                        char lastline[1000], lastrec[1000];
                        featurefile.getline(lastline, 1000);
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
                          records.push_back(lastrec);
                          return records;
                        } 
                        begin -= (recsize +1);
                    } 
              }                          
            featurefile.close();        
            return records;
   }

  void Loadbinfile(const char* binfilename, size_t start, size_t stop)
  {
       ifstream in(binfilename);
       if(!in.is_open())
       {
          cout<<"could not open bin file\n";
          exit(1);
       }
       size_t llimit,ulimit;
       size_t initialindex;
       size_t max;
       int count;
       bool flag = false;
       bool readrecflag = false;
       string mapquery;
       vector <size_t> vecmax,vecindex;
       while(in>>llimit>>ulimit>>initialindex >> count >> max )
       {
           size_t llimitindex = llimit;
           int temprange = start/50000;
	   int maxrange = stop/50000;
	   size_t minrange = temprange*50000;
           vecmax.push_back(max);
           vecindex.push_back(initialindex);
           while(flag == false && llimit > minrange )
	   {
			 minrange += 50000;
	   }
           if(minrange == llimit)
           {
               SetInitialindex(initialindex);
               for (int i =vecindex.size()-1 ; i > 0; i--)
               { 
                   if(vecmax[i] > start)
                   {
                      SetInitialindex(vecindex[i]);
                   }
               }
               flag = true;
               readrecflag = true;
           }         
           int adjslimit = temprange;
           int elimit = maxrange + 1;
           int adjelimit = elimit*50000;           
           if(readrecflag == true && llimit >= adjslimit && ulimit <= adjelimit )
           {
                  reccount += count;
           } 
       }
       in.close();
   }

 void Setparameter(const char* species, const char* sampletype,const char* feature, const char* chr,size_t startpos, size_t stoppos)
 {
     char flname[500],binfilename[500];
     strcpy(flname,"../project/");
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
     Loadbinfile(binfilename, startpos, stoppos);
     Setrecfile(flname);
  }
};

int main(int argc, char* argv[] )
{
    RangeSearch r; 
    r.Setparameter(argv[1],argv[2],argv[3],argv[4],(size_t) atol(argv[5]), (size_t)atol(argv[6]));
    vector<string> feature = r.Getrecord(r.Getrecfile(),(size_t)atol(argv[5]),(size_t)atol(argv[6]),atoi(argv[7]));
    vector <string> featuree = r.scaling(feature,atoi(argv[7]));
    for(int i =0 ; i< featuree.size(); i++)
    {
       cout<<featuree[i]<<endl;
    }
    feature.clear();
    featuree.clear();
    return 0;
}
