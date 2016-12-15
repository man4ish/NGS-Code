#include <iostream>
#include <map>
#include <fstream>
#include <vector>
using namespace std;
class KeywordSearch
{
    string mapquery;
    char ffile[100], fbinfile[100];
    char featurefile[100], featurebinfile[100];
    char gifile[100],gibinfile[100];
    size_t startpos, endpos;
    public:
    void SetMapquery(string query)
    {
         mapquery= query;
    }
    string GetMapQuery()
    {
         return mapquery;
    }
    void Setgifile(char* gifilename)
    {
          strcpy(gifile,gifilename );
    }
    const char* Getgifile()
    {
          return gifile;
    }
    void Setgibinfile(char* gibinfilename)
    {
          strcpy(gibinfile,gibinfilename );
    }
    const char* Getgibinfile()
    {
          return gibinfile;
    }

    void Setfeaturefile(char* featurefile)
    {
       strcpy(ffile,featurefile);
    }
    void Setfeaturebinfile(char* featurebinfile)
    {
       strcpy(fbinfile,featurebinfile);
    }
    const char* Getfeaturefile ()
    {
       return ffile;
    }
    const char* Getfeaturebinfile()
    {
       return fbinfile;
    }
    void Setstartpos(const char* pos)
    {
       startpos = (size_t)atol(pos);
    }
    
    void Setendpos(const char* pos)
    {
       endpos = (size_t)atol(pos);
    }

    size_t Getstartpos ()
    {
       return startpos;
    }
    size_t  Getendpos()
    {
       return endpos;
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
        string gid;
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

    vector<string> Get_Record(const char* filename, size_t start, size_t stop)
    {
        ifstream snpfile(filename, ios::in | ios::binary);
        if(!snpfile.is_open())
        {
            cout<<"could not open feature file\n";
            exit(1);
        }
        vector<string> Records;
             size_t inde;
             inde = (size_t) atol (GetMapQuery().c_str()) ;
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

    void SetParameter( char* species , char* sample, char* feature, char* search_type, char* id )
    {
         char flname[500],gifilename[500],giindex[500],featurefile[500],featurebinfile[500],featureindexfile[500];
         strcpy(flname,"../species/");
         strcat(flname,species);
         strcat(flname,"/");
         strcat(flname,sample);
         strcpy(featurefile,flname);
         strcat(featurefile,"/");
         strcat(featurefile,search_type);
         strcat(featurefile,"/");
         strcat(featurefile,species);
         strcat(featurefile,".");
         strcat(featurefile,sample);
         strcat(featurefile,".");
         strcat(featurefile,search_type);
         strcat(featurefile,".");
         strcat(featurefile,"chr");
         Setfeaturefile(featurefile);
         strcpy(gifilename,flname);
         strcat(gifilename,"/keyword/");
         strcat(gifilename,species);
         strcat(gifilename,".");
         strcat(gifilename,sample);
         strcat(gifilename,".");
         strcat(gifilename,feature);
         strcat(gifilename,".txt");
         strcpy(giindex,gifilename);
         strcat(giindex,".index");
         Setgifile(gifilename);
         Setgibinfile(giindex);
     }

   void GetRec(char* ginum)
   {
         map<string ,size_t > df;
         ifstream GI(Getgibinfile());
         string idd;
         size_t inde, maxrange;
         ifstream sfile(Getgifile(), ios::in | ios::binary);
         char lin[1000],reco[1000];
         while(GI >> idd >>inde)
         {
          if(strstr(idd.c_str(),ginum))
          {
             sfile.seekg(inde,ios::beg);
             sfile.getline(lin, 1000);
             if(strlen(lin)== 0)
             exit(1); 
             strcpy(reco,lin);
             char* phch;
             phch =  strtok(lin,"\t");
             vector <const char* > vech;
             while(phch != NULL)
             {
                vech.push_back(phch);
                phch = strtok(NULL,"\t");
             }
            strcpy(featurefile, Getfeaturefile());
             strcat(featurefile,vech[1]);
             strcat(featurefile,".txt");
             strcpy(featurebinfile,featurefile);
             strcat(featurebinfile,".bin");        
         load_binfile(featurebinfile,(size_t)atol(vech[2]));
         vector<string> ved =  Get_Record(featurefile, (size_t) atol(vech[2]), (size_t) atol(vech[3]));
         for(int i =0; i< ved.size(); i++)
           cout<<ved[i]<<endl; 
       }
      }
   sfile.close();
  }
};

int main(int argc, char* argv[] )
{
    KeywordSearch k ;
    if(argc != 6)
     {
        cout<<"Check the number of parametrs\n";
        exit(1);
     }
     k.SetParameter( argv[1] , argv[2], argv[3] , argv[4], argv[5]);
     k.GetRec(argv[5]);
     return 0;
}
