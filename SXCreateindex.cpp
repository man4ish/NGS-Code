#include <cstdlib>
#include <iostream>
#include <string.h>
#include <fstream>
#include <dirent.h>
#include <stdio.h>
#include <sstream>
#include <map>
#include <vector>
using namespace std;
class SXCreateFileIndex
{
      multimap <string,string> m;
      multimap <string,string> patternmap;
      public:
      bool patternflag;

      void manual()
      {
           cout<<"\n\nSXCreateindex\n";
           cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
           cout<<"Usage:\n";
           cout<<"\tSXCreateindex option <input_directory> <indexfile_prefix> \n\n";
      }

      bool endswithtxt(const char* filename)
      {
           int len = strlen(filename);
           if(filename[len-1] == 't' && filename[len-2] == 'x' && filename [len-3] == 't' && filename[len-4] == '.')
              return 1;
           else 
              return 0;
      }

      void GetPattern(char* str, char* index)
      {
            string strr = str;
            if(strr.size() > 2 )
            {
                for (int i = 0; i< strr.size() -2;i++)
	        {
                     string pattern = strr.substr(i,3);   
                     if(strstr(pattern.c_str(),".") == NULL && strstr(pattern.c_str(),"-") == NULL && strstr(pattern.c_str(),"_") == NULL && pattern != "" )
                     {
                            patternmap.insert(multimap<string,string>::value_type(ConvertTolowercase(pattern),index));
                      }
                 }
              }
      }

     char *trim(char *str)
     {
         char ptr[strlen(str)+1];
         int i,j=0;
         for(i=0;str[i]!='\0';i++)
         {
             if(str[i] != ' ' && str[i] != '\t')
                ptr[j++] = str[i];
         }
         ptr[j]= '\0';
         return ptr;
      }
  
     vector<string> split(const char* kw, const char* token)
     {
           char* pch;
           char str[200] = "\0";
           strcpy(str,kw);
           vector<string> vec;
           pch = strtok (str,token);
           while (pch != NULL)
           {
                  vec.push_back(pch);
                  pch = strtok (NULL, token);
           }
           return vec;
     }

      string ConvertTolowercase(string str)
      {
      for (int i=0;i<strlen(str.c_str());i++)
      if (str[i] >= 0x41 && str[i] <= 0x5A)
      str[i] = str[i] + 0x20;
      return str;
      }
 
      void Indexfile(char* filename,char* fname)
      {
           char line[1000];
           ifstream in(filename);
           int index = 0;
           if(!in.is_open())
           {
              cout<<"could not open the file\n";
           }
           while(in)
          {
                char line[1000],rec[1000];
                in.getline(line,1000);
                if(in)
                {
                  char* pch; 
                  strcpy(rec,line); 
                  pch = strtok(line,"\t");
                  char filepath[200] = "\0";
                  vector<const char*> kword, findex; 
                  while(pch != NULL)
                  {
                        char value[100];
                        string j = fname;
                        int h1 = j.find("bi.");
                        int k1 =  j.find(".txt");
                        string hk1 =  j.substr(h1+3, k1-h1-3 );
                        strcpy(value,hk1.c_str());
                       strcat(value,":");
                        char* buffer;
                        std::stringstream ss;
                        ss<< index;
                        strcat(value,(ss.str()).c_str());
                        kword.push_back(pch);
                        findex.push_back(value);  
                        pch = strtok(NULL,"\t");
                  }
                  char sf[300] = "\0";
                  strcpy(sf,findex[findex.size()-1]); 
                  vector<string> Kkeyword = split(kword[kword.size()-1],",");
                  for ( int ii = 0; ii < Kkeyword.size(); ii++)
                  {
                       if(strcmp(Kkeyword[ii].c_str(),"-") != 0 &&  strcmp(Kkeyword[ii].c_str(),"0") != 0 )
                       {
                           char keywordindex[500] = "\0";
                           strcpy(keywordindex,Kkeyword[ii].c_str() );
                           if(patternflag == 1)     
                             GetPattern(trim(keywordindex), sf); 
                           m.insert(multimap<string,string>::value_type(ConvertTolowercase(trim(keywordindex)),sf));
                       }
                  }
                  findex.clear();
                  Kkeyword.clear();
                  index += strlen(rec);
                  index += 1;
               }
          }
      }
      void ListFiles(char* drname)
      {
           DIR *pDIR;
           struct dirent *entry;
           if( pDIR=opendir(drname) )
           {
               while(entry = readdir(pDIR))
               {
                     if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 )
                     {
                         DIR *pDIR1;
                         struct dirent *entry1;
                         char flname[100]; 
                         strcpy(flname,drname);
                         strcat(flname,"/");
                         strcat(flname,entry->d_name);
                         strcat(flname,"/");
                         if( pDIR1=opendir(flname) )
                         {
                             while(entry1 = readdir(pDIR1))
                             {
                                   if(strcmp(entry1->d_name, ".") != 0 && strcmp(entry1->d_name, "..") != 0 )
                                   {
                                       char filname[100]; 
                                       strcpy(filname,drname);
                                       strcat(filname,"/");
                                       strcat(filname,entry->d_name);
                                       strcat(filname,"/");
                                       strcat(filname,entry1->d_name);
                                       if(endswithtxt(filname))
                                       Indexfile(filname,entry1->d_name); 
                                    }
                               }
                               closedir(pDIR1);
                          }
                          else 
                          { 
                               char filname[100];
                               strcpy(filname,drname);
                               strcat(filname,"/");
                               strcat(filname,entry->d_name); 
                               if(endswithtxt(filname))
                               Indexfile(filname, entry->d_name);
                          }
                      }
                }
                closedir(pDIR);
        }
  }
  multimap<string,string> GetIndexmap() const
  {
       return m;
  } 
  multimap<string,string> GetPatternmap() const
  {
      return patternmap;
  }
};
void CreateIndex(int paramcount, char* param[])
{
    bool patflag;
     SXCreateFileIndex I;
     int i =1;

     if(paramcount != 3 && paramcount != 4 )
     {
        I.manual();
        exit(1);
     }
     if(strcmp(param[i],"-p") == 0)
     {
        i++;
        I.patternflag = 1; 
     }

     I.ListFiles(param[i]);
     char patternindexfile[200], keywordindexfile[200];
     strcpy(patternindexfile, param[i]);
     strcpy(keywordindexfile, param[i]);
     i++;
     strcat(patternindexfile, "/");
     strcat(keywordindexfile, "/");
     strcat(patternindexfile, param[i]);
     strcat(keywordindexfile, param[i]);
     strcat(patternindexfile,".kw.pt");
     strcat(keywordindexfile,".kw");
     ofstream outfile1(keywordindexfile);
     multimap<string,string> IndexMap = I.GetIndexmap();
     for (multimap<string, string>::iterator it = IndexMap.begin();it != IndexMap.end(); it++)
     {
        outfile1 <<(it-> first) << "\t" << it->second <<endl;
     }
     outfile1.close();
      if(I.patternflag == 1)
     {
        multimap<string,string> PatternMap = I.GetPatternmap();
        ofstream outfile2(patternindexfile);
        for (multimap<string, string>::iterator it = PatternMap.begin();it != PatternMap.end(); it++)
        {
            outfile2 <<it->first << "\t" << it->second <<endl;
        }
        outfile2.close();
     }
}

int main(int argc, char* argv[])
{
   CreateIndex(argc, argv);
/*     bool patflag;
     SXCreateFileIndex I;
     int i =1; 
     
     if(argc != 3 && argc != 4 )
     {
        I.manual();
        exit(1);
     }   
     if(strcmp(argv[i],"-p") == 0)
     {
        i++;
        I.patternflag = 1; 
     }
 
     I.ListFiles(argv[i]);
     char patternindexfile[200], keywordindexfile[200];
     strcpy(patternindexfile, argv[i]);
     strcpy(keywordindexfile, argv[i]);
     i++;
     strcat(patternindexfile, "/");
     strcat(keywordindexfile, "/");
     strcat(patternindexfile, argv[i]);
     strcat(keywordindexfile, argv[i]);
     strcat(patternindexfile,".kw.pt");
     strcat(keywordindexfile,".kw");
     ofstream outfile1(keywordindexfile);
     multimap<string,string> IndexMap = I.GetIndexmap();
     for (multimap<string, string>::iterator it = IndexMap.begin();it != IndexMap.end(); it++)
     {
        outfile1 <<(it-> first) << "\t" << it->second <<endl; 
     }
     outfile1.close();
     if(I.patternflag == 1)
     {
        multimap<string,string> PatternMap = I.GetPatternmap();
        ofstream outfile2(patternindexfile);
        for (multimap<string, string>::iterator it = PatternMap.begin();it != PatternMap.end(); it++)
        {
            outfile2 <<it->first << "\t" << it->second <<endl;
        }
        outfile2.close();
     }
*/
     return 0;
}



