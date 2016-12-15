/*
 * File:   SXCreateFileIndex.3.0.0.cpp
 * Author: manish
 *
 * Created on October 21, 2010, 9:54 AM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <stdio.h>
#include <sstream>
#include <map>
#include <vector>
#include <ctype.h>
using namespace std;
class SXCreateFileIndex
{
      struct ltstr
      {
        bool operator()(const char* _lhs, const char* _rhs) const
        {
          return strcmp(_lhs, _rhs) < 0;
        }
      };
      typedef multimap <const char*,const char*,ltstr> PMap;
      PMap m,patternmap;
      public:
      bool patternflag;
      void manual()
      {
           cout<<"\n\nSXCreateindex\n";
           cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
           cout<<"Usage:\n";
           cout<<"\tSXCreateindex option <input_directory> <indexfile_prefix> \n\n";
      }

      const char*  ConvertTolowercase(const char* str)
      {
         char* temp = new char[strlen(str)+1];
         for (int i=0;i<strlen(str);i++)
         {
            temp[i] = tolower(str[i]);
         }
         temp[strlen(str)] = '\0';
         return temp;
      }

      bool Endswithtxt(const char* filename)
      {
           int len = strlen(filename);
           if(filename[len-1] == 't' && filename[len-2] == 'x' && filename [len-3] == 't' && filename[len-4] == '.')
              return 1;
           else
              return 0;
      }

      void GetPattern(const char* str, const char* index)
      {
            char strr[4000];
            strcpy(strr,str);
            char *pattern;
            if(strlen(strr) > 2 )
            {
                for (int i = 0; i< strlen(strr)-2;i++)
	        {
                     pattern = new char[4];
                     strncpy(pattern,strr+i,3); pattern[3]=0;
                     if(isalpha(pattern[0]) && isalpha(pattern[1]) && isalpha(pattern[2]))
                     {
                            patternmap.insert(PMap::value_type(ConvertTolowercase(pattern),index));
                     }
                 }
              }
       }

      char *trim(char *str)
      {
       char *ibuf = str, *obuf = str;
       int i = 0, cnt = 0;
       if (str)
       {
            for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf);
            if (str != ibuf)
                  memmove(str, ibuf, ibuf - str);
            while (*ibuf)
            {
                  if (isspace(*ibuf) && cnt)
                        ibuf++;
                  else
                  {
                        if (!isspace(*ibuf))
                              cnt = 0;
                        else
                        {
                              *ibuf = ' ';
                              cnt = 1;
                        }
                        obuf[i++] = *ibuf++;
                  }
            }
            obuf[i] = NULL;
            while (--i >= 0)
            {
                  if (!isspace(obuf[i]))
                        break;
            }
            obuf[++i] = NULL;
       }
       return str;
     }

     vector<const char*> split(const char* kw, const char* token)
     {
           char* pch;
           char str[4000];
           strcpy(str,kw);
           vector<const char*> vec;
           pch = strtok (str,token);
           while (pch != NULL)
           {
                  if(pch != NULL && strlen(pch)> 2 )
                  vec.push_back(pch);
                  pch = strtok (NULL, token);

           }
           return vec;
     }

      void Indexfile(char* filename)
      {
           char line[5000];
           ifstream in(filename); char *indexflname;
           int index = 0;
           if(!in.is_open())
           {
              cout<<"could not open "<<filename<<"\n";
           }
           while(in)
           {
                char line[5000],rec[5000];
                in.getline(line,5000);
                if(in)
                {
                  char* pch;
                  strcpy(rec,line);
                  pch = strtok(line,"\t");
                  char filepath[500] = "\0";
                  vector<const char*> kword, findex;
                  while(pch != NULL)
                  {
                        char value[4000];
                        char * chr_pos;
                        char *chr_ptr = strstr(filename, "bi.")+3;
                        chr_pos = strchr(chr_ptr,'.');
                        char chrnum[3];
                        strncpy(chrnum, chr_ptr, chr_pos-chr_ptr);
                        chrnum[chr_pos-chr_ptr] = '\0';
                        strcpy(value,chrnum);
                        strcat(value,":");
                        char str_index[100];
                        sprintf(str_index, "%d", index);
                        strcat(value,str_index);
                        kword.push_back(trim(pch));
                        findex.push_back(value);
                        pch = strtok(NULL,"\t");
                  }
                  indexflname = new char[4000];
                  strcpy(indexflname,findex[findex.size()-1]);
                  vector<const char*> Kkeyword = split(kword[kword.size()-1],"`|");

                  for ( int ii = 0; ii < Kkeyword.size(); ii++)
                  {
                       if(strcmp(Kkeyword[ii],"-") != 0 &&  strcmp(Kkeyword[ii],"0") != 0 )
                       {
                           char keywordindex[4000];
                           strcpy(keywordindex,Kkeyword[ii] );
                           if(patternflag == 1)
                             GetPattern(keywordindex, indexflname);
                           m.insert(PMap::value_type(ConvertTolowercase(keywordindex),indexflname));
                       }
                  }
                  findex.clear();
                  Kkeyword.clear();
                  index += strlen(rec);
                  index += 1;
               }
          }
      }

      void ListFiles( const char* drname)
      {
           DIR *pDIR;
           struct dirent *entry;
           if( pDIR=opendir(drname) )
           {
               while(entry = readdir(pDIR))
               {
                    if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 )
                     {
                         char flname[1000];
                         strcpy(flname,drname);
                         strcat(flname,"/");
                         strcat(flname,entry->d_name);
                         if(Endswithtxt(flname))
                         {
                            Indexfile(flname);
                         }
                      }
                }
                closedir(pDIR);
        }
  }

  PMap GetIndexmap() const
  {
       return m;
  }

  PMap GetPatternmap() const
  {
      return patternmap;
  }

  void CreateIndex(int paramcount, const char* param[])
  {
     bool patflag;
     int i =1;
     if(paramcount != 3 && paramcount != 4 )
     {
        manual();
        exit(1);
     }
     if(strcmp(param[i],"-p") == 0)
     {
        i++;
        patternflag = 1;
     }
     ListFiles(param[i]);
     char patternindexfile[2000], keywordindexfile[2000];
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
     for (PMap::iterator it = m.begin();it != m.end(); it++)
     {
        outfile1 <<(it-> first) << "\t" << it->second <<endl;
     }
     outfile1.close();
     if(patternflag == 1)
     {
        ofstream outfile2(patternindexfile);
        for (PMap::iterator it = patternmap.begin();it != patternmap.end(); it++)
        {
            outfile2 <<it->first << "\t" << it->second <<endl;
        }
        outfile2.close();
     }
   }
   void FreeMemory()
   {
       m.clear();
       patternmap.clear();
   }
};

int main(int argc,const char* argv[])
{
   SXCreateFileIndex I;
   I.CreateIndex(argc, argv);
   I.FreeMemory();
   return 0;
}




