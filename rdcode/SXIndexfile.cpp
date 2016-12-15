
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <fstream>
#include <dirent.h>
#include <stdio.h>
#include <sstream>
#include <map>
using namespace std;
class SXIndexFile
{
      multimap <string,string> m;
      public:
      bool endswithtxt(const char* filename)
      {
           int len = strlen(filename);
           if(filename[len-1] == 't' && filename[len-2] == 'x' && filename [len-3] == 't' && filename[len-4] == '.')
              return 1;
           else 
              return 0;
      }
      void indexfile(char* filename)
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
                  while(pch != NULL)
                  {
                        char value[100];
                        strcpy(value,filename);
                        strcat(value,":");
                        char* buffer;
                        std::stringstream ss;
                        ss<< index;
                        strcat(value,(ss.str()).c_str());          
                        m.insert(multimap<string,string>::value_type(pch,value));
                        pch = strtok(NULL,"\t");
                  }
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
                                       indexfile(filname); 
                                    }
                               }
                               closedir(pDIR1);
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
};
int main(int argc, char* argv[])
{
    SXIndexFile I;
    I.ListFiles(argv[1]);
    multimap<string,string> IndexMap = I.GetIndexmap();
    for (multimap<string, string>::iterator it = IndexMap.begin();it != IndexMap.end(); it++)
    {
       cout <<(*it).first << "\t" << (*it).second << endl;
    }
    return 0;
}



