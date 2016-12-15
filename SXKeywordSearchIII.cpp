#include <iostream>
#include <fstream>
#include <map>
#include <vector>
using namespace std;
multimap<string,string> m;
vector<string> split(char* str, const char* token)
{
     char* pch;
     vector<string> vec;
     pch = strtok (str,token);
     while (pch != NULL)
     {
          vec.push_back(pch);  
          pch = strtok (NULL, token);
     }     
     return vec;
}

vector<string>  SearchKeyword(const char* keyword)
{ 
   vector<string> recvec;
   if(m.count(keyword) != 0)          // exact keyword search
   { 
      multimap<string,string>::iterator iter;
      pair<multimap<string,string>::iterator,multimap<string,string>::iterator> ret;
      ret = m.equal_range(keyword);
      for (iter =ret.first; iter != ret.second; iter++)
      {   
           recvec.push_back(iter->second);
      }
   }
   else                               // partial keyword search
   {
      for (multimap<string, string>::iterator it = m.begin(); it != m.end(); ++it)
      {
           if(strstr((it->first).c_str(),keyword)) 
           {
              recvec.push_back(it->second);
           }
      }
   }
   return recvec;
}
const char* GetRecord(string filename)
{
    char flname[1000];
    strcpy(flname,filename.c_str());
    vector<string> vecb = split(flname,":");
    ifstream in(vecb[0].c_str());
    char line[1000];
    in.seekg(atoi((vecb[1]).c_str()), ios::beg);
    in.getline(line, 1000);
    return line;
}
void LoadIndexMap(const char* filename)
{  
     ifstream in(filename);
    char line[1000];
    while(in)
    {
         in.getline(line,1000);
         if(in)
         {
            char* pch;
            vector<const char* > vec;
            pch = strtok(line,"\t");
            while(pch != NULL)
            { 
                 vec.push_back(pch);
                 pch = strtok(NULL,"\t");
            }
            m.insert(pair<string, string>(vec[0],vec[1]));
         }
    }
    in.close();

}
int main(int argc, char* argv[])
{
    LoadIndexMap(argv[1]);
    vector<string> vecc = SearchKeyword(argv[2]);
    for(int i =0; i < vecc.size(); i++ )
    { 
       cout<< GetRecord(vecc[i]) <<endl;;
    }
     m.clear();
}

