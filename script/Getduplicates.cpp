#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
using namespace std;
map<string,string> m;
vector<string> vec;
 struct ltstr
{
    bool operator()(const char * _lhs, const char * _rhs) const
    {
        return strcmp(_lhs, _rhs) < 0;
    }
};
typedef std::map<string,string> QueryMap;
typedef QueryMap::iterator QueryMapIT;
typedef std::pair<string, string> QueryMapPAIR;
std::pair<QueryMapIT, bool> ret;
QueryMap qryMap;

void LoadSnpFile1(string line1)
{
  if(line1 != "") 
  {
     m[line1] =line1;
  }	
}

void LoadSnpFile2(string line2)
{
  if(line2 != "") 
  {
     vec.push_back(line2);
  }	
}

void GetDuplicates() 
{
  for(int i =0; i < vec.size(); i++)
  {
     ret = qryMap.insert(QueryMapPAIR(vec[i], vec[i]));
        if (ret.second == false)
             cout<<vec[i]<<endl;
           // dupQuery.push_back(new Query(qryVector[i]->queryHeader, qryVector[i]->querySequence));
  }
}


int main(int argc ,char* argv[]) 
{ 
   ifstream myfile1 (argv[1]);
    string line1,line2;
    while(!myfile1.eof() ) 
    {
          getline (myfile1,line1);
          LoadSnpFile1(line1);
    }

   ifstream myfile2 (argv[2]);
   while (!myfile2.eof() ) 
   {
           getline (myfile2,line2);
           LoadSnpFile2(line2);
   } 
    GetDuplicates();
    return 0;
}
