#include <iostream>
#include <map>
#include <utility>
#include <fstream>
#include <vector>
#include <set>
using namespace std;
int main(int argc,char* argv[])
{
  multimap <string,string> m;
  ifstream in(argv[1]);
  vector<string> vec;
  set<string> s;
  char line[1000],recline[1000];
  while(in)
  {
    in.getline(line,1000);
    strcpy(recline,line);
    if(in)
    {
       char* pch;
       pch = strtok(line,"\t");
       vector<const char*> vec;
       while(pch != NULL)
       {
         vec.push_back(pch);
         pch = strtok(NULL,"\t");
       }
       m.insert(pair<string,string>(vec[8],recline));
       s.insert(vec[8]);
       //cout << vec[8]<< endl;
    }
  }
  for (set<string>:: iterator it = s.begin(); it != s.end(); it++)
  {
        if(m.count(*it)> 1)
        cout << *it<<endl;
  }
  s.clear();
  vec.clear();
  in.close();
  return 0;
}
