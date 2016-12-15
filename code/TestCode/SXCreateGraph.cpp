#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <iterator>
using namespace std;
map<string,int> m;
int Getmax(vector <int> v)
{
   std::ostream_iterator<int>output(cout," ");
   return *(std::max_element(v.begin(),v.end()));
}
typedef vector <int> Inputvec;
Inputvec v1,v2; 
int main(int argc, char* argv[])
{
  char line[500];
  ifstream in(argv[1]);
  while(in)
  {
     in.getline(line,500);
     if(in)
     {
       char* pch;
       vector <const char*> vec;
       pch = strtok(line,"\t");
       while(pch != NULL)
       {
          vec.push_back(pch);  
          pch = strtok(NULL,"\t");
       }
       v1.push_back(atoi(vec[5]));
       v2.push_back(atoi(vec[14]));
       string s;
       s += vec[5];
       s += "-";
       s += vec[14];
       if(m.count(s) == 0)
          m[s] = 1;
       else 
          m[s]++;  
      } 
  }
  in.close();
  int row  = Getmax(v1);
  int col =  Getmax(v2);

  char flname[1000];
  strcpy(flname,argv[1]);
  strcat(flname,"_graph");
  ofstream outgraph(flname); 
  outgraph<< "SR-RS\t";
  for (int ii = 0;ii <= col; ii++)
  {
     outgraph<< ii << "\t";
  }
  outgraph << "\n"<< endl;
  for (int i = 1;i <= row; i++)
  {
     outgraph <<i<<"\t";
     for (int j = 0;j <= col; j++)
     {
        std::stringstream ss1,ss2;
        ss1 << i;
        ss2 <<j;
        string ss;
        ss += ss1.str();
        ss += "-";
        ss += ss2.str();
        string out;
        if(m.count(ss) == 0)
         outgraph <<"0\t";
        else 
         outgraph <<m[ss]<<"\t"; 
     }
     outgraph <<"\n";
  }
  m.clear();
  v1.clear();
  v2.clear();
  return 0;
}
