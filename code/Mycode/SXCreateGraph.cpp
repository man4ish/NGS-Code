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
       v1.push_back(atoi(vec[0]));
       v2.push_back(atoi(vec[1]));
       string s;
       s += vec[0];
       s += "-";
       s += vec[1];
       if(m.count(s) == 0)
          m[s] = 1;
       else 
          m[s]++;  
      } 
  }
  in.close();
  int row  = Getmax(v1);
  int col =  Getmax(v2); 
  cout<< "SR-RS\t";
  for (int ii = 0;ii <= col; ii++)
  {
     cout << ii << "\t";
  }
  cout << "\n"<< endl;
  for (int i = 1;i <= row; i++)
  {
     cout <<i<<"\t";
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
         cout <<"0\t";
        else 
         cout <<m[ss]<<"\t"; 
     }
     cout <<"\n";
  }
  m.clear();
  v1.clear();
  v2.clear();
  return 0;
}
