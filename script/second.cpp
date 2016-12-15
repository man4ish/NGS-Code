#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <set>
using namespace std;
int main(int argc, char* argv[])
{
  ifstream in(argv[1]);
  if(!in.is_open())
  {
     cout<< "could not open the file\n";
     exit(1);
  }  
  multimap <string,vector<string> > m;
  set<string> vect;
  char line[1000];
  while(in)
  {
     in.getline(line,1000);
     if(in)
     {
        char* pch; 
        pch = strtok(line,"\t" );
        vector<string> vec;
        while(pch != NULL)
        {
            vec.push_back(pch);
            pch = strtok(NULL, "\t");
        }
        vect.insert(vec[3]);
        m.insert(pair<string, vector<string> >(vec[3],vec));
     }
  }
     set<string>::iterator theIterator;
     multimap <string, double> mn;
     for( theIterator = vect.begin(); theIterator != vect.end(); theIterator++ )
     {
         pair<multimap<string, vector<string> >::iterator, multimap<string, vector<string> >::iterator> ppp;
         ppp = m.equal_range(*theIterator);
         vector <string> vecbg;
         double sum = 0; 
         string a ; 
         set<string> vect2;
         for (multimap<string, vector<string> >::iterator it2 = ppp.first;it2 != ppp.second; ++it2)
         {
              a = (*it2).first;
              sum += (double) atof((((*it2).second)[15]).c_str());
          }
         
           mn.insert(pair<string,double>(a,(double)sum/(double)m.count(*theIterator) ));
     }

    for (multimap<string, double>::iterator it = mn.begin();it != mn.end(); ++it)
    {
         pair<multimap<string, vector<string> >::iterator, multimap<string, vector<string> >::iterator> ppp;
         ppp = m.equal_range((*it).first);
        // cout<< (*it).first <<endl;
         for (multimap<string, vector<string> >::iterator it2 = ppp.first;it2 != ppp.second; ++it2)
         { 
            for(int i = 0; i < (*it2).second.size(); i++)
            {
               cout<< (*it2).second[i]<< "\t";
            }
              cout <<  (*it).second<<endl;
              //  cout << (*it).first <<endl;
         }
    }
}

