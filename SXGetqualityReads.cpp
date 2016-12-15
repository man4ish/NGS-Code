#include  <map>
#include  <vector>
#include  <algorithm>
#include  <string>
#include  <iostream>
#include  <stdio.h>
#include  <time.h>
#include <fstream>
using namespace std;

typedef vector<string> V;
typedef multimap <string,vector<string> > M;
map <int,string> ma;
map <string, string> ha;
M m;
V v,vec;

void manual()
{
    cout<<"\n\nSXGetQualityreads\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXGetQualityreads <Latest_row_format>\n\n";
}



vector<string> splitv(char* s,const char* token)
{
   char* pch;
   vector<string> v;
   pch = strtok(s,token);
   while(pch!= NULL)
   {
       v.push_back(pch);
       pch = strtok (NULL,token);
   }
   return v;
}

string getbiggest( string a)
{
   pair<multimap<string,vector<string> >::iterator, multimap<string, vector<string> >::iterator> ppp;
   string read;
   if(m.count(a)>0)
   {
      ppp = m.equal_range(a);
      for (multimap<string,vector<string> >::iterator it2 = ppp.first;  it2 != ppp.second; ++it2)
      {
           string h;
           if(((*it2).second).size()==14||((*it2).second).size()==15)
           {
                for( int i =0; i< ((*it2).second).size()-1;i++)
                {
                      h += ((*it2).second)[i];
                      h +="\t";
                }
           
           h += ((*it2).second)[((*it2).second).size()-1];
           ma[atoi(((*it2).second)[4].c_str())] = h;
           }
       }

      map<int,string>::reverse_iterator  it = ma.rbegin();
      read = it->second;
      ma.clear();
   }
   return read;
}

int main(unsigned argc, char* argv[])
{
    if ( argc != 2)
    {
         manual();
         exit(1);
    }
    FILE * queryFile = fopen(argv[1], "r");
    char qline[1024] = "\0";
    char qHeader[1024] = "\0";
    char qSequence[65536] = "\0";
    char qSequence1[65536] = "\0";

    while (!feof(queryFile))
    {
        fgets(qline, 1024, queryFile);
        v = splitv(qline,"\t");
        string s ;
        s += v[0];
        s += "-";
        s += v[1];
        m.insert(std::pair<string,vector<string> > (s,v));   
        ha[s] = s;     
    }
    fclose(queryFile);

   char flname[500];
   strcpy(flname,argv[1]);
   strcat(flname,".qualityreads"); 
   fstream out (flname,ios::out);

   for (  map<string,string>::iterator ii= ha.begin(); ii!= ha.end(); ++ii)
   {
        out << getbiggest(ii->first);
        ha.erase(ii->first);
   }

   ha.clear();
   m.clear();
   return 0;
}
