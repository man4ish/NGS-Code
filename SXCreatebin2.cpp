#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
using namespace std;
struct ltstr
{
  bool operator()(string s1, string s2) const
  {
    int mark1 = s1.find("-");
    int mark2 = s2.find("-");
    return atoi(s1.substr(0,mark1-1).c_str()) < atoi(s2.substr(0,mark2-1).c_str());
  }
};

map <string,vector<int>,ltstr> IndexbinMap;
void createbin(const char* filename )
{
   int index = 0;
   int reccount = 0;
   ifstream infile(filename);
   if(!infile.is_open())
   {
      cout<<"could not open the feature file\n";
      exit(1);
   }
   char line[500],rec[500];
   while(infile)
   {
        infile.getline(line,500);
        strcpy(rec,line);
        if(infile)
        {
           char *pch ;
           pch = strtok(line,"\t ");
           vector<const char* > vec;
           while(pch != NULL )
           {
                vec.push_back(pch);
                 pch = strtok(NULL,"\t ");
           }
              size_t offset = (size_t)atol(vec[0]);
              int i = offset/50000;
              size_t min = i*50000;
              size_t max = i*50000 + 50000;
              std::ostringstream minval, maxval;
              minval << min;
              maxval << max;
              string bin;
              bin += minval.str();
              bin += "\t";
              bin += maxval.str();
              vector <int> vecint;
              if(IndexbinMap.count(bin) == 0)
              {
                 vecint.push_back(index);
                 vecint.push_back(1);
                 IndexbinMap[bin] = vecint;
              }
              else 
              {
                  
                  IndexbinMap[bin][1]++;
              }              
              index += strlen(rec);
              index += 1;
              vec.clear();
        }
    }
    infile.close();
}
int main(int argc, char* argv[])
{
   createbin(argv[1]);
   char flname[1000];
   strcpy(flname, argv[1]);
   strcat(flname,".bin");
   ofstream outfile(flname);
   if(!outfile.is_open())
   {
      cout<< "could not open the output bin file\n";
      exit(1);
   }

   for( map<string, vector<int> >::iterator ii= IndexbinMap.begin(); ii!= IndexbinMap.end(); ++ii)
   {
       outfile << (*ii).first << "\t" << (*ii).second[0]<<"\t"<<(*ii).second[1] << endl;
   }
   return 0;
}
