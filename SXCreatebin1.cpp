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
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};

map <string,size_t,ltstr> IndexbinMap;
void createbin(const char* filename )
{
   int index = 0;
   ifstream infile(filename);
   if(!infile.is_open())
   {
      cout<<"could not open the feature file\n";
      exit(1);
   }
   char flname[200];
   strcpy(flname, filename);
   strcat(flname,".bin");
   ofstream outfile(flname);
   if(!outfile.is_open())
   {
      cout<< "could not open the output bin file\n";
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
              if(IndexbinMap.count(bin) == 0)
              {
                 IndexbinMap[bin] = index;
                 outfile << bin<< "\t"<<index <<endl;
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
   return 0;
}
