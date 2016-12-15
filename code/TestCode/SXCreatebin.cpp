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
   float index = 0;
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
        if(infile && line[0] != '\0')
        {
           char *pch ;
           pch = strtok(line,"\t ");
           vector<const char* > vec;
           while(pch != NULL )
           {
                vec.push_back(pch);
                 pch = strtok(NULL,"\t ");
           }
           size_t offset = (size_t)(atof(vec[14])*1000);
           if(offset > 6000)
           cout << rec << endl;
           int i = offset/250;
           int min = i*250;
           int max = i*250 + 250;
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
   outfile << "Min Total-Score\tMax Total-Score\tCount\n";
   for( map<string, vector<int> >::iterator ii= IndexbinMap.begin(); ii!= IndexbinMap.end(); ++ii)
   {
       outfile << (*ii).first  <<"\t"<<(*ii).second[1] << endl;
   }
   outfile.close();
   return 0;
}

