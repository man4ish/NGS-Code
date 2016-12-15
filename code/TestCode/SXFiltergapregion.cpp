#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;
bool FilterGapsnps(int i, size_t offset)
{
    string filename = "../../syhwah/Manish/gapdata/chr";
    std::stringstream out;
    out << i;       
    string s = filename;
    s  += out.str();
    s  += ".txt";
    size_t start,stop; 
    ifstream in(s.c_str());
    char* pch;
    char line[1000];
    pch = strtok(line,"\t");
    while(in)
    {
      in.getline(line,1000);
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
         size_t start = (size_t) atol(vec[0]);
         size_t stop = (size_t) atol(vec[1]);

         if(offset >= start && offset <= stop)
         {
           in.close();
           return 1;
         }
      } 
    }
    in.close();
    return 0;
} 


int main(int argc, char* argv[])
{
   ifstream inn(argv[1]);
   char lline[5000],recline[5000];
   while(inn)
    {
      inn.getline(lline,5000);
      if(inn)
      {
         strcpy(recline,lline);
         char* pchh;
         vector <const char*> vecc;
         pchh = strtok(lline,"\t");
         while(pchh != NULL)
         {
            vecc.push_back(pchh);
            pchh = strtok(NULL,"\t");
         } 
         if(FilterGapsnps(atoi(vecc[0]), (size_t) atol(vecc[1])) == 0)
         {
           cout << recline <<endl;
         }
      }
    }
    inn.close();
  return 0;
}

