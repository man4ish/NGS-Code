#include<iostream>
#include <fstream>
#include <vector>

using namespace std;
size_t CreateBin(size_t a ,size_t b,size_t c)
{
  if(b <= c && b >= a)
  {
     return b;
  }
  return 0;
}
int main(int argc, char* argv[])
{
   ifstream infile(argv[1]);
   if(!infile.is_open())
   {
      // print message
      exit(1);
   }
   char flname[200];
   strcpy(flname, argv[1]);
   strcat(flname,".bin");
   ofstream outfile(flname);
   if(!outfile.is_open())
   {
      // print message
      exit(1);
   }
 
   size_t min = 1;
   size_t max = 500000	;
   char line[500],rec[500];
   while(infile)
   {
        infile.getline(line,500);
        strcpy(rec,line); 
        if(infile)
        {
           char *pch ;
           pch = strtok(line,"\t");
           vector<const char* > vec;
           while(pch != NULL )
           {
                vec.push_back(pch);
                 pch = strtok(NULL,"\t");
           }
           if(CreateBin(min,(size_t)atol(vec[1]),max)!= 0)
           {
              outfile<< min << "\t" <<max<< "\t"<<vec[0]<<"\n";
              min += 500000;
              max += 500000; 
            }         
 
           vec.clear();
        }
   }
   infile.close();
   return 0;
}
