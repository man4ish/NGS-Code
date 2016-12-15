#include <iostream>
using namespace std;
std::string GetMapIndexFile(std::string filename)
{
   int j =0;
   for (int i = filename.size()-1; i >= 0 ; i--  )
   {  if(filename[i] == '.')
       j++;
     if(j == 3)
       return filename.substr(0,i) + ".MapIndex.txt";
   }
}
int main()
{
  string s = "hs_mgb_1.0.ref_genome.37.1.gene.chr1.txt";
  std::cout << GetMapIndexFile(s)<< endl;
  return 0;
}
