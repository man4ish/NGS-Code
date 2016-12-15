#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
vector <char> GetBases(const char* filename,size_t start, size_t stop)
{
      vector <char> charvec; 
      if(start > stop)
      {
          cout << "min_range is greater than max_range\n";
          return charvec;
      }
      else if(start == stop)
      {
            cout << "min_range and max_range are same\n";
            return charvec;
      }
      ifstream is;
      char base;
      is.open(filename, ios::binary );
      size_t index = start -1;
      while(1)
      {
          is.seekg(index, ios::beg);
          base = is.get();          
          if(index >= stop)
          {
             return charvec;
          }
          charvec.push_back(base);
          index++;
      }
      is.close();
}

int main (int argc, char* argv[])
{
     vector<char> charvec = GetBases(argv[1], (size_t)atol(argv[2]), (size_t)atol(argv[3]));
     for(int y = 0; y < charvec.size(); y++)
     {
        cout << charvec[y];
     }
     cout << "\n";
     return 0;
}

