#include <iostream>
#include <fstream>
using namespace std;

char* GetBases(const char* filename,size_t start, size_t stop)
{     
      if(start > stop)
      {
            cout << "min_range is greater than max_range\n";
            return 0;
      }
      else if(start == stop)
      {
            cout << "min_range and max_range are same\n";
            return 0;
      }

      ifstream seqfile;
      size_t filesize = 0;
      seqfile.open(filename, ios::binary );
      {
          if(seqfile.is_open())
          {
             seqfile.seekg(0,ios::end);
             filesize = seqfile.tellg();
          }  
      }
      if(filesize < stop )
      {
          cout << "Out of Range\n";
          return 0;
      }

      char base;
      int len = stop - start;
      char* extseq = new char[len+1];
      size_t irdx = 0;
      size_t index = start -1;

      while(index < stop)
      {
          seqfile.seekg(index, ios::beg);
          base = seqfile.get();          
          extseq[irdx++] = base;
          index++;
      }
      seqfile.close();

      extseq[len] = '\0';  
      return extseq;
}

int main (int argc, char* argv[])
{
      cout << GetBases(argv[1], (size_t)atol(argv[2]), (size_t)atol(argv[3]))<< endl;
      return 0;
}

