#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <set>
#include <fstream>
#include <vector>
using namespace std;
vector <const char*> records;
void manual()
{
    cout<< "\n\nSXRandomRecordsGenerator\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXRandomRecordsGenerator <query file> <number of records> \n\n";
    exit(1);
}

int main (int argc, char* argv[])
{
  if(argc != 3)
  { 
     manual();
  }
  ifstream in(argv[1]);
  char line[10000];
  while(in)
  {
     in.getline(line,10000); 
     if(in)
     {
        records.push_back(line);  
     }     
  }
  in.close();
  int randnum;
  srand ( time(NULL) );
  set <int> s;
  for (int i =0; i < (atoi(argv[2])+100); i++)
  {
      randnum = rand() % records.size()+ 1;    
      s.insert(randnum);
  }
  set<int>::iterator theIterator;
  vector <size_t> vec;
  ofstream out("randomreads");
  for( theIterator = s.begin(); theIterator != s.end(); theIterator++ ) 
  {
       vec.push_back(*theIterator);
       if(vec.size() == atoi(argv[2]))
       {         
          for (int j=0; j < vec.size(); j++)
          out << records[vec[j]]<< endl;
          out.close();
          s.clear();
          vec.clear();
          records.clear();
          exit(1);
       }
  }
  return 0;
}

