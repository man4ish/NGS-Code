#include <iostream>
#include <queue>
#include <fstream>
#include <set>
using namespace std;
int Getoffset(string line)
{
    char record[10000];
    strcpy(record,line.c_str());
    char* pch;
    pch = strtok(record,"\t");
    vector<const char*> vec;
    while(pch != NULL)
    {
       vec.push_back(pch);
       pch = strtok(NULL,"\t");  
    }  
   return (size_t)atol(vec[1]);
}

int main (int argc, char* argv[])
{
  ofstream out1("data_exceptionfile.txt");
  ofstream out2("Normalfile.txt");
  char line[10000];
  queue<string> myqueue;
  set<string> s1;
  set<string> s2;
  ifstream in(argv[1]);
  while(in)
  {
    in.getline(line,10000);
    if(in)
    {
       myqueue.push(line);
       if(myqueue.size() == 2)
       {
         string record1 = myqueue.front();
         string record2 = myqueue.back();
         int previousindex = Getoffset(record1);
         int currentindex = Getoffset(record2);
         if(currentindex - previousindex <= 5)
         {
           s1.insert(record1);
           s1.insert(record2);
          // out1<< record1 << endl;
          // out1 << record2<<endl;
          // cout << previousindex << "\t"<< currentindex << endl;
         }
         else 
         {
             s2.insert(record1);
             s2.insert(record2);
           //  out2<< record1 << endl;
           //  out2 << record2<<endl;
         }  
         myqueue.pop();
      //   cout << myqueue.front()<<endl;
       //  cout << myqueue.size()<<endl;
       }
    }
  }

     set<string>::iterator theIterator1;
     for( theIterator1 = s1.begin(); theIterator1 != s1.end(); theIterator1++ ) {
       out1<< *theIterator1<<endl;
     }
    set<string>::iterator theIterator2;
     for( theIterator2 = s2.begin(); theIterator2 != s2.end(); theIterator2++ ) {
       out2 << *theIterator2<<endl;
     }

  out1.close();
  out2.close();
  return 0;
}
