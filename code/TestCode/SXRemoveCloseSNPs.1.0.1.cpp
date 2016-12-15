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
   return (size_t)atol(vec[2]);
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
         int diff = currentindex - previousindex ;

        // if(diff <= 5 && diff > 0)
          if(diff == 0) 
         {
           s1.insert(record1);
           s1.insert(record2);
         }
         else 
         {
             s2.insert(record1);
             s2.insert(record2);
         }  
         myqueue.pop();
       }
    }
  }
    char filename[10000];
    strcpy(filename,argv[1]);
    strcat(filename,"_5bpCloseSNPSummary");
    ofstream summary(filename);
    summary <<s1.size()<<endl;
    summary.close(); 
    set<string>::iterator theIterator1;
    for( theIterator1 = s1.begin(); theIterator1 != s1.end(); theIterator1++ ) 
    {
       out1<< *theIterator1<<endl;
    }
 /*   set<string>::iterator theIterator2;
     for( theIterator2 = s2.begin(); theIterator2 != s2.end(); theIterator2++ ) {
       out2 << *theIterator2<<endl;
     }
*/
  out1.close();
  out2.close();
  return 0;
}
