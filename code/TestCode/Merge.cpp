#include <iostream>
#include <cstring>
#include <queue>
#include <cctype>
#include <fstream>
#include <stdlib.h>
using namespace std;

vector<const char*> Getoffset(string line)
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
    return vec;
}

int main(int argc, char* argv[])
{    
  char line[10000];
  queue<string> myqueue;
  ifstream in(argv[1]);
  while(in)
  {
    in.getline(line,10000);
    if(in)
    {
         if(myqueue.size() > 0 )
         {
             string record1 = myqueue.back();
             int previousindex = atol(Getoffset(record1)[1]);
             int currentindex = atol(Getoffset(line)[1]);
             int diff = currentindex - previousindex;
             if(diff == 1)
             {
                myqueue.push(line);
             }
             else
             {
                string allele;
                size_t offset = atol(Getoffset(myqueue.front())[1]);
                while(!myqueue.empty())
                {
                      allele += Getoffset(myqueue.front())[4];
                      myqueue.pop();
                }
                cout << offset << "\t"<< allele << endl;
                myqueue.push(line); 
             }
        }
        else 
        {
           myqueue.push(line);
        }
   }
 }  
 return 0;
}
