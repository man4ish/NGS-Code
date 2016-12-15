#include <iostream>
#include <map>
#include <fstream>

using namespace std;

struct ltstr
{
    bool operator()(string _lhs, string _rhs) const
    {
        return strcmp(_lhs.c_str(), _rhs.c_str()) < 0;
    }
};


void manual()
{
    cout<< "\n\nSXRemoveShortReads\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXRemoveShortReads <query file>\n\n";
}

int main(int argc, char* argv[])
{
        if(argc != 2)
        {
           manual();
           exit(1);
        }

               char qline[1024] = "\0";
               char qHeader[1024] = "\0";
               char qSequence[65536] = "\0";
               char flname[400] = "\0";
               strcpy(flname,argv[1]);
               strcat(flname,".filtdShort" );
               ofstream out (flname);     
               FILE * queryFile = fopen(argv[1], "r");
               while (!feof(queryFile))
               {
                      fgets(qline, 1024, queryFile);
                      if((qline[0] == '>') || (feof(queryFile)))
                      {        
                          if (qSequence[0] == '\0')
                          {
                              strcpy(qHeader, qline);
                          }
                          else
                          {
                              if(qSequence[strlen(qSequence)-1]== '\n')
                                 qSequence[strlen(qSequence)-1]= '\0'; 
                              if ( strlen(qSequence) >= 18)
                              {
                                     out<<qHeader<<qSequence<<endl;  
                              }
                              qSequence[0] = '\0';
                              strcpy(qHeader, qline);
                          }
                     }
                     else
                     {
                         strcat(qSequence, qline);
                     }
              }
       fclose(queryFile);
return 0; 
}

