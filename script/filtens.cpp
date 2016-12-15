#include <iostream>
#include <fstream>
using namespace std;
const char*  PLOYA1 = "AATAAA";
const char*  PLOYA2 = "ATTAAA";
const char*  PLOYA3 = "AGTAAA";
const char* ADAP_seq = "TCGTATGCCGTCTTCTGCTTGT";

void manual()
{
    cout<< "\n\nSXFilterPolyAs\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXFilterPolyAs <query file>\n\n";
}

bool IsallAs(const char* seq)
{
    while(*seq)
    {
        if(*seq!='A')
          return 0;  
        *seq++;
    }
    return 1;
}

bool IstenAs(const char* seq )
{
   for(int i =0; i<= 10;i++)
   {
       if(seq[i] != 'A')
          return 0;  
   }
   return 1;
}


int ThreePrimeAscount(char* seq)
{
    int count = 0;
    int len = strlen(seq);
    for(int ii =len-1; ii >= len -29 ;ii--)
    {
          if(seq[ii] !='A')
             return count;
            count++;
    }
   return count; 
} 


int main(int argc ,char* argv[])
{
         if(argc != 2)
         {
           manual();
           exit(1);
         }
         char qline[1024] = "\0";
         char qHeader[1024] = "\0";
         char qSequence[65536] = "\0";
         FILE * queryFile = fopen(argv[1], "r");
         char flname[200]; 
         strcpy(flname,argv[1]);
         strcat(flname,".filtTens");
         ofstream out(flname);
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
                       if(qSequence[strlen(qSequence)-1]=='\n')
                           qSequence[strlen(qSequence)-1]='\0';
 
                        if (IstenAs(qSequence)==0)
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
