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

bool IstwntyetAs(const char* seq)
{
   int len = strlen(seq);
   if(len < 29)
   return 0;
   for(int ii =len-1; ii >= len -29 ;ii--)
   {
          if(seq[ii] !='A')
             return 0;        
   }
   return 1;
}

bool IstwelveAs(const char* seq)
{
   int len = strlen(seq);
   if(len < 13)
   return 0;
   for(int ii =len-1; ii >= len -13 ;ii--)
   {
          if(seq[ii] !='A')
             return 0;
   }
   return 1; 
}

bool IsfiveAs(const char* seq)
{
   int len = strlen(seq);
   if(len < 5)
   return 0;
   for(int ii =len-1; ii >= len -5 ;ii--)
   {
          if(seq[ii] !='A')
             return 0;
   }
   return 1;       
}

bool IsAdapter(const char* seq)
{
     string sequi = seq;
     string ada = ADAP_seq;
     reverse(sequi.begin(), sequi.end());
     reverse(ada.begin(), ada.end());
     for( int i = 0; i< ada.length(); i++)
     {
          if(ada.at(i) != sequi.at(i))
             return 0;
     }
     return 1;
}

bool IsPoly(const char* seq)
{
     if((strstr(seq,PLOYA1)!= NULL)||(strstr(seq,PLOYA2)!=NULL)||(strstr(seq,PLOYA3)!=NULL)) 
     {
        return 1;
     }
     return 0;
}

const char* Getsequence(const char* seq)
{
      int i =0;
      int seqend;
      for(;i<strlen(seq);i++)
      {
          if(*(seq+i)=='\n' && *(seq+i+1)== '+')
          {
             seqend = i;
          }
      }
      char* sequence = new char[seqend+1];
      memcpy(sequence, seq,seqend);
      sequence[seqend]='\0'; 
      return sequence;     
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
         strcat(flname,".mask");
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
 
                       if (ThreePrimeAscount(qSequence) > 12 )
                        {
                              int len = strlen(qSequence);
                              int range = len - 12;
                              char seqw[range+1]; 
                              strncpy(seqw, qSequence,range);
                              seqw[range] = '\0';  
                                 out<<qHeader<< seqw <<endl;
                        }
                        else if (  ThreePrimeAscount(qSequence) < 12 && ThreePrimeAscount(qSequence) > 5)
                        {
                            
                            
                             char ss[100];
                             memcpy(ss,qSequence,strlen(qSequence)-6);
                                   out<<qHeader<<ss<<endl;
                             
                        }
                        else 
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
