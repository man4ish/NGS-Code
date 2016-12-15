#include <iostream>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <set>
#include <time.h>
#include <string.h>
using namespace std;
map<string,string> qual;
vector<string> reads; 

void manual()
{
    fprintf(stderr, "\n\nSXQualSort\n");
    fprintf(stderr, "~~~~~~~~~~~~~~~~~~\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\tSXQualSort  <read file> <quality file>\n\n");
}

char* Getfile(char* filename)
 {
    char* pfile;
    int leninputfile = strlen(filename);
    int x, fi = 0;
     for(;fi < leninputfile;fi++)
         if(*(filename + fi)=='.')
             x= fi;
    pfile = new char[x+1];
     memcpy(pfile,filename,x);
     pfile[x]='\0';
     return pfile;
 }

void get_sortedqual(char* filename) {
char flname[200];
strcpy(flname,filename);
strcat(flname,"_sorted.qual");
ofstream out(flname);
for( int ii=0;ii < reads.size(); ii++)   {
      out << reads[ii]<<endl;
       string h = reads[ii];
       out<<qual[h];
   }
}


int main(int argc, char*argv[]) {       
         if(argc== 1){ 
         manual();
         exit(9);
       }
               clock_t start= clock(); 
               char qline[1024] = "\0";
               FILE * queryFile = fopen(argv[1], "r");
               while (!feof(queryFile)) {
                  fgets(qline, 1024, queryFile);
                  if(qline[0] == '>')  {        
                         qline[strlen(qline) - 1] = '\0';
                         reads.push_back(qline);
                     
                        }
                 }
               fclose(queryFile); 

               char qline2[1024] = "\0";
               char qHeader2[1024] = "\0";
               char qSequence2[655366] = "\0";
               FILE * queryFile2 = fopen(argv[2], "r");
               while (!feof(queryFile2))
                 {
                  fgets(qline2, 1024, queryFile2);
                  if((qline2[0] == '>') || (feof(queryFile2)))  {
                     if (qSequence2[0] == '\0')  {
                         qline2[strlen(qline2) - 1] = '\0';
                         strcpy(qHeader2, qline2);
                       }
                     else
                       {
                         qual[qHeader2] = qSequence2;
                         qSequence2[0] = '\0';
                         qline2[strlen(qline2) - 1] = '\0';
                         strcpy(qHeader2, qline2);
                       }
                   }
                 else {
                   strcat(qSequence2, qline2);
                  }
                }
    fclose(queryFile2);
    char * str = Getfile(argv[2]);
    get_sortedqual(str);
    cout<< "Elaspsed: "<<(clock() - start)/CLOCKS_PER_SEC<<" seconds\n";
    return 0;
}

