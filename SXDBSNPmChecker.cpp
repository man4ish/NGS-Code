#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>

using namespace std;
typedef map<string,string> M; 
M m,n;
clock_t start = clock();
int novcount =0;
int knowncount =0;
int snpcount=0;

void manual() {
  cout<<"\n\nSXDBSNPmChecker\n";
  cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  cout<<"Usage:\n\n";
  cout<<"\tSXDBSNPmChecker <Snpfile> <SnpRnaDB>\n\n";
 }

void GetStatistics(char* filename) {
  char flname[200] ="\0";
  strcpy(flname,filename);
  strcat(flname,".log");
  ofstream out(flname);
  out<<"Total number of Snps : "<<snpcount<<endl;
  out<<"Novel Snps : "<<novcount<<"	"<<float (novcount*100)/float(snpcount)<<"%"<<endl;
  out<<"Known Snps : "<<knowncount<<"	"<<float (knowncount*100)/float (snpcount)<<"%"<<endl;
  out<<"Elapsed:     "<<(clock() - start)*1000/CLOCKS_PER_SEC<<" "<<"ms\n";
}

void GetNovelSnps(char* filename,string lo ) {
 if(lo != "") {
     char* ff = &lo[0],*pchh;
     vector<string> ssdf;
     string ss,snprec;
     snpcount++;
     snprec = lo;
     pchh = strtok (ff," ");
     while (pchh != NULL) {
       ssdf.push_back(pchh);
       pchh = strtok (NULL, " ");
     }
  ss += ssdf[0]; ss += "-"; ss += ssdf[1]; ss += "-"; ss += ssdf[2]; ss += "-";ss += ssdf[4];
  char flname[500]="\0";
  strcpy(flname,filename);
  strcat(flname,".mapresult");
  ofstream outc(flname);
      if(n.count(ss)==0) {
           novcount++;
           outc <<snprec<<"No"<<endl; 
      }
      else {
            knowncount++; 
            outc <<snprec<<"Yes"<<endl; 
      }
   }

}

void LoadSnpDB(string lo) {
  if(lo != "")  {
    char* f = &lo[0],*pch;
    vector<string> sdf;
    string s;
    pch = strtok (f,"\t");
    while (pch != NULL) {
	sdf.push_back(pch);
        pch = strtok (NULL,"\t");
    }
    s += sdf[0]; s += "-"; s += sdf[1]; s += "-"; s += sdf[2]; s += "-";s += sdf[3];
    n[s] =s;
  }
}

int main(int argc, char* argv[]) 
{
   if(argc != 3) {
   manual();
   exit(1);
   }
   string line1,line2;
   ifstream rnasnpdb (argv[2]);
   if (rnasnpdb.is_open()) {
   while (! rnasnpdb.eof() ) {
     getline (rnasnpdb,line2);
     LoadSnpDB(line2);
    }
    rnasnpdb.close();
  }
 else cout<<"Unable to load the database\n";
   ifstream snpfile (argv[1]);
   if (snpfile.is_open()) {
    while (! snpfile.eof() ) {
      getline (snpfile,line1);
      GetNovelSnps(argv[1],line1);
     }
    snpfile.close();
  }
  else cout << "Unable to load snp file\n";
  GetStatistics(argv[1]);
  return 0;
}
