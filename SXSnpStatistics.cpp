#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
#include <string>


using namespace std;
typedef map<string, string > m;
m m1,m2,m3;
int commcount =0;
int fcount1 =0;
int fcount2 =0;
clock_t start = clock();

void manual() 
{
     cout<< "\n\nSXSnpstatistics\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXSnpstatistics <Normal Sample> <Cancerous Sample> \n\n";
     exit(1);
}

void readsnpfile1(string line1)
{
     if(line1 != "") 
     {
        char* f = &line1[0],*pch;
        vector<string> sdf;
        string s;
        pch = strtok (f," ");
        while(pch != NULL) 
        {
              sdf.push_back(pch);
              pch = strtok (NULL," ");
        }
      s += sdf[0]; s += "-"; s += sdf[1]; s += "-"; s += sdf[2];s += "-"; s += sdf[4];
      m1[s] =line1;
     }
}

void readsnpfile2(string line2)
{
     if(line2 != "") 
     {
        char* f = &line2[0],*pch;
        vector<string> sdf;
        string s;
        pch = strtok (f," ");
        while(pch != NULL) 
        {
              sdf.push_back(pch);
              pch = strtok (NULL," ");
        }
        s += sdf[0]; s += "-"; s += sdf[1]; s += "-"; s += sdf[2];s += "-"; s += sdf[4];
        m2[s] =line2;
     }
}

void getcommon_snps(char* filename1,char* filename2) 
{
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_");
     strcat(flname,filename2);
     strcat(flname,"_common");
     ofstream fh(flname);
     string h,h1,h2;
     vector<string> vec;
	for(m::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
        {
	    m::const_iterator it = m2.find(ii->first);
            if(it != m2.end()) 
            {
	       commcount++;
	       m3[ii->first]= ii->second;
               fh << ii->second <<endl;
            }
         }
}

void getuniquetofile1(char* filename1)
{
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_unique");
     ofstream fh(flname);
     for( m::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
     {
          m::const_iterator it = m3.find(ii->first);
          if(it == m3.end()) 
          {
             fcount1++;
             fh << ii->second<<endl;
          }
     }
}

void getuniquetofile2(char* filename2)
{ 
     char flname[200] = "\0";
     strcpy(flname,filename2);
     strcat(flname,"_unique");
     ofstream fh(flname);
     for( m::iterator ii=m2.begin(); ii!=m2.end(); ++ii) 
     {
          m::const_iterator it = m3.find(ii->first);
	  if(it == m3.end()) 
          {
             fcount2++;
             fh << ii->second <<endl;
          }
     }
}

void get_statistics(char* filename1,char* filename2) 
{
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,filename2);
     strcat(flname,".log");
     ofstream out(flname);
     out<<filename1<<"_"<<filename2<<"_common"<< ":"<<commcount<<endl;
     out <<filename1<<"_only"<< ":"<<fcount1<<endl;
     out <<filename2<<"_only"<<":"<<fcount2<<endl;
     out <<"Elapsed"<<":"<<1000*(clock() - start)/CLOCKS_PER_SEC<<" "<<"ms\n";
}

int main(int argc ,char* argv[]) 
{
  if(argc < 2||argc > 3) 
  {
    manual();
    exit(1);
  }
  string line1,line2;

   ifstream snpfile1 (argv[1]);
   if(!snpfile1.is_open())
   {
      cout<<"Unable to open the snpfile1\n";
      exit(1);
   }
   while(!snpfile1.eof()) 
   {
          getline (snpfile1,line1);
          readsnpfile1(line1);
   }

   ifstream snpfile2 (argv[2]);
   if(!snpfile2.is_open())
   {
      cout<<"Unable to open the snpfile2\n";
      exit(1);
   }
   while(!snpfile2.eof()) 
   {
          getline (snpfile2,line2);
          readsnpfile2(line2);
   }
   getcommon_snps(argv[1],argv[2]);
   getuniquetofile1(argv[1]);
   getuniquetofile2(argv[2]);
   get_statistics(argv[1],argv[2]);
   return 0;
}
