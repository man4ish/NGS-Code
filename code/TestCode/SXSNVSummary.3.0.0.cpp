#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
#include <string>
using namespace std;

struct Key
{
    int gi;
    size_t offset;
    string allele;
    Key(int _gi, size_t  _offset, string _allele)
    {
      gi = _gi ;
      offset = _offset;
      allele = _allele; 
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi)
         return ( _lhs->gi < _rhs->gi);
        else if(_lhs->offset != _rhs->offset)
           return  ( _lhs->offset < _rhs->offset );
          return (strcmp((_lhs->allele).c_str(),(_rhs->allele).c_str())< 0);
    }
};

typedef map<Key*,string,ltstr> m;
m m1,m2,m3;
int commcount =0;
int fcount1 =0;
int fcount2 =0;
int reccount2  =0;
clock_t start = clock();
char Title[200];
int RecordSample1 = 0; 
int RecordSample2 = 0;
void manual() 
{
     cout<< "\n\nSXSNVSummary\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXSNVSummary <filter limit> <Normal Sample> <filter limit> <Cancerous Sample> \n\n";
}

void SetTitle(const char* title)
{
   strcpy(Title,title);
}

const char* GetTitle()
{
   return Title;
}

void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && line1[0]!= '#') 
     {
           if(line1[0] == 'C' && line1[1] == 'h' && line1[2] == 'r')
           SetTitle(line1);
           char recline1[40000];
           strcpy(recline1,line1);
           char* pch;
           vector<string> sdf;
           string s;
           pch = strtok (recline1,"\t");
           while(pch != NULL) 
           {
              sdf.push_back(pch);
              pch = strtok (NULL,"\t");
           }
           int chrnum;
           size_t offset;
            if(sdf[0] == "X")
           chrnum = 23;
            else if(sdf[0]== "Y")
           chrnum = 24;
           else
           chrnum = atoi(sdf[0].c_str());
           Key* key= new Key(abs(chrnum),atoi(sdf[2].c_str()),sdf[3].c_str());
           m1[key] = line1;
     }
}

void getcommon_snps(char* filename1,char* filename2) 
{
     char Uniquefile2[1000];
     strcpy(Uniquefile2,filename2);
     strcat(Uniquefile2,"_unique");
     ofstream uniquefile2(Uniquefile2);
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_");
     strcat(flname,filename2);
     strcat(flname,"_common");
     ofstream common_file(flname);
     string h,h1,h2;
     vector<string> vec;
     ifstream snpfile2(filename2);
     uniquefile2 << GetTitle() << endl;
     while(snpfile2)
     {
       char line2[40000];
       snpfile2.getline(line2,40000);
       if(snpfile2)
       {
        if( line2[0] != '\0' && line2[0] != '#')
        {
              reccount2++;
              char recline2[40000];
              strcpy(recline2,line2);
              char* pch;
              vector<string> sdf;
              string s;
              pch = strtok (recline2,"\t");
              while(pch != NULL)
              {
                  sdf.push_back(pch);
                  pch = strtok (NULL,"\t");
              }
              int chrnum;
              size_t offset;
              if(sdf[0] == "X")
              chrnum = 23;
              else if(sdf[0]== "Y")
              chrnum = 24;
              else
              chrnum = atoi(sdf[0].c_str());
              Key* key= new Key(abs(chrnum),atoi(sdf[2].c_str()),sdf[3].c_str());
              if(m1.count(key) != 0 )
              {
                 commcount++;
                 m3[key]= line2;
                 common_file << line2 <<"\t"<< m1[key]<< endl;
              }
              else 
              {
                     fcount2++;
                     uniquefile2 << line2 <<endl;
              } 
              RecordSample2++;
         }
      }
    }
     uniquefile2.close();
     common_file.close();
}

void getuniquetofile1(char* filename1)
{
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_unique");
     ofstream uniquefile1(flname);
     uniquefile1 <<GetTitle()<<endl;
     for( m::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
     {
          m::const_iterator it = m3.find(ii->first);
          if(it == m3.end()) 
          {
             fcount1++;
             uniquefile1 << ii->second<<endl;
          }
     }
     uniquefile1.close();
}

void get_statistics(char* filename1,char* filename2) 
{
     if(RecordSample1 != (commcount + fcount1))
       cout << "Duplicate Record Found in "<< filename1 << endl;
     else  if(RecordSample2 != (commcount + fcount2))
       cout << "Duplicate Record Found in "<< filename2 << endl;
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,".");
     strcat(flname,filename2);
     strcat(flname,".log");
     ofstream summaryfile(flname);
     summaryfile<<filename1<<"_"<<filename2<<"_common"<<"\t"<<commcount <<"\t"<<commcount<<endl;
     summaryfile<<filename1<<"_unique\t"<<fcount1<<"\t"<<(float)(commcount*100)/(float)RecordSample1<<endl;
     summaryfile<<filename2<<"_unique\t"<<fcount2<<"\t"<<(float)(commcount*100)/(float)RecordSample2<<endl;
     summaryfile<<"Elapsed"<<"\t"<<1000*(clock() - start)/CLOCKS_PER_SEC<<" "<<"ms\n";
     summaryfile.close();
}

int main(int argc ,char* argv[]) 
{
  if(argc != 3) 
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
   while(snpfile1) 
   {
          char line1[40000];
          snpfile1.getline(line1,40000);
          if(snpfile1)
            readsnpfile1(line1);
   }   
   snpfile1.close();
   cout << "snpfile1 loaded\n"; 

   ifstream snpfile2 (argv[2]);
   if(!snpfile2.is_open())
   {
      cout<<"Unable to open the snpfile2\n";
      exit(1);
   }
   getcommon_snps(argv[1],argv[2]);
   getuniquetofile1(argv[1]);
   get_statistics(argv[1],argv[2]);
   m1.clear();
   m3.clear();
   return 0;
}
