#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <time.h>
using namespace std;

struct Key
{
    int gi;
    size_t offset;
    char allele[5000];
    Key(int _gi, size_t  _offset, const char* _allele)
    {
      gi = _gi ;
      offset = _offset;
      strcpy(allele, _allele); 
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
          return (strcmp(_lhs->allele,_rhs->allele)< 0);
    }
};

typedef map<Key*,char*,ltstr> SNVMap;
SNVMap m1,m2,m3;

int commcount =0;
int fcount1 =0;
int fcount2 =0;
int reccount2  =0;

clock_t start = clock();
char Title[2000];

int RecordSample1 = 0; 
int RecordSample2 = 0;

bool SNP = false;
bool DEL = false;
bool INS = false;


void manual() 
{
     cout<< "\n\nSXMapperAnalysis\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXMapperAnalysis <SNV Type> <Introduced SNV> <Mapper SNV> \n\n";
}

void SetTitle(const char* title)
{
   strcpy(Title,title);
}

const char* GetTitle()
{
   return Title;
}

void readIntroducedsnvfile(char* line)
{
     if(line[0] != '\0' && line[0]!= '#') 
     {
           if(line[0] == 'C' && line[1] == 'h' && line[2] == 'r')
           SetTitle(line);
           char recline[40000];
           strcpy(recline,line);
           char* pch;
           vector<const char*> IntroSNVvec;
           pch = strtok (recline,"\t");
           while(pch != NULL) 
           {
              IntroSNVvec.push_back(pch);
              pch = strtok (NULL,"\t");
           }
           
           int chrnum;
           
           if(IntroSNVvec[0] == "X")
              chrnum = 23;
           else if(IntroSNVvec[0]== "Y")
              chrnum = 24;
           else
              chrnum = atoi(IntroSNVvec[0]);
           
           size_t offset = (size_t) atol(IntroSNVvec[1]);

           char allele[4];
           if(SNP)           
           {
              allele[0] = IntroSNVvec[3][0];
              allele[1] = '>';
              allele[2] = IntroSNVvec[4][0];
              allele[3] = '\0';                  
           }
           else if(INS)
           {
              allele[0] = IntroSNVvec[4][0];
              allele[1] = '\0';   
           }
           else if(DEL)
           { 
              allele[0] = IntroSNVvec[3][0];
              allele[1] = '\0'; 
           } 
           //cout << chrnum <<"\t"<< offset <<"\t"<< allele <<"\t"<< endl;
           Key* key = new Key(abs(chrnum),offset,allele);
           int len = strlen(line);
           char* inprec = new char[len+1];
           strcpy(inprec,line);
           inprec[len] = '\0';
           m1[key]  = inprec;
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
     common_file << GetTitle() <<"\t"<< GetTitle()<< endl;      
     vector<const char*> vec;

     ifstream snpfile2(filename2);
     uniquefile2 << GetTitle() << endl;
     while(snpfile2)
     {
        char line[40000];
        snpfile2.getline(line,40000);
        if(snpfile2)
        {
           if( line[0] != '\0' && line[0] != '#')
           {
               reccount2++;
               char recline[40000];
               strcpy(recline,line);
               char* pch;
               vector<const char*> MapperSNVvec;
               pch = strtok (recline,"\t");
               while(pch != NULL)
               {
                  MapperSNVvec.push_back(pch);
                  pch = strtok (NULL,"\t");
               }
               int chrnum;
               size_t offset;

               if(strcmp(MapperSNVvec[0],"X") == 0)
                  chrnum = 23;
               else if(strcmp(MapperSNVvec[0],"Y") == 0)
                  chrnum = 24;
               else
                  chrnum = atoi(MapperSNVvec[0]);
            
               size_t Mapper_offset = (size_t) atol(MapperSNVvec[2]);
               Key* key = new Key(chrnum,Mapper_offset,MapperSNVvec[3]);
                
               //cout << chrnum << "\t" << Mapper_offset << "\t" << MapperSNVvec[3].c_str()<<endl;
              
               if(m1.count(key) != 0 )
               {
                  commcount++;
                  int len = strlen(line);
                  char* unirec = new char[len+1];
                  strcpy(unirec,line);
                  unirec[len] = '\0';
                  m3[key]= unirec;
                  common_file << line <<"\t"<< m1[key]<< endl;
               }
               else 
               {
                     fcount2++;
                     uniquefile2 << line <<endl;
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
     char flname[1000] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_unique");
     ofstream uniquefile1(flname);
     uniquefile1 <<GetTitle()<<endl;
     for( SNVMap::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
     {
          SNVMap::const_iterator it = m3.find(ii->first);
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
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,".");
     strcat(flname,filename2);
     strcat(flname,".log");
     ofstream summaryfile(flname);
     summaryfile << filename1 << "_" << filename2 <<"_common(True Positive)\t"<<commcount<< "\t"<< (float)(100*commcount)/(float)m1.size()<<endl;
     summaryfile << filename1<<"_unique(False Negative)\t" << fcount1 << "\t" <<(float)(fcount1*100)/(float)m1.size() <<endl;
     summaryfile << filename2<<"_unique(False Positive)\t" << fcount2 <<  "\t" <<(float)(fcount2*100)/(float)m1.size() <<endl;
     summaryfile << "Elapsed"<<"\t"<<1000*(clock() - start)/CLOCKS_PER_SEC<<" "<<"ms\n";
     summaryfile.close();
}

int main(int argc ,char* argv[]) 
{
     if(argc != 4) 
     {
        manual();
        exit(1);
     }

     if(strcmp(argv[1],"SNP") == 0)
        SNP = true;
     else if(strcmp(argv[1],"INS") == 0)
        INS = true;
     else if(strcmp(argv[1],"DEL") == 0)   
        DEL = true;

    char line1[1000],line2[40000];
    ifstream snpfile1 (argv[2]);
    if(!snpfile1.is_open())
    {
       cout<<"Unable to open the Introduced SNV file\n";
       exit(1);
    }
    while(snpfile1) 
    {
          char line1[40000];
          snpfile1.getline(line1,40000);
          if(snpfile1)
            readIntroducedsnvfile(line1);
    }   
    snpfile1.close();
    cout << "Introduced SNV File Loaded\n"; 

    ifstream snpfile2 (argv[3]);
    if(!snpfile2.is_open())
    {
       cout<<"Unable to open the snpfile2\n";
       exit(1);
    }
    getcommon_snps(argv[2],argv[3]);
    getuniquetofile1(argv[2]);
    get_statistics(argv[2],argv[3]);

    m1.clear();
    m3.clear();
    return 0;
}
