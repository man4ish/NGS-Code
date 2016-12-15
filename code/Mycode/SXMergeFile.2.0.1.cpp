#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <string.h>
#include <time.h>
#include <sstream>
#include <string>
using namespace std;

struct Key
{
    int gi;
    size_t offset;
    string allele;
    Key(int _gi, size_t  _offset,string _allele):gi(_gi),offset(_offset), allele(_allele) {}
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

unsigned int GetReadDensity(int i, unsigned int uioffset)
{
    string filename;
    filename = "raw_read/Chr";
    std::stringstream out;
    out << i;
    string s = filename;
    s  += out.str();
    s += "_read_den";
    FILE* pf = fopen(s.c_str(),"rb");
    if (!pf) return 0;
    unsigned int g_unReads;
    unsigned char g_acBuf[2];
    fseek(pf,uioffset*2,SEEK_SET);
    fread(g_acBuf,2,1,pf);
    g_unReads = g_acBuf[1];g_unReads<<=8;
    g_unReads|= g_acBuf[0];
    fclose(pf);
    return g_unReads;
}

vector <string> split(const char* line)
{
   char record[5000];
   vector <string> vec;
   strcpy(record,line);
   char* pch;
   pch = strtok(record,"\t");
   while(pch != NULL)
   {
      vec.push_back(pch);
      pch = strtok(NULL,"\t");
   }
   return vec;
}
typedef map<int,float> rdmap;
rdmap nrdmap,trdmap;
map <int, int > mread;
typedef map<Key*,string,ltstr> m;
m m1,m2,m3;
int commcount =0;
int fcount1 =0;
int fcount2 =0;
int reccount2  =0;
clock_t start = clock();
char Title[200];

void manual() 
{
     cout<< "\n\nSXMergeFile\n";
     cout<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
     cout<< "Usage:\n\n";
     cout<< "\tSXMergeFile <Normal_mean_ReadDensity File><Normal Sample><Tumor_mean_ReadDensity File> <Cancerous Sample> \n\n";
}

void SetTitle(const char* title)
{
   strcpy(Title,title);
}

const char* GetTitle()
{
   return Title;
}

void readsnpfile1(char* line1, int suppreads)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chr"))
           SetTitle(line1);
        else
        {
           char recline1[10000];
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
           offset = atoi(sdf[1].c_str());
           Key* key= new Key(abs(chrnum),atoi(sdf[1].c_str()),sdf[2].c_str());
           m1[key] = line1;
         }
     }
}

void getcommon_snps(char* filename1,char* filename2) 
{
     char Uniquefile2[1000];
     strcpy(Uniquefile2,filename2);
     strcat(Uniquefile2,"_unique");
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,"_");
     strcat(flname,filename2);
     strcat(flname,"_common");
     string h,h1,h2;
     vector<string> vec;
     ofstream common("file1.file2");
     ofstream uniquefile1("file2_unique");
     uniquefile1 <<"Chr\tOffset\tAllele\tCG1.SNP_Density\tCG2.SNP_Density\tTotal SNP Desnity\tCG1.Read_Density\tCG2.Read_Density\tTotal_Read_Desity\tCG1.AVQS\tCG2.AVQS\tAVQS\tdbSNP\tHomo-Het\tRS\n";
     common <<"Chr\tOffset\tAllele\tCG1.SNP_Density\tCG2.SNP_Density\tTotal SNP Desnity\tCG1.Read_Density\tCG2.Read_Density\tTotal_Read_Desity\tCG1.AVQS\tCG2.AVQS\tAVQS\tdbSNP\tHomo-Het\tRS\n";
     ifstream snpfile2(filename2);
     while(snpfile2)
     {
       char line2[1000];
       snpfile2.getline(line2,1000);
       if(snpfile2)
       {
        if( line2[0] != '\0' && !strstr(line2,"#"))
        {
           reccount2++;
           if(strstr(line2,"Chromosome"))
           SetTitle(line2);
           else
           {
              char recline2[1000];
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
              offset = atoi(sdf[1].c_str());
              Key* key= new Key(abs(chrnum),atoi(sdf[1].c_str()),sdf[2].c_str());
              if(m1.count(key) !=0)
              {
                 m3[key] = line2;
                 vector <string> records1 = split(line2);
                 string s = m1[key].c_str();
                 vector <string> records2 = split(s.c_str());
                 int TotalSd = 0;
                 int sr1 = atoi(records1[3].c_str());
                 int sr2 = atoi(records2[3].c_str());
                 TotalSd += sr1;
                 TotalSd += sr2;
                 int TotalRd = 0;
                 TotalRd += atoi(records1[4].c_str());
                 TotalRd += atoi(records2[4].c_str());
                 float avqs1 = atof(records1[5].c_str());
                 float avqs2 = atof(records2[5].c_str());
                 float avqs = (float)(avqs1*sr1+ avqs2*sr2)/(float)(sr1+sr2);
                 int rs = (int)(( (float)TotalSd*100 / (float)TotalRd ) + 0.5);
                 float hethomcount =  (float)TotalSd/(float)(TotalRd-TotalSd);
                 string hethom;
                 if(hethomcount>= 2.333) 
                 hethom = "homo";
                 else 
                 hethom = "het";
             /*    int tot_rd = 0;
                 int chrnum;
                 if(records1[0]== "X")
                 chrnum = 23;
                 else if(records1[0]== "Y")
                 chrnum = 24;
                 else 
                 chrnum = atoi(records1[0].c_str());  
                 size_t offset = (size_t)atol(records1[1].c_str());
                 for (int index = 0; index <= 200;index++)
                 tot_rd += GetReadDensity(chrnum,offset-100+index);
                 float avg_rd = ((float)tot_rd/(float)201);
                 float cutoff = (float) avg_rd/(float)nrdmap[chrnum]; */
                // cout << "cutoff\t"<< cutoff << endl;
common << records1[0]<<"\t"<<records1[1]<<"\t"<<records1[2]<<"\t"<<records1[3]<<"\t"<<records2[3]<<"\t"<<TotalSd<<"\t"<<records1[4]<<"\t"<< records2[4]<<"\t"<<TotalRd<<"\t"<<avqs1<<"\t"<<avqs2<<"\t"<<avqs<<"\t"<<records1[6]<<"\t"<<hethom<<"\t"<<rs<< endl;
                  commcount++;
              }
              else 
              {
                 vector <string > records1 = split(line2);
                 float hethomcount =  (float)atof(records1[3].c_str())/(float)(atof(records1[4].c_str())-atof(records1[3].c_str()));
                 string hethom;
                 if(hethomcount>= 2.333)
                 hethom = "homo";
                 else
                 hethom = "het";
               /*
                 int tot_rd = 0;
                 int chrnum;
                 if(records1[0]== "X")
                 chrnum = 23;
                 else if(records1[0]== "Y")
                 chrnum = 24;
                 else
                 chrnum = atoi(records1[0].c_str());
                 size_t offset = (size_t)atol(records1[1].c_str());
                 for (int index = 0; index <= 200;index++)
                 tot_rd += GetReadDensity(chrnum,offset-100+index);
                 float avg_rd_unique1 = ((float)tot_rd/(float)201);
                 float cutoff = (float) avg_rd_unique1/(float)nrdmap[chrnum];*/
                // cout << "cutoff\t"<< cutoff << endl;
                 int rs_s1 =  (int)( atoi(records1[3].c_str() )*100/atoi( records1[4].c_str() ) +0.5);
                 uniquefile1 << records1[0]<<"\t"<<records1[1]<<"\t"<<records1[2]<<"\t"<<records1[3]<<"\t-\t"<<records1[3]<<"\t"<<records1[4]<<"\t-\t"<<records1[4]<<"\t"<<records1[5]<<"\t-\t"<<records1[5]<<"\t"<<records1[6]<<"\t"<<hethom<<"\t"<<rs_s1<<endl;
                 fcount2++;
              } 
           }
         }
      }
    }
    uniquefile1.close();
    common.close();
}

void getuniquetofile1(char* filename1)
{   
     char line2[5000];
     ofstream uniquefile2("file1_unique");
     uniquefile2 <<"Chr\tOffset\tAllele\tCG1.SNP_Density\tCG2.SNP_Density\tTotal_SNP_Desnity\tCG1.Read_Density\tCG2.Read_Density\tTotal_Read_Desity\tCG1.AVQS\tCG2.AVQS\tAVQS\tdbSNP\tHomo-Het\tRS\n";
     for( m::iterator ii=m1.begin(); ii!=m1.end(); ++ii) 
     {
          m::const_iterator it = m3.find(ii->first);
          if(it == m3.end()) 
          {
             fcount1++;
             strcpy(line2,ii->second.c_str());
             vector <string > records = split(line2);
             float hethomcount =  (float)atof(records[3].c_str())/(float)(atof(records[4].c_str())-atof(records[3].c_str()));
             string hethom;
             if(hethomcount>= 2.333)
             hethom = "homo";
             else
             hethom = "het";

             int tot_rd = 0;
             int chrnum;
             if(records[0]== "X")
             chrnum = 23;
             else if(records[0]== "Y")
             chrnum = 24;
             else
             chrnum = atoi(records[0].c_str());
            /* size_t offset = (size_t)atol(records[1].c_str());
             for (int index = 0; index <= 200;index++)
             tot_rd += GetReadDensity(chrnum,offset-100+index);
             float avg_rd_unique2 = ((float)tot_rd/(float)201);
             float cutoff = (float) avg_rd_unique2/(float)trdmap[chrnum];
          //   cout << "cutoff\t"<< cutoff << endl; */ 
             int rs_s2 =  (int)(atoi(records[3].c_str())*100/atoi(records[4].c_str())+0.5);
             uniquefile2 << records[0]<<"\t"<< records[1]<<"\t"<< records[2]<<"\t-\t"<<records[3]<<"\t"<<records[3]<<"\t-\t"<< records[4]<<"\t"<<records[4]<<"\t-\t"<<records[5]<<"\t"<<records[5]<<"\t"<<records[6]<<"\t"<<hethom<<"\t"<<rs_s2<<endl;
          }
     }
     uniquefile2.close();
}

void get_statistics(char* filename1,char* filename2) 
{
     char flname[200] = "\0";
     strcpy(flname,filename1);
     strcat(flname,".");
     strcat(flname,filename2);
     strcat(flname,".log");
     ofstream summaryfile(flname);
     summaryfile<<filename1<<"_"<<filename2<<"_common"<< "\t" <<commcount <<"\t"<<endl;
     summaryfile<<filename1<<"_unique"<<"\t"<<fcount1<< endl;
     summaryfile<<filename2<<"_unique"<<"\t"<<fcount2<< endl;
     summaryfile<<"Elapsed"<<"\t"<<1000*(clock() - start)/CLOCKS_PER_SEC<<" "<<"ms\n";
     summaryfile.close();
}
void LoadmeanRDNormal(const char* flname)
{
   ifstream rdfile(flname);
   int chrn = 1;
   char nrdensity[500];
   while(rdfile)
   {
     rdfile.getline(nrdensity,500);
     if(rdfile)
      nrdmap[chrn]= atof(nrdensity);
     chrn++;
   }
   rdfile.close();
}

void LoadmeanRDTumor(const char* flname)
{
   ifstream rdfile(flname);
   int chrn = 1;
   char trdensity[500];
   while(rdfile)
   {
     rdfile.getline(trdensity,500);
     if(rdfile)
      trdmap[chrn]= atof(trdensity);
      chrn++;
   }
   rdfile.close();
}


int main(int argc ,char* argv[]) 
{
  if(argc != 5) 
  {
    manual();
    exit(1);
  }
  LoadmeanRDNormal(argv[1]);
  LoadmeanRDTumor(argv[3]);
  string line1,line2;
   ifstream snpfile1 (argv[2]);
   if(!snpfile1.is_open())
   {
      cout<<"Unable to open the snpfile1\n";
      exit(1);
   }
   while(snpfile1) 
   {
          char line1[10000];
          snpfile1.getline(line1,10000);
          if(snpfile1)
            readsnpfile1(line1,atoi(argv[2]));
   }   
   snpfile1.close();

   getcommon_snps(argv[2],argv[4]);
   getuniquetofile1(argv[2]);
   get_statistics(argv[2],argv[4]);
   mread.clear();
   m1.clear();
   m3.clear();
   return 0;
}
