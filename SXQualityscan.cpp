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
#include <set> 
using namespace std;

struct Query
{
    char * queryHeader;
    char * querySequence;
    Query(char * _qHeader, char * _qSequence)
    {
        queryHeader = new char [strlen(_qHeader) + 1];
        strcpy(queryHeader, _qHeader);
        querySequence = new char [strlen(_qSequence) + 1];
        strcpy(querySequence, _qSequence);
    }
    ~Query()
    {
        if (queryHeader != NULL)
        {
            delete [] queryHeader;
            queryHeader = NULL;
        }
        if (querySequence != NULL)
        {
            delete [] querySequence;
            querySequence = NULL;
        }
    }
};

struct ltstr
{
    bool operator()(const char * _lhs, const char * _rhs) const
    {
        return strcmp(_lhs, _rhs) < 0;
    }
};

unsigned long long bp=0;
unsigned long long length=0;
unsigned long long nonagctcount=0;
unsigned long long gccount=0;
unsigned long long sxtnqs = 0;
unsigned long long twntyqs = 0;
unsigned long long avtenqs = 0;
unsigned long long avsxtnqs = 0;
unsigned long long avtwntyqs = 0;
int maxlength=0;
int minlength=5000;  
bool flag;
int readcount=0;
int xmin= 100;
int xmax = -100;

vector<string> veccon,panvec,seqd;
vector<unsigned int> AVQSset;
map<int, int>m, mlen;
std::vector <Query *> qryVector;
std::vector <Query * > dupQuery;
typedef std::map<char *, char *, ltstr> QueryMap;
typedef QueryMap::iterator QueryMapIT;
typedef std::pair<char *, char *> QueryMapPAIR;
std::pair<QueryMapIT, bool> ret;
QueryMap qryMap;

void manual()
{
    fprintf(stderr, "\n\nSXQualityScan\n");
    fprintf(stderr, "~~~~~~~~~~~~~~~~~~\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\tSXQualityScan -v <query file> <quality file> -d duplicates -u unique >output_file\n\n");
}

char* chomp(char* s)
{
      int l = strlen(s)-1;
      if(s[l]=='\n')
      s[l]='\0';
}

void SetLength(int i)
{
	if(mlen.count(i)== 0)
	mlen[i]=1;
	else 
	mlen[i] += 1;
}

void Setaveragelength(int a)
{
  	length +=a; 
}
  
void Setmaxlength(int sc)
{
     if(sc>maxlength)
     {
        int templength= sc;
        sc=maxlength;
        maxlength=templength;                   
     }   
}

void Setminlength(int sc)
{
     if(sc<minlength)
     {
        int templength=sc;
        sc=minlength;
        minlength=templength;
     }     
}

void Getlength(char*s2)
{
     chomp(s2);
     int lenread=0;
     while(*s2 !='\0' &&*s2!='\n')
     {
           lenread++;
           *s2++;
     }
     SetLength(lenread);
     Setmaxlength(lenread);
     Setminlength(lenread);
     Setaveragelength(lenread);
}
  

int  Getmaxlength()
{
     return maxlength;
}

int Getminlength()
{
     return minlength;
}

void Setreadcount()
{
     readcount ++;
}


unsigned long long Getreadcount()
{
     return readcount;
}

long double Getaveragelength()
{
      long double f = double(length)/double(Getreadcount());
      return f;
}

void SetNonAGCTcount(char* a)
{
     while(*a)
  	{      
       	if(*a!='A'&& *a!='C'&& *a!='G'&& *a!='T'&& *a!='\n')
          {
            nonagctcount++;           
          }
   	*a++;
  	}
}

unsigned long long GetNonAGCTbases()
{
    return nonagctcount;
}

unsigned long long GetNonAGCTcount()
{
  	return nonagctcount;
}

void SetMinimumQScore(int a )
{
  	int temp;
  	if(a<xmin)
  	{
  		temp = a;
  		a=xmin;
  		xmin=temp;
  	}
}

void  SetMaximumQScore(int a)
{
  	int temp;
  	if(a >xmax)
  	{
  		temp = a;
  		a=xmax;
  		xmax=temp;
  	}
}

long double GetMorethansxtenqs()
{
     long double persxtnqs = double ((sxtnqs)*100)/double (bp);
     return persxtnqs; 
}


long double GetMorethantwntyqs()
{
     long double pertwntyqs = double ((twntyqs)*100)/double (bp);
     return pertwntyqs;
}

long double GetMorethantenavqs()
{
     long double peravtenqs = double ((avtenqs)*100)/double (readcount);
      return peravtenqs;
}

long double GetMorethansxtenavqs()
{
     long double peravsxtnqs = double ((avsxtnqs)*100)/double (readcount);
     return peravsxtnqs;
}

long double GetMorethatwntyavqs()
{
     long double peravtwntyqs = double ((avtwntyqs)*100)/double (readcount);
     return peravtwntyqs;
}

int GetMinimumQScore()
{
  	return xmin;
}

int GetMaximumQScore()
{
  	return xmax;
}

void Settotalbp(char* a)
{
     int end = strlen(a)-1;
    if(a[end]== '\n')
       a[end]='\0';
  	while(*a)
  	{
  	bp++;
  	*a++;
  	}	
}

unsigned long long Gettotalbp()
{
  	return bp;
}

long double GetpercentageNONAGCT()
{
       
  	long double f =double (nonagctcount*100)/double (Gettotalbp());
        return f;
}

unsigned int GetMaxAVQS()
{
         return *max_element(AVQSset.begin(), AVQSset.end());
}

unsigned int GetMinAVQS()
{
         return *min_element(AVQSset.begin(), AVQSset.end());
}

void SetGCcount(char* a)
{   
        int end = strlen(a)-1;
        if(a[end]=='\n')
        a[strlen(a)-1]='\0';
        int i=0;
   	while(*a)
   	{            
   		if(*a=='G'||*a=='C')
   		gccount++;
   		*a++;
   	}
}

long double GetGCcount()
{
      long double f = double(gccount*100)/double(Gettotalbp());
    	return f;
}

void Setscoredistribution(int i)
{
	if(m.count(i) == 0)
	m[i]=1;
	else 
	m[i] += 1;
}

void split(char* array )
{
    int lenseq = 0;
    int totalqs =0;
    const char* delimiter = " ";
    char * pch;
    
    pch = strtok (array,delimiter);
    while (pch != NULL && *pch!='\n')
    {
       string buffer = pch;
       int integer=atoi(buffer.c_str());
       if ( integer > 16) 
       {
            sxtnqs++;        
       }
       if ( integer > 20)
       {
            twntyqs++;
       } 
       Setscoredistribution(integer);
       totalqs += integer;
       lenseq++;
       SetMinimumQScore(integer);
       SetMaximumQScore(integer);
       pch = strtok (NULL, delimiter);
    }   
       unsigned int avqs = int (totalqs/lenseq);

       if ( avqs > 10)
       {
            avtenqs++;
       }
       if ( avqs > 16)
       {
            avsxtnqs++;
       }
       if ( avqs > 20)
       {
            avtwntyqs++;
       } 
       AVQSset.push_back(avqs);
}


void GetQualityStatistics(char* Qualityfile)
{ 
     ofstream outfile;
     outfile.open(Qualityfile);
     if (outfile.is_open())
     {
         outfile<<"Quality Score:frequency"<<endl;
         for( map<int,int>::iterator ii=m.begin(); ii!=m.end(); ++ii)
         {
              outfile<<(*ii).first<<" : "<<(*ii).second<<endl;
         }
        outfile.close();
    }
    else
    {
        cout << "Error opening file" << Qualityfile<<endl;
    }
}

void Removeduplicates()
{
    unsigned totalQueries = qryVector.size();   
    for (unsigned i = 0; i < totalQueries; i++)
    { 
         ret = qryMap.insert(QueryMapPAIR(qryVector[i]->querySequence, qryVector[i]->queryHeader));           
         if (ret.second == false)          
             dupQuery.push_back(new Query(qryVector[i]->queryHeader, qryVector[i]->querySequence));
    }
    QueryMapIT qryMapIterator;
}

void GetUnique(char* Fileuniqe)
{
     QueryMapIT qryMapIterator;
     unsigned totalUnique = qryMap.size();
     unsigned uniqueCnt = 0;
     ofstream outfile;
     outfile.open(Fileuniqe);
     if (outfile.is_open())
     {
        for(qryMapIterator = qryMap.begin(); qryMapIterator != qryMap.end(); qryMapIterator++)
        {
            if(flag==1)
            {
               string tempmodh =  (*qryMapIterator).second;
               string tempmods =  (*qryMapIterator).first;
               int y = tempmodh.find('\n');
               tempmodh.insert(y+1,tempmods);
               outfile<<tempmodh;  
            }
            else 
           {
              outfile<<(*qryMapIterator).second<<endl;
              outfile<< (*qryMapIterator).first<<endl;
           }
        }    
        outfile.close();
     }
     else
     {
        cout << "Error opening file";
     } 
}

void Getduplicates(char* duplicatesFileName)
 {
     FILE * duplicatesFile = fopen(duplicatesFileName, "w");
     if (!dupQuery.empty())
     {
        for (unsigned i = 0; i < dupQuery.size(); i++)
        {
            if(flag ==1)
            { 
               string tempmodh =  dupQuery[i]->queryHeader;
               string tempmods = dupQuery[i]->querySequence;
               int y = tempmodh.find('\n');
               tempmodh.insert(y+1,tempmods);
               fprintf(duplicatesFile, "%s",&tempmodh[0]);
            }
            else
            {
            fprintf(duplicatesFile, "%s\n%s\n", dupQuery[i]->queryHeader, dupQuery[i]->querySequence);
            }
            delete dupQuery[i];
        }
        dupQuery.clear();
    }
    fclose(duplicatesFile);
 }

void FreeMemory()
{
   if (!qryVector.empty())
    {
        for (unsigned i = 0; i < qryVector.size(); i++)
        {
            delete qryVector[i];
        }
        qryVector.clear();
    } 
    qryMap.clear();
    m.clear();
    mlen.clear();
}

void GetLengthStatistics(char* Qualityfile)
{ 
     ofstream outfile;
     outfile.open(Qualityfile);
     if (outfile.is_open())
     {
         outfile<<"Read length : percentage of occurence"<<endl;
         for( map<int,int>::iterator ii=mlen.begin(); ii!=mlen.end(); ++ii)
         {
              double percentage = float(((*ii).second)*(100))/float (readcount);
              outfile<<(*ii).first<<" : "<<percentage<<endl;
         }
        outfile.close();
     }
    else
    {
        cout << "Error opening file" << Qualityfile<<endl;
    }
}

bool Matchextenstion(char* s, const char*q)
{
     int len = strlen(s);
     int qlen= strlen(q);
     int i=0;
     while(i<qlen)
     {
           if(q[qlen-i]!=s[len-i])
           {
              return false;
           }         
           i++;
     }     
     return true;
 }

int main(int argc, char*argv[])
{ 
    bool uniqueflag,dupflag;
    char inputfile[20]="\0";
    char inputfile2[40]="\0";
    char qlength[40]="\0";
    char qscore[40]="\0";
    char summary[40]="\0";
    char qHeader[40]="\0";
    char qSequence[40]="\0";
    char unique[40]="\0";
    char duplicate[40]="\0";
  
    if(argc > 8||argc==1)
    {
       manual();
       exit(1);
    }   
      
   clock_t start= clock(); 
   for (int i =0;i<argc;i++)
   {  
       if( strstr(argv[i],"-v")) 
       {
           argv[++i];
           if (Matchextenstion(argv[i],".fna") && Matchextenstion(argv[i+1],".qual"))
           {
               strcpy(inputfile,argv[i]);
               strcpy(inputfile2,argv[i+1]);
               strcat(qlength,inputfile);
               strcat(qlength,".qscan.lengthstat");
               char qline[1024] = "\0";
               char qHeader[1024] = "\0";
               char qSequence[65536] = "\0";
               FILE * queryFile = fopen(inputfile, "r");
               while (!feof(queryFile))
               {
                      fgets(qline, 1024, queryFile);
                      if((qline[0] == '>') || (feof(queryFile)))
                      {        
                          if (qSequence[0] == '\0')
                          {
                              qline[strlen(qline) - 1] = '\0';
                              strcpy(qHeader, qline);
                          }
                          else
                          {
                              Setreadcount();
                              qryVector.push_back(new Query(qHeader, qSequence));                  
                              Getlength( qSequence);
                              SetNonAGCTcount(qSequence);
                              Settotalbp(qSequence);
                              SetGCcount(qSequence);
                              qSequence[0] = '\0';
                              qline[strlen(qline) - 1] = '\0';
                              strcpy(qHeader, qline);
                          }
                     }
                     else
                     {
                         qline[strlen(qline) - 1] = '\0';
                         strcat(qSequence, qline);
                     }
       }
       fclose(queryFile);
 
       FILE* fp;
       char qSequence4[65536] = "\0";
       char qSequence5[65536] = "\0";
       char sss[1000];
       if ((fp = fopen(inputfile2,"r"))==NULL)
       {
         printf("Cannot open Quality file.\n");
         exit(1);
       }
       fp = fopen(inputfile2,"r");  
       while (!feof(fp))
          {
            fgets(sss,999,fp);
             if(sss[0]=='>' || (feof(fp)))
             {
                if (qSequence4[0] != '\0')
                {
                    strcpy(qSequence5,qSequence4);
                    qSequence5[strlen(qSequence5)-1]= '\0';
                    split(qSequence5);
                    qSequence4[0] ='\0';
                } 
             }
            else
            {
               if(sss[strlen(sss)-1]=='\n')
               {
                sss[strlen(sss)-1]=' ';
               }
               strcat(qSequence4, sss);
            }
      }
      fclose(fp);
      strcat(qscore,inputfile);
      strcat(qscore,".qscan.scorestat"); 
      strcat(summary,inputfile);
      strcat(summary,".qscan.summary");
      argv[i+2];
   }
   else if(Matchextenstion(argv[i],".fastq"))
   {
      flag =1;
      strcpy(inputfile,argv[i]);
      strcat(qlength,inputfile); 
      strcat(qlength,".qscan.lengthstat");
      strcat(qscore,inputfile);
      FILE * queryFile = fopen(inputfile, "r");
      if (queryFile == NULL)
      {
          manual();
          return -1;
      }
      char qline[1024] = "\0";
      char qHeader[1024] = "\0";
      char qSequence[65536] = "\0";
      while(!feof(queryFile))
      {
            fgets(qline, 1024, queryFile);

            if(isalpha(qline[0]))
            {        
               string adf = qline;
               veccon.push_back(adf);
               Setreadcount(); 
               Getlength(qline);
               SetNonAGCTcount(qline);
               Settotalbp(qline);
               SetGCcount(qline);                
           }
           else
           {
               panvec.push_back(qline);
               if(isdigit(qline[0]))
               {
                  split(qline);
               }
           }
      }
            
      strcat(qscore,".qscan.scorestat"); 
      strcat(summary,inputfile);
      strcat(summary,".qscan.summary");
      argv[++i]; 
  }
  else
  {
      manual();
      exit(1);
  }
 } 

 if(Matchextenstion(argv[i],"-d") && Matchextenstion(argv[i+1],"duplicates"))
 {
    strcat(duplicate,inputfile);
    strcat(duplicate,".qscan.duplicates");
    dupflag=1;
  } 

  if(Matchextenstion(argv[i],"-d") && Matchextenstion(argv[i+1],"duplicates")=='\0')
  { 
     manual();
     exit(1);
  }
  
  if(Matchextenstion(argv[i],"-u") && Matchextenstion(argv[i+1],"unique")=='\0')
  { 
     manual(); 
     exit(1);
  }

  if(Matchextenstion(argv[i],"-u") && Matchextenstion(argv[i+1],"unique"))
  {
     char* pfile,*pff;
     int leninputfile = strlen(inputfile);
     int x, fi = 0;
     for(;fi < leninputfile;fi++)
     { 
         if(*(inputfile + fi)=='.')
         {
            x= fi;
         }
     }
     pfile = new char[x+1];
     memcpy(pfile,inputfile,x);
     pfile[x]='\0'; 
     strcat(unique,pfile);     
     strcat(unique,".unique");
     pff = new char[leninputfile-x];
     memcpy(pff,inputfile+x,leninputfile-x);
     pff[leninputfile-x]='\0';
     strcat(unique,pff);
     uniqueflag=1;
   }        
 }
   
   if ( flag ==1)
   {
        for ( int ind =0; ind <panvec.size();ind++)
        {
              if ( panvec[ind].at(0)=='@')
              {
                   string sequene;  
                   sequene += panvec[ind];                          
                   sequene +=panvec[ind+1];
                   sequene += panvec[ind+2];
                  seqd.push_back(sequene);
              }
        }
  
        for ( int index =0;index<veccon.size();index++)
        {
              string hd = seqd[index];
              string se = veccon[index];
              qryVector.push_back(new Query(&hd[0], &se[0]));
        }   
   }

   GetQualityStatistics(qscore);
   GetLengthStatistics(qlength);
   FILE * summaryFile = fopen(summary, "w");
   if (summaryFile == NULL)
   {
       fprintf(stderr, "\nOutput file error\n");
       return -1;
   }
   if(dupflag==1||uniqueflag==1)
   {
       Removeduplicates();
   }
   if( dupflag==1)
   {
       Getduplicates(duplicate);
   }
   if( uniqueflag==1)
   {
       GetUnique(unique);
   }

   time_t rawtime;
   struct tm * timeinfo;
   char curtime [80];
   char curdate [80];
   time ( &rawtime );
   timeinfo = localtime ( &rawtime );
   strftime (curdate,80,"%x",timeinfo);
   strftime (curtime,80,"%X",timeinfo);
   fprintf(summaryFile,"File Name:		                   %s\n",inputfile);
   fprintf(summaryFile,"Date/Time:		                   %s %s\n",curdate,curtime);
   fprintf(summaryFile, "Total Elaspse:         	                   %.2f seconds\n", (float)(clock() - start)/CLOCKS_PER_SEC);
   fprintf(summaryFile, "Total Base Pair:	                   %llu\n", Gettotalbp());
   fprintf(summaryFile, "Number of Reads:	                   %d\n", qryVector.size());
   if(uniqueflag==1)
   {
   fprintf(summaryFile, "Number of Unique Reads:	                   %d\n", qryMap.size());
   }
   if(dupflag==1)
   {
   fprintf(summaryFile, "Number of Duplicate Reads:	           %d (%.2f%%)\n", qryVector.size() - qryMap.size(), (float)  (qryVector.size() - qryMap.size())/(float)qryVector.size() * 100.0);
   }
   fprintf(summaryFile,"Average Length:                            %Lf\n",Getaveragelength());
   fprintf(summaryFile, "Shortest Read:	                           %d\n",Getminlength());
   fprintf(summaryFile, "Longest Read:	                           %d\n",Getmaxlength());
   fprintf(summaryFile, "Minimum Quality Score:	                   %d\n",GetMinimumQScore());
   fprintf(summaryFile, "Maximum Quality Score:            	   %d\n", GetMaximumQScore());
   fprintf(summaryFile, "Non AGCT bases (%):             	   %Lf\n", GetpercentageNONAGCT());
   fprintf(summaryFile, "Bases QS >16   (%):                        %Lf\n", GetMorethansxtenqs());
   fprintf(summaryFile, "Bases QS >20   (%):                        %Lf\n", GetMorethantwntyqs());
   fprintf(summaryFile, "Reads with Average QS > 10 (%):            %Lf\n", GetMorethantenavqs());
   fprintf(summaryFile, "Reads with Average QS > 16 (%):            %Lf\n", GetMorethansxtenavqs());
   fprintf(summaryFile, "Reads with Average QS > 20 (%):            %Lf\n", GetMorethatwntyavqs());
   fprintf(summaryFile, "NON AGCT Bases (bp):           	           %llu\n", GetNonAGCTbases());
   fprintf(summaryFile, "Minimum AVQS:                              %llu\n",  GetMinAVQS());
   fprintf(summaryFile, "Maximum AVQS:                              %llu\n",  GetMaxAVQS());  
   fprintf(summaryFile, "Average GC%:                    	   %Lf\n",GetGCcount());   
   FreeMemory();
   return 0;
}


