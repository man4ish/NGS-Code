#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>

using namespace std;
map<string,string> ChrtoGi;
typedef multimap<string,vector<string> > M;
pair< multimap <string,vector<string> > ::iterator,multimap <string,vector<string> > ::iterator> ret;
M m,n;

clock_t start = clock();
int novcount =0;
int knowncount =0;
int snpcount =0;

void manual()
{
    cout<<"\n\nSXDBChecker\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXDBChecker <Chromosome_list> <dbSNP> <Snp_file>\n\n";
}

char* ReverseComplement( char* s)  //function for caluclating the reverse complement of thre allele
{    
      int l = strlen(s);
      string rc ;
      for(int i = l-1; i >= 0; i--) 
      {
          if(s[i]== 'A')
             rc += 'T';
          else  if(s[i]== 'T')
                   rc += 'A';
          else  if(s[i]== 'C')
                   rc += 'G';
          else  if(s[i]=='G')
                   rc += 'C';
       }
 return &rc[0];
}

bool match(string snpa ,string orient,string db) 
{
     vector <string> snpvec,dbvec;
     if(snpa != "")
     {
       char* f = &snpa[0],*pch;
       pch = strtok(f," > ");
       while(pch != NULL) 
       {
             snpvec.push_back(pch);
             pch = strtok(NULL," > ");
       }
    }

  if(db != "") 
    {
     char* f = &db[0];
       if(orient == "-") 
       {
          char* pch1;
          pch1 = strtok(f,"/");
          while(pch1 != NULL) {
          char* v =  ReverseComplement(pch1);
          dbvec.push_back(v);
          pch1 = strtok(NULL,"/");
       }
    }
 else {
         char* pch2;
         pch2 = strtok(f,"/");
         while(pch2 != NULL){
         dbvec.push_back(pch2);
         pch2 = strtok(NULL,"/");
      }
   }
}

  for(int ii =0; ii< dbvec.size();ii++)
  {  
      if(dbvec[ii] == snpvec[1])    
         return 1;   
  }
  dbvec.clear();
  return 0;
}

bool MapAllele(string key,vector<string>snpvector)
{
     M::iterator it;
     ret = n.equal_range(key);
     for(it=ret.first; it!=ret.second; ++it)
     {
         if(match(snpvector[2],(it->second)[3],(it->second)[4])== 1)
         return 1;
     }
     return 0;
}


 
void GetStatistics(char* filename)
{
     char flname[200] ="\0";
     char support_read[10];
     strcpy(flname,filename);
     strcat(flname,".summary");
     ofstream out(flname);
     out<<"Total number of Snps : "<< snpcount<<endl;
     out<<"Novel Snps           : "<<novcount<<"        "<<float (novcount*100)/float(snpcount)<<"%"<<endl;
     out<<"Known Snps           : "<<knowncount<<"      "<<float (knowncount*100)/float (snpcount)<<"%"<<endl;
     out<< "Elapsed: "<<(clock() - start)/CLOCKS_PER_SEC<<" "<<"seconds\n";
}


void GetNovelSnps(char* filename,string snpval) 
{
     char flname[500]="\0";
     char support_read[10];
     strcpy(flname,filename);
     strcat(flname,".mapresult");
     ofstream outc(flname);

     if ( snpval != "")
     {
          string snpline = snpval;
          char* ff = &snpval[0],*pchh;
          vector<string> snpvect;
          string ss;
          pchh = strtok (ff,"\t");
          while ( pchh != NULL)
          {
                  snpvect.push_back(pchh);
                  pchh = strtok (NULL, "\t");
          }
          snpcount++;
          string gi,chr;
          gi = snpvect[0];
          ss += ChrtoGi[gi];
          ss += "-";
          ss += snpvect[1];
          string chrno = ChrtoGi[snpvect[0]];
          if ( chrno == "X" )
          {
               chrno = "23";
          }
          if ( chrno == "Y" )
          {
               chrno = "24";
          }
          int chrint = atoi(chrno.c_str());
          if ( chrint <= 9)
          {
               chrno = "0" + chrno;
          }
          if ( n.count(ss) > 0)
          {
 
              if(MapAllele(ss,snpvect)==1)
              {
                 knowncount++; 
                 outc <<"chr"<<chrno<<"\t"<<snpline<<"\tYes"<<endl;
              }
              else 
              { 
                 novcount++; 
                 outc<<"chr"<<chrno<<"\t"<<snpline<<"\tNo"<<endl;
              } 
          }
          else 
          {
              novcount++; 
              outc <<"chr"<<chrno<<"\t"<<snpline<<"\tNo"<<endl; 
          }
    }
}

void SetChrtoGi(string li) 
{
   if(li != "")
   {
      char* f = &li[0],*pch;
      vector <string> sdf;
      string s;
      pch = strtok(f,"\t");
      while(pch != NULL) 
      {
            sdf.push_back(pch);
            pch = strtok(NULL,"\t"); 
      }
      ChrtoGi[sdf[0]]= sdf[1];
   }
}

void loadsnpdb(string lo) 
{
     if(lo != "")  
     {
        char* f = &lo[0],*pch;
        vector<string> dblist;
        string s,val;
        pch = strtok (f,"\t");
        while(pch != NULL) 
        {
	      dblist.push_back(pch);
              pch = strtok (NULL,"\t");
        }
      s += dblist[1]; s += "-"; s += dblist[2];
      n.insert( pair<string,vector<string> >(s,dblist));
     } 
}



int main(int argc, char* argv[])
{
   if(argc !=4 && argc != 6) 
   {
      manual();
      exit(1);
   }
        
  int i = 1;
  string dbline,snpline,chrlist;
  ifstream chrfile (argv[i]);
  if(chrfile.is_open()) 
  {
     while(!chrfile.eof()) 
     {
            getline(chrfile,chrlist);
            SetChrtoGi(chrlist);
     }
    chrfile.close();
  }
  else
  {
        cout << "Unable to open chromosome file\n";
  }
  i++;
  
   ifstream myfile1 (argv[i]);

  ifstream snpfile (argv[i]);
    if ( myfile1.is_open())
    {
         while ( !myfile1.eof())
         {
                 getline (myfile1,dbline);
                 loadsnpdb(dbline);
         }
         myfile1.close();
   }
   else
   {
         cout << "Unable to open the database\n";
   }

   i++;
   char* filename = argv[i];
   ifstream myfile (argv[i]);

   if ( myfile.is_open())
   {
        while ( !myfile.eof() )
        {
                getline (myfile,snpline);
                GetNovelSnps(filename,snpline);
        }
        myfile.close();
   }
   else
   {
       cout << "Unable to open snp file\n";
   }

   GetStatistics(filename);
   n.clear();
   return 0;
}   
