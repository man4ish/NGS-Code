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
M n;
 

clock_t start = clock();
int novcount =0;
int knowncount =0;
int snpcount = 0;


void manual() 
{
    cout<<"\n\nSXDBIndelChecker\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXDBIndelChecker <Chromosome_list> <indelfile> <dbIndel>\n\n";
} 


bool match(string snpa ,string orient,string db) 
{
     string rc; 
     if(orient == "-")
     {
       if(db.at(2)== 'A')
          rc += 'T';
        else  if(db.at(2)== 'T')
                 rc += 'A';
        else  if(db.at(2)== 'C')
                   rc += 'G';
        else  if(db.at(2)=='G')
                   rc += 'C';
     }
    else
    { 
       rc += db.at(2);
    }
  
    if(rc == snpa)
    {    
       return 1;
    } 
    return 0;
}

bool MapAllele(string key,vector<string>indelvector)
{   
     M::iterator it;
     ret = n.equal_range(key);
     for(it=ret.first; it!=ret.second; ++it)
     {
         if(match(indelvector[2],(it->second)[3],(it->second)[4])== 1)
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
     out<<"Total number of indel : "<< snpcount<<endl;
     out<<"Novel indels          : "<<novcount<<"	"<<float (novcount*100)/float(snpcount)<<"%"<<endl;
     out<<"Known indels          : "<<knowncount<<"	"<<float (knowncount*100)/float (snpcount)<<"%"<<endl;
     out<< "Elapsed: "<<(clock() - start)/CLOCKS_PER_SEC<<" "<<"seconds\n";
}
 


void GetNovelIndels(char* filename, string snprec) 
{
     int aa =0;
     char flname[500]="\0";
     char support_read[10];
     strcpy(flname,filename);
     strcat(flname,".mapresult");
     ofstream outc(flname);
     if ( snprec != "")
     {
          string snpl = snprec;
          char* refsnprec = &snprec[0],*pchh;
          vector<string> indellist;
          string ss;
          pchh = strtok (refsnprec,"\t");
          while ( pchh != NULL)
          {
                  indellist.push_back(pchh);
                  pchh = strtok (NULL, "\t");
          }
          string gi,chr;
          gi = indellist[0];
          ss += ChrtoGi[gi];
          ss += "-";
          ss += indellist[1]; 
          snpcount++;
          string chrno = ChrtoGi[indellist[0]];
          if ( chrno == "X" ) 
          {  
               chrno = "23";
          }
          if(  chrno == "Y" )  
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
              aa++;
              if ( MapAllele(ss,indellist)==1)
              {
                   knowncount++; 
                   outc <<"chr"<<chrno<<"\t"<<snpl<<"Yes"<<endl;
              }
              else 
              { 
                  novcount++; 
                  outc<<"chr"<<chrno<<"\t"<<snpl<<"No"<<endl;
              } 
         }
        else 
        { 
             novcount++; 
             outc <<"chr"<<chrno<<"\t"<<snpl<<"No"<<endl; 
        }
    }
}

void SetChrtoGi(string li) 
{
     if ( li != "")
     {
          char* f = &li[0],*pch;
          vector <string> sdf;
          string s;
          pch = strtok(f,"\t");
          while ( pch != NULL) 
          {
                  sdf.push_back(pch);
                  pch = strtok(NULL,"\t"); 
          }
          ChrtoGi[sdf[0]]= sdf[1];
     }
}

void loadindeldb(string lo) 
{
     if ( lo != "")  
     {
          char* f = &lo[0],*pch;
          vector<string> dblist;
          string s,val;
          pch = strtok (f,"\t");
          while ( pch != NULL) 
          {
	          dblist.push_back(pch);
                  pch = strtok (NULL,"\t");
          }
      s += dblist[1]; 
      s += "-"; 
      s += dblist[2];
      n.insert( pair<string,vector<string> >(s,dblist));
     } 
}


int main(int argc, char* argv[])
{
   if(argc !=4) 
   {
      manual();
      exit(1);
   }
        
  int i = 1;
  string dbline,indelline,chrlist;
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

  ifstream indeldb (argv[i]);
  if(indeldb.is_open()) 
  {
     while (!indeldb.eof()) 
     {
      getline (indeldb,dbline);
      loadindeldb(dbline);
      }
     indeldb.close();
  }
  else 
  {
     cout << "Unable to open the database\n";
  }

  i++;
  char* filename = argv[i];
  ifstream indelfile (argv[i]);
  if (indelfile.is_open())
  {
      while(!indelfile.eof() )
      {
             getline (indelfile,indelline);
             GetNovelIndels(filename,indelline); 
      }
      indelfile.close();
  }
  else
  {
      cout << "Unable to open indelfile\n";
  }
  GetStatistics(filename);
  n.clear();
  return 0;
}
