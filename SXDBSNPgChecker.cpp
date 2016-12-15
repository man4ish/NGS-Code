#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <time.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

struct Key
{
    int gi;
    size_t offset;
    Key(int  _gi, size_t  _offset):gi(_gi),offset(_offset) {}
};

struct Value
{
    char * orient;
    char * allele ;
    char* rsid;
    Value(const char* _rsid, const char * _orient, const char * _allele)
    {
        rsid = new char [strlen(_rsid) + 1];
        strcpy(rsid, _rsid); 
        orient = new char [strlen(_orient) + 1];
        strcpy(orient, _orient);
        allele = new char [strlen(_allele) + 1];
        strcpy(allele, _allele);
    }
    ~Value()
    {
        if (rsid != NULL)
        {
            delete [] rsid;
            rsid = NULL;
        }
 
        if (orient != NULL)
        {
            delete [] orient;
            orient = NULL;
        }
        if (allele != NULL)
        {
            delete [] allele;
            allele = NULL;
        }
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi)        
         return ( _lhs->gi < _rhs->gi);
           return  ( _lhs->offset < _rhs->offset );
    }
};

map<int ,long> ChrtoGi;
typedef multimap< Key*,Value*,ltstr> dbmap;
pair< multimap<Key*,Value*, ltstr> ::iterator,multimap<Key*,Value*, ltstr> ::iterator> ret;
dbmap n;

clock_t start = clock();
int novcount =0;
int knowncount =0;
int snpcount =0;

void manual()
{
    cout<<"\n\nSXDBSNPgChecker\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXDBSNPgChecker <Chromosome_list> <dbSNP> <Snp_file> > <outfile>\n\n";
}

const char* ReverseComplement(char* seq)
{
            string revseq;
            int l = strlen(seq);
            for( int i = l-1; i >= 0; i--)
            {
                 switch (seq[i])
                {
                         case 'A':
                         revseq += "T";
                         break;
                         case 'C':
                         revseq += "G";
                         break;
                         case 'G':
                         revseq += "C";
                         break;
                         case 'T':
                         revseq += "A";
                         break;
                }
            }
            return revseq.c_str();
}

bool match(const char* allele ,const char* orient,const char* db) 
{
     vector <string> snpvec,dbvec;
     if(allele)
     {
        char *pch;
        int lenallele = strlen(allele);
        char* nallele = new char[lenallele+1];
        strcpy(nallele,allele);
        nallele[lenallele] = '\0';
        pch = strtok(nallele," > ");
        while(pch != NULL) 
        {
             snpvec.push_back(pch);
             pch = strtok(NULL," > ");
        }
         delete [] nallele;
      }
     if(db) 
     {          
        int lendb = strlen(db);
        char* ndb = new char[lendb+1];
        strcpy(ndb,db);
        ndb[lendb] = '\0';
        char* pch1;
        pch1 = strtok(ndb,"/");
        while(pch1 != NULL) 
        {
              if(strcmp(orient,"-")==0)
              {
                 dbvec.push_back(ReverseComplement(pch1));
              }
              else 
              {
                 dbvec.push_back(pch1);
              }
              pch1 = strtok(NULL,"/");
        }
        delete[] ndb;  
     }


  for(int ii =0; ii< dbvec.size();ii++)
  {  
      if(dbvec[ii] == snpvec[1])
      {
         return 1;
      }   
  }
  dbvec.clear();
  return 0;
}

bool MapAllele( Key* location ,const char* variant)
{
     dbmap::iterator it;
     ret = n.equal_range(location);     
     for(it=ret.first; it!=ret.second; ++it)
     {
         if(match(variant,(it->second)->orient,(it->second)->allele))
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
     out.close();
}


void GetNovelSnps(char* filename,char* snpval) 
{
     char flname[500]="\0";
     char support_read[10];
     strcpy(flname,filename);
     strcat(flname,".mapresult");
     ofstream outc(flname);
     if (!outc.is_open())
     {
         cout<<"could not open the wrtie file\n";
         manual();
     }

     if (snpval)
     {
          char snprec[1000]="\0";
          strcpy(snprec,snpval);
          char* pchh;
          vector<const char*> snpvect;
          string ss;
          pchh = strtok (snpval,"\t");
          while ( pchh != NULL)
          {
                  snpvect.push_back(pchh);
                  pchh = strtok (NULL, "\t");
          }
          snpcount++;
          int chrno = ChrtoGi[atoi(snpvect[0])];
  
          Key * snpkey = new Key(chrno, (size_t)atol(snpvect[1]));
         dbmap::const_iterator elefound; 
         elefound = n.find (snpkey);
        if( chrno > 9)
        {
           cout<<"chr"<<chrno<<"\t"<<snprec<<"\t"; 
        }
        else 
        {
          cout<<"chr"<<setfill('0')<<setw(2)<<chrno<<"\t"<<snprec<<"\t";
        }
        if (elefound != n.end ())
        {
            if(MapAllele(snpkey,snpvect[2])==1)
              {
                 knowncount++; 
                 cout << elefound->second->rsid<<endl;
              }
              else 
              { 
                 novcount++; 
                 cout<<"No"<<endl;
              }
  
        } 
        else 
          {
              novcount++; 
              cout <<"No"<<endl; 
          }
    }
  outc.close();
}

void SetChrtoGi(char* chromlist) 
{
   if(chromlist)
   {
      char *pch;
      vector <const char*> chrgi;
      pch = strtok(chromlist,"\t");
      while(pch != NULL) 
      {
            chrgi.push_back(pch);
            pch = strtok(NULL,"\t"); 
      }
      ChrtoGi[atoi(chrgi[0])]= atol(chrgi[1]);
   }
}

void loadsnpdb(char* db) 
{
     if(db)  
     {
        char* pch;
        vector< const char* > dblist;
        pch = strtok (db,"\t");
        while(pch != NULL) 
        {
	      dblist.push_back(pch);
              pch = strtok (NULL,"\t");
        }
        int chr_num;
        unsigned int offset;
        if(strcmp(dblist[1],"X")==0)
        {
           chr_num = 23;
        }
        else if(strcmp(dblist[1],"Y")==0)
        {           
           chr_num = 24;
        }
        else 
        {
          chr_num = atoi(dblist[1]);
        }
        Key * key = new Key(chr_num, (size_t)atol(dblist[2]));
	Value * value = new Value(dblist[0],dblist[3], dblist[4]);
	n.insert(pair<Key *, Value *>(key, value));
     } 
}

int main(int argc, char* argv[])
{

        if (argc !=4) 
        {
            manual();
            exit(1);
        }
        int i =0;
        i++;
        char snpline[1000];
        char chrlist[1000];

        ifstream chrfile (argv[i]);
        if (!chrfile.is_open())
        {
             manual();
             exit(1);
        } 
        while(chrfile) 
        {
             chrfile.getline(chrlist,1000);
             if(chrfile)
             {
                SetChrtoGi(chrlist);
             }
        }
        chrfile.close();
     //   cout<<"chromosome file loaded\n";
        i++;
  
        ifstream dbfile (argv[i]);
        if (!dbfile.is_open())
        {
             manual();
             exit(1);
        } 
        char dbline[1000];
        while (dbfile)
        {
               dbfile.getline (dbline,1000);
               if(dbfile)
               {
                  loadsnpdb(dbline);
               }
         }
         dbfile.close();
        // cout<<"dbSNP file loaded\n";
         i++;
 
         char* filename = argv[i];
         ifstream snpfile (argv[i]);
         if (!snpfile.is_open())
         {
              manual();
              exit(1);
         } 
         while (snpfile )
         {
                snpfile.getline (snpline,1000);
                if( snpfile)
                { 
                  GetNovelSnps(filename,snpline);
                }
         }
         snpfile.close();
         GetStatistics(filename);
         n.clear();
     return 0;
}   
