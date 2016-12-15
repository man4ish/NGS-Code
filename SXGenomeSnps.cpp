#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <fstream>

using namespace std;
typedef map<string,string> M;
M m;

void manual()
{
    fprintf(stderr, "\n\nSXGenomeSnps\n");
    fprintf(stderr, "~~~~~~~~~~~~~~~~~~\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\tSXGenomeSnps <query file>\n\n");
}


vector<string> splitv(char* s,const char* token)
{
   char* pch;
   vector<string> v;
   pch = strtok(s,token);
   while(pch!= NULL)
   {
	 v.push_back(pch);
     pch = strtok (NULL,token);
   }
   return v;
}

bool matchquery(const char* a,const char* b)
{
   return (strstr(a,b))?1:0;
}

bool matchchromosome(const char* ch1,const char* ch2)
{
	char s[] = "Homo sapiens chromosome ";
	strcat(s,ch2);
	return (strstr(ch1,s))? 1:0;
}

bool charmatch(string q,string t,int pos)
{
     return (q.at(101-pos)==t.at(101-pos))?1:0;
}


char* getid(char* header)
{
	int pos;
	int i =0;
	char* id;
	int l = strlen(header)-1;
    for(; i< l; i++)
       if(header[i]=='=' && header[i+1]==' ')
          pos=i+2;
          id = new char[l-pos+1];
          memcpy(id,header+pos,l-pos+1);
          id[l-pos+1]='\0';
    return id;
}

unsigned int GetGenomePos(string refid, string chrno)
{
	     if ( m.count(refid)>0)
	     {
	          string a = m[refid];
   	          vector<int> vindex;
		      vindex.push_back(0);
		      int pos = a.find("*>");
		      while (pos!= string::npos)
		      {
		     	     vindex.push_back(pos);
			         pos = a.find("*>",pos+1);
		      }
		      int size = a.size();
		      vindex.push_back(size);
		      for ( int i = 0; i<vindex.size()-1;i++)
		      {
		            string h = a.substr(vindex[i],vindex[i+1]);
		            vector<string>t =  splitv(&h[0],"*");
		            if ( matchchromosome(t[0].c_str(),chrno.c_str()) )
		            {
		       	         for ( int ii =0;ii < t.size();ii++)
			             {
		                       string alignment = t[ii];
		                       vector<string> query = splitv (&alignment[0]," ");
		       	               string a = query[1];
			                   string b = query[3];
			                   int llim = atoi(&a[0]);
                                           int ulim = atoi(&b[0]);
		                       if ( strstr(query[0].c_str(),"Query:") && (llim <= 101 && ulim >=101))
		                       {
				                    string sbjct = t[ii+1];
				                    vector<string> ref  = splitv (&sbjct[0]," ");
				                    if ( charmatch(query[2],ref[2],llim)==1)
				                       {
					                     char* p;
                                                             char* q;
					                     unsigned int te = strtoul ((ref[3]).c_str(), &p, 0);
                                                             unsigned int ts = strtoul ((ref[1]).c_str(), &q, 0);
					                     unsigned int genome_location;
                                                             if(te > ts )
                                                             genome_location    = te - (ulim - 101);
                                                             else
                                                             genome_location = te + (ulim - 101);
					                     return genome_location;
				                       }
			                    }
		                 }
		            }
              }
	 }
return 0;
}


int main(unsigned argc, char* argv[])
{
	if(argc != 3)
	{
		manual();
		exit(1);
	}

    FILE * queryFile = fopen(argv[1], "r");
    char qline[1024] = "\0";
    char qHeader[1024] = "\0";
    char qSequence[65536] = "\0";
    char qSequence1[65536] = "\0";
    clock_t start = clock();
    while(!feof(queryFile))
    {
        fgets(qline, 1024, queryFile);
        if ((qline[0] == 'Q' && qline[5]=='=') || (feof(queryFile)))
        {
            if (qSequence[0] == '\0' )
            {
                qline[strlen(qline) - 1] = '\0';
                strcpy(qHeader, qline);
            }
            else
            {
				char* refseqid = getid(qHeader);
				m[refseqid] = qSequence;
                qSequence[0] = '\0';
                qline[strlen(qline) - 1] = '\0';
                strcpy(qHeader, qline);
            }
        }
        else if(qline[0]=='>'||(qline[0]=='Q' && qline[5]==':')||(qline[0]=='S' && qline[5]==':'))
        {
                        qline[strlen(qline) - 1] = '*';
			            strcat(qSequence, qline);
        }
    }

    fclose(queryFile);

    ifstream infile(argv[2]);
    char filename[500]= "\0";
    strcpy(filename,argv[2]);
    strcat(filename,"_genome.location");
    ofstream out(filename);
    char refseqid[100],chrmno[20],snploc[100],base[10];
    while ( infile>>refseqid>>snploc>>chrmno>>base)
    {
            if ( strcmp(refseqid,"Acc_number")==0)
            {
                  out<<refseqid<<"\t"<<snploc<<"\t"<<chrmno<<"\t"<<base<<"\t"<<"GenomePosition"<<endl;
            }
            else
            {
                 if ( GetGenomePos(refseqid,chrmno)!=0)
                    out<<refseqid<<"\t"<<snploc<<"\t"<<chrmno<<"\t"<<base<<"\t"<<GetGenomePos(refseqid,chrmno)<<endl;       
          else
                     out<<refseqid<<"\t"<<snploc<<"\t"<<chrmno<<"\t"<<base<<"\t"<<"XXXX"<<endl;
            }
	}

    printf("Elapsed: %.4f ms\n\n", (float)(clock() - start)*1000/CLOCKS_PER_SEC);
    m.clear();
    return 0;
}
 
