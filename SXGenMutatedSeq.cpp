#include <iostream>
#include <vector>
#include <stdio.h>
#include <time.h>
#include <map>
#include <fstream>
#include <time.h>
using namespace std;

clock_t start = clock();
typedef map<string, string> M;
M m;

void manual()
{
    cout<< "\n\n SXMutatedSequence\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXMutatedSequence <rna file> <transfile>\n\n";
}

char getbase(const char* Refseqid , int pos)
{
     if ( m.count(Refseqid)>0 )
     {
          string sequence = m[Refseqid];
          return sequence.at(pos-1);
      }
}

string  getid(char* header)
{
  	string s = header;
  	int y =s.find("|ref|");
	int z= s.size();
	int len = z-y-6;
	string gid = s.substr(y+5,len) ;
        int p = gid.find("|");
        string geneid =  gid.substr(0,p);
	return geneid;
}

char* getsequence(const char* transid, int pos)
{
       if(m.count(transid)>0) 
       {
       string seqi = m[transid];
       char* seq = &seqi[0];
       int start = pos-101;
       char* sequence = new char[201];
       memcpy(sequence,seq+start,201);
       sequence[201]= '\0';
       return sequence;
       }
       else 
       {
        return "";
       }
}

int main(int argc, char* argv[])
{
       if(argc != 3)
       {
          manual();
          exit(1);
       }
    char qline[1024] = "\0";
    char qHeader[1024] = "\0";
    char qSequence[655366] = "\0";
    string RefseqAcc;
    FILE * queryFile = fopen(argv[1], "r");
    if (queryFile == NULL)
    {
        manual();
        return -1;
    }

    while (!feof(queryFile))
    {
        fgets(qline, 1024, queryFile);
        if ((qline[0] == '>') || (feof(queryFile)))
        {
            if (qSequence[0] == '\0')
            {
                qline[strlen(qline) - 1] = '\0';
                strcpy(qHeader, qline);
                RefseqAcc=getid(qHeader);
            }
            else
            {
                m[RefseqAcc]= qSequence;
                qSequence[0] = '\0';
                qline[strlen(qline) - 1] = '\0';
                strcpy(qHeader, qline);
                RefseqAcc = getid(qHeader);
            }
        }
        else
        {
            qline[strlen(qline) - 1] = '\0';
            strcat(qSequence, qline);
        }
    }
 fclose(queryFile);
  
    char refid[100],offset[100],chrno[100];
    ifstream inFile (argv[2]);
    char filename[500]="\0";
    char outfile[500]="\0";
    strcpy(filename,argv[2]);
    strcat(filename, ".mutatedseq.fna");
    strcpy(outfile,argv[2]);
    strcat(outfile,".base");
    ofstream basefile (outfile);
    ofstream out(filename);

    if(!inFile.is_open())
    {
        manual();
        return -1;
    }
    while (inFile >> refid >> offset >>chrno ) 
    {
           if(strcmp(refid,"Acc_number")!=0)
           {
              int l = strlen(offset);
              char* off = new char[l+1];
              memcpy(off,offset+1,l-2);
              off[l-2]='\0';
              if(getsequence(refid,atoi(off))!="")
              out<<">"<<refid<<"\n"<<getsequence(refid,atoi(off))<<endl;
              basefile << refid<<"\t"<<offset<<"\t"<<chrno<<"\t"<<getbase(refid,atoi(off))<<endl;
              delete[] off;
          }
    }
   inFile.close();
   cout<<"Elapsed:  "<< (clock() - start)*1000/CLOCKS_PER_SEC<<" ms\n\n";
return 0;
}    
   
