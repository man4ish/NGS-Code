#include <iostream>
#include <map>
#include <vector>
#include <time.h>
#include <string.h>
#include <fstream>
using namespace std;
typedef map<string,string> M;
typedef vector<char*> V;
M m,qmap;

void manual()
{
    cout<<"\n\nSXSNPparser\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXSNPParser <contaminated_readfile> <Quality file> <Row_Format>\n\n";
}

bool Ismutation(char* s)
{
     while(*s)
     {
          if(*s=='.')
              return 1;
          *s++;
     }
     return 0;
}

V splitv(char* s,const char* token)
{
         char* pch;
         vector<char*> v;
         pch = strtok(s,token);
         while(pch!= NULL)
         {
               v.push_back(pch);
               pch = strtok (NULL,token);
         }
         return v;
}

unsigned int getAVQS(char qscore[])
{
             int tsq = 0;
             int qcount = 0;
             char * pch;
             pch = strtok (qscore," ");
             while (pch != NULL)
             {
                    tsq  +=  atoi(pch);
                    qcount++ ;
                   pch = strtok (NULL, " ");
             }
             unsigned int avsq  = tsq/qcount;
             return avsq;
}

string getid(char* s,const char* token)
{
            char* pch;
            vector<string> v;
            pch = strtok(s,token);
            while(pch!= NULL)
            {
                  v.push_back(pch);
                 pch = strtok (NULL,token);
            }
       return ">"+ v[3];
}


bool match (char* header)
{
     while(*header)
     {
         if(*header == 'a' && *(header+1)=='l' && *(header+2)=='i'&& *(header+3) == 'g' && *(header+4)=='n' && *(header+5)=='m' )
            return 1;
            *header++;
     }
   return 0;
}

bool Isinsert(char* ins)
{
     while (*ins)
     {
            if( *ins== '.' && *(ins+2)== '-')
                return 1;
            else if(*ins == '.' && *(ins+3)== '-')
                return 1; 
            *ins++;
     }
     return 0;
}

bool Isdelete(char* del)
{
     while (*del)
     {
            if( *del== '.' && *(del+1)== '-')
                return 1;
            *del++;
     }
     return 0;
}


int getindex(char* st,int ts,int qs )
{
    int l = strlen(st);
    string s,s1; 
    for(int i =0; i<l;i++)
    {
        if(st[i]!='.')
        {
           if ( st[i]!= ' ')
                s += st[i];
         }
         else
           break;
     }

     int y = atoi(s.c_str());
     int diff = y- ts;
     int index = diff + qs;
     return index;
}

int main(int argc, char* argv[])
{
       if(argc != 4)
       {
          manual();
          exit(1);
       }

    clock_t start = clock();
    char qline[1024] = "\0";
    char qualline[1024] = "\0";
    char snpline[1024] = "\0";
    char qHeader[1024] = "\0";
    char qHdr[1024] = "\0";
    char qSequence[655366] = "\0";
    char qScore[655366] = "\0"; 

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
            }
            else
            {
                m[qHeader]= qSequence;
                qSequence[0] = '\0';
                qline[strlen(qline) - 1] = '\0';
                strcpy(qHeader, qline);
            }
        }
        else if(match(qline)==1) 
        { 
             qline[strlen(qline) - 1] = '\0';
             strcat(qSequence, qline);          
        }
    }
   fclose(queryFile);

 FILE * qualFile = fopen(argv[2], "r");
    if (qualFile == NULL)
    {
        manual();
        return -1;
    }

    while (!feof(qualFile))
    {
        fgets(qualline, 1024, qualFile);
        if ((qualline[0] == '>') || (feof(qualFile)))
        {
            if (qScore[0] == '\0')
            {
                qualline[strlen(qualline) - 1] = '\0';
                strcpy(qHdr, qualline);
            }
            else
            {
                qmap[qHdr]= qScore;
                qScore[0] = '\0';
                qualline[strlen(qualline) - 1] = '\0';
                strcpy(qHdr, qualline);
            }
        }
        else 
        {
             qualline[strlen(qualline) - 1] = '\0';
             strcat(qScore, qualline);
        }
    }
    fclose(qualFile);

    char flname[500];
    strcpy(flname,argv[3]);
    strcat(flname,".elim_contmin_reads");
    ofstream out(flname);

  FILE * snpFile = fopen(argv[3], "r");
    if (snpFile == NULL)
    {
        manual();
        return -1;
    }

    while (!feof(snpFile))
    {
            fgets ( snpline, 1024, snpFile);
                    if (snpline[0] == '>' && (snpline != '\0' )) 
                    { 
                        snpline[strlen(snpline)-1]= '\0';
                        V rowline  = splitv(snpline,"\t");
                        if(rowline.size()==14|| rowline.size()==13) 
                        {
                            string a;
                            a +=  ">";
                            a += rowline[3];           
                            if ( (m.count(a) == 0)||((m[a]).at(0)== ' ' && m[a].at(1)=='0' ))
                            {  
                                 
                                 if ( qmap.count(a) > 0 )
                                 {                                      
                                      string ds = qmap[a];
                                      string df = ds; 
                                      V qvec = splitv(&ds[0]," ");
                                      unsigned int avqs = getAVQS(&df[0]);
                                      int t_start = atoi(rowline[1]);
                                      int q1 = atoi(rowline[4]);
                                      int q2 = atoi(rowline[5]);
                                      int q_start;
                                      if(q1 < q2)
                                         q_start = q1;
                                      else 
                                         q_start = q2;
                                      float qvalue,qval;
                                      int vc = 0; 
                                      for ( ;vc < 4;vc++)
                                      {
                                           out<<rowline[vc]<<"\t";
                                      }
                                       out<<avqs<<"\t";
                                      int ins_count = 0;
                                      int del_count = 0;
                                      for ( ;vc<11;vc++)
                                       {
                                            out<<rowline[vc]<<"\t";
                                       }
                                      bool flag = false; 
                                      if ( Ismutation(rowline[11])==1)
                                      {
                                           vc++;
                                           V vb = splitv(rowline[11],",");
 
                                           
                                           for (int cc =0 ; cc < vb.size()-1; cc++)
                                           {     
                                                                                                                    
                                                int in =  getindex(vb[cc],t_start,q_start);
                                               
                                                if ( Isinsert(vb[cc])==1)
                                                {
                                                     ins_count++;
                                                
                                                     qvalue =  atoi(qvec[in + ins_count- del_count -1]);             
                                                     out<<vb[cc]<<qvalue<<",";  
                                                }
                                                else if(Isdelete(vb[cc])==1)
                                                {  
                                                             del_count++;
                                                     qval = (atoi(qvec[in+ins_count-del_count]) + atoi(qvec[in+ins_count-del_count-2]))/2;
                                                     out<<vb[vb.size()-1]<<qval;
                                                }                                                                                                                            else  
                                                {                
                                                     qvalue = atoi(qvec[in + ins_count -del_count -1]);
                                                     out<<vb[cc]<<qvalue<<",";   
                                                }  
                                            }                                                                                                                            int ind =  getindex(vb[vb.size()-1],t_start,q_start);

                                        
                                            if ( Isinsert(vb[vb.size()-1])==1)
                                            {
                                                 ins_count++;
                                                 qval = atoi(qvec[ind+ins_count-del_count-1]);
                                                 out<<vb[vb.size()-1]<<qval;
                                            }
                                            else if ( Isdelete(vb[vb.size()-1])==1)
                                            { 
                                                    del_count++;
                                                     qval = (atoi(qvec[ind+ins_count-del_count]) + atoi(qvec[ind+ins_count-del_count-2]))/2;
                                                     out<<vb[vb.size()-1]<<qval;
                                            }   
                                            else
                                            {
                                              out<<vb[vb.size()-1]<<qvec[ind+ins_count-del_count-1];
                                            }
                                           out<<"\t";
                                       }
                                       for ( ;vc<rowline.size()-1;vc++)
                                       {
                                             out<<rowline[vc]<<"\t";
                                       }
                                      out<<rowline[rowline.size()-1]<<endl;  
                                                    
                                  }          
                            }
                  } 
        }
  }
    fclose(snpFile);
    cout<< "Elapsed: "<<(clock() - start)/CLOCKS_PER_SEC<<" "<<"seconds\n"; 
}
