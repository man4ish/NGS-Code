#include <iostream>
#include <SXSynabase.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

void manual()
{
    cout<<"\n\nSXGetSequence\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXGetSequence -gi  <gi>  -i <index> -r <range> -s <side>\n\n";
}

void Getsequence(char* ginum,char* index, char* _range, char* side)
{    
     char flname[500];
     strcpy(flname,ginum);
     strcat(flname,"_");
     strcat(flname,index);
     strcat(flname,"_"); 
     strcat(flname,_range);
     strcat(flname,".seq_side");
     strcat(flname,side);
     ofstream out(flname);
     t_sindex sequenceidx;
     const char * seq;
     int size;
     int range = atoi(_range); 
     if( strcmp(side,"1")==0)
     {
         size = atoi(index)-atoi(_range);
     }
     else 
     {
         size = atoi(index);
     }
 
     SXSynabase* sb = new SXSynabase( "hsr_gnom_ncbi-36.3_sb3.0.9", NXRead);
     if (sb->SXGetGiIdx(ginum, &sequenceidx))
     {
         seq = sb->SXGetSequence(sequenceidx);
     }
     out<<">gi|"<<ginum<<"|index|"<<index<<"|\tHomo Sapiens Chromosome 17 range: "<<_range<<" side: "<<side<<"\n";
     char* s = new char[range + 1];
     memcpy(s,seq+size,range);
     s[range] = '\0';
     for(int i=1;i<=strlen(s); i++)
     {
         out<<s[i-1];
         if ( i%70 == 0) 
         {
              out<<"\n";
         }
     }
     out<<"\n";
}

int main(int argc, char* argv[])
{
    if ( argc != 9)
    {
         manual();
         return -1;
    }
    if ( (strcmp(argv[1],"-gi")== 0) && (strcmp(argv[3],"-i")==0) && (strcmp(argv[5],"-r")==0) && (strcmp(argv[7],"-s")==0))
    {
         Getsequence(argv[2],argv[4],argv[6],argv[8]);
    }
    else 
    {
         manual();
         return -1;
    }     
    return 0;
}


