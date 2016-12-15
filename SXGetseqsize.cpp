#include <iostream>
#include <SXSynabase.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

void manual()
{
    cout<<"\n\nSXGetSeqsize\n";
    cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    cout<<"Usage:\n";
    cout<<"\tSXGetSeqsize <gi_file>\n\n";
}

unsigned int getseqsize(char* ginum)
{
         unsigned int seqlen;
         t_sindex sequenceidx;
         const char * seq;
         SXSynabase* sb = new SXSynabase( "hsr_gnom_ncbi-36.3_sb3.0.9", NXRead);
         if(sb->SXGetGiIdx(ginum, &sequenceidx)) 
         {     
            seq = sb->SXGetSequence(sequenceidx);   
         }

         seqlen = strlen(seq);
        return seqlen;
}



int main(int argc, char* argv[])
{
    ifstream in(argv[1]);
    char qline[100] = "\0";
   
    while (in)
    {
        in.getline(qline,100);
        if(in)
        cout<<qline<<"\t"<<getseqsize(qline)<<endl;   
    }
    in.close();
 return 0;
}


