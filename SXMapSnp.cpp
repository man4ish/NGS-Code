#include <iostream>
#include <map>
#include <vector>
#include <fstream>
using namespace std;
map<string,string> ChrtoGi;
map <string,string> snp; 

void manual()
{
    cout<< "\n\n SXMapSnps\n";
    cout<< "~~~~~~~~~~~~~~~~~~\n\n";
    cout<< "Usage:\n";
    cout<< "\tSXMapSnps <Chromosome_list> <translated_Snp_file> <snpfile>\n\n";
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

void SetChrtoGi(string li) {
   if(li != ""){
   char* f = &li[0],*pch;
   vector <string> sdf;
   string s;
   pch = strtok(f,"\t");
   while(pch != NULL) {
         sdf.push_back(pch);
         pch = strtok(NULL,"\t"); 
       }
   ChrtoGi[sdf[0]]= sdf[1];
   }
}

int main(int argc, char* argv[])
{
 string line3;
 ifstream myfile3 (argv[1]);
  if (myfile3.is_open()) {
    while (! myfile3.eof() ) {
      getline (myfile3,line3);
       {
        SetChrtoGi(line3);
       }
    }
    myfile3.close();
  }
  else cout << "Unable to open chromosome file\n";

 ifstream snpfile(argv[2]);
 if(!snpfile.is_open())
    {
        manual();
        return -1;
    }
char refid[100],offset[100],chrno[100],base[10],Genome_location[100];

    while (snpfile >> refid >> offset >>chrno>>Genome_location) 
    {
           if(strcmp(refid,"Acc_number")!=0)
           {
               string s;
               s += chrno;
               s += "-";
               s += Genome_location;
               snp[s] = s;   
                                        
           }
    }

char flname[1000];
strcpy(flname,argv[3]);
strcat(flname,"_mapped");
ofstream outfile(flname);

    FILE * snpFile = fopen(argv[3], "r");
    if (snpFile == NULL)
    {
        manual();
        return -1;
    }
    char snpline[1000000] ="\0";
    while (!feof(snpFile))
    {
        fgets(snpline, 1024, snpFile);
        if ((snpline[0] == '>') || (snpFile))
        {
           snpline[strlen(snpline)-1]=='\0';
           string snpread = snpline;
           vector<string> vec;  
           char* pch;
           pch = strtok (&snpread[0],"\t");
           while (pch != NULL)
           {
                  vec.push_back(pch);
                  pch = strtok (NULL, "\t");
           }
         
           string a; 
           a += ChrtoGi[vec[0]];
           a += "-";
           a += vec[1]; 
           cout<<a<<endl; 
            vec.clear();
           if(snp.count(a)!= 0)
            outfile<<snpline;
               
        }
    }
    fclose(snpFile);
return 0;
}
