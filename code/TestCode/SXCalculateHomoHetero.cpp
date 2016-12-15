#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
int homocount =0;
int heterocount =0;
void readsnpfile1(char* line1)
{
     if(line1[0] != '\0' && !strstr(line1,"#")) 
     {
        if(strstr(line1,"Chromosome"))
           return;
        else
        {
           char recline1[10000];
           strcpy(recline1,line1);
           char* pch;
           vector<string> sdf;
           string s;
           pch = strtok (recline1,"\t");
           while(pch != NULL) 
           {
              sdf.push_back(pch);
              pch = strtok (NULL,"\t");
           }
            int diff;
            float x;
            diff  =  atoi(sdf[8].c_str()) - atoi(sdf[5].c_str());
            int  snp =   atoi(sdf[5].c_str()); 
            if(diff  > snp )
            {
              int num = diff;
              x = (float)num/(float)snp;
             // cout << x <<"\t"<< snp<<"\t"<< num<<"\t"<< atoi(sdf[5].c_str())<< endl;
            }
            else 
              x= (float)snp/(float)diff;      
           
          // if(x >= 5.6  )
           if(x >= 6.14)
           {
             homocount++;
            // cout <<line1<<"\t"<< "homo" << endl;
           }
           else 
           { 
             heterocount++; 
           }
              //cout << line1 <<"\t"<<"hetero"<< endl;
          // else if(x >2.33 && x< 5.6)
            //          cout <<line1<<"\t"<<"fuzzy" << endl;
          // else if(x <= 2.33)
            //          cout <<line1<<"\t"<< "hetero"<< endl;

         }
     }
}

int main(int argc, char* argv[])
{
  ifstream in(argv[1]);
  char line[5000],recline[5000];
  while(in)
  {
     in.getline(line,5000);
     if(in)
     {
        readsnpfile1(line);
     }
  }
  in.close();
  cout << "homocount=\t"<< homocount<<endl;
    cout << "heterocount=\t"<< heterocount<<endl;
}

