#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;
map<string,string> genekwds;


string GetKeyword(const char* genecol)
{
   string kwd;
   char line[1000];
   strcpy(line,genecol);
   vector <string> vec;
   char* pch;
   pch = strtok(line,", ");
   while(pch != NULL)
   {
       vec.push_back(pch);
       pch = strtok(NULL,", ");
   }
  
   for(int i =0; i< vec.size()-1; i++)
   {
           map<string,string>::const_iterator it = genekwds.find(vec[i]);
           if(it != genekwds.end())
           {
             kwd += it->second;
             kwd += ", ";
           }
   }
   map<string,string>::const_iterator it1 = genekwds.find(vec[vec.size()-1]);
           if(it1 != genekwds.end())
           {
             kwd +=  it1->second;
           }
   return kwd;
}
int main(int argc, char* argv[])
{
  ifstream in1(argv[1]);
  char line1[1000], lstline1[5000], lstline2[2000];  
  while(in1)
  {
     in1.getline(line1,1000);
     if(in1)
     {
       vector <const char*> vec;
       in1.getline(line1,1000);
       char* pch;
       pch = strtok(line1,"\t");
       while(pch != NULL)
       {
          vec.push_back(pch);
          pch = strtok(NULL,"\t");
       }
        genekwds[vec[0]] = vec[1];                 
     }   
  }
  in1.close();


// for (map <string,string>:: iterator it  = genekwds.begin(); it!= genekwds.end(); it++)
//  cout <<"******"<<it->first<<"\t"<<it->second << endl;

  ifstream in2(argv[2]);
  while(in2)
  {
     in2.getline(lstline1,5000);
    // cout << lstline1 << endl;
     if(in2)
     { 
        if( lstline1[0] != '\0' && !strstr(lstline1,"#"))
        {
        strcpy(lstline2,lstline1);
        vector <const char*> vec2;
        strcpy(lstline2,lstline1);
        char* pch2;
        pch2 = strtok(lstline1,"\t");
        while(pch2 != NULL)
        {
          vec2.push_back(pch2);
          pch2 = strtok(NULL,"\t");
        }
        if(strcmp(vec2[9],"-") != 0)
        {
           cout << vec2[0]<<"\t"<<vec2[1]<<"\t"<<vec2[2]<<"\t"<<vec2[3]<<"\t"<<vec2[4]<<"\t"<<vec2[5]<<"\t"<<vec2[6]<<"\t"<<vec2[7]<<"\t"<<vec2[8]<<"\t"<<vec2[9]<<"\t";
           if(GetKeyword(vec2[9]) ==  "")
           cout <<"-";
           else 
                      cout <<GetKeyword(vec2[9]);
           cout <<"\t"<<vec2[10]<<"\t"<<vec2[11]<<"\t"<<vec2[12]<<"\t"<<vec2[13]<<"\t"<<vec2[14]<<"\t"<<vec2[15]<<"\t"<<vec2[16]<<"\t"<<vec2[17]<<"\t"<<vec2[18]<<"\t"<<vec2[19]<<endl;
       }
        else
         cout << vec2[0]<<"\t"<<vec2[1]<<"\t"<<vec2[2]<<"\t"<<vec2[3]<<"\t"<<vec2[4]<<"\t"<<vec2[5]<<"\t"<<vec2[6]<<"\t"<<vec2[7]<<"\t"<<vec2[8]<<"\t"<<vec2[9]<<"\t-\t"<<vec2[10]<<"\t"<<vec2[11]<<"\t"<<vec2[12]<<"\t"<<vec2[13]<<"\t"<<vec2[14]<<"\t"<<vec2[15]<<"\t"<<vec2[16]<<"\t"<<vec2[17]<<"\t"<<vec2[18]<<"\t"<<vec2[19]<<endl; 
        }
     }
    
  }
  in2.close(); 
  return 0;
}
