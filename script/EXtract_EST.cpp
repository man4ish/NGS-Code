#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
vector <string > vec;

vector <size_t > veccount; 


string getseq(char* seq)
{
      string s ;
      while(*seq)
      {
          if(isalpha(*seq))
          s += *seq;
           *seq++;
      }
      s += "\n";
      return s;
}

string getid(char* l)
{
   vector<string> vec;
   char *pch;
   pch = strtok(l, " ");
   while(pch != NULL)
   {
      vec.push_back(pch);
      pch = strtok(NULL," ");
   }
  return vec[2];
}
char *trim(char *s)
{

     size_t len = strlen(s);
     while (len > 0 && isspace(s[len - 1]))
     s[--len] = '\0'; return s;
}

void getdata(int a , int b)
{
  string s ;
   bool seqflag = 0,entryflag= 0; 
  for (int i = a ; i < b-1 ; i++ )
  {

   if(strstr(vec[i].c_str(),"SEQUENCE"))
      {
         seqflag =1;
         entryflag = 0;
      }

      if(strstr(vec[i].c_str(),"dbEST Id:"))
      {
         s += ">dbEST Id: ";
         s += getid(&(vec[i])[0]);
      }

      if(strstr(vec[i].c_str(),"EST name:"))
      {
         s += ", EST name: ";
         s += getid(&(vec[i])[0]);
         s += ", Organaism Name: ";
         s += "Elaeis guineensis";
      }

      if(strstr(vec[i].c_str(),"GenBank gi:"))
      {
         s += ", gi: ";
         s += getid(&(vec[i])[0]);
         s += "\n";
      }
      if(strstr(vec[i].c_str(),"Entry Created:"))
           entryflag =1;
      if(seqflag == 1 && entryflag == 0 )
      {
         if(strlen(vec[i].c_str()) != 0 && !strstr(vec[i].c_str(), "SEQUENCE") && !strstr(vec[i].c_str(), "Quality") )
         {
            s += getseq(&(vec[i])[0]);
         }
      }
      char ser[200];
      strcpy(ser,"Organism:       ");
      strcat(ser,"Elaeis guineensis");
      if(strstr(vec[i].c_str(),ser))
      {
         cout<<s<<"\n";
         s = "";
      }  
   }
}

int main(int argc, char* argv[])
{
   int vecc = 0;
   ifstream in(argv[1]);
   char line[1000];
   while(in)
   {
      
      in.getline(line,1000);
      vecc++;
      if(strstr(line,"IDENTIFIERS"))
       veccount.push_back(vecc);   
      vec.push_back(line);
   }
   
    veccount.push_back(vecc);
   for (int y = 0; y < veccount.size()-1; y++)
    {
    //    cout<<veccount[y]<<"\t"<<veccount[y+1]<<endl;
       getdata( veccount[y],veccount[y+1]);
    }
//   cout<<vec.size()<<endl; 
    vec.clear();
    veccount.clear();   
   return 0;
}
