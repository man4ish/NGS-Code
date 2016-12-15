#include <iostream>
#include <map>
#include <fstream>
#include <string>
using namespace std;
typedef multimap<string , size_t> CDSmap;
typedef multimap <string,string> AssociationMap;
AssociationMap AssMap;
CDSmap cdsindexmap;
void GetRecord( string keyword ,const char* filename)
{
   pair<multimap<string, size_t>::iterator, multimap<string, size_t>::iterator> ppp;
   ppp = cdsindexmap.equal_range(keyword);
   ifstream in(filename, ios::in | ios::binary); 
   for (CDSmap::iterator it2 = ppp.first;it2 != ppp.second; ++it2)
   { 
        in.seekg(it2->second, ios::beg);
        char ch;
        while(in.get(ch))
        {
              if(ch== '\n')
              break;
              cout << ch;
         }
     cout<<"\n";
   }
   in.close(); 
}

int main(int argc , char* argv[])
{
   if(argc != 6)
   {
      cout<<"Check the parameters\n";
      exit(1);
   }
   bool flag;
   char flname[400],secfile[400],indexfile[400];
   strcpy(flname,"../Species/");
   strcat(flname,argv[1]);
   strcat(flname,"/");
   strcat(flname,argv[2]);
   strcat(flname,"/KeyWord/");

   if(strcmp(argv[4],"Gene_Name")== 0||strcmp(argv[4],"Gene_ID") == 0||strcmp(argv[4],"GI")==0) 
   {  
       strcpy(secfile,flname); 
       strcat(secfile,argv[4]);
       strcat(secfile,"_");
       strcat(secfile,argv[3]);
       strcat(secfile,".txt");
       flag = 1;
       ifstream Assfile(secfile);
       if(!Assfile.is_open())
       {
          cout<<"NO Association file\n";
          exit(1);
       }
       string Geneinfo, ID;
       while(Assfile>>Geneinfo>>ID) 
       {
             AssMap.insert(pair<string, string>(Geneinfo,ID)); 
       }
   }
   strcat(flname,argv[3]);
   strcat(flname,".txt");
   strcpy(indexfile,flname);
   strcat(indexfile,".index");
      ifstream infile(indexfile);
      if(!infile.is_open())
      {
          cout<<"cannot open the input file\n";
          exit(1);
      }
      char id[100];
      size_t index;
      while(infile >> id >>index)
      {
           cdsindexmap.insert(pair<string, size_t>(id, index)); 
      }
      if(AssMap.size()== 0)
      {
         GetRecord(argv[5],flname); 
      }
      else 
      {
            pair<multimap<string, string>::iterator, multimap<string, string>::iterator> ppp;
            ppp =  AssMap.equal_range(argv[5]);
            for ( AssociationMap::iterator it2 = ppp.first;it2 != ppp.second; ++it2)
             {
                  GetRecord(it2->second,flname);                      
             }
      }
     return 0;
}

