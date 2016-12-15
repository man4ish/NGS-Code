#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <iostream>
#include <vector>
#include <fstream>
using namespace boost::interprocess;
typedef managed_shared_memory::segment_manager segment_manager_t;
typedef allocator<void, segment_manager_t> void_allocator;
typedef allocator<char, segment_manager_t> char_allocator;
typedef basic_string<char, std::char_traits<char>, char_allocator> char_string;
typedef std::pair<const char_string, char_string> map_value_type;
typedef allocator<map_value_type, segment_manager_t>  map_value_type_allocator;
typedef multimap< char_string, char_string, std::less<char_string>, map_value_type_allocator>  complex_map_type;

std::string trim(const char *str)
{
   std::string ptr;
   int i,j=0;
   for(i=0;str[i]!='\0';i++)
   {
      if (str[i] != ' ' && str[i] != '\t')
         ptr += str[i];
   }
   return ptr;
} 

const char* GetFeatureName (char * filename)
{
     char* pch;
     pch = strtok(filename,".");
     std::vector <const char* > vec;
     while(pch != NULL)
     {
        vec.push_back(pch);
        pch = strtok(NULL,".");
     }
     return vec[vec.size()-3];
}

const char* Getsharedmemoryregion(char* filename)
{
       int pos = ((std::string) filename).find("chr");
       std::string region_name =  ((std::string) filename).substr(0,pos-1); 
       return region_name.c_str();
 }

 const char* GetPartialsharedregion(char* filename)
 {
       int pos = ((std::string) filename).find("chr");
       std::string region_name =  ((std::string) filename).substr(0,pos-1);         
       return (region_name + ".pi").c_str();
 }

 const char* Getchrnum(char* str)
 {
   int j = ((std::string)str).find("chr");
   int k = ((std::string)str).find("txt"); 
   return ((std::string)str).substr(j+3,k-1-j-3).c_str();   
 }

 std::string ConvertTolowercase(std::string str)
 {
      for (int i=0;i<strlen(str.c_str());i++)
      if (str[i] >= 0x41 && str[i] <= 0x5A)
      str[i] = str[i] + 0x20;
      return str;
 }

 std::vector <std::string>split (std::string strr, const char* token)
  {
       char str[300];
       strcpy(str,strr.c_str());
       char* pch;
       std::vector<std::string> vec;
       pch = strtok (str,token);
       while (pch != NULL)
       {
          vec.push_back(pch);
          pch = strtok (NULL, token);
       }
       return vec;
  }

 void UpdateIndexMap(std::string str, std::string index, const char* SharedMemory)
 {
    managed_shared_memory segment(open_only,SharedMemory);
    void_allocator alloc_inst (segment.get_segment_manager());
    complex_map_type *mymap = segment.find<complex_map_type>("MyMap").first;
    std::string h = ConvertTolowercase(str);
    char_string patternvalue(h.c_str(), alloc_inst);
    char_string indexvalue((index).c_str(), alloc_inst);
    std::cout<< patternvalue<<"\t"<<indexvalue<<std::endl;    
    mymap->insert(std::pair<char_string, char_string>(patternvalue,indexvalue));
 }

void UpdatePatternMap(std::string str, std::string index, const char* SharedMemory)
 {
    std::cout<< "yes\n";
    std::cout << SharedMemory << "\n";
    managed_shared_memory segment(open_only,SharedMemory);
     std::cout<< "yes\n";
    void_allocator alloc_inst (segment.get_segment_manager());
    complex_map_type *mymap = segment.find<complex_map_type>("MyMap").first; 
         std::string strr = str;
         if(strr.size() > 2 )
         {
                for (int j = 0; j< strr.size() -2;j++)
                {
                     std::string pattern = strr.substr(j,3);
                     if(strstr(pattern.c_str(),".") == NULL && strstr(pattern.c_str(),"-") == NULL && strstr(pattern.c_str(),"_") == NULL && pattern != "" ) 
                     {
                           std::string h = ConvertTolowercase(pattern);
                           char_string patternvalue(h.c_str(), alloc_inst);
                           char_string indexvalue((index).c_str(), alloc_inst);
                           std::cout<< patternvalue<<"\t"<<indexvalue<<"\n";
                           mymap->insert(std::pair<char_string, char_string>(patternvalue,indexvalue));
                     }
                 }
          }
 }

void IndexFile(char* params[])
{
   bool pflag = 0;
   int argcount = 1;
   if(strstr(params[argcount],"-p"))
   {
     argcount++;
     pflag = 1;
   }
      
   using namespace std;
   char filename[500];
   char inputfile[500];
   strcpy(inputfile, params[argcount]);
    ifstream in(inputfile);
   strcpy(filename,inputfile);

   char line[1000];
   char recline[100];
   int index = 0;
   bool patternflag = 0;
   bool indexflag = 0;
   while(in)
   {
     in.getline(line,1000);
     int len = strlen(line);
     strcpy(recline,line);
     if(in)
     {
        std::stringstream ss;
        ss<< index;
        char indexname[1000];
        strcpy(indexname, Getchrnum(inputfile));
        strcat(indexname, ":");
        char* pch;
        pch = strtok(line,"\t");
        vector <const char*> vec;
        while(pch != NULL)
        {
          vec.push_back(pch);
          pch = strtok(NULL,"\t");
        }
        std::string Keywordcol  = trim(vec[vec.size()-1]);
        std::vector <std::string> vecin = split(Keywordcol,",");
        strcat(indexname, (ss.str()).c_str());
        for (int i =0; i< vecin.size(); i++ )
        {
              UpdateIndexMap(ConvertTolowercase(vecin[i]), indexname, Getsharedmemoryregion(inputfile));
             // if(pflag == 1)
                 UpdatePatternMap(vecin[i],indexname,GetPartialsharedregion(inputfile));
        }
        index += len;
        index += 1;
     }
  }
    in.close();
}


int main (int argc, char* argv[])
{
   IndexFile(argv);
   return 0;
}
