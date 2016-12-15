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
char filename[300];
char feature[100];
void SetFilename(const char* flname)
{
   strcpy(filename,flname);
}

const char* Getsharedmemoryregion()
{
       int pos = ((std::string) filename).find("chr");
       std::string region_name =  ((std::string) filename).substr(0,pos-1); 
       return region_name.c_str();
}

const char* GetPartialsharedregion()
{
       int pos = ((std::string) filename).find("chr");
       std::string region_name =  ((std::string) filename).substr(0,pos-1);         
       return (region_name + ".pi").c_str();
}

void SetFeaturename (char * filename)
{
     char* pch;
     pch = strtok(filename,".");
     std::vector <const char* > vec;
     while(pch != NULL)
     {
        vec.push_back(pch);
        pch = strtok(NULL,".");
     }
     strcpy(feature,vec[vec.size()-3]);
}

const char* GetFeatureName()
{
   return feature;
}

std::string GetIndexMapFile(std::string filename)
{
   int j =0;
   for (int i = filename.size()-1; i >= 0 ; i--  )
   {  if(filename[i] == '.')
       j++;
     if(j == 3)
       return filename.substr(0,i) + ".MapIndex.txt";
   }
}

std::string GetMapIndexFile(std::string filename)
{
   int j =0;
   for (int i = filename.size()-1; i >= 0 ; i--  )
   {  if(filename[i] == '.')
       j++;
     if(j == 3)
       return filename.substr(0,i) + ".MapIndex.txt";
   }
}

std::string GetMapIndexFile(std::string filename)
{
   int j =0;
   for (int i = filename.size()-1; i >= 0 ; i--  )
   {  if(filename[i] == '.')
       j++;
     if(j == 3)
       return filename.substr(0,i) + ".MapIndex.txt";
   }
}
void CreateIndexMap()
{
    shared_memory_object::remove(Getsharedmemoryregion());
    managed_shared_memory segment(create_only,Getsharedmemoryregion(), 10000000);
    void_allocator alloc_inst (segment.get_segment_manager());
    complex_map_type *mymap = segment.construct<complex_map_type>("MyMap")(std::less<char_string>(), alloc_inst);
    std::ofstream outfile("IndexMap.txt");
    outfile << GetFeatureName() << "\t" << Getsharedmemoryregion() << "\n";
    outfile.close();              
}

void CreatePatternMap()
{
    std::cout << GetPartialsharedregion()<<"\n";
    shared_memory_object::remove( GetPartialsharedregion());
    managed_shared_memory segment(create_only, GetPartialsharedregion(), 400000000);
    void_allocator alloc_inst (segment.get_segment_manager());
    complex_map_type *mymap = segment.construct<complex_map_type>("MyMap")(std::less<char_string>(), alloc_inst);
}

void SetParameter(char* params[],int paramcount)
{
    char filename[500];
    int i = 1;
    bool pflag = 0;
    if(strstr(params[i],"-p"))
    {
      pflag = 1;
      i++;
    }  
    SetFilename(params[i]);
    strcpy(filename,params[i]);
    SetFeaturename(filename);
  
    CreateIndexMap();
    CreatePatternMap(); 
    /*if(pflag == 1)
    {
     std::cout<< "yes\n";
     CreatePatternMap(); 
    }   */
}
int main (int argc, char* argv[])
{
    SetParameter(argv, argc);
    return 0;
}
