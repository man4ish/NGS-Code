#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <iostream>
#include <vector>
#include <map>
#include <boost/interprocess/containers/vector.hpp>
#include <algorithm>
using namespace boost::interprocess;
typedef managed_shared_memory::segment_manager                       segment_manager_t;
typedef allocator<void, segment_manager_t>                           void_allocator;
typedef allocator<int, segment_manager_t>                            int_allocator;
typedef vector<int, int_allocator>                                   int_vector;
typedef allocator<int_vector, segment_manager_t>                     int_vector_allocator;
typedef vector<int_vector, int_vector_allocator>                     int_vector_vector;
typedef allocator<char, segment_manager_t>                           char_allocator;
typedef basic_string<char, std::char_traits<char>, char_allocator>   char_string;

class complex_data
{
   int id_;
   char_string char_string_;
   int_vector_vector int_vector_vector_;
   public:
   complex_data(int id, const char *name, const void_allocator &void_alloc): id_(id), char_string_(name, void_alloc), int_vector_vector_(void_alloc)
   {}
};

typedef std::pair<const char_string, int> map_value_type;
typedef std::pair<char_string, int> movable_to_map_value_type;
typedef allocator<map_value_type, segment_manager_t>  map_value_type_allocator;
typedef map< char_string, int,std::less<char_string>, map_value_type_allocator> complex_map_type;
void DeleteKeywdmap(const char* sharedmemory)
{
  using namespace boost::interprocess;
  try
  {
      managed_shared_memory segment(open_only  ,sharedmemory);
      std::cout << "yes"<< "\n";
      typedef map< char_string, int, std::less<char_string>, map_value_type_allocator> complex_map_type;
           std::cout << "yes"<< "\n"; 
      complex_map_type *MyMap = segment.find<complex_map_type>("MyMap").first;
            std::cout << "yes"<< "\n"; 
      std::cout<< MyMap->size() << "\n";
      segment.destroy<complex_map_type>("MyMap");
   }
   catch(...)
   {
      shared_memory_object::remove(sharedmemory);
      throw;
   }
   shared_memory_object::remove(sharedmemory); 
}

void DeletePatternmap(const char* sharedmemory)
{
  using namespace boost::interprocess;
  try
  {
      managed_shared_memory segment(open_only,sharedmemory);
      typedef map< char_string, int, std::less<char_string>, map_value_type_allocator> complex_map_type;
      complex_map_type *MyMap = segment.find<complex_map_type>("MyMap").first;
      segment.destroy<complex_map_type>("MyMap");
   }
   catch(...)
   {
      shared_memory_object::remove(sharedmemory);
      throw;
   }
   shared_memory_object::remove(sharedmemory);
}

void Updateindexfile(char* params[])
{
       int i = 1;
       bool pflag = 0;
       if(strstr(params[i],"-p")) 
       {
          i++;
          pflag = 1;
       }                              
       char MemoryIndexFile[400];
       char featurename[100];
       char partialsharedregion[500];
       strcpy(MemoryIndexFile,params[i]);
       i++;
       strcpy(featurename,params[i]);                          
       std::map<std::string,std::string> m1;
       vector <std::string> vec;
       std::ifstream in(MemoryIndexFile);
 
       char line[1000];
       while(in)
       {
          in.getline(line,1000);
          if(in)
          {
              char* pch; 
              pch = strtok(line,"\t");
              while(pch != NULL)
              {
                 vec.push_back(pch);
                 pch = strtok(NULL,"\t");
              }     
             m1[vec[0]] = vec[1];
          }
       } 
       in.close();
       strcpy(partialsharedregion,(m1[featurename]).c_str());
       strcat(partialsharedregion, ".pi");
       DeleteKeywdmap((m1[featurename]).c_str());
      // if(pflag == 1)
          DeletePatternmap(partialsharedregion);
       std::ofstream out(MemoryIndexFile);        
       for( std::map<std::string,std::string>::iterator ii= m1.begin(); ii != m1.end(); ii++)
       {
           if(strstr((ii->first).c_str(), featurename))
            out<< ii->first <<"\t"<< ii->second;
       } 
       out.close();
}


int main (int argc, char* argv[])
{
   DeleteKeywdmap("hs_mgb_1.0.ref_genome.37.1.gene");
   std::cout<< "done\n";
   DeletePatternmap("hs_mgb_1.0.ref_genome.37.1.gene.pi");
   std::cout<< "done\n";
  // Updateindexfile(argv);
   return 0;
}

