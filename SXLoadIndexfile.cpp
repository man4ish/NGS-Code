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


int main (int argc, char* argv[])
{
   using namespace std;
   char MySharedMemory[400];
   strcpy(MySharedMemory, argv[1]);
   shared_memory_object::remove(MySharedMemory);
   managed_shared_memory segment(create_only,MySharedMemory, 1999999999 ); // 1999999999);
   void_allocator alloc_inst (segment.get_segment_manager());
   complex_map_type *mymap = segment.construct<complex_map_type>("MyMap")(std::less<char_string>(), alloc_inst);
   ifstream infile(argv[1]); 
   char line[1000];
   int i = 0;
   while(infile)
   {
     infile.getline(line,1000);
     if(infile)
     {
        if(infile)
        {
           char* pch; 
           pch = strtok(line,"\t");
           vector <const char*> vec;
           while(pch != NULL)
           {
               vec.push_back(pch);
               pch = strtok(NULL,"\t");
           }
           char_string cs(vec[0], alloc_inst);
           char_string indexvalue(vec[1], alloc_inst);
           if( i%10000 == 0 ) { 
             std::cout << "["  << i << "]" <<  cs <<"\t"<<indexvalue<< " ,";  //endl;
           }
           mymap->insert(std::pair<char_string, char_string>(cs,indexvalue));      
         }
      }
     i++;
   }
   infile.close();
   return 0;
}
