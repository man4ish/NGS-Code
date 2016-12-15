#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <iostream>
#include <boost/interprocess/containers/vector.hpp>
#include <algorithm>
#include <fstream>
#include <sys/types.h>
#include <regex.h>
using namespace boost::interprocess;
typedef managed_shared_memory::segment_manager                       segment_manager_t;
typedef allocator<void, segment_manager_t>                           void_allocator;
typedef allocator<int, segment_manager_t>                            int_allocator;
typedef vector<int, int_allocator>                                   int_vector;
typedef allocator<int_vector, segment_manager_t>                     int_vector_allocator;
typedef vector<int_vector, int_vector_allocator>                     int_vector_vector;
typedef allocator<char, segment_manager_t>                           char_allocator;
typedef basic_string<char, std::char_traits<char>, char_allocator>   char_string;
typedef std::pair<const char_string, char_string> map_value_type;
typedef std::pair<char_string, char_string> movable_to_map_value_type;
typedef allocator<map_value_type, segment_manager_t>  map_value_type_allocator;
typedef multimap< char_string, char_string,std::less<char_string>, map_value_type_allocator> complex_map_type;
string patt = "abcdefghijklmnopqrstuvwxyz0123456789";
char dirpath[500];

vector <const char* > Records;

  vector<const char*> DisplayRecords()
  {
    return Records;
  }

  vector<std::string> split(std::string strr, const char* token)
  {
       char str[300];
       strcpy(str,strr.c_str());
       char* pch;
       vector<std::string> vec;
       pch = strtok (str,token);
       while (pch != NULL)
       {
          vec.push_back(pch);  
          pch = strtok (NULL, token);
       }     
       return vec;
  }

  void Setfilepath(const char* filepath)
  {
      strcpy(dirpath,filepath);                
  }
 
  const char* Getfilepath()
  {
     return dirpath;
  }

  const char* GetRecord(std::string filename)
  {
     char flname[1000];
     strcpy(flname,filename.c_str());
     vector<std::string> vecb = split(flname,":");
     char str[300];
     strcpy(str,Getfilepath()); 
     strcat(str,vecb[0].c_str());
     strcat(str,".txt");
     std::ifstream in(str);
     char line[1000];
     in.seekg((size_t)atol((vecb[1]).c_str()), std::ios::beg);
     in.getline(line, 1000);
     return line;
  }
 
  std::string toLowerCase(std::string str)
  {
       for(int i=0;i<strlen(str.c_str());i++)
           if(str[i] >= 0x41 && str[i] <= 0x5A)
              str[i] = str[i] + 0x20;
          return str;
  }

  bool match(const char *string, char *pattern)
  {
        int status;
        regex_t re;
        if(regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0)
        {
           return 0;
        }
        status = regexec(&re, string, (size_t)0, NULL, 0);
        regfree(&re);
        if(status != 0)
        {
           return 0;
        }
          return 1;
  }

  bool matchrecord(std::string record,std::string pattern)
  {         
       vector<std::string> keywd = split(record,"\t");
       vector<std::string> keywdid = split(keywd[keywd.size()-1],",");
       char Quermatch[500];
       strcpy(Quermatch,(toLowerCase(pattern)).c_str()); 
       for (int rec =0; rec < keywdid.size(); rec++)
       {
            std::string kwwd = keywdid[rec]; 
            if(match(kwwd.c_str(), Quermatch) == 1)
            {
               return true;
            }
       } 
       return false;
  }

  vector <string> splitv(std::string inrec, const char* token)
  {
         char line[1000];
         strcpy(line, inrec.c_str());
         char* pch; 
         pch = strtok(line,token);
         vector<string > splitedvec;
         while(pch != NULL)
         {
           splitedvec.push_back(pch);
           pch = strtok(line,token);
         }
         return splitedvec;
    }

   void GetPatternrecord(std::string Query, std::string ori_query, const char* MySharedMemory)
   {
        std::string pattern = toLowerCase(Query.substr(0,3));
        using namespace boost::interprocess;
        try
        {
            managed_shared_memory segment(open_only,MySharedMemory);
            complex_map_type *MyMap = segment.find<complex_map_type>("PatternMap").first; 
            void_allocator alloc_inst (segment.get_segment_manager());
            std::pair<complex_map_type::iterator, complex_map_type::iterator> ppp;               
            char_string cs(pattern.c_str(), alloc_inst);
            if(MyMap->count(cs) != 0)
            {
               ppp = MyMap->equal_range(cs);
               for (complex_map_type::iterator it2 = ppp.first; it2 != ppp.second; it2++)
               {
                   std::string y(((*it2).second).c_str());
                   char line[500]; 
                   strcpy(line, GetRecord(y));
                   if(matchrecord((std::string)(line),ori_query)== 1) 
                   {
                       Records.push_back(line);
                   }
               }
            } 
         }
         catch(...)
         {
            shared_memory_object::remove(MySharedMemory);
            throw;
         }         
     }

   void Addcharright1(std::string query, std::string ori_query,const char* MySharedMemory  )
   {
        for (int i = 0; i < patt.size(); i++)
        {
             std::string outquery ;
	     outquery += patt[i];
	     outquery += query;
             GetPatternrecord(outquery, ori_query, MySharedMemory);
       	}
    }

    void Addchar2(std::string query, std::string ori_query, const char* MySharedMemory)
    {
		for (int i = 0; i < patt.size(); i++)
		{
		    for (int j =0;j < patt.size(); j++)
                     {
 	     		std::string outquery = query;
		        outquery += patt[i];
		        outquery += patt[j];
                        GetPatternrecord(outquery, ori_query, MySharedMemory);
                     }
		}
     }

     void Addchar1(std::string query, std::string ori_query, const char* MySharedMemory)
     {
          for (int ii =0;ii < patt.size(); ii++)
          {
               std::string outquery = query;
	       outquery += patt[ii];
               GetPatternrecord(outquery, ori_query, MySharedMemory);
	  }
     }

     void Getwildcardrecord(std::string query,std::string Query,const char* MySharedMemory)
     {           
           int len = query.size();
           if(query.find("?") != -1)
           {
		int pos = query.find("?");
                if(len == 3)
                {
                   if(pos ==1)
                   {
                      std::string Query1;
                      Query1 += query[0];
                      Query1 += query[2];
                      Addchar1(Query1,Query,MySharedMemory);
                      std::string Query2;
                      Query2 += query[2];
                      Addchar2(Query2, Query, MySharedMemory);
                   }
                   else if(pos == 2)
                   {
                       std::string Query1;
                       Query1 += query.substr(0,2);
                       Query1 += query[1];
                       GetPatternrecord(Query1, Query, MySharedMemory);
                       std::string Query2;
                       Query2 +=  query.substr(0,2);
                       Addchar1(Query2, Query, MySharedMemory);
	           }
		}
		else if(len > 3)
		{
                       if(pos == 1)
                       {
			      std::string Query1;
                              Query1 += query[0];
	                      Query1 += query.substr(2,2);
                              GetPatternrecord(Query1, query, MySharedMemory);
		              std::string Query2 = query.substr(2,2);
		              Addcharright1(Query2, Query, MySharedMemory);
			}
			else if(pos == 2)
			{
                             std::string Query1 = query.substr(0,3);
                             GetPatternrecord(Query1, query,MySharedMemory);
                             std::string Query2;
                             Query2 += query[0];
                             Query2 += query.substr(0,2);
                             Addchar1(Query2, Query, MySharedMemory);
			 }
			 else if(pos ==3)
			 {
                             std::string Query1 = query.substr(0,3);
                             std::string Query2 = query.substr(0,2);
                             Addchar1(Query2, Query, MySharedMemory);
                        }
	         }
    }
    else if(query.find("*") != -1)
    {
        std::string queryh = query;
        int aspos = query.find("*");
        if(len == aspos + 1)
        {
           query = query.substr(aspos-2, 3);
        }
        else
        {
            query = query.substr(aspos-1,3);
        }
        int apos = query.find("*");
            if(apos == 1)
	    {
                 std::string Query1;
                 Query1 += query[0];
                 Query1 += query[0];
                 Query1 += query[1];
                 GetPatternrecord(Query1, Query, MySharedMemory);
                 std::string Query2;
                 Query2 += query[0];
		 Query2 += query[0];
                 Query2 += query[0];
                 GetPatternrecord(Query2, Query, MySharedMemory);
                 std::string Query3;
                 Query3 += query[2];
                 Addchar2(Query3, Query, MySharedMemory);
             }
	     else if(apos == 2)
	     {
                  std::string Query1;
                  Query1 += query[0];
                  Addchar2(Query1, Query,MySharedMemory);
                  std::string Query2;
                  Query2 += query.substr(0,2);
                  Addchar1(Query2, Query, MySharedMemory);
                  std::string Query3;
                  Query3 += query.substr(0,2);
                  Query3 += query[2];
                  GetPatternrecord(Query3, Query, MySharedMemory);
              }
           }
}

void PartialSearch(std::string str, const char* MySharedMemory)
{
    if(str.find("*")!= -1 || str.find("?")!= -1)
         Getwildcardrecord((str).substr(0,3),str,MySharedMemory);    // take always * and ?
    else 
         GetPatternrecord((str).substr(0,3),str,MySharedMemory);        
}
// Need to implement logical operations
// Need to solve substr pblm
int main (int argc, char* argv[])
{  
    const char* MySharedMemory= "PatternIndex"; 
    Setfilepath("/archive/mgrc_testbed/working_directory/manish/MGBSearchEngine/programs/planning/hs_mgb_1.0.ref_genome.37.1.snp.chr");
    PartialSearch(argv[1],MySharedMemory); 
    vector<const char* > vec = DisplayRecords() ;     
    for (int i =0; i < vec.size(); i++ )
      std::cout<< vec[i]<<std::endl;
    return 0; 
}
