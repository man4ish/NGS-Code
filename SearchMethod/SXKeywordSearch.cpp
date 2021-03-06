#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <iostream>
#include <boost/interprocess/containers/vector.hpp>
#include <algorithm>
#include <fstream>
#include <set>
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
char dirpath[500];
char SharedMemory[200];
char PartialSharedMemory[200];
std::set<std::string> recset;
string patt = "abcdefghijklmnopqrstuvwxyz0123456789";
vector <const char* > Records;

 vector<std::string> split(char* str, const char* token)
 {
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

 std::string ConvertTolowercase(std::string str)
 {
      for (int i=0;i<strlen(str.c_str());i++)
      if (str[i] >= 0x41 && str[i] <= 0x5A)
      str[i] = str[i] + 0x20;
      return str;
 }

 void Setfilepath(const char* filepath)
 {
     strcpy(dirpath,filepath);                
 } 

 const char* Getfilepath()
 {
     return dirpath;
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
    std::cout << line << "\n";
    return line;
 }

 const char* GetPartialSharedMemoryName()
 {
    return  PartialSharedMemory;
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

 void GetPatternrecord(std::string Query, std::string ori_query)
 {
        std::string pattern = ConvertTolowercase(Query.substr(0,3));
        using namespace boost::interprocess;
        try
        {
            managed_shared_memory segment(open_only,GetPartialSharedMemoryName());
            complex_map_type *MyMap = segment.find<complex_map_type>("MyMap").first; 
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
                   recset.insert(GetRecord(y));        
               }
            } 
         }
         catch(...)
         {
            shared_memory_object::remove(GetPartialSharedMemoryName());
            throw;
         }         
  }

 void Addcharright1(std::string query, std::string ori_query  )
 {
        for (int i = 0; i < patt.size(); i++)
        {
             std::string outquery ;
	     outquery += patt[i];
	     outquery += query;
             GetPatternrecord(outquery, ori_query);
       	}
  }

  void Addchar2(std::string query, std::string ori_query)
  {
		for (int i = 0; i < patt.size(); i++)
		{
		    for (int j =0;j < patt.size(); j++)
                     {
 	     		std::string outquery = query;
		        outquery += patt[i];
		        outquery += patt[j];
                        GetPatternrecord(outquery, ori_query);
                     }
		}
  }

  void Addchar1(std::string query, std::string ori_query)
  {
          for (int ii =0;ii < patt.size(); ii++)
          {
               std::string outquery = query;
	       outquery += patt[ii];
               GetPatternrecord(outquery, ori_query);
	  }
  }

  void Getwildcardrecord(std::string query,std::string Query)
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
                      Addchar1(Query1,Query);
                      std::string Query2;
                      Query2 += query[2];
                      Addchar2(Query2, Query);
                   }
                   else if(pos == 2)
                   {
                       std::string Query1;
                       Query1 += query.substr(0,2);
                       Query1 += query[1];
                       GetPatternrecord(Query1, Query);
                       std::string Query2;
                       Query2 +=  query.substr(0,2);
                       Addchar1(Query2, Query);
	           }
		}
		else if(len > 3)
		{
                       if(pos == 1)
                       {
			      std::string Query1;
                              Query1 += query[0];
	                      Query1 += query.substr(2,2);
                              GetPatternrecord(Query1, query);
		              std::string Query2 = query.substr(2,2);
		              Addcharright1(Query2, Query);
			}
			else if(pos == 2)
			{
                             std::string Query1 = query.substr(0,3);
                             GetPatternrecord(Query1, query);
                             std::string Query2;
                             Query2 += query[0];
                             Query2 += query.substr(0,2);
                             Addchar1(Query2, Query);
			 }
			 else if(pos ==3)
			 {
                             std::string Query1 = query.substr(0,3);
                             std::string Query2 = query.substr(0,2);
                             Addchar1(Query2, Query);
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
                 GetPatternrecord(Query1, Query);
                 std::string Query2;
                 Query2 += query[0];
		 Query2 += query[0];
                 Query2 += query[0];
                 GetPatternrecord(Query2, Query);
                 std::string Query3;
                 Query3 += query[2];
                 Addchar2(Query3, Query);
             }
	     else if(apos == 2)
	     {
                  std::string Query1;
                  Query1 += query[0];
                  Addchar2(Query1, Query);
                  std::string Query2;
                  Query2 += query.substr(0,2);
                  Addchar1(Query2, Query);
                  std::string Query3;
                  Query3 += query.substr(0,2);
                  Query3 += query[2];
                  GetPatternrecord(Query3, Query);
              }
           }
  }

  void PartialSearch(std::string str)
  {
  
         GetPatternrecord((str).substr(0,3),str);        
  }

  bool FilterAND(std::string Rec, const char* recvalue)
  {
      char Keywd[500];
      strcpy(Keywd,recvalue);
      char Recline[500];
      strcpy(Recline,Rec.c_str());
      std::set<std::string>::iterator theIterator;
      vector <std::string> Keywordcol = split(Recline,",");
      for(int i = 0; i < Keywordcol.size(); i++)
      {                 
                  if (match((Keywordcol[i]).c_str(),Keywd) == 0)
                  {
                      return true;
                  }
       }   
       return false;
  }

  bool FilterOR(std::string Rec, const char* recvalue)
  {
      char Keywd[500];
      char Keyrecord[500];
      strcpy(Keyrecord,Rec.c_str());
      strcpy(Keywd,recvalue);
      std::set<std::string>::iterator theIterator;
      vector <std::string> Keywordcol = split(Keyrecord,",");
      for(int i = 0; i < Keywordcol.size(); i++)
      {                 
                  if (match((Keywordcol[i]).c_str(),Keywd) == 1)
                  {
                      return true;
                  }
       }   
       return false;
  }  

  void SetSharedMemoryName(const char* shm)
  {
     strcpy(SharedMemory,shm);
     strcpy(PartialSharedMemory,shm);
     strcat(PartialSharedMemory,".pi");
  }

const char* GetSharedMemoryName()
{
  return  SharedMemory;
}



void SearchKeyword( const char* Keyword,bool pflag )
{   
   using namespace boost::interprocess;
   try
   {
      managed_shared_memory segment(open_only,GetSharedMemoryName());
      typedef map< char_string,char_string, std::less<char_string>, map_value_type_allocator> complex_map_type;
      complex_map_type *MyMap = segment.find<complex_map_type>("MyMap").first;
      void_allocator alloc_inst (segment.get_segment_manager());
      std::pair<complex_map_type::iterator, complex_map_type::iterator> ppp;
      char_string cs(ConvertTolowercase(Keyword).c_str(), alloc_inst);
      if(MyMap->count(cs) != 0)
      {
         ppp = MyMap->equal_range(cs);
         for (complex_map_type::iterator it2 = ppp.first; it2 != ppp.second; it2++)
         {
           std::string fileindex  ((*it2).second.c_str());
           recset.insert(GetRecord(fileindex));
           std::cout << fileindex <<"\n";
           std::cout << GetRecord(fileindex) << std::endl;
         }
      }
      else
      {
         if(pflag == 1)
         PartialSearch(ConvertTolowercase(Keyword));
      }
   }
   catch(...)
   {
      shared_memory_object::remove(GetSharedMemoryName());
      throw;
   }
}


std::set<std::string> FilterRecord(char* params[], int paramcount)
{
     bool Andflag = false;
     bool pflag = 0;
     int argcount = 1;
     if(strstr(params[argcount],"-p"))
     {
        pflag = 1;
        argcount++;
     }
     
     for(int i = argcount; i < paramcount; i++)
     {
        
          if(!strstr(params[i],"+"))
          {
             if(strstr(params[i],"*") || strstr(params[i],"?"))
                Getwildcardrecord(((std::string)params[i]).substr(0,3),params[i]);
             else
              SearchKeyword(params[i],pflag);
          }
         else
         {
           Andflag = true;
         }
      }

     std::set<std::string> filteredrec;
      std::set<std::string>::iterator theIterator; 
     for (int j =0; j < paramcount; j++)
     {
      for( theIterator = recset.begin(); theIterator != recset.end(); theIterator++ )
      {
              if(Andflag == true)
              {
                 if(FilterAND(*theIterator,params[j]) == true)
                  filteredrec.insert(*theIterator);
              }
              else
              {
                 if(FilterOR(*theIterator,params[j]) == true)
                  filteredrec.insert(*theIterator);
              }
        }
      }
      return filteredrec;
  }
int main (int argc, char* argv[])
{
 Setfilepath("hs_mgb_1.0.ref_genome.37.1.gene.chr");
 SetSharedMemoryName("hs_mgb_1.0.ref_genome.37.1.gene");
 std::cout<<GetSharedMemoryName()<<"\n";
 std::cout<<GetPartialSharedMemoryName()<<"\n";
 FilterRecord(argv,argc);  
 return 0;
}

