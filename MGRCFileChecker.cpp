/*
 * File:   MGRCFileChecker.cpp
 * Author: manish
 *
 * Created on September 15, 2011, 10:42 AM
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

struct stat_desc
{
     size_t filesize;
     size_t readcount;
     size_t basecount;
     stat_desc(size_t _filesize, size_t _readcount, size_t _basecount):filesize(_filesize),readcount(_readcount),basecount(_basecount) {}
};

class SXFileChecker
{
    bool pair_end;
    bool long_reads;
    bool solid_reads;
    bool cg_reads;
    bool ill_reads;
    char dirname[1000];
    char readtype[1000];
    vector <const char*> filevec;
    typedef map <const char*, stat_desc*> Statmap;
    Statmap outmap;

    public:
    SXFileChecker()
    {
         pair_end    = false;
         ill_reads   = false;
         long_reads  = false;
         solid_reads = false;
         cg_reads    = false;
    }

    void manual()
    {
         fprintf(stderr, "\nUsage: MGRCFilechecker -pe -r <Read Type> -d <dirname>\n");
         fprintf(stderr, "\noptions:\n");
         fprintf(stderr, "\n\t-pe\t\t\tpair_end\n");
         fprintf(stderr, "\n\t-r  <Read Type>\t\till/454/abi/cg\n\n");
    }

    void generateStatistics(const char* flname)
    {
          size_t data_size = 0;
          size_t numreads  = 0;
          size_t numbases  = 0;

          if(strcmp(readtype,"ill") == 0)
             ill_reads = true;
          else if(strcmp(readtype,"454") == 0)
             long_reads = true;
          else if(strcmp(readtype,"abi") == 0)
             solid_reads = true;
          else if(strcmp(readtype,"cg") == 0)
             cg_reads = true;

          char header;

          if(ill_reads )
             header = '@';
          else
             header = '>';

          bool headflag = false;
          char qline[2000];
          char qSequence[2000];
          char qHeader[2000];

          FILE * queryFile = fopen(flname, "r");
          if (queryFile == NULL)
          {
              fprintf(stderr,"Could not open input file\n");
              exit(1);
          }

          fseek(queryFile, 0L, SEEK_END);
          data_size = ftell(queryFile);
          fseek(queryFile, 0L, SEEK_SET);

          while (!feof(queryFile))
          {
               fgets(qline, 1024, queryFile);

               if(cg_reads)
               {
                   qline[strlen(qline) - 1] = '\0';
                   if(strlen(qline) != 0 && qline[0] != '#' && qline[0] != '>' )
                   {
                          char* pch;
                          pch   = strtok(qline,"\t");
                          pch   = strtok(NULL,"\t");         //reads
                          numreads++;
                          numbases +=  strlen(pch);
                   }
               }
               else
               {
                   if ((qline[0] == header) || (feof(queryFile)))
                   {
                       if(ill_reads)
                       {
                          if(qline[0] == header && strstr(qline,"length"))
                             numreads++;
                       }
                       else if (qline[0] == header)
                            numreads++;

                       headflag = true;

                       if(qSequence[0] == '\0')
                       {
                            qline[strlen(qline) - 1] = '\0';
                            strcpy(qHeader, qline);
                       }
                       else
                       {
                             if(long_reads || solid_reads)
                                numbases +=  strlen(qSequence);

                             qSequence[0] = '\0';
                             qline[strlen(qline) - 1] = '\0';
                             strcpy(qHeader, qline);
                       }
                   }                                                                                                                                                            else
                   {
                          bool flag = true;
                          int ncount = 0;
                          qline[strlen(qline) - 1] = '\0';

                          for(int i = 0; i < strlen(qline); i++)
                          {
                                if(!isalpha(qline[i]))
                                {
                                    flag = false;
                                }
                          }

                          if(flag == true && headflag == true)
                          {
                                if(ill_reads)
                                    numbases += strlen(qline);
                                headflag = false;
                          }
                          strcat(qSequence, qline);
                   }
               }
         }
         fclose(queryFile);
         outmap[flname] = new stat_desc (data_size,numreads,numbases);
    }

   void loadfiles( const char* drname)
   {
	   DIR *pDIR;
	   struct dirent *entry;
	   if( pDIR=opendir(drname) )
	   {
	       while(entry = readdir(pDIR))
	       {
	            if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 )
	            {
	                   char flname[1000];
	                   strcpy(flname,drname);
	                   strcat(flname,"/");
	                   strcat(flname,entry->d_name);
	                   int len = strlen(flname);
                           char* name = new char[strlen(flname)+1];
                           strcpy(name,flname);
                           filevec.push_back(name);
	             }
	        }
	        closedir(pDIR);
	  }
          else
             fprintf(stderr, "could not open %s\n", drname);
   }

   bool integritycheck()
   {
         char line1[1000], line2[1000];

         for (int i =0; i < filevec.size(); i++)
         {
              bool integretiyflag = false;

              FILE * fp1 = fopen(filevec[i], "r");
              if (fp1 == NULL)
              {
                  fprintf(stderr,"Could not open input file\n");
                  exit(1);
              }

              fseek(fp1, 0, SEEK_SET);
              fgets(line1,1000, fp1);
              fclose(fp1);


              for (int j =0; j < filevec.size(); j++)
              {
                  if(strcmp(filevec[i], filevec[j]) != 0)
                  {
                     FILE * fp2 = fopen(filevec[i], "r");
                     if (fp2 == NULL)
                     {
                        fprintf(stderr,"Could not open input file\n");
                        exit(1);
                     }

                     fseek(fp2, 0, SEEK_SET);
                     fgets(line2,1000, fp2);
                     fclose(fp2);


                     if(strcmp(line1,line2) == 0)
                        integretiyflag = true;

                     if(integretiyflag == false)
                         return 0;
                   }
             }
         }
         return 1;
   }

   void setdirname(char* _dirname)
   {
        strcpy(dirname,_dirname);
        loadfiles(dirname);
   }

   void setreadtype(char* _readtype)
   {
        strcpy(readtype,_readtype);
   }

   void setintegrityflag(char dirflag)
   {
        if(dirflag == 'Y')
           pair_end = true;
   }

   void parseCmdLine(int argc, char *argv[])
    {
	 for (int i = 1; i < argc; i++)
         {
                if (strcmp(argv[i], "-r") == 0)
                {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -r requires a readtype.\n");
				exit(9);
			}
                        setreadtype(argv[++i]);
			continue;
		 }

		 if (strcmp(argv[i], "-d") == 0)
                 {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -d requires a filename.\n");
				exit(9);
			}

			setdirname(argv[++i]);
			continue;
		 }

                 if (strcmp(argv[i], "-pe") == 0)
                 {
                        pair_end = true;
                        continue;
                 }
         }
     }

   void run()
   {
            if(pair_end)
               printf("Data Integrity: %s\n", integritycheck() ?"OK":"Not OK");

            for (int j = 0; j < filevec.size(); j++)
                generateStatistics(filevec[j]);
   }

   Statmap getstatistics()
   {
           cout <<"File Name\tFile Size\tNumber of reads\tNumber of bases\n";

           for( Statmap::iterator ii=outmap.begin(); ii!= outmap.end(); ++ii)
               printf("%s\t%d\t%d\t%d\n",(*ii).first,((*ii).second)->filesize,((*ii).second)->readcount,((*ii).second)->basecount);

           return outmap;
   }

   void freememory()
   {
         outmap.clear();
         filevec.clear();
   }
};

int main(int argc, char* argv[])
{
    SXFileChecker* sb = new SXFileChecker();

    char dirname[1000];
    char readtype[100];
    char dirflag;

    #ifdef _WIN32
    cout << "Pls enter directory name\n";
    cin >> dirname;
    sb->setdirname(dirname);

    cout << "Pls enter readtype\n";
    cin >> readtype;
    sb->setreadtype(readtype);

    cout << "Data integrity test needed [Y/N]\n";
    cin >> dirflag;
    sb->setintegrityflag(dirflag);
    
    #elif __unix__
    if(argc != 5 &&  argc != 6)
    {
       sb->manual();
       delete sb;
       exit(1);
    }

    sb->parseCmdLine(argc,argv);
    #endif
    
    sb->run();
    sb->getstatistics();
    sb->freememory();

    delete sb;
}


