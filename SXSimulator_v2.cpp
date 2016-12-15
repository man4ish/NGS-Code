/*
 * File:   SXSimulator.cpp
 * Author: manish
 *
 * Created on March 10, 2011, 8:55 AM
 */


#include <stdio.h>
#include <ctype.h>
#include <cstring>
#include <limits>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <set>
#include <sys/types.h>
#include <dirent.h>
#include <queue>
#include <vector>
#include <time.h>
#include <fstream>
#include <map>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;

class Sequencer
{
    double inspercent;
    double delpercent;
    double snppercent;
    double cnvpercent;
    int transcount;
    double invpercent;
    double SRmean;
    double RDmean;
    double PEGAPmean;
    double Bases_apartmean;
    double SRstdev;
    double RDstdev;
    double PEGAPstdev;
    double Bases_apartstdev;
    int minins;
    int mindel;
    int maxins;
    int maxdel;
    int minmnv;
    int maxmnv;
    size_t maxcnv;
    size_t mincnv;
    size_t mintrans;
    size_t maxtrans;
    size_t mininv;
    size_t maxinv;
    size_t maxloh;
    size_t minloh;
    double serr;
    double sd;
    double Nerr;
    double pbaseerr;
    double lohpercent;
    bool cnvflag;
    bool invflag;
    bool transflag;
    bool rc;
    bool seed;
    bool perfect;
    bool complete_genomics;
    bool hetero;
    bool pair_end;
    bool random;
    bool LOH;
    int lohcount;
    float mutpercent;
    int coverage;
    int size;
    size_t used;
    size_t NCount;
    size_t NucleotideCount;
    size_t rdcount;
    size_t seednum;
    int meanpd;
    char inputFileName[256];
    char outFileName[256];
    char dirName[500];
    char outfileprefix[500];
    bool dist;
    char* bases;
    int sdsize;
    size_t genomesize;
    int margin;
    int cnvnum;
    typedef map <int,int> Distmap;
    vector <const char*> Genomevec;
    Distmap pedistmap,mutatiogapdist,Lengthmap;
    vector <int> chrvec;

    public:
    Sequencer()
    {
        rc = false;
        genomesize = 0;
        mutpercent = 0.00;
        coverage = 30;
        size = 30;
        perfect = false;
        complete_genomics = false;
        used = 0;
        hetero = false;
        pair_end = false;
        cnvflag = false;
        invflag = false;
        transflag = false;
        serr = 0.00;
        Nerr = 0.00;
        sd = 0.00;
        meanpd = 100;
        random = false;
        NCount = 0;
        lohcount =0;
        NucleotideCount = 0;
        dist = false;
        SRmean = 0;
        RDmean = 0;
        PEGAPmean = 0;
        Bases_apartmean = 0;
        SRstdev = 0.0;
        RDstdev = 0.0;
        PEGAPstdev = 0;
        Bases_apartstdev = 0.0;
        inspercent = 0.0;
        delpercent = 0.0;
        snppercent = 0.0;
        cnvpercent = 0.0;
        transcount = 0;
        invpercent = 0.0;
        seed = false;
        rdcount = 0;
        minins = 1;
        mindel = 1;
        maxins = 1;
        maxdel = 1;
        minmnv = 1;
        maxmnv = 1;
        mininv = 1;
        maxinv = 1;
        maxcnv = 1000;
        mincnv = 1000000;
        maxinv = 1000;
        mininv = 1000000;
        mintrans = 1000;
        maxtrans = 1000000;
        seednum = 1;
        pbaseerr = 0.00;
        sdsize = 0;
        LOH = false;
        maxloh = 1000000;
        minloh = 1000;
    }

    ~Sequencer()
    {

    }

    struct Key
    {
        int chr;
        size_t offset;
        char snv_type;
        Key(int _chr, size_t _offset, char _snv_type):chr(_chr),offset(_offset),snv_type(_snv_type) {}
    };

    struct SNV
    {
        char* original_allele;
        char* mutated_allele;
        char snv_type;
        int suppreads;

        SNV(char* _original_allele, char* _mutated_allele)
        {
            original_allele = new char [strlen(_original_allele) + 1];
            strcpy(original_allele, _original_allele);
            mutated_allele = new char [strlen(_mutated_allele) + 1];
            strcpy(mutated_allele, _mutated_allele);
        }
        ~SNV()
        {
            if (original_allele != NULL)
            {
                delete [] original_allele;
                original_allele = NULL;
            }
            if (mutated_allele != NULL)
            {
                delete [] mutated_allele;
                mutated_allele = NULL;
            }
        }
    };

    struct CNV
    {
        char* ref_allele;
        char* dupallele;
        int copynum;
        CNV(char* _ref_allele, char* _dupallele,int _copynum)
        {
            ref_allele = new char [strlen(_ref_allele) + 1];
            strcpy(ref_allele, _ref_allele);
            dupallele = new char [strlen(_dupallele) + 1];
            strcpy(dupallele, _dupallele);
            copynum = _copynum;
        }
        ~CNV()
        {
            if (ref_allele != NULL)
            {
                delete [] ref_allele;
                ref_allele = NULL;
            }
            if (dupallele != NULL)
            {
                delete [] dupallele;
                dupallele = NULL;
            }
        }
    };

    struct ltstr
    {
	   inline bool operator()(const Key* _lhs, const Key* _rhs) const
	   {
	            if(_lhs->chr != _rhs->chr)
	                return ( _lhs->chr < _rhs->chr);
	            return  ( _lhs->offset < _rhs->offset );
	   }
    };

    struct SV
    {
		int chr;
		size_t start;
		size_t stop;
		SV(int _chr, size_t _start, size_t _stop)
		{
			chr   = _chr;
			start = _start;
			stop  = _stop;
		}
    };

    struct comparator
    {
		inline bool operator()(const SV* _lhs, const SV* _rhs) const
		{
			  if(_lhs->chr != _rhs->chr)
			       return ( _lhs->chr < _rhs->chr);
			  return  ( _lhs->start < _rhs->start);
	         }
    };

    typedef multimap< Key*,SNV*,ltstr> SNVmap;
    SNVmap snvoffmap;
    typedef map <int,char*> Genome;
    Genome hap1, hap2;
    typedef map<Key*, CNV*, ltstr> CNVmap;
    CNVmap cmap, lohmap;
    typedef map <SV*,size_t,comparator> INVmap;
    INVmap imap;
    typedef map <SV*, SV*,comparator> Transmap;
    Transmap tmap;

    void manual()
    {
         fprintf(stderr, "\nUsage: SXSimulator [-options]\n");
         fprintf(stderr, "\nInput options\n");
         fprintf(stderr, "\n\t-d <input directory>                            Input directory\n");
         fprintf(stderr, "\t-f <file name>                                  Input file\n");
         fprintf(stderr, "\nMutation options\n\n");
         fprintf(stderr, "\t-snp <double> <INT> <INT>                       <SNP percentage> <min size> <max size>\n");
         fprintf(stderr, "\t-ins <double> <INT> <INT>                       <Insert percentage> <min size> <max size>\n");
         fprintf(stderr, "\t-del <double> <INT> <INT>                       <Delete percentage> <min size> <max size>\n");
         fprintf(stderr, "\t-cnv <double> <INT> <INT>                       <CNV percentage> <min size> <max size>\n");
         fprintf(stderr, "\t-loh <double> <INT> <INT>                       <LOH percentage> <min size> <max size>\n");
         fprintf(stderr, "\t-trans <INT> <INT> <INT>                        <Translocation Count> <min size> <max size>\n");
         fprintf(stderr, "\t-inv <double> <INT> <INT>                       <Inversion> <min size> <max size>\n");
         fprintf(stderr, "\t-r                                              Reverse complement\n");
         fprintf(stderr, "\t-rand                                           Random forward/reverse, left/right\n");
         fprintf(stderr, "\t-seed                                           Seed for random mutation position\n");
         fprintf(stderr, "\t-c <INT>                                        Coverage\n");
         fprintf(stderr, "\t-s <INT> <INT>                                  <Mean read length> <Standard Deviation>\n");
         fprintf(stderr, "\t-cg                                             Read type\n");
         fprintf(stderr, "\t-pe_gap <INT> <double>                          Mean pair end distance & standard deviation\n");
         fprintf(stderr, "\t-het                                            Hetero\n");
         fprintf(stderr, "\t-serr <double>                                  Sequencing Error\n");
         fprintf(stderr, "\t-N <double>                                     Ambiguous Base\n");
         fprintf(stderr, "\t-pbaseerr <double>                              PolyBase Error\n");
         fprintf(stderr, "\nOutput options\n");
         fprintf(stderr, "\n\t-o <Output File Prefix>                         Output file prefix\n");
         fprintf(stderr, "\t-dist                                           Output distribution file for SNV, CNV, INV, Translocation, base apart & Pair end\n\n");
    }

    void parseCmdLine(int argc, char *argv[])
    {
	 for (int i = 1; i < argc; i++)
         {
                if (strcmp(argv[i], "-o") == 0)
                {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -o requires a output file prefix.\n");
				exit(9);
			}

			strcpy(outfileprefix, argv[++i]);
			continue;
		 }

		 if (strcmp(argv[i], "-f") == 0)
                 {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -f requires a filename.\n");
				exit(9);
			}

			strcpy(inputFileName, argv[++i]);
                        Genomevec.push_back(inputFileName);
			continue;
		  }

                  if(strcmp(argv[i], "-d") == 0)
                  {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -d input a directory.\n");
				exit(9);
			}
			strcpy(dirName, argv[++i]);
                        loadGenome(dirName);
			continue;
		  }

                 

		  if (strcmp(argv[i], "-c") == 0)
                  {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -c requires a coverage value.\n");
				exit(9);
			}

			coverage = (atoi(argv[++i]));
			continue;
		  }

		  if (strcmp(argv[i], "-p") == 0)
                  {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -p requires a mutation percent for SNPs.\n");
				exit(9);
			}

			mutpercent = atof(argv[++i]);
			continue;
		   }

		   if (strcmp(argv[i], "-r") == 0)
                   {
                        rc = true;
			continue;
		   }


                   if (strcmp(argv[i], "-seed") == 0)
                   {
			seed = true;
			continue;
		   }

                   if (strcmp(argv[i], "-s") == 0)
                   {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -s requires a read lengh.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -s requires sd of read length\n");
				exit(9);
			}
			size = atoi(argv[++i]);
                        sdsize = atoi(argv[++i]);
			continue;
                   }

                   if (strcmp(argv[i], "-cnv") == 0)
                   {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -cnv requires a CNV %.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -cnv requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -cnv requires maximum insert size\n");
				exit(9);
			}

                        cnvpercent = atof(argv[++i]);
                        cnvflag = true;
                        mincnv = (size_t)atol(argv[++i]);
                        maxcnv = (size_t)atol(argv[++i]);
			continue;
                   }

                 if (strcmp(argv[i], "-loh") == 0)
                   {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -loh requires a LOH %.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -loh requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -loh requires maximum insert size\n");
				exit(9);
			}

                        lohpercent = atof(argv[++i]);
                        LOH = true;
                        minloh = (size_t)atol(argv[++i]);
                        maxloh = (size_t)atol(argv[++i]);
			continue;
                   }

                   if (strcmp(argv[i], "-trans") == 0)
                   {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -trans requires a CNV %.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -trans requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -trans requires maximum insert size\n");
				exit(9);
			}

                        transcount = atoi(argv[++i]);
                        if(transcount > 300)
                        {
                             fprintf(stderr, "Maximum translocation limited to 300.\n");
                             exit(9);
                        }
                        transflag = true;
                        mintrans = (size_t)atol(argv[++i]);
                        maxtrans = (size_t)atol(argv[++i]);
			continue;
                   }

                   if (strcmp(argv[i], "-inv") == 0)
                   {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -inv requires a CNV %.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -inv requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -inv requires maximum insert size\n");
				exit(9);
			}
                        invflag = true;
                        invpercent = atof (argv[++i]);
                        mininv = (size_t) atol(argv[++i]);
                        maxinv = (size_t) atol(argv[++i]);
			continue;
                   }

                   if (strcmp(argv[i], "-cg") == 0)
                   {
			complete_genomics = true;
			continue;
		   }

                   if (strcmp(argv[i], "-het") == 0)
                   {
			hetero = true;
			continue;
		   }

                   if (strcmp(argv[i], "-rand") == 0)
                   {
			random = true;
			continue;
		   }

                   if (strcmp(argv[i], "-dist") == 0)
                   {
			dist = true;
			continue;
		   }

                   if (strcmp(argv[i], "-pe_gap") == 0)
                   {

                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -pe_gap requires a pair end distance\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -pe_gap requires standard deviation\n");
				exit(9);
			}

			pair_end = true;
                        PEGAPmean = atoi(argv[++i]);
                        PEGAPstdev = atof(argv[++i]);
			continue;
		   }

                   if (strcmp(argv[i], "-snp") == 0)
                   {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -snp requires a SNP mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -snp requires minimum MNV size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -snp requires maximum MNV size\n");
				exit(9);
			}

                        snppercent = atof(argv[++i]);
                        minmnv = atoi(argv[++i]);
                        maxmnv = atoi(argv[++i]);
			continue;
		  }

                  if (strcmp(argv[i], "-ins") == 0)
                  {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -ins requires a INS mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -ins requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -ins requires maximum insert size\n");
				exit(9);
			}

                        inspercent = atof(argv[++i]);
                        minins = atoi(argv[++i]);
                        maxins = atoi(argv[++i]);
			continue;
		  }

                  if (strcmp(argv[i], "-del") == 0)
                  {
                        if (i >= argc - 1)
                        {
				fprintf(stderr,  "Parameter -del requires a DEL mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -del requires minimum delete size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -del requires maximum delete size \n");
				exit(9);
			}

                        delpercent = atof(argv[++i]);
                        mindel = atoi(argv[++i]);
                        maxdel = atoi(argv[++i]);
			continue;
		  }

                  if (strcmp(argv[i], "-serr") == 0)
                  {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -serr requires a sequencing error percentage\n");
				exit(9);
			}

                        serr = atof(argv[++i]);
			continue;
		  }

                  if (strcmp(argv[i], "-N") == 0)
                  {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -N requires a N(wildcard) percentage\n");
				exit(9);
			}

                        Nerr = atof(argv[++i]);
			continue;
		  }

                 if (strcmp(argv[i], "-pbaseerr") == 0)
                 {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -N requires a N(wildcard) percentage\n");
				exit(9);
			}

                        pbaseerr = atof(argv[++i]);
			continue;
		 }
	 }
  }

  void loadGenome( const char* drname)
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
		           if(flname[len-1] == 'a' && flname[len-2] == 'f' && flname [len-3] == '.')
                           {
                              char* name = new char[strlen(flname)+1];
                              strcpy(name,flname);
                              Genomevec.push_back(name);
                           }
	             }
	        }
	        closedir(pDIR);
	  }
  }

  void loadChrFile(const char* inputFileName )
  {
        size_t used = 0;
        char qString[1024] = {0};
        size_t allocated = 0;
        char header[1024] = {0};
        int chr_num;
        FILE* inFile = fopen(inputFileName, "r");
        if (inFile ==NULL)
        {
            printf("Could not open file\n");
            exit(1);
        }

        while (fgets(qString, 1024, inFile))
        {
            char *p = strchr(qString, '\n');
            if (p != NULL)
                *p = '\0';
            if (qString[0] == '>')
            {
                strcpy(header, qString);
                char *chr_ptr = strstr(qString, "NC_") + 3;
                chr_num = (unsigned)atoi(chr_ptr);
                continue;
            }
            else
            {
                unsigned sLen = strlen(qString);
                if (sLen + used > allocated)
                {
                    fprintf(stderr, "Used = %lu\r", used);
                    bases = (char*)realloc(bases, allocated + 2000);
                    allocated += 2000;
                }
                for (size_t i =0; i < sLen; i++)
                {
                      if(qString[i] == 'N')
                         NCount++;
                      else
                         NucleotideCount++;
                }
                memmove(bases + used, qString, sLen);
                used += sLen;
            }
        }

        fclose(inFile);
        bases[used] = '\0';
        char* chrbase = new char[used+1];
        strcpy(chrbase,bases);
        chrbase[used]= '\0';
        chrvec.push_back(chr_num);
        hap1[chr_num] = chrbase;
        hap2[chr_num] = chrbase;
        genomesize += used;
        free(bases);
        bases = NULL;
    }

    int genUniformRandnumber(size_t min, size_t max)
    {
        if(min == max)
            return max;
        return (min + rand()%(max-min));
    }

    bool makeReverse(char *dest, const char *src)
    {
        int sLen = strlen(src);
        src += sLen - 1;
        while (sLen-- > 0)
        {
            *dest = *src;
            dest++;
            src--;
        }
        return true;
    }

    bool makeReverseComplement(char *dest, const char *src)
    {
        int srcLen = strlen(src);
        src += srcLen - 1;
        while (srcLen-- > 0)
        {
            switch(*src)
            {
                case 'a':
                case 'A':
                    *dest = 'T';
                    break;
                case 'c':
                case 'C':
                    *dest = 'G';
                    break;
                case 'g':
                case 'G':
                    *dest = 'C';
                    break;
                case 't':
                case 'T':
                    *dest = 'A';
                    break;
                default:
                    *dest = toupper(*src);
                    break;
            }
            src--;
            dest++;
        }
        *dest = '\0';
        return true;
    }

  void genSNV(int chrnum)     // need to implement homo-hetero
  {
       typedef unsigned int UINT;
       //bool stop = false;
       size_t nsize1,nsize2;
       bool mutation = false;
       //bool strand = true;

       char* mutation_bases;
       char* pmutation_bases;

       if(seed)
       {
            seednum =   time(NULL);
            srand ( seednum );
       }
       char nucleotides[4] = {'A','C','G','T'};

       size_t slen1 = strlen(hap1[chrnum]);
       size_t slen2 = strlen(hap2[chrnum]);

       size_t used = slen1;

       size_t snpcount = size_t(used * snppercent)/(size_t)(100);
       size_t inscount = size_t(used * inspercent)/(size_t)(100);
       size_t delcount = size_t(used * delpercent)/(size_t)(100);

       if(snpcount > 0 ||  inscount > 0 || delcount > 0 )
         mutation = true;

       if (mutation)
       {
           mutation_bases = new char[used];
           pmutation_bases = new char[used];
           for (size_t i = 0; i < used; i++)
           {
                pmutation_bases[i] = 0;
                mutation_bases[i] = 0;
           }

           char* chrnome1 = new char[slen1+1];
           char* chrnome2 = new char[slen2+1];
           strcpy(chrnome1,hap1[chrnum]);
           strcpy(chrnome2, hap2[chrnum]);
           chrnome1[slen1] = '\0';
           chrnome2[slen2] = '\0';

           fprintf(stderr, "\n");
           nsize1 = slen1;
           nsize2 = slen2;

           while (inscount > 0)
           {
                char insallele[2000];                             //max 2000 bases
                bool mut_ok = true;
                fprintf(stderr, "Mutating Insert %lu        \r", inscount);
                size_t mut_pos = rand() % used;
                while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

                 for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
                 {
                      if (chrnome1[i] == 'N')
                      {
                                mut_ok = false;
                                break;
                      }
                  }
                  if (mut_ok == false)
                     continue;

                int select = rand()%2;

                size_t randinssize = genUniformRandnumber(minins,maxins);

                char* refallele = new char[2];
	        refallele[0] = '-';
                refallele[1] = '\0';

                int irdx =0;
                for (size_t sz =0; sz < randinssize; sz++)
                {
                    insallele[irdx++] = nucleotides[rand()%4];
                }
                insallele[irdx++] = '\0';

                snvoffmap.insert(pair<Key*, SNV*>(new Key(chrnum,mut_pos+1,'I'), new SNV(refallele,insallele)));

                if(hetero)
                {
                    if(select == 0)
                    {
                        mutation_bases[mut_pos] = '+';
                        pmutation_bases[mut_pos] = '+';
                        nsize1 += randinssize;
                        nsize2 += randinssize;
                    }
                    else
                    {
                        mutation_bases[mut_pos] = '+';
                        nsize1 += randinssize;
                    }
                }
                else
                {
                        mutation_bases[mut_pos] = '+';
                        pmutation_bases[mut_pos] = '+';
                        nsize1 += randinssize;
                        nsize2 += randinssize;
                }
                inscount--;
         }
         fprintf(stderr, "\n");

         while (snpcount > 0)
	 {
	    	 fprintf(stderr, "Introducing SNP %lu        \r", snpcount);

                 bool mut_ok = true;
	    	 size_t mut_pos = rand() % used;
	    	 while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
	  	       mut_pos = rand() % used;
	    	 for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
	    	 {
	    	      if (chrnome1[i] == 'N')
	    	      {
		                mut_ok = false;
		                break;
	    	      }
	    	  }
	  	  if (mut_ok == false)
		     continue;
		  int select = rand()%2;

		  size_t randsnpsize = genUniformRandnumber(minmnv,maxmnv);
		  for (size_t sz =0; sz < randsnpsize; sz++)
		  {
			    int ipos = rand()%4;
		            while(nucleotides[ipos] == hap1[chrnum][mut_pos+sz])
				  ipos = rand()%4;

		            char* refallele = new char[2];
	                    refallele[0] = chrnome1[mut_pos+sz];
			    refallele[1] = '\0';

			    char* variant = new char[2];
			    variant[0] = nucleotides[ipos];
			    variant[1] = '\0';
		            snvoffmap.insert(pair<Key*, SNV*>(new Key(chrnum,mut_pos+sz+1,'S'), new SNV(refallele,variant)));

                            if(hetero)
                            {
                                if(select == 0)
			        {
				    mutation_bases[mut_pos+sz] = 'S';
				    pmutation_bases[mut_pos+sz] = 'S';
			        }
			        else
		                     mutation_bases[mut_pos+sz] = 'S';
                             }
                             else
                             {
                                   mutation_bases[mut_pos+sz] = 'S';
                                   pmutation_bases[mut_pos+sz] = 'S';
                             }
		  }
		  snpcount--;
         }

         while (delcount > 0)
	 {
		           bool mut_ok = true;
		           fprintf(stderr, "Introducing delete %lu        \r", delcount);
		           size_t mut_pos = rand() % used;
		           while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
		               mut_pos = rand() % used;

		           for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
		           {
		                if (chrnome1[i] == 'N')
		                {
		                      mut_ok = false;
		                      break;
		                }
		           }
		           if (mut_ok == false)
		              continue;
		           int select = rand()%2;

		           size_t randelsize = genUniformRandnumber(mindel,maxdel);
		           for (size_t sz =0; sz < randelsize; sz++)
		           {
	                        char* refallele = new char[2];
			        refallele[0] = chrnome1[mut_pos+sz];
			        refallele[1] = '\0';

                                char* variant = new char[2];
			        variant[0] = '-';
			        variant[1] = '\0';

		                snvoffmap.insert(pair<Key*, SNV*>(new Key(chrnum,mut_pos+sz+1,'D'), new SNV(refallele,variant)));
                                if(hetero)
                                {
		                    if(select == 0)
		                    {
				        mutation_bases[mut_pos+sz] = '-';
				        pmutation_bases[mut_pos+sz] = '-';
                                        nsize1--;
                                        nsize2--;
			            }
			            else
                                    {
		                        mutation_bases[mut_pos+sz] = '-';
                                        nsize1--;
                                    }
                                 }
                                 else
                                 {
                                        mutation_bases[mut_pos+sz] = '-';
                                        pmutation_bases[mut_pos+sz] = '-';
                                        nsize1--;
                                        nsize2--;
                                 }
		           }
		           delcount--;
        }

       char* chrgenome1 = new char[nsize1+1];
       char* chrgenome2 = new char[nsize2+1];
       size_t ind1 = 0;
       size_t ind2 = 0;

       for (size_t i = 0; i < slen1;i++)                                //hap1
       {
	     if ((mutation_bases == NULL) || (mutation_bases[i] == 0))
	     {
	             switch (chrnome1[i])
	             {
	                       case 'A':
	                       case 'C':
	                       case 'G':
	                       case 'T':
	                       case 'a':
	                       case 'c':
	                       case 'g':
	                       case 't':
	                       case 'N':
	                          chrgenome1[ind1++] = chrnome1[i];
	                           break;
	                       default:
	                           break;
	               }
	        }
	        else if(mutation_bases[i] == '+')
	        {
                       SNVmap::const_iterator insfound;
                       Key* inskey = new Key(chrnum, i+1, 'I');
                       insfound = snvoffmap.find (inskey);
	               chrgenome1[ind1++] = chrnome1[i];
                       if (insfound != snvoffmap.end ())
                       {
                           int inslen = strlen(snvoffmap.find(inskey)->second->mutated_allele);
                           char* varallele = new char[inslen+1];
                           strcpy(varallele,snvoffmap.find(inskey)->second->mutated_allele);
                           varallele[inslen] = '\0';

	                   for (int sz = 0; sz <inslen; sz++)
	                   {
	                       chrgenome1[ind1++] = varallele[sz];
	                   }
                       }
	        }
                else if(mutation_bases[i] == 'S')
                {
                       SNVmap::const_iterator snpfound;
                       Key* snpkey = new Key(chrnum, i+1, 'S');
                       snpfound = snvoffmap.find (snpkey);
                       if (snpfound != snvoffmap.end ())
                       {
                           int mutlen = strlen(snvoffmap.find(snpkey)->second->mutated_allele);
                           char* varallele = new char[mutlen+1];
                           strcpy(varallele,snvoffmap.find(snpkey)->second->mutated_allele);
                           varallele[mutlen] = '\0';

                           for (int sz = 0; sz < mutlen; sz++)
                           {
                               chrgenome1[ind1++] = varallele[sz];
                           }
                       }
                }
                else if(mutation_bases[i] == 'D')
                {

                }
	     }
	     chrgenome1[ind1++] = '\0';
             hap1[chrnum] = chrgenome1;


      for (size_t i = 0; i < slen2;i++)                                //hap2
      {
	     if ((pmutation_bases == NULL) || (pmutation_bases[i] == 0))
	     {
	             switch (chrnome2[i])
	             {
	                       case 'A':
	                       case 'C':
	                       case 'G':
	                       case 'T':
	                       case 'a':
	                       case 'c':
	                       case 'g':
	                       case 't':
	                       case 'N':
	                          chrgenome2[ind2++] = chrnome2[i];
	                           break;
	                       default:
	                           break;
	               }
	        }
	        else if(pmutation_bases[i] == '+')
	        {
	               chrgenome2[ind2++] = chrnome2[i];
                       int inslen = strlen(snvoffmap.find(new Key(chrnum, i+1, 'I'))->second->mutated_allele);
                       char* varallele = new char[inslen+1];
                       strcpy(varallele,snvoffmap.find(new Key(chrnum, i+1, 'I'))->second->mutated_allele);
                       varallele[inslen] = '\0';

	               for (int sz = 0; sz <inslen; sz++)
	               {
	                    chrgenome2[ind2++] = varallele[sz];
	               }
	        }
                else if(pmutation_bases[i] == 'S')
                {
                       int mutlen = strlen(snvoffmap.find(new Key(chrnum, i+1, 'S'))->second->mutated_allele);
                       char* varallele = new char[mutlen+1];
                       strcpy(varallele,snvoffmap.find(new Key(chrnum, i+1, 'S'))->second->mutated_allele);
                       varallele[mutlen] = '\0';

                       for (int sz = 0; sz < mutlen; sz++)
                       {
                            chrgenome2[ind2++] = varallele[sz];
                       }
                }
                else if(pmutation_bases[i] == 'D')
                {

                }
	     }
	     chrgenome2[ind2++] = '\0';
             hap2[chrnum] = chrgenome2;

            delete [] chrnome1;
            delete [] chrnome2;
            delete [] mutation_bases;
            delete [] pmutation_bases;
      }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
void genLOH(int chrnum)    // chromosome wise
{
    size_t slen1 = strlen(hap1[chrnum]);
    size_t slen2 = strlen(hap2[chrnum]);
    
    /* Getting hap1 */
    char* chrhapnome1 = new char[slen1+1];
    strcpy(chrhapnome1,hap1[chrnum]);
    chrhapnome1[slen1] = '\0';

    /* Getting hap2 */
    char* chrhapnome2 = new char[slen2+1];
    strcpy(chrhapnome2,hap2[chrnum]);
    chrhapnome2[slen2] = '\0';

    char* mutation_bases = new char[slen1];
    for (size_t i = 0; i < slen1; i++)
	mutation_bases[i] = 0;

    char* pmutation_bases = new char[slen2];
    for (size_t i = 0; i < slen2; i++)
	pmutation_bases[i] = 0;

    size_t nsize1 = slen1;
    size_t nsize2 = slen2;

    lohcount = (lohpercent*slen1)/(100*maxloh);
    
    while (lohcount > 0)
    {
         bool mut_ok = true;
         fprintf(stderr, "Introducing LOH %lu        \r", lohcount);         
         
         size_t mut_pos = rand() % slen1;
         while ((mut_pos <= size_t(size)) || (mut_pos >= (slen1-size)) || (mutation_bases[mut_pos] != 0))
                mut_pos = rand() % slen1;
         
         int select = rand()%2;

         size_t lohsize = genUniformRandnumber(minloh,maxloh);
         
         for (size_t sz = 0; sz < lohsize; sz++)
         {
               if(mutation_bases[mut_pos+sz] != 0 || chrhapnome1[mut_pos+sz] == 'N')
               {
                      mut_ok = false;
                      break;
               }
         }

         if (mut_ok == false)
            continue;

         char* refallele = new char[2];
         refallele[0] = '-';
         refallele[1] = '\0';

         char* dupallele = new char[lohsize+1];
         memcpy(dupallele,chrhapnome1+mut_pos,lohsize);
         dupallele[lohsize] = '\0';

         lohmap.insert(pair<Key*, CNV*>(new Key(chrnum,mut_pos+1,'I'), new CNV(refallele,dupallele,1)));

         mutation_bases[mut_pos] = 'l';             // need to update the method for both chromosome
         nsize1 += lohsize;

         for (size_t sz = 0; sz < lohsize; sz++)
         {
               if(mutation_bases[mut_pos+sz] == 0)
               {
                      mutation_bases[mut_pos+sz] = 1;
               }
         }

         for (size_t sz = 0; sz < lohsize; sz++)
         {
               if(pmutation_bases[mut_pos+sz] == 0)
               {
                      pmutation_bases[mut_pos+sz] = 'l';
               }
         }
         lohcount--;
    }
    fprintf(stderr, "\n");

    char* chrgenome1 = new char[nsize1];
    size_t ind1 =0;

    for (size_t i = 0; i < slen1;i++)                                //for Hap1
    {
               if ((mutation_bases == NULL) || (mutation_bases[i] == 0))
               {
                  switch (chrhapnome1[i])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome1[ind1++] = chrhapnome1[i];
                        break;
                    default:
                        break;
                  }
                }
                else if(mutation_bases[i] == 'l')
                {
                     Key* cnvkey = new Key(chrnum,i+1,'I');
                     CNVmap::const_iterator dupinsert;
                     dupinsert = cmap.find (cnvkey);
                     int cpnum;
                     char* dupall;
                     chrgenome1[ind1++] = hap1[chrnum][i];
                     if (dupinsert != cmap.end ())
                     {
                         cpnum =  (*dupinsert->second).copynum;
                         size_t duplen = strlen((*dupinsert->second).dupallele);
                         dupall = new char[duplen+1];
                         strcpy(dupall,(*dupinsert->second).dupallele);
                         dupall[duplen] ='\0';
                         for (int cp = 0; cp < cpnum; cp++)
                         {
                            for (size_t sz =0; sz < duplen; sz++)
                            {
                               chrgenome1[ind1++] = dupall[sz];
                            }
                         }
                     }
                 }
      }
      chrgenome1[ind1++] = '\0';
      hap1[chrnum] = chrgenome1;

      size_t ind2 =0;                                    /* Processing Hap2 */
      char* chrgenome2 = new char[nsize2];
      for (size_t i = 0; i < slen2;i++)
       {
               if((pmutation_bases == NULL) || (pmutation_bases[i] == 0))
               {
                  switch (chrhapnome2[i])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome2[ind2++] = chrhapnome2[i];
                        break;
                    default:
                        break;
                   }
               }
               else if(pmutation_bases[i] == 'l')
               {

               }
      }
      chrgenome2[ind2++] = '\0';
      hap2[chrnum] = chrgenome2;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genCNV(int chrnum)    // chromosome wise
{

    size_t slen1 = strlen(hap1[chrnum]);
    size_t slen2 = strlen(hap2[chrnum]);

    /* Getting hap1 */
    char* chrhapnome1 = new char[slen1+1];
    strcpy(chrhapnome1,hap1[chrnum]);
    chrhapnome1[slen1] = '\0';

    /* Getting hap2 */
    char* chrhapnome2 = new char[slen2+1];
    strcpy(chrhapnome2,hap2[chrnum]);
    chrhapnome1[slen2] = '\0';

    char* mutation_bases = new char[slen1];
    char* pmutation_bases = new char[slen2];

    for (size_t i = 0; i < slen1; i++)
	mutation_bases[i] = 0;

    for (size_t i = 0; i < slen2; i++)
        pmutation_bases[i] = 0;

    size_t nsize1 = slen1;
    size_t nsize2 = slen2;

    int cnvcount = (int)(cnvpercent*(slen1/100));

    int dupcount = cnvcount/maxcnv;
    int delcount = 4*dupcount;

    while (dupcount > 0)
    {
         bool mut_ok = true;
         fprintf(stderr, "Introducing duplicates (cnv) %lu        \r", dupcount);

         size_t mut_pos = rand() % slen1;
         while ((mut_pos <= size_t(size)) || (mut_pos >= (slen1-size)) || (mutation_bases[mut_pos] != 0))
                mut_pos = rand() % slen1;

         int copynum = 1+rand()%4;
         int select = rand()%2;

         

         size_t randupsize =genUniformRandnumber(mincnv,maxcnv);


         for (size_t sz = 0; sz < randupsize; sz++)
         {
               if(mutation_bases[mut_pos+sz] != 0 || chrhapnome1[mut_pos+sz] == 'N')
               {
                      mut_ok = false;
                      break;
               }
         }

         if (mut_ok == false)
            continue;

         char* refallele = new char[2];
         refallele[0] = '-';
         refallele[1] = '\0';

         char* dupallele = new char[randupsize+1];
         memcpy(dupallele,chrhapnome1+mut_pos,randupsize);
         dupallele[randupsize] = '\0';

         cmap.insert(pair<Key*, CNV*>(new Key(chrnum,mut_pos+1,'I'), new CNV(refallele,dupallele,copynum)));

         if(select == 0)
         {
              mutation_bases[mut_pos] = '+';             // need to update the method for both chromosome
              //pmutation_bases[mut_pos] = '+';
              //nsize2 += copynum*randupsize;
              nsize1 += copynum*randupsize;
         }
         else
         {
              mutation_bases[mut_pos] = '+';
              nsize1 += copynum*randupsize;
         }

         for (size_t sz = 0; sz < randupsize; sz++)
         {
               if(mutation_bases[mut_pos+sz] == 0)
               {
                      mutation_bases[mut_pos+sz] = 1;
               }
         }
         dupcount--;
    }
    fprintf(stderr, "\n");


    while (delcount > 0)
    {
         bool mut_ok = true;
         fprintf(stderr, "Introducing deletes (cnv) %lu        \r", delcount);
         size_t mut_pos = rand() % slen1;

         while ((mut_pos <= size_t(size)) || (mut_pos >= (slen1 - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % slen1;

         int select = rand()%2;

         for (size_t i = mut_pos - size; i <= mut_pos + size; i++)
         {
                  if (chrhapnome1[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
         }

         if (mut_ok == false)
            continue;

         size_t randelsize = genUniformRandnumber(mincnv,maxcnv);
         for (size_t sz =0; sz < randelsize; sz++)
         {
               if(mutation_bases[mut_pos+sz] != 0 || chrhapnome1[mut_pos+sz] == 'N')
               {
                      mut_ok = false;
                      break;
               }
         }

         if (mut_ok == false)
            continue;

         char delallele[1000000];
         memcpy(delallele,chrhapnome1+mut_pos,randelsize);
         delallele[randelsize] = '\0';

         char varallele[2];
         varallele[0] = '-';
         varallele[1] = '\0';


         cmap.insert(pair<Key*, CNV*>(new Key(chrnum,mut_pos+1,'D'), new CNV(delallele,varallele,0)));
         for (size_t sz =0; sz < randelsize; sz++)
         {
              if(select == 0)
              {
                 mutation_bases[mut_pos+sz] = '-';
                 pmutation_bases[mut_pos+sz] = '-';
                 nsize2--;
                 nsize1--;
              }
              else
              {
                mutation_bases[mut_pos+sz] = '-';
                nsize1--;
              }
         }
         delcount--;
    }
    fprintf(stderr, "\n");

    char* chrgenome1 = new char[nsize1];
    size_t ind1 =0;

    for (size_t i = 0; i < slen1;i++)                                //for Hap1
    {
               if ((mutation_bases == NULL) || (mutation_bases[i] == 0))
               {
                  switch (chrhapnome1[i])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome1[ind1++] = chrhapnome1[i];
                        break;
                    default:
                        break;
                  }
                }
                else if(mutation_bases[i] == '+')
                {
                     Key* cnvkey = new Key(chrnum,i+1,'I');
                     CNVmap::const_iterator dupinsert;
                     dupinsert = cmap.find (cnvkey);
                     int cpnum;
                     char* dupall;
                     chrgenome1[ind1++] = hap1[chrnum][i];
                     if (dupinsert != cmap.end ())
                     {
                         cpnum =  (*dupinsert->second).copynum;
                         size_t duplen = strlen((*dupinsert->second).dupallele);
                         dupall = new char[duplen+1];
                         strcpy(dupall,(*dupinsert->second).dupallele);
                         dupall[duplen] ='\0';
                         for (int cp = 0; cp < cpnum; cp++)
                         {
                            for (size_t sz =0; sz < duplen; sz++)
                            {
                               chrgenome1[ind1++] = dupall[sz];
                            }
                         }
                     }
                 }
                 else if(mutation_bases[i] == '-')
                 {

                 }
      }
      chrgenome1[ind1++] = '\0';
      hap1[chrnum] = chrgenome1;

      size_t ind2 =0;                                    /* Processing Hap2 */
      char* chrgenome2 = new char[nsize2];
      for (size_t i = 0; i < slen2;i++)
       {
               if((pmutation_bases == NULL) || (pmutation_bases[i] == 0))
               {
                  switch (chrhapnome2[i])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome2[ind2++] = chrhapnome2[i];
                        break;
                    default:
                        break;
                   }
               }
               else if(pmutation_bases[i] == '+')
               {
                     Key* cnvkey = new Key(chrnum,i+1,'I');
                     CNVmap::const_iterator dupinsert;
                     dupinsert = cmap.find (cnvkey);
                     int cpnum;
                     char* dupall;
                     chrgenome1[ind1++] = hap1[chrnum][i];
                     if (dupinsert != cmap.end ())
                     {
                         cpnum =  (*dupinsert->second).copynum;
                         size_t duplen = strlen((*dupinsert->second).dupallele);
                         dupall = new char[duplen+1];
                         strcpy(dupall,(*dupinsert->second).dupallele);
                         dupall[duplen] ='\0';
                         for (int cp = 0; cp < cpnum; cp++)
                         {
                            for (size_t sz =0; sz < duplen; sz++)
                            {
                               chrgenome1[ind1++] = dupall[sz];
                            }
                         }
                     }
                }
                else if(pmutation_bases[i] == '-')
                {

                }
      }
      chrgenome2[ind2++] = '\0';
      hap2[chrnum] = chrgenome2;
}

void genInvesrion(int nchr)
{
     size_t chrsize1 = strlen(hap1[nchr]);
     char* chrnome1 = new char[chrsize1+1];
     strcpy(chrnome1,hap1[nchr]);
     chrnome1[chrsize1]= '\0';

     size_t chrsize2 = strlen(hap2[nchr]);
     char* chrnome2 = new char[chrsize2+1];
     strcpy(chrnome2,hap2[nchr]);
     chrnome2[chrsize2]= '\0';

     size_t invcount = (size_t)((invpercent*chrsize1)/(maxinv*100));

     char* mutation_bases =  new char[chrsize1];
     for (size_t i = 0; i < chrsize1; i++)
        mutation_bases[i] = 0;

     char* pmutation_bases = new char[chrsize2];
     for (size_t i = 0; i < chrsize2; i++)
        pmutation_bases[i] = 0;

     while (invcount > 0)
     {
         fprintf(stderr, "Introducing Inversion %lu        \r", invcount);
         size_t chrsize;
         int haptype = rand()%2;
         if(haptype == 1)
          chrsize = chrsize1;
         else
          chrsize = chrsize2;

         size_t randpos1 = rand()% chrsize;
         size_t randlen = genUniformRandnumber(mininv,maxinv);
         size_t randpos2 = randpos1+randlen;

         if(randpos2 > chrsize)
           randpos2 -= randpos1;

         if(randpos1 > randpos2)
         {
            size_t temp = randpos2;
            randpos2 = randpos1;
            randpos1 = temp ;
         }

         bool ok = true;
         char* seq = new char[randpos2-randpos1+1];
         size_t idx =0;
         for (size_t i = randpos1; i < randpos2 ; i++ )
         {
                if (haptype == 1)
                  seq[idx++] = chrnome1[i];
                else
                  seq[idx++] = chrnome2[i];
         }
         seq[randpos2-randpos1] = '\0';
         for (size_t i = randpos1; i < randpos2 ; i++ )
         {
           if(mutation_bases[i] != 0 || chrnome1[i] == 'N')
           {
               ok = false;
           }
         }

        if( ok == false )
           continue;

       bool hetero = true;

       if(ok)
       {
           imap.insert(pair<SV*, size_t>(new SV(nchr,randpos1+1,randpos2+1), randpos2+1));
           char* revseq = new char[randpos2-randpos1+1];
           revseq[randpos2 - randpos1] = '\0';
           makeReverse(revseq,seq);
           size_t range = randpos2-randpos1;

           for (size_t i = 0; i < range ; i++ )
           {
                if(hetero)
                {
                   mutation_bases[randpos1+i] = revseq[i];
                }
                else
                {
                   mutation_bases[randpos1+i] = revseq[i];
                   pmutation_bases[randpos1+i] = revseq[i];
                }
           }
           delete [] revseq;
        }
        delete [] seq;
        invcount--;
    }

    char* chrgenome1 = new char[chrsize1+1];
    char* chrgenome2 = new char[chrsize2+1];
    size_t ind1 = 0;
    size_t ind2 = 0;

    for (size_t gsize =0; gsize < chrsize1; gsize++)         //hap1
    {
         if((mutation_bases == NULL) || (mutation_bases[gsize] == 0))
         {
                  switch (chrnome1[gsize])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome1[ind1++] = chrnome1[gsize];
                        break;
                    default:
                        break;
                   }
          }
          else
          {
                 chrgenome1[ind1++] = mutation_bases[gsize];
          }
    }
    chrgenome1[ind1++] = '\0';
    hap1[nchr] = chrgenome1;

    for (size_t gsize =0; gsize < chrsize2; gsize++)         //hap2
    {
         if((pmutation_bases == NULL) || (pmutation_bases[gsize] == 0))
         {
                  switch (chrnome2[gsize])
                  {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                    case 'N':
                        chrgenome2[ind2++] = chrnome2[gsize];
                        break;
                    default:
                        break;
                   }
          }
          else
          {
                 chrgenome2[ind2++] = pmutation_bases[gsize];
          }
    }
    chrgenome2[ind2++] = '\0';
    hap2[nchr] = chrgenome2;
 }


 void genTranslocation()
 {
     typedef map <int, char*> Mgenome;
     Mgenome mutgenome, pmutgenome;

     for( map<int,char*>::iterator ii= hap1.begin(); ii!= hap1.end(); ii++)
     {
          size_t chromlen1 = strlen(hap1[(*ii).first]);
          char* mutationbases = new char[chromlen1+1];
          for (size_t i =0; i < chromlen1; i++)
             mutationbases[i] = 0;

          size_t chromlen2 = strlen(hap2[(*ii).first]);
          char* pmutationbases = new char[chromlen2+1];
          for (size_t i =0; i < chromlen2; i++)
             pmutationbases[i] = 0;

          mutgenome[(*ii).first]  = mutationbases;
          pmutgenome[(*ii).first] = pmutationbases;
     }

     int randcnt = hap1.size();
     int trcount = transcount;

     while(trcount > 0)         // Need Enhancement and include Intrachromosomal translocation
     {
         fprintf(stderr, "Introducing Translocation %lu        \r", trcount);

         bool ok = true;
         int randchr1 = chrvec[rand()%randcnt];
         int randchr2 = chrvec[rand()%randcnt];

         size_t chrlen1 = strlen(hap1[randchr1]);
         size_t chrlen2 = strlen(hap2[randchr2]);

         size_t randpos1 = genUniformRandnumber(1,chrlen1-maxtrans);
         size_t randlen1 = genUniformRandnumber(mintrans,maxtrans);
         size_t randpos2 = (randlen1+randpos1);

         size_t randpos3 = genUniformRandnumber(1,chrlen2-maxtrans);
         size_t randlen2 = genUniformRandnumber(mintrans,maxtrans);
         size_t randpos4 = (randlen2+randpos3);

         char* markgenome1 = new char[randlen2+1];
         memcpy(markgenome1,mutgenome[randchr1]+randpos3,randlen2);
         markgenome1[randlen2] = '\0';

              for (size_t j = 0; j < randlen2; j++ )
              {
                  if(markgenome1[j] != 0 || hap1[randchr1][j] == 'N')
                  {
                     ok = false;
                     break;
                  }
              }

              delete [] markgenome1;

              if(ok == false)
                continue;

              char* markgenome2 = new char[randlen1+1];
              memcpy(markgenome2,mutgenome[randchr2]+randpos1,randlen1);
              markgenome2[randlen1] = '\0';

             for (size_t j = 0; j < randlen1; j++ )
             {
                if(markgenome2[j] != 0 || hap2[randchr2][j] == 'N')
                {
                   ok = false;
                   break;
                }
             }

             delete [] markgenome2;

            if(ok == false)
              continue;
             
          tmap.insert(pair<SV*, SV*>(new SV(randchr1,randpos1,randpos2), new SV(randchr2,randpos3,randpos4)));

          if(randchr1 == randchr2) // Interatranslocation
          {            
              char* mutchrom1 = new char[chrlen1+1];
              char* genchrom1 = new char[chrlen1+1];
            
              char* pmutchr1 = new char[randpos1+1];
              char* pmutchr2 = new char[randpos2-randpos1+1];
              char* pmutchr3 = new char[randpos3-randpos2+1];
              char* pmutchr4 = new char[randpos4-randpos3+1];
              char* pmutchr5 = new char[chrlen1-randpos4+1];

              char* pgenchr1 = new char[randpos1+1];
              char* pgenchr2 = new char[randpos2-randpos1+1];
              char* pgenchr3 = new char[randpos3-randpos2+1];
              char* pgenchr4 = new char[randpos4-randpos3+1];
              char* pgenchr5 = new char[chrlen1-randpos4+1];

              memcpy(pmutchr1,mutgenome[randchr1],randpos1);
              memcpy(pmutchr2,mutgenome[randchr1]+randpos1,randpos2-randpos1);
              memcpy(pmutchr3,mutgenome[randchr1]+randpos2,randpos3-randpos2);
              memcpy(pmutchr4,mutgenome[randchr1]+randpos3,randpos4-randpos3);
              memcpy(pmutchr5,mutgenome[randchr1]+randpos4,chrlen1-randpos4);

              memcpy(pgenchr1,hap1[randchr1],randpos1);
              memcpy(pgenchr2,hap1[randchr1]+randpos1,randpos2-randpos1);
              memcpy(pgenchr3,hap1[randchr1]+randpos2,randpos3-randpos2);
              memcpy(pgenchr4,hap1[randchr1]+randpos3,randpos4-randpos3);
              memcpy(pgenchr5,hap1[randchr1]+randpos4,chrlen1-randpos4);
           
              pmutchr1[randpos1] = '\0';
              pmutchr3[randpos3-randpos2] = '\0';        
              pmutchr5[chrlen1-randpos4] = '\0';
            
              pgenchr1[randpos1] = '\0';
              pgenchr2[randpos2-randpos1] = '\0';
              pgenchr3[randpos3-randpos2] = '\0';
              pgenchr4[randpos4-randpos3] = '\0';
              pgenchr5[chrlen1-randpos4] = '\0';
         
              strcpy(mutchrom1,pmutchr1);
              strcat(mutchrom1,pgenchr2);
              strcat(mutchrom1,pmutchr3);
              strcat(mutchrom1,pgenchr4);
              strcat(mutchrom1,pmutchr5);

              strcpy(genchrom1,pgenchr1);
              strcat(genchrom1,pgenchr4);
              strcat(genchrom1,pgenchr3);
              strcat(genchrom1,pgenchr2);
              strcat(genchrom1,pgenchr5);
            
            
              delete [] pmutchr1;
              delete [] pmutchr3;
              delete [] pmutchr5;

              delete [] pgenchr1;
              delete [] pgenchr2;
              delete [] pgenchr3;
              delete [] pgenchr4;
              delete [] pgenchr5;

              if(hetero)
                  hap1[randchr1] = genchrom1;
              else
              {
                 hap1[randchr1] = genchrom1;
                 hap2[randchr1] = genchrom1;
              }
          }
          else  // InterTranslocation
          {
              size_t mutchromlen1 = randpos1 + (randpos4-randpos3) + (chrlen1-randpos2);
              size_t mutchromlen2 = randpos3 + (randpos2-randpos1) + (chrlen2-randpos4);

              char* mutchrom1 = new char[mutchromlen1+1];
              char* mutchrom2 = new char[mutchromlen2+1];
              char* genchrom1 = new char[mutchromlen1+1];
              char* genchrom2 = new char[mutchromlen2+1];

              char* pmutchr1 = new char[randpos1+1];
              char* pmutchr2 = new char[randpos2-randpos1+1];
              char* pmutchr3 = new char[chrlen1-randpos2+1];
              char* pmutchr4 = new char[randpos3+1];
              char* pmutchr5 = new char[randpos4-randpos3+1];
              char* pmutchr6 = new char[chrlen2-randpos4+1];

              memcpy(pmutchr1,mutgenome[randchr1],randpos1);
              memcpy(pmutchr2,hap1[randchr1]+randpos1,randpos2-randpos1);
              memcpy(pmutchr3,mutgenome[randchr1]+randpos2,chrlen1-randpos2);
              memcpy(pmutchr4,mutgenome[randchr2],randpos3);
              memcpy(pmutchr5,hap1[randchr2]+randpos3,randpos4-randpos3);
              memcpy(pmutchr6,mutgenome[randchr2]+randpos4,chrlen2-randpos4);

            char* pgenchr1 = new char[randpos1+1];
            char* pgenchr3 = new char[chrlen1-randpos2+1];
            char* pgenchr4 = new char[randpos3+1];
            char* pgenchr6 = new char[chrlen2-randpos4+1];

            memcpy(pgenchr1,hap1[randchr1],randpos1);
            memcpy(pgenchr3,hap1[randchr1]+randpos2,chrlen1-randpos2);
            memcpy(pgenchr4,hap1[randchr2],randpos3);
            memcpy(pgenchr6,hap1[randchr2]+randpos4,chrlen2-randpos4);

            pgenchr1[randpos1] ='\0' ;
            pgenchr3[chrlen1-randpos2] = '\0' ;
            pgenchr4[randpos3] = '\0';
            pgenchr6[chrlen2-randpos4] = '\0';

            pmutchr1[randpos1] ='\0' ;
            pmutchr2[randpos2-randpos1] = '\0';
            pmutchr3[chrlen1-randpos2] = '\0' ;
            pmutchr4[randpos3] = '\0';
            pmutchr5[randpos4-randpos3] = '\0';
            pmutchr6[chrlen2-randpos4] = '\0';

            strcpy(mutchrom1,pmutchr1);
            strcat(mutchrom1,pmutchr5);
            strcat(mutchrom1,pmutchr3);

            strcpy(mutchrom2,pmutchr4);
            strcat(mutchrom2,pmutchr2);
            strcat(mutchrom2,pmutchr6);

            strcpy(genchrom1,pgenchr1);
            strcat(genchrom1,pmutchr5);
            strcat(genchrom1,pgenchr3);

            strcpy(genchrom2,pgenchr4);
            strcat(genchrom2,pmutchr2);
            strcat(genchrom2,pgenchr6);

            mutchrom1[mutchromlen1] = '\0';
            mutchrom2[mutchromlen2] = '\0';
            genchrom1[mutchromlen1] = '\0';
            genchrom2[mutchromlen2] = '\0';

            mutgenome[randchr1] = mutchrom1;
            mutgenome[randchr2] = mutchrom2;

            if(hetero)
            {
               hap1[randchr1] = genchrom1;
               hap1[randchr2] = genchrom2;
            }
            else
            {
               hap1[randchr1] = genchrom1;
               hap1[randchr2] = genchrom2;
               hap2[randchr1] = genchrom1;
               hap2[randchr2] = genchrom2;
            }

           delete [] pmutchr1;
           delete [] pmutchr2;
           delete [] pmutchr3;
           delete [] pmutchr4;
           delete [] pmutchr5;
           delete [] pmutchr6;

           delete [] pgenchr1;
           delete [] pgenchr3;
           delete [] pgenchr4;
           delete [] pgenchr6;
       }
       trcount--;
   }
     for (Mgenome::iterator it = mutgenome.begin(); it != mutgenome.end(); it++)
     {
        mutgenome.erase(it->first);
        pmutgenome.erase(it->first);
     }
     mutgenome.clear();
     pmutgenome.clear();
 }

void genSeqErr ( char* seq)
{
        int len = strlen(seq);
        int pos = rand() % len;
        char nucleotide[4] = {'A','C','G','T'};
        int randbase = rand() % 4;
        while(seq[pos] == nucleotide[randbase])
            randbase = rand() % 4;
        seq[pos] = nucleotide[randbase];
}

void genWildcardErr ( char* seq)
{
        int len = strlen(seq);
        int pos = rand() % len;
        seq[pos] = 'N';
}

void genPolybaseError(char* seq)
{
        int len = strlen(seq);
        int pos =  rand() % (len-1);
        if(seq[pos] != seq[pos+1])
           seq[pos] = seq[pos+1];
}

void SequenceReads(int chrnum)
{
      boost::normal_distribution<> nd(PEGAPmean, PEGAPstdev);
      boost::mt19937 rng;
      boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);

      size_t pos;
      int pairend_gap;
      int paired_size;
      int haptype;

      int perrcount =0;
      int Nerrcount=0;
      int serrcount=0;
      if(pbaseerr !=0 )
         perrcount = (int)(100/pbaseerr);
      if(Nerr !=0 )
         Nerrcount = (int)(100/Nerr);
      if (serr !=0)
         serrcount = (int)(100/serr);

      char* qual_score = new char[size+1];

      for(int i =0; i < size; i++)
       qual_score[i] = '[';
      qual_score[size] = '\0';

      size_t chrsize1 = strlen(hap1[chrnum]);
      char* chrnome1 = new char[chrsize1+1];
      strcpy(chrnome1,hap1[chrnum]);
      chrnome1[chrsize1] = '\0';

      size_t chrsize2 = strlen(hap2[chrnum]);
      char* chrnome2 = new char[chrsize2+1];
      strcpy(chrnome2,hap2[chrnum]);
      chrnome2[chrsize2] = '\0';

      char* fwdbaseshap1 = new char[chrsize1];
      char* revbaseshap1 = new char[chrsize1];
      for (size_t gs =0; gs < chrsize1; gs++)
      {
          fwdbaseshap1[gs] = 0;
          revbaseshap1[gs] = 0;
      }

      char* fwdbaseshap2 = new char[chrsize2];
      char* revbaseshap2 = new char[chrsize2];
      for (size_t gs =0; gs < chrsize2; gs++)
      {
          fwdbaseshap2[gs] = 0;
          revbaseshap2[gs] = 0;
      }

      size_t total_reads = size_t (chrsize1*coverage)/(size_t)(2*size);
      size_t numreads = total_reads;

      int strand;
      char chr[3];
      sprintf(chr, "%d", chrnum);
      strcpy(outFileName,outfileprefix);
      strcat(outFileName,"CHR_");
      strcat(outFileName,chr);

      if(complete_genomics)
          strcat(outFileName,".fqx");
      else
          strcat(outFileName,".fastq");

      FILE * OutFile = fopen(outFileName, "w");
         assert(OutFile != NULL);

      size_t read_count = 0;

      //bool ok = true;
      char seq[10000];
      char seq_rc[10000];
      int gap1 = 0;
      int gap2 = 0;

      while(numreads > 0)
      {
            haptype = rand() % 2;
            strand = rand() % 2;
            pairend_gap = (int)var_nor();

                 if(complete_genomics)
                 {
                          switch(rand() % 10)
                          {
                            case 0:
                            gap1 = 5;
                                break;
                            case 1:
                            gap1 = 5;
                                break;
                            case 2:
                            gap1 = 7;
                                break;
                            default:
                            gap1 = 6;
                                break;
                           }

                           switch(rand() % 10)
                           {
                              case 0:
                              gap2 = 5;
                                break;
                                case 1:
                              gap2 = 5;
                                break;
                                case 2:
                              gap2 = 7;
                                break;
                                default:
                              gap2 = 6;
                                break;
                            }
                 }

                 paired_size = 2*size + gap1 + gap2 + pairend_gap;

                 if(haptype)
                 {
                     pos = rand() % chrsize1;
                     if (strand)
                     {
                         while ((pos >= (chrsize1 - paired_size - 1))||(fwdbaseshap1[pos] == true))
                         pos = rand() % chrsize1;
                     }
                     else
                     {
                         while ((pos >= (chrsize1 - paired_size - 1))||(revbaseshap1[pos] == true))
                         pos = rand() % chrsize1;
                     }
                     memcpy(seq,chrnome1+pos,paired_size);
                 }
                 else
                 {
                     pos = rand() % chrsize2;
                     if (strand)
                     {
                         while ((pos >= (chrsize2 - paired_size - 1))||(fwdbaseshap2[pos] == true))
                         pos = rand() % chrsize2;
                     }
                     else
                     {
                         while ((pos >= (chrsize2 - paired_size - 1))||(revbaseshap2[pos] == true))
                         pos = rand() % chrsize2;
                     }
                     memcpy(seq,chrnome2+pos,paired_size);
                 }

            seq[paired_size] = '\0';
            makeReverseComplement(seq_rc, seq);

            ///////////////////////////CG//////////////////////////////////
            char output_leftseq[200];
            char output_rightseq[200];
            char output_leftseq_rc[200];
            char output_rightseq_rc[200];

            int seq_idx_lt = 0;
            int seq_idx_rt = 0;
            int seq_idx_lt_rc = 0;
            int seq_idx_rt_rc = 0;

            if(complete_genomics)
            {
                  int readIdx = 0;

                  for (; readIdx < 30 + gap1; readIdx++)
                  {
                        output_leftseq[seq_idx_lt++] = seq[readIdx];
                        output_leftseq_rc[seq_idx_lt_rc++] = seq_rc[readIdx];
                        if (readIdx == 19)
                            readIdx += gap1;
                  }

                  readIdx += pairend_gap;

                  for (; readIdx < 2*size + pairend_gap + gap1 + gap2; readIdx++)
                  {
                        output_rightseq[seq_idx_rt++] = seq[readIdx];
                        output_rightseq_rc[seq_idx_rt_rc++] = seq_rc[readIdx];
                        if (readIdx == (size+gap1+pairend_gap+9))
                            readIdx += gap2;
                  }

                  output_leftseq[30]    = '\0';
                  output_rightseq[30]   = '\0';
                  output_leftseq_rc[30] = '\0';
                  output_rightseq_rc[30]= '\0';
            }
            else
            {
                 memcpy(output_leftseq,seq,size+gap1);
                 output_leftseq[size+gap1] = '\0';

                 memcpy(output_rightseq,seq+size+gap1+pairend_gap,size+gap2);
                 output_rightseq[size+gap2] = '\0';

                 memcpy(output_leftseq_rc,seq_rc,size+gap1);
                 output_leftseq_rc[size+gap1] = '\0';

                 memcpy(output_rightseq_rc,seq_rc+size+gap1+pairend_gap,size+gap2);
                 output_rightseq_rc[size+gap2] = '\0';
            }


            if(pedistmap.count(pairend_gap) == 0)
                pedistmap[pairend_gap] = 1;
            else
                pedistmap[pairend_gap]++;

            int lrflag = rand()%2;

            if (strand == false)    // rev
            {
                  if(!strstr(output_leftseq_rc,"N") && !strstr(output_rightseq_rc,"N"))
                  {
                       if(haptype)
                          revbaseshap1[pos] = true;
                       else
                          revbaseshap2[pos] = true;
                       fprintf(OutFile,"@RD_CHR:%d_MateGap:%lu_STRAND:REV\n", chrnum,pairend_gap);

                       // Keeping on any one of the side
                       if( Nerr !=0 && (rand()%Nerrcount == 0))
                       {
                             if(lrflag)
                                genWildcardErr(output_leftseq_rc);
                             else
                                genWildcardErr(output_rightseq_rc);
                       }
                       if( serr !=0 && (rand()%serrcount) == 0)
                       {
                           if(lrflag)
                             genSeqErr(output_leftseq_rc);
                           else
                             genSeqErr(output_rightseq_rc);
                       }
                       if( pbaseerr !=0 && (rand() % perrcount == 0))
                       {
                            if(lrflag)
                               genPolybaseError(output_leftseq_rc);
                            else
                               genPolybaseError(output_rightseq_rc);
                       }

                       if(complete_genomics)
                             fprintf(OutFile,"%s\t%s\n", output_leftseq_rc,output_rightseq_rc);
                       else
                       {
                           char revcom_output_rightseq_rc[200];
                           makeReverseComplement(revcom_output_rightseq_rc,output_rightseq_rc);
                           fprintf(OutFile,"%s\t%s\n+\n", output_leftseq_rc,revcom_output_rightseq_rc);
                       }
                       fprintf(OutFile,"%s\t%s\n", qual_score,qual_score);
                  }
            }
            else    //fwd
            {
                  if(!strstr(output_leftseq,"N") && !strstr(output_rightseq,"N"))
                  {
                       if(haptype)
                          fwdbaseshap1[pos] = true;
                       else
                          fwdbaseshap2[pos] = true;

                       fprintf(OutFile,"@RD_CHR:%d_MateGap:%lu_STRAND:FWD\n", chrnum,pairend_gap);


                       if( Nerr !=0 && (rand()%Nerrcount == 0))
                       {
                          if(lrflag)
                             genWildcardErr(output_leftseq);
                          else
                             genWildcardErr(output_rightseq);
                       }
                       if( serr !=0 && (rand()%serrcount) == 0)
                       {
                          if(lrflag)
                             genSeqErr(output_leftseq);
                          else
                               genSeqErr(output_rightseq);
                       }
                       if( pbaseerr !=0 && (rand() % perrcount == 0))
                       {
                           if(lrflag)
                             genPolybaseError(output_leftseq);
                           else
                               genPolybaseError(output_rightseq);
                       }
                       if(complete_genomics)
                             fprintf(OutFile,"%s\t%s\n", output_leftseq,output_rightseq);
                       else
                       {
                          char revcom_output_rightseq[200];
                          makeReverseComplement(revcom_output_rightseq,output_rightseq);
                          fprintf(OutFile,"%s\t%s\n+\n", output_leftseq,revcom_output_rightseq);
                       }
                       fprintf(OutFile,"%s\t%s\n", qual_score,qual_score);
                  }
            }

            read_count++;
            fprintf(stderr, "Read: %lu (%.2f%%)\r", read_count, double(read_count)/double(total_reads)*100);
            numreads--;
      }

      rdcount += read_count;
      delete [] chrnome1;
      delete [] chrnome2;
      delete [] fwdbaseshap1;
      delete [] revbaseshap1;
      delete [] fwdbaseshap2;
      delete [] revbaseshap2;
}


void SequenceSingleEndReads(int chrnum)
{
      boost::normal_distribution<> nd(size, sdsize);
      boost::mt19937 rng;
      boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);

      size_t pos;
      int haptype;

      int perrcount=0;
      int Nerrcount=0;
      int serrcount =0;
      if(pbaseerr !=0 )
         perrcount = (int)(100/pbaseerr);
      if(Nerr !=0 )
         Nerrcount = (int)(100/Nerr);
      if (serr !=0)
         serrcount = (int)(100/serr);

      size_t chrsize1 = strlen(hap1[chrnum]);
      char* chrnome1 = new char[chrsize1+1];
      strcpy(chrnome1,hap1[chrnum]);
      chrnome1[chrsize1] = '\0';

      size_t chrsize2 = strlen(hap2[chrnum]);
      char* chrnome2 = new char[chrsize2+1];
      strcpy(chrnome2,hap2[chrnum]);
      chrnome2[chrsize2] = '\0';

      char* fwdbaseshap1 = new char[chrsize1];
      char* revbaseshap1 = new char[chrsize1];
      for (size_t gs =0; gs < chrsize1; gs++)
      {
          fwdbaseshap1[gs] = 0;
          revbaseshap1[gs] = 0;
      }

      char* fwdbaseshap2 = new char[chrsize2];
      char* revbaseshap2 = new char[chrsize2];
      for (size_t gs =0; gs < chrsize2; gs++)
      {
          fwdbaseshap2[gs] = 0;
          revbaseshap2[gs] = 0;
      }

      size_t total_reads = size_t (chrsize1*coverage)/(size_t)(size);
      size_t numreads = total_reads;

      char qualFileName[200];
      int strand;
      char chr[3];
      sprintf(chr, "%d", chrnum);
      strcpy(outFileName,outfileprefix);
      strcat(outFileName,"CHR_");
      strcat(outFileName,chr);
      strcat(outFileName,".fasta");                     // may include option for Pac Bio

      strcpy(qualFileName,outfileprefix);
      strcat(qualFileName,"CHR_");
      strcat(qualFileName,chr);
      strcat(qualFileName,".qual");

      FILE * seqFile = fopen(outFileName, "w");
         assert(seqFile != NULL);

      FILE * qualFile = fopen(qualFileName, "w");
         assert(qualFile != NULL);

      size_t read_count = 0;

      //bool ok = true;
      char seq[10000];
      char seq_rc[10000];
      int varsize;

      while(numreads > 0)
      {
            haptype = rand() % 2;
            strand = rand() % 2;
            varsize = (int)var_nor();
            char* qual_score = new char[varsize+1];
            for(int i =0; i < varsize; i++)
              qual_score[i] = '[';
            qual_score[varsize] = '\0';

            if(haptype)
            {
                     pos = rand() % chrsize1;
                     if (strand)
                     {
                         while ((pos >= (chrsize1 - varsize - 1))||(fwdbaseshap1[pos] == true))
                         pos = rand() % chrsize1;
                     }
                     else
                     {
                         while ((pos >= (chrsize1 - varsize - 1))||(revbaseshap1[pos] == true))
                         pos = rand() % chrsize1;
                     }
                     memcpy(seq,chrnome1+pos,varsize);
             }
             else
             {
                     pos = rand() % chrsize2;
                     if (strand)
                     {
                         while ((pos >= (chrsize2 - varsize - 1))||(fwdbaseshap2[pos] == true))
                         pos = rand() % chrsize2;
                     }
                     else
                     {
                         while ((pos >= (chrsize2 - varsize - 1))||(revbaseshap2[pos] == true))
                         pos = rand() % chrsize2;
                     }
                     memcpy(seq,chrnome2+pos,varsize);
             }

             seq[varsize] = '\0';

            makeReverseComplement(seq_rc, seq);

            if(Lengthmap.count(varsize) == 0)
                Lengthmap[varsize] = 1;
            else
                Lengthmap[varsize]++;


            if(random)
            {
                 if (strand == false)    // rev
                 {
                     if(!strstr(seq_rc,"N"))
                     {
                         if(haptype)
                           revbaseshap1[pos] = true;
                         else
                           revbaseshap2[pos] = true;

                         if( Nerr !=0 && (rand()%Nerrcount == 0))
                             genWildcardErr(seq_rc);
                         if( serr !=0 && (rand()%serrcount == 0))
                             genSeqErr(seq_rc);
                         if( pbaseerr !=0 && (rand()%perrcount == 0))
                             genPolybaseError(seq_rc);

                         fprintf(seqFile,">RD_STRAND:REV_Chr_%d_Length_%d\n",chrnum,varsize);
                         fprintf(seqFile,"%s\n", seq_rc);
                         fprintf(qualFile,">RD_STRAND:REV_Chr_%d_Length_%d\n",chrnum,varsize);
                         fprintf(qualFile,"%s\n", qual_score);
                     }
                 }
                 else    //fwd
                 {
                     if(!strstr(seq,"N"))
                     {
                          if(haptype)
                             fwdbaseshap1[pos] = true;
                          else
                             fwdbaseshap2[pos] = true;

                          fprintf(seqFile,">RD_STRAND:FWD_Chr_%d_Length_%d\n",chrnum,varsize);
                          if( Nerr !=0 && (rand()%Nerrcount == 0))
                             genWildcardErr(seq);
                          if( serr !=0 && (rand()%serrcount == 0))
                             genSeqErr(seq);
                          if( pbaseerr !=0 && (rand()%perrcount == 0))
                             genPolybaseError(seq);

                          fprintf(seqFile,"%s\n", seq);
                          fprintf(qualFile,">RD_STRAND:FWD_Chr_%d_Length_%d\n",chrnum,varsize);
                          fprintf(qualFile,"%s\n", qual_score);
                      }
                 }
            }
            else if(rc == true)
            {
                     if(!strstr(seq_rc,"N"))
                     {
                         if(haptype)
                           revbaseshap1[pos] = true;
                         else
                           revbaseshap2[pos] = true;

                         if( Nerr !=0 && (rand()%Nerrcount == 0))
                             genWildcardErr(seq_rc);
                         if( serr !=0 && (rand()%serrcount == 0))
                             genSeqErr(seq_rc);
                         if( pbaseerr !=0 && (rand()%perrcount == 0))
                             genPolybaseError(seq_rc);

                         fprintf(seqFile,">RD_STRAND:REV_Chr_%d_Length_%d\n",chrnum,varsize);
                         fprintf(seqFile,"%s\n", seq_rc);

                         fprintf(qualFile,">RD_STRAND:REV_Chr_%d_Length_%d\n",chrnum,varsize);
                         fprintf(qualFile,"%s\n", qual_score);
                     }
            }
            else if(rc == false)
            {
                  if(!strstr(seq,"N"))
                  {
                         if(haptype)
                            fwdbaseshap1[pos] = true;
                         else
                            fwdbaseshap2[pos] = true;

                         if( Nerr !=0 && (rand()%Nerrcount == 0))
                             genWildcardErr(seq);
                         if( serr !=0 && (rand()%serrcount == 0))
                             genSeqErr(seq);
                         if( pbaseerr !=0 && (rand()%perrcount == 0))
                             genPolybaseError(seq);

                         fprintf(seqFile,">RD_STRAND:FWD_Chr_%dLength_%d\n",chrnum,varsize);
                         fprintf(seqFile,"%s\n", seq);

                         fprintf(qualFile,">RD_STRAND:FWD_Chr_%dLength_%d\n",chrnum,varsize);
                         fprintf(qualFile,"%s\n", qual_score);
                   }
            }
            delete [] qual_score;
            read_count++;
            fprintf(stderr, "Read: %lu (%.2f%%)\r", read_count, double(read_count)/double(total_reads)*100);
            numreads--;
      }

      rdcount += read_count;
      fclose(seqFile);
      fclose(qualFile);
      delete [] chrnome1;
      delete [] chrnome2;
      delete [] fwdbaseshap1;
      delete [] revbaseshap1;
      delete [] fwdbaseshap2;
      delete [] revbaseshap2;
}

void gen_pairend_distribution()
{
      char pdistpairendflname[500];
      strcpy(pdistpairendflname,outfileprefix);
      strcat(pdistpairendflname,"_Pairend.dist");
      FILE * pdistpairendFile = fopen (pdistpairendflname,"w");
      if (pdistpairendFile == NULL)
          perror ("Error opening in paired end distribution file");
      size_t pegaptotal =0;
      size_t pegapcount =0;
      double pegapsdsum =0.0;

      fprintf (pdistpairendFile, "PEGap\tFrequency\n");

      for( map<int, int>::iterator ii=pedistmap.begin(); ii!=pedistmap.end(); ++ii)
      {
           pegaptotal += (*ii).second;
           pegapcount++;
           fprintf (pdistpairendFile, "%d\t%d\n",(*ii).first,(*ii).second);
      }

      PEGAPmean = (double) pegaptotal/(double) pegapcount;

      for( map<int, int>::iterator ii=pedistmap.begin(); ii!=pedistmap.end(); ++ii)
      {
           pegapsdsum = (double)(PEGAPmean -(*ii).second)*(PEGAPmean -(*ii).second);
      }

      PEGAPstdev = (double)sqrt(pegapsdsum/pegaptotal);

      fclose (pdistpairendFile);
}

void gen_length_distribution()
{
      char lengthflname[500];
      strcpy(lengthflname,outfileprefix);
      strcat(lengthflname,"_Length.dist");
      FILE * lengthFile = fopen (lengthflname,"w");
      if (lengthFile == NULL)
          perror ("Error opening in length distribution file");

      fprintf (lengthFile, "Length\tFrequency\n");

      for( map<int, int>::iterator ii=Lengthmap.begin(); ii!=Lengthmap.end(); ++ii)
      {
           fprintf (lengthFile, "%d\t%d\n",(*ii).first,(*ii).second);
      }

      fclose (lengthFile);
}

void genSNVList()
{
       char SNVflname[500];
       strcpy(SNVflname,outfileprefix);
       strcat(SNVflname,"_SNV.dist");
       FILE * SNVFile = fopen (SNVflname,"w");
       if (SNVFile == NULL)
          perror ("Error opening in SNV distribution file");

       fprintf(SNVFile,"Chromosome\tOffset\tType\tReference\tVariant\n");
       queue <size_t> diffqueue;

       for (SNVmap::iterator it = snvoffmap.begin(); it != snvoffmap.end(); it++)
       {
           diffqueue.push((size_t)it->first->offset);
           if(diffqueue.size() == 2)
           {
                int previousoffset = diffqueue.front();
                int currentoffset = diffqueue.back();
                int diff  = currentoffset - previousoffset;
                if(mutatiogapdist.count(diff) == 0)
                    mutatiogapdist[diff]  = 1;
                else
                    mutatiogapdist[diff]++;
                diffqueue.pop();
           }

           fprintf(SNVFile,"%d\t%lu\t%c\t%s\t%s\n",it->first->chr,it->first->offset,it->first->snv_type,it->second->original_allele,it->second->mutated_allele);
       }

       fclose(SNVFile);
}

void gen_gapbase_distribution()
{
       char gapbaseflname[500];
       strcpy(gapbaseflname,outfileprefix);
       strcat(gapbaseflname,"_gapbase.dist");
       FILE * mutgapFile = fopen (gapbaseflname,"w");
       if (mutgapFile == NULL)
          perror ("Error opening in gap distribution file");

       size_t gaptotal =0;
       size_t gapcount =0;
       double gapsdsum =0.0;

       fprintf (mutgapFile, "Gap Size\tFrequency\n");
       for( map<int, int>::iterator ii=mutatiogapdist.begin(); ii!=mutatiogapdist.end(); ++ii)
       {
           gaptotal += (*ii).second;
           gapcount++;
           if((*ii).first > 0)
             fprintf (mutgapFile, "%d\t%d\n",(*ii).first,(*ii).second);
       }


       Bases_apartmean = double (gaptotal)/ double (gapcount);
       for( map<int, int>::iterator ii=mutatiogapdist.begin(); ii!=mutatiogapdist.end(); ++ii)
       {
           gapsdsum = (double)(Bases_apartmean -(*ii).second)*(Bases_apartmean -(*ii).second);
       }

       Bases_apartstdev = (double)sqrt(gapsdsum/gaptotal);
       fclose (mutgapFile);
}

void genCNVList()
{
       char CNVflname[500];
       strcpy(CNVflname,outfileprefix);
       strcat(CNVflname,"_CNV.dist");

       FILE * CNVFile = fopen (CNVflname,"w");
       if (CNVFile == NULL)
          perror ("Error opening in CNV distribution file");

       fprintf(CNVFile,"Chromosome\tOffset\tType\tDelete_Length\tDuplicate_Length\tCopy_Number\n");
       CNVmap::iterator it;

       for ( it=cmap.begin(); it != cmap.end(); it++ )
       {
           if((*it).first->snv_type == 'I')
               fprintf (CNVFile,"%d\t%lu\t%c\t%d\t%d\t%d\n", (*it).first->chr,(*it).first->offset, (*it).first->snv_type,0,strlen((*it).second->dupallele),(*it).second->copynum+1);
           else
              fprintf (CNVFile,"%d\t%lu\t%c\t%d\t%d\t%s\n", (*it).first->chr,(*it).first->offset, (*it).first->snv_type,strlen((*it).second->ref_allele),0,"0/1");
       }

       fclose (CNVFile);
}

void genLOHList()
{
       char LOHflname[500];
       strcpy(LOHflname,outfileprefix);
       strcat(LOHflname,"_LOH.dist");

       FILE * LOHFile = fopen (LOHflname,"w");
       if (LOHFile == NULL)
          perror ("Error opening in CNV distribution file");

       fprintf(LOHFile,"Chromosome\tOffset\tType\tDelete_Length\tDuplicate_Length\tCopy_Number\n");
       CNVmap::iterator it;

       for ( it=lohmap.begin(); it != lohmap.end(); it++ )
       {
           if((*it).first->snv_type == 'I')
               fprintf (LOHFile,"%d\t%lu\t%c\t%d\t%d\t%d\n", (*it).first->chr,(*it).first->offset, (*it).first->snv_type,0,strlen((*it).second->dupallele),(*it).second->copynum+1);
           else
              fprintf (LOHFile,"%d\t%lu\t%c\t%d\t%d\t%s\n", (*it).first->chr,(*it).first->offset, (*it).first->snv_type,strlen((*it).second->ref_allele),0,"0/1");
       }

       fclose (LOHFile);
}

void genINVList()
{
       char INVflname[500];
       strcpy(INVflname,outfileprefix);
       strcat(INVflname,"_INV.dist");

       FILE * INVFile = fopen (INVflname,"w");
       if (INVFile == NULL)
          perror ("Error opening in inversion distribution file");

       INVmap::iterator it;
       fprintf(INVFile,"Chromosome\tStart\tEnd\n");

       for ( it=imap.begin(); it != imap.end(); it++ )
           fprintf (INVFile,"%d\t%lu\t%lu\n", (*it).first->chr,(*it).first->start, (*it).first->stop);

       fclose (INVFile);
 }

void genTrasnList()
{
       char Transflname[500];
       strcpy(Transflname,outfileprefix);
       strcat(Transflname,"_Trans.dist");

       FILE * TRANSFile = fopen (Transflname,"w");
       if (TRANSFile == NULL)
          perror ("Error opening in translocation distribution file");

       Transmap::iterator it;

       fprintf(TRANSFile,"Chromosome1\tStart\tStop\tChromsome2\tStart\tStop\n");

       for ( it=tmap.begin(); it != tmap.end(); it++ )
           fprintf (TRANSFile,"%d\t%lu\t%lu\t%d\t%lu\t%lu\n", (*it).first->chr,(*it).first->start, (*it).first->stop,(*it).second->chr,(*it).second->start, (*it).second->stop);

       fclose (TRANSFile);
}

void genSummary()
{
       time_t rawtime;
       time ( &rawtime );

       char SummaryFileName[500];
       strcpy(SummaryFileName,outfileprefix);
       strcat(SummaryFileName,".Summary");
       FILE *SummaryFile = fopen (SummaryFileName,"w");
       if (SummaryFile == NULL)
          perror ("Error opening in summary file");

       fprintf (SummaryFile, "Input Data\n");

       fprintf (SummaryFile, "Input File Name\t%s\n",inputFileName);
       fprintf (SummaryFile,"Current Local Time\t%s",ctime (&rawtime));
       fprintf (SummaryFile, "N Bases\t%lu\n",(unsigned long)NCount );
       fprintf (SummaryFile, "ACGT Bases\t%lu\n",(unsigned long)NucleotideCount);
       fprintf (SummaryFile, "Total Bases\t%lu\n",(unsigned long)(NCount+NucleotideCount));
       fprintf (SummaryFile, "Output File Name\t%s\n",outFileName);

       if(random)
           fprintf (SummaryFile, "Random\t%s\n",(random)?"true":"false" );
       else
       {
           fprintf (SummaryFile, "Reverse Complement\t%s\n",(rc)?"true":"false" );
           fprintf (SummaryFile, "Forward\t%s\n",(left)?"true":"false" );
       }
       fprintf (SummaryFile, "Mutation Data\n");

       fprintf (SummaryFile, "Zygosity\t%s\n",(hetero)?"hetero":"homo" );
       fprintf (SummaryFile, "Seed\t%lu\n",(seed)?(unsigned long)seednum:1 );
       fprintf (SummaryFile, "Read Type\t%s\n",(pair_end)?"Pair End":"Single End" );

       if(pair_end)
          fprintf (SummaryFile, "Sample Type\t%s\n",(complete_genomics)?"CG":"Illumina" );
       else
          fprintf (SummaryFile, "Sample Type\t454/PacBio");
       
       fprintf (SummaryFile, "SNP Percentage\t%4.2f\n",snppercent);
       fprintf (SummaryFile, "Delete Precentage\t%4.2f\n",delpercent);
       fprintf (SummaryFile, "Insert Percentage\t%4.2f\n",inspercent);
       fprintf (SummaryFile, "CNV Percentage\t%4.2f\n",cnvpercent);
       fprintf (SummaryFile, "Number of Tranlocations\t%d\n",transcount);
       fprintf (SummaryFile, "Inversion Percentage\t%4.2f\n",invpercent);
       fprintf (SummaryFile, "Minimum Delete Size \t%d\n",mindel);
       fprintf (SummaryFile, "Maximum Delete Size \t%d\n",maxdel);
       fprintf (SummaryFile, "Minimum Insert Size \t%d\n",minins);
       fprintf (SummaryFile, "Maximum Insert Size \t%d\n",maxins);
       fprintf (SummaryFile, "Minimum MNV Size \t%d\n",minmnv);
       fprintf (SummaryFile, "Maximum MNV Size \t%d\n",maxmnv);
       fprintf (SummaryFile, "Minimum CNV Size \t%d\n",mincnv);
       fprintf (SummaryFile, "Maximum CNV Size \t%d\n",maxcnv);
       fprintf (SummaryFile, "Minimum Inversion Size \t%d\n",mininv);
       fprintf (SummaryFile, "Maximum Inversion Size \t%d\n",maxinv);
       fprintf (SummaryFile, "Minimum Translocation Size \t%d\n",mintrans);
       fprintf (SummaryFile, "Maximum Translocation Size \t%d\n",maxtrans);


       fprintf (SummaryFile, "Sequencing Stat\n");

       fprintf (SummaryFile, "Sequencing Error (illumina)\t%4.2f\n",serr);
       fprintf (SummaryFile, "Ambiguous Bases (N %)\t%4.2f\n",Nerr);
       fprintf (SummaryFile, "PolyBase Error(454)\t%4.2f\n",pbaseerr);
       fprintf (SummaryFile, "Read Length\t%d\n",size);
       fprintf (SummaryFile, "Coverage\t%d\n",coverage);
       fprintf (SummaryFile, "Number of Reads\t%lu\n",(unsigned long) rdcount);

       fprintf (SummaryFile, "SIM Density\n");

       if(dist)
       {
           fprintf (SummaryFile, "Mean Pair End Gap\t%4.2f\n",PEGAPmean);
           fprintf (SummaryFile, "Mean Distance between SNV\t%4.2f\n",Bases_apartmean);
           fprintf (SummaryFile, "Standard Deviation Pair End Gap\t%4.2f\n",PEGAPstdev);
           fprintf (SummaryFile, "Standard Deviation SNV Distance\t%4.2f\n",Bases_apartstdev);
       }

       fprintf (SummaryFile, "\n");
       fclose (SummaryFile);
}

void freeMemory()
{
      for (Distmap::iterator it=pedistmap.begin();it!=pedistmap.end();it++)
         pedistmap.erase(it->first);
      pedistmap.clear();

      for (Distmap::iterator it=mutatiogapdist.begin();it!=mutatiogapdist.end();it++)
         mutatiogapdist.erase(it->first);
      mutatiogapdist.clear();

      for (Distmap::iterator it=Lengthmap.begin();it!=Lengthmap.end();it++)
         Lengthmap.erase(it->first);
      Lengthmap.clear();

      for (Genome::iterator it= hap1.begin();it != hap1.end();it++)
      {
         hap1.erase(it->first);
         hap2.erase(it->first);
      }

      hap1.clear();
      hap2.clear();

      Genomevec.clear();

      for (CNVmap::iterator it= cmap.begin();it != cmap.end();it++)
        cmap.erase(it->first);
      cmap.clear();

      for (CNVmap::iterator it= lohmap.begin();it != lohmap.end();it++)
        lohmap.erase(it->first);
      lohmap.clear();


      for (INVmap::iterator it= imap.begin();it != imap.end();it++)
        imap.erase(it->first);
      imap.clear();

      for (Transmap::iterator it= tmap.begin();it != tmap.end();it++)
        tmap.erase(it->first);
      tmap.clear();

      for (SNVmap::iterator it= snvoffmap.begin();it !=snvoffmap.end();it++)
        snvoffmap.erase(it->first);
      snvoffmap.clear();
 }

   void run()
   {
          for (size_t i =0; i< Genomevec.size(); i++)
             loadChrFile(Genomevec[i]);

           for( map<int,char*>::iterator ii=hap1.begin(); ii!=hap1.end(); ii++)
           {
              genSNV(ii->first);
              if(cnvflag)
                genCNV(ii->first);
              if(invflag)
                genInvesrion(ii->first);            
              if(LOH)
                    genLOH(ii->first);
           }

            if(transflag)
      	       genTranslocation();

      	    for( map<int,char*>::iterator ii=hap1.begin(); ii!=hap1.end(); ii++) //read generation
            {
                       if(pair_end)
         		 SequenceReads(ii->first);
                       else
                         SequenceSingleEndReads(ii->first);
            }

      	    genSummary();

      	    if(dist)
      	    {
                 if(cnvflag)
          	    genCNVList();

                 if(invflag)
          	    genINVList();
                 
                 if(LOH)
                    genLOHList();

                 if(pair_end)
          	    gen_pairend_distribution();
                 else
                    gen_length_distribution();

          	 genSNVList();
          	 gen_gapbase_distribution();

                 if(transflag)
          	   genTrasnList();
       	    }

       	    freeMemory();
   }
};

int main(int argc, char** argv)
{
     Sequencer *gs = new Sequencer();

     if (argc == 1)
     {
        gs->manual();
        return 0;
     }

     gs->parseCmdLine(argc,argv);
     gs->run();
     delete gs;
     return (EXIT_SUCCESS);
}
