
/*
 * File:   Sequencer.cpp
 * Author: manish
 *
 * Created on November 15, 2010, 11:45 AM
 */

#include <stdlib.h>
#include <map>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
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
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;

class Sequencer
{
    double inspercent;
    double delpercent;
    double snppercent;
    double cnvpercent;
    double transpercent;
    double invpercent;
    double SRmean;
    double RDmean;
    double PEGAPmean;
    double Bases_apartmean;
    double SRstdev;
    double RDstdev;
    double PEGAPstdev;
    double Bases_apartstdev;
    int insmin;
    int delmin;
    int insmax;
    int delmax;
    int mnvmin;
    int mnvmax;
    size_t maxcnv;
    size_t mincnv;
    size_t mintrans;
    size_t maxtrans;
    size_t mininv;
    size_t maxinv;
    double serr;
    double sd;
    double N;
    double pbaseerr;
    bool rc;
    bool left;
    bool all;
    bool seed;
    bool perfect;
    bool illumina;
    bool hetero;
    bool pair_end;
    bool random;
    float mutpercent;
    int coverage;
    int size;
    size_t used;
    size_t NCount;
    size_t NucleotideCount;
    size_t read_count;
    size_t seednum;
    int meanpd;
    unsigned int chr_num;
    char inputFileName[256];
    char outFileName[256];
    char dirName[500];
    char outfileprefix[500];
    char chr[3];
    char *bases;
    char* pairbases;
    bool dist;
    typedef map <int,int> Distmap;
    vector <const char*> Filevec;
    Distmap pedistmap,mutatiogapdist;

    public:
    Sequencer()
    {
        rc = false;
        left = false;
        all = false;
        mutpercent = 0.00;
        coverage = 30;
        size = 30;
        perfect = false;
        illumina = false;
        used = 0;
        bases = NULL;
        hetero = false;
        pair_end = false;
        serr = 0.00;
        N = 0.00;
        sd = 0.00;
        meanpd = 100;
        random = false;
        NCount = 0;
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
        transpercent = 0.0;
        invpercent = 0.0;
        seed = false;
        read_count = 0;
        insmin = 1;
        delmin = 1;
        insmax = 1;
        delmax = 1;
        mnvmin = 1;
        mnvmax = 1;
        mininv = 1;
        maxinv = 1;
        maxcnv = 1000;
        mincnv = 1000000;
        maxinv = 1000;
        mininv = 1000000;
        mintrans = 1000;
        maxtrans = 1000000;
        seednum = 1;
        double pbaseerr = 0.00;
    }

    ~Sequencer()
    {

    }

    struct Key
    {
       int chr;
       size_t offset;
       char snv_type;
       Key(int  _chr, size_t  _offset, char _snv_type):chr(_chr),offset(_offset),snv_type(_snv_type) {}
    };

    struct SNV
    {
        char* original_allele;
        char * mutated_allele;
        char snv_type;
        int suppreads;

        SNV(char* _original_allele, char* _mutated_allele)
        {
            original_allele = new char [strlen(_original_allele) + 1];
            strcpy(original_allele, _original_allele);
            mutated_allele = new char [strlen(_mutated_allele) + 1];
            strcpy(mutated_allele, _mutated_allele);
            //uppreads = _suppreads;
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

    struct ltstr
    {
        inline bool operator()(const Key* _lhs, const Key* _rhs) const
        {
            if(_lhs->chr != _rhs->chr)
                return ( _lhs->chr < _rhs->chr);
            return  ( _lhs->offset < _rhs->offset );
        }
    };

    typedef multimap< Key*,SNV*,ltstr> SNVmap;
    SNVmap snvoffmap;

    void manual()
    {
        fprintf(stderr, "\nUsage: SXSimulator [-options]\n");
        fprintf(stderr, "\nInput options\n");
        fprintf(stderr, "\n\t-d <input directory>                            Input directory\n");
        fprintf(stderr, "\t-f <file name>                                  Input file\n");
        fprintf(stderr, "\nMutation options\n\n");
        fprintf(stderr, "\t-SNP <double> <INT> <INT>                       <SNP percentage> <min size> <max size>\n");
        fprintf(stderr, "\t-INS <double> <INT> <INT>                       <Insert percentage> <min size> <max size>\n");
        fprintf(stderr, "\t-DEL <double> <INT> <INT>                       <Delete percentage> <min size> <max size>\n");
        fprintf(stderr, "\t-CNV <double> <INT> <INT>                       CNV %\n");
        fprintf(stderr, "\t-Tranlocation <double> <INT> <INT>              Translocation %\n");
        fprintf(stderr, "\t-Inversion <double> <INT> <INT>                 Inversion %\n");
        fprintf(stderr, "\t-r                                              Reverse complement\n");
        fprintf(stderr, "\t-rand                                           Random forward/reverse, left/right\n");
        fprintf(stderr, "\t-l                                              Left read\n");
        fprintf(stderr, "\t-seed                                           Seed for random mutation position\n");
        fprintf(stderr, "\t-c <INT>                                        Coverage\n");
        fprintf(stderr, "\t-s <INT>                                        Read length\n");
        fprintf(stderr, "\t-ill                                            Read type\n");
        fprintf(stderr, "\t-pe <INT> <double>                              Mean pair end distance & standard deviation\n");
        fprintf(stderr, "\t-het                                            Hetero\n");
        fprintf(stderr, "\t-serr <double>                                  Sequencing Error\n");
        fprintf(stderr, "\t-N <double>                                     Ambiguous Base\n");
        fprintf(stderr, "\t-pbaseerr <double>                              PolyBase Error\n");
        fprintf(stderr, "\nOutput options\n");
        fprintf(stderr, "\n\t-o <Output File Prefix>                         Output file prefix\n");
        fprintf(stderr, "\t-dist                                           Output distribution file for SR, RD, SNP base apart & Pair end\n\n");
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
                        Filevec.push_back(inputFileName);
			continue;
		}

                if (strcmp(argv[i], "-d") == 0)
                {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -d input a directory.\n");
				exit(9);
			}
			strcpy(dirName, argv[++i]);
                        LoadFiles(dirName);
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

		if (strcmp(argv[i], "-l") == 0)
                {
			left = true;
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

			size = atoi(argv[++i]);
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
                        mincnv = (size_t)atol(argv[++i]);
                        maxcnv = (size_t)atol(argv[++i]);
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

                        transpercent = atof(argv[++i]);
                        mintrans = (size_t)atol(argv[++i]);
                        maxtrans = (size_t)atol(argv[++i]);
			continue;
                }

                if (strcmp(argv[i], "-INV") == 0)
                {
			if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -INV requires a CNV %.\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -INV requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -INV requires maximum insert size\n");
				exit(9);
			}

                        invpercent = atof (argv[++i]);
                        mininv = (size_t) atol(argv[++i]);
                        maxinv = (size_t) atol(argv[++i]);
			continue;
                }

                if (strcmp(argv[i], "-ill") == 0)
                {
			illumina = true;
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

                if (strcmp(argv[i], "-pe") == 0)
                {

                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -pe requires a pair end distance\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -pe requires standard deviation\n");
				exit(9);
			}

			pair_end = true;
                        PEGAPmean = atoi(argv[++i]);
                        PEGAPstdev = atof(argv[++i]);
			continue;
		}

                if (strcmp(argv[i], "-SNP") == 0)
                {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -SNP requires a SNP mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -SNP requires minimum MNV size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -SNP requires maximum MNV size\n");
				exit(9);
			}

                        snppercent = atof(argv[++i]);
                        mnvmin = atoi(argv[++i]);
                        mnvmax = atoi(argv[++i]);
			continue;
		}

                if (strcmp(argv[i], "-INS") == 0)
                {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -INS requires a INS mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -INS requires minimum insert size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -INS requires maximum insert size\n");
				exit(9);
			}

                        inspercent = atof(argv[++i]);
                        insmin = atoi(argv[++i]);
                        insmax = atoi(argv[++i]);
			continue;
		}

                if (strcmp(argv[i], "-DEL") == 0)
                {
                        if (i >= argc - 1)
                        {
				fprintf(stderr,  "Parameter -DEL requires a DEL mutation percentage\n");
				exit(9);
			}

                        if (i >= argc - 2)
                        {
				fprintf(stderr, "Parameter -DEL requires minimum delete size\n");
				exit(9);
			}

                        if (i >= argc - 3)
                        {
				fprintf(stderr, "Parameter -DEL requires maximum delete size \n");
				exit(9);
			}

                        delpercent = atof(argv[++i]);
                        delmin = atoi(argv[++i]);
                        delmax = atoi(argv[++i]);
			continue;
		}

                if (strcmp(argv[i], "-serr") == 0)
                {
                        if (i >= argc - 1)
                        {
				fprintf(stderr, "Parameter -INS requires a INS mutation percentage\n");
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

                        N = atof(argv[++i]);
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

    void LoadFiles( const char* drname)
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
                                     Filevec.push_back(name);
                                 }
	             }
	        }
	        closedir(pDIR);
	  }
    }


    void loadChrFile()
    {
        char qString[1024] = {0};
        size_t allocated = 0;
        char header[1024] = {0};

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
                sprintf(chr, "%d", chr_num);
                continue;
            }
            else
            {
                unsigned sLen = strlen(qString);
                for (int i =0; i < strlen(qString); i++)
                {
                    if(qString[i] == 'N')
                    {
                        NCount++;
                    }
                    else if(qString[i] == 'A' || qString[i] == 'C' ||qString[i] == 'G' ||qString[i] == 'T')
                    {
                        NucleotideCount++;
                    }
                }
                if (sLen + used > allocated)
                {
                    fprintf(stderr, "Used = %lu\r", used);
                    bases = (char*)realloc(bases, allocated + 2000);
                    allocated += 2000;
                }
                memmove(bases + used, qString, sLen);
                used += sLen;
            }
        }
        fclose(inFile);
        bases[used] = '\0';
    }

    char* substr( size_t start, size_t stop, const char* src)
    {
        if(stop <= start || strlen(src) == 0 || stop > strlen(src)|| start < 0 )
            return 0;
        size_t len= stop-start;
        char* dest = new char[len+1];
        for (size_t idx =0; idx < len; idx++)
            dest[idx] = src[idx+start];
        dest[len] = '\0';
        return dest;
    }

    int genUniformRandnumber(size_t min, size_t max)
    {
        if(min == max)
            return max;
        return (min + rand()%(max-min));
    }
/*
   void createbin(const char* filename )
   {
       float index = 0;
       int reccount = 0;
       ifstream infile(filename);
       if(!infile.is_open())
       {
          cout<<"could not open the feature file\n";
          exit(1);
       }
       char line[500],rec[500];
       while(infile)
       {
          infile.getline(line,500);
          strcpy(rec,line);
          if(infile && line[0] != '\0')
          {
            char *pch ;
            pch = strtok(line,"\t ");
            vector<const char* > vec;
            while(pch != NULL )
            {
                vec.push_back(pch);
                 pch = strtok(NULL,"\t ");
            }
            size_t offset = (size_t)(atof(vec[14])*1000);
            if(offset > 6000)
            cout << rec << endl;
            int i = offset/250;
            int min = i*250;
            int max = i*250 + 250;
            std::ostringstream minval, maxval;
            minval << min;
            maxval << max;
            string bin;
            bin += minval.str();
            bin += "\t";
            bin += maxval.str();
            vector <int> vecint;
            if(IndexbinMap.count(bin) == 0)
            {
              vecint.push_back(index);
              vecint.push_back(1);
              IndexbinMap[bin] = vecint;
           }
           else
           {
              IndexbinMap[bin][1]++;
           }
           index += strlen(rec);
           index += 1;
           vec.clear();
        }
    }
    infile.close();
}
*/
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

   void gen_pairend_distribution()
   {
       char pdistpairendflname[500];
       strcpy(pdistpairendflname,outfileprefix);
       strcat(pdistpairendflname,"_CHR");
       strcat(pdistpairendflname,chr);
       strcat(pdistpairendflname,"Pairend.dist");
       FILE * pdistpairendFile = fopen (pdistpairendflname,"w");
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

   void gen_gapbase_distribution()
   {
       char gapbaseflname[500];
       strcpy(gapbaseflname,outfileprefix);
       strcat(gapbaseflname,"_CHR");
       strcat(gapbaseflname,chr);
       strcat(gapbaseflname,"gapbase.dist");
       FILE * mutgapFile = fopen (gapbaseflname,"w");
       size_t gaptotal =0;
       size_t gapcount =0;
       double gapsdsum =0.0;

       for( map<int, int>::iterator ii=mutatiogapdist.begin(); ii!=mutatiogapdist.end(); ++ii)
       {
           gaptotal += (*ii).second;
           gapcount++;
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

    const char* genSeqErr ( char* seq)
    {
        int len = strlen(seq);
        char* sread = new char[len+1];
        strcpy(sread,seq);

        int pos = rand() % 36;
        char nucleotide[4] = {'A','C','G','T'};
        int randbase = rand() % 4;
        while(seq[pos] == nucleotide[randbase])
            randbase = rand() % 4;

        sread[pos] = nucleotide[randbase];
        sread[len]= '\0';
        return sread;
    }

    const char* genWildcardErr ( char* seq)
    {
        int len = strlen(seq);
        char* sread = new char[len+1];
        strcpy(sread,seq);
        int pos = rand() % 36;
        sread[pos] = 'N';
        sread[len]= '\0';
        return sread;
    }

    char* PolyBaseerror(char* seq)
    {
        char mutBases[4] = {'A', 'C', 'G', 'T'};
        char* sread = new char[strlen(seq)+1];
        strcpy(sread,seq);

        for (int i = 0; i < strlen(seq)-5; i++)
        {
             if(seq[i] == 'A' && seq[i+1] == 'A' && seq[i+2] == 'A' && seq[i+3] == 'A' && seq[i+4]== 'A')
            {
                char base =  mutBases[rand() % 4];
                while(base == 'A')
                    base = mutBases[rand() % 4];
                sread[i+4] = base;
            }
            else if(seq[i] == 'A' && seq[i+1] == 'A' && seq[i+2] == 'A' && seq[i+3] == 'A' && seq[i+4] != 'A')
            {
                sread[i+4] = 'A';
            }
       }
       sread[strlen(seq)] = '\0';
       return sread;
    }



   void genSimCGPairEndSNV()
   {
        boost::normal_distribution<> nd(PEGAPmean, PEGAPstdev);
        boost::mt19937 rng;
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);

        typedef unsigned int UINT;
        bool stop = false;
        size_t SRCount = 0;
        size_t RDCount = 0;


        bool mutation = true;
        bool strand = true;
        char* mutation_bases;
        char* pmutation_bases;
        map <size_t, string> insmap;
        map <size_t, char> snpmap;
        int *delarray;

        int *SR_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            SR_bases[i] = 0;

        int *RD_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            RD_bases[i] = 0;

        bool * fwd_random_bases = new bool[used];
        for (size_t i = 0; i < used; i++)
            fwd_random_bases[i] = false;

        bool * rev_random_bases = new bool[used];
        for (size_t i = 0; i < used; i++)
            rev_random_bases[i] = false;

        strcpy(outFileName,outfileprefix);
        strcat(outFileName,"CHR_");
        strcat(outFileName,chr);
        if(illumina)
             strcat(outFileName,".fastq");
        else
             strcat(outFileName,".fqx");

         FILE * OutFile = fopen(outFileName, "w");
           assert(OutFile != NULL);

        if(seed)
        {
            seednum =   time(NULL);
            srand ( seednum );
        }
        char nucleotides[4] = {'A','C','G','T'};

        if (mutation)
        {
            size_t snpcount = size_t(used * snppercent)/(size_t)(100);
            size_t inscount = size_t(used * inspercent)/(size_t)(100);
            size_t delcount = size_t(used * delpercent)/(size_t)(100);
            mutation_bases = new char[used];
            pmutation_bases = new char[used];
            char snpmutallele[256] = {0}, insmutallele[256] = {0},delmutallele[256] = {0};
            char snprefallele[256] = {0}, insrefallele[256] = {0},delrefallele[256] = {0};
            for (size_t i = 0; i < used; i++)
            {
                pmutation_bases[i] = 0;
                mutation_bases[i] = 0;
            }
           fprintf(stderr, "\n");

           while (snpcount > 0)
           {
               int irdx = 0;
               bool mut_ok = true;
               fprintf(stderr, "Mutating SNP %lu        \r", snpcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - (2*size+14+PEGAPmean); i <= mut_pos + (2*size+14+PEGAPmean); i++)
               {
                  if (bases[i] == 'N' || bases[i] == 'a' || bases[i] == 'c' || bases[i] == 'g' || bases[i] == 't')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int select = rand()%2;
               size_t randmnvsize = genUniformRandnumber(mnvmin,mnvmax);

               
               char* snpmutallele= new char[randmnvsize+1];
               char* snprefallele = new char[randmnvsize+1];
               int irefdx =0;
               int imutdx =0;


               // 
               
               if(select == 0)
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = 'S';
                     pmutation_bases[mut_pos+kk] = 'S';
                     SR_bases[mut_pos] = 1;
                     // snprefallele[irefdx++]= bases[mut_pos+kk];

                     int rspos = rand()%4;
                     while(nucleotides[rspos] == bases[mut_pos+kk])
                      rspos = rand()%4;
                      snpmap[mut_pos+kk]= nucleotides[rspos];
                     // isnpmutallele[imutdx++]= nucleotides[rspos];
                  }
               }
               else
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                    mutation_bases[mut_pos+kk] = 'S';
                    SR_bases[mut_pos] = 1;
                    int rspos = rand()%4; 
                    while(nucleotides[rspos] == bases[mut_pos+kk])
                      rspos = rand()%4;
                     snpmap[mut_pos+kk]= nucleotides[rspos];
                    // snpmutallele[imutdx++]= nucleotides[rspos];
                  }
               }
               //snpmutallele[randmnvsize] = '\0';
              // snprefallele[randmnvsize] = '\0';
               
               
               // snvoffmap.insert(pair<Key*, SNV*>(new Key(chr_num,mut_pos+1,'S'), new SNV(snprefallele,snpmutallele)));
               
               snpcount--;
          }
         fprintf(stderr, "\n");

         while (inscount > 0)
         {
               string insallele;
               bool mut_ok = true;
               fprintf(stderr, "Mutating Insert %lu        \r", inscount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;
               int select = rand()%2;


               size_t randinssize = genUniformRandnumber(insmin,insmax);
               for (int sz =0; sz < randinssize; sz++)
               {
                    insallele += nucleotides[rand()%4];
               }

               insmap[mut_pos] = insallele;


               if(select == 0)
               {
                  mutation_bases[mut_pos] = '+';
                  pmutation_bases[mut_pos] = '+';
               }
               else
                   mutation_bases[mut_pos] = '+';
               // SR_bases[mut_pos] = 1;
               inscount--;
        }
        fprintf(stderr, "\n");

        while (delcount > 0)
        {
               bool mut_ok = true;
               fprintf(stderr, "Mutating Delete %lu        \r", delcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int dsize =5;
               int select = rand()%2;

               size_t randdelsize = genUniformRandnumber(delmin,delmax);

               if(select == 0)
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                     pmutation_bases[mut_pos+kk] = '-';
                  }
               }
               else
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                  }
               }
               //SR_bases[mut_pos] = 1;
               delcount--;
        }
        fprintf(stderr, "\n");
      }
       ///////////////////////////////////////////////////////////////////////////////////
       size_t pos;
       int PEgap; 
       int psize = 2*size + PEGAPmean + 12;   
       size_t Totreads = size_t (used*coverage)/(size_t)(2*size);
       size_t numreads = Totreads;
       int gap1, gap2;
       while(numreads > 0)
       {
          PEgap = var_nor();
          if(pedistmap.count(PEgap) == 0)
              pedistmap[PEgap] = 1;
          else
              pedistmap[PEgap]++;

          strand = rand() % 2;
          int lf = rand() % 2;
          if(lf)
            left = true;
          else 
            left = false;
         
          if (random)
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

             psize = 2*size+ gap1+gap2+PEgap;
             
             pos = rand() % used;
             if (strand)
             {
                while ((pos >= (used - psize - 1)) || (fwd_random_bases[pos] == true))
                    pos = rand() % used;
             }
             else
             {
                while ((pos >= (used - psize - 1)) || (rev_random_bases[pos] == true))
                    pos = rand() % used;
             }
          }
         bool ok = true;
         char illqual[1024] = {0}, illread[1024] = {0};
         char illqual_rc[1024] = {0}, illread_rc[1024] = {0};
         char testread[1024] = {0};
         string sread, sillqual;
         int dels = 0;
         char snpallele[200];
         char delallele[200];
         char insallele[200];

         for (long j = 0; j < (psize + dels) && ok; j++)
         {
             if (pos + j >= (used - 1) || (islower(bases[pos+j])))
             {
                ok = false;
                break;
             }
             if ((mutation_bases == NULL) || (mutation_bases[pos + j] == 0))
             {
                switch (bases[pos + j])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    //case 'a':
                    //case 'c':
                    //case 'g':
                    //case 't':
                        illread[j] = bases[pos + j];
                        sread += bases[pos + j];
                        illqual[j] = '[';
                        sillqual += '[';
                        break;
                    default:
                        ok = false;
                        break;
                }
             }                                         // nucleotides
             else if(mutation_bases[pos + j] == 'S')
             {
                 mutchar = snpmap[pos+j];
                 //Key* key = new Key(chr_num,mut_pos,'S');
                 SR_bases[pos+j]++;
                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] = mutchar;
                        sread += mutchar;
                        illqual[j] = '#';
                        sillqual += '#';
                        // SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == 'S')
                        {
                            illread[j] =  mutchar;
                            sread += mutchar;
                            illqual[j] = '#';
                            sillqual += '#';
                            //SR_bases[pos+j]++;
                        }
                        else
                        {
                            illread[j] =  bases[pos+j];
                            sread += bases[pos+j];
                            illqual[j] = '[';
                            sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] = mutchar;
                     sread += mutchar;
                     illqual[j] = '#';
                     sillqual += '#';
                     // SR_bases[pos+j]++;
                 }
             }
             else if(mutation_bases[pos + j] == '-')
             {
                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] = mutation_bases[pos + j];
                        //test += mutation_bases[pos + j];
                        illqual[j] = '-';
                        sillqual += '-';
                        dels++;
                        // SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == '-')
                        {
                            illread[j] =  pmutation_bases[pos + j];
                            //test += pmutation_bases[pos + j];
                            illqual[j] = '-';
                            sillqual += '-';
                            dels++;
                        }
                        else
                        {
                           illread[j] =  bases[pos + j];
                           sread += bases[pos + j];
                           illqual[j] = '[';
                           sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] = mutation_bases[pos + j];
                     //test += mutation_bases[pos + j];
                     illqual[j] = '-';
                     sillqual += '-';
                     dels++;
                     //SR_bases[pos+j]++;
                 }
              }
              else if(mutation_bases[pos + j] == '+')
              {
                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] =  bases[pos + j];
                        sillqual += '[';
                        illread[j] = mutation_bases[pos + j];
                        sread += bases[pos + j];
                        sread += insmap[pos + j];
                        for (int sz =0; sz < insmap[pos + j].size(); sz++)
                        {
                            sillqual += '+';
                        }
                        illqual[j] = '+';
                        //sillqual += '+';
                        //SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == '+')
                        {
                            illread[j] =  bases[pos + j];
                            sillqual += '[';
                            illread[j] =  pmutation_bases[pos + j];
                            sread += bases[pos + j];
                            sread += insmap[pos + j];
                           for (int sz =0; sz < insmap[pos + j].size(); sz++)
                            {
                               sillqual += '+';
                            }
                            illqual[j] = '+';

                        }
                        else
                        {
                           illread[j] =  bases[pos + j];
                           sread += bases[pos + j];
                           illqual[j] = '[';
                           sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] =  bases[pos + j];
                     sread += bases[pos + j];
                     sillqual += '[';
                     illread[j] = mutation_bases[pos + j];
                     sread += insmap[pos + j];
                     for (int sz =0; sz < insmap[pos + j].size(); sz++)
                     {
                          sillqual += '+';
                     }
                     illqual[j] = '+';

                     //SR_bases[pos+j]++;
                 }
              }
        }


       ////////////////////////////////////////////////////////////////////////////

        if (ok)
        {
            // for (int k =0; k < psize; k++)
                // RD_bases[pos+k]++;

            int rlen = sread.size();
            // string tempqual;             need to modify quality score for delete
            
            strcpy(illqual,sillqual.c_str());
            
            for (int kk =0; kk< (rlen-psize+1); kk++)
            {
                 int dcount = 0;
                 int icount = 0;
                 int id = 0;
                 int inpos =0;
                 while(inpos < kk+1 )
                 {
                     if(illqual[id] == '[')
                        inpos++;
                     else if(illqual[id] == '-')
                        dcount++;
                     else if(illqual[id] == '+')
                        {icount--;inpos++;}
                     id++;
                 }

                 
            
                 strcpy(illread,substr(kk,psize+kk,sread.c_str()));
                 size_t npos =  pos + dcount+kk+icount;

                 makeReverse(illqual_rc, illqual);
                 makeReverseComplement(illread_rc, illread);

                 char* outread1 = new char[size+gap1+1];
                 memcpy(outread1,illread,size+gap1);
                 outread1[size+gap1] = '\0';

                 char* outread2 = new char[size+gap2+1];
                 memcpy(outread2,illread+size+gap1+PEgap,size+gap2);
                 outread2[size+gap2] = '\0';

                 char* outqual1 = new char[size+gap1+1];
                 memcpy(outqual1,illqual,size+gap1);
                 outqual1[size+gap1] = '\0';

                 char* outqual2 = new char[size+gap2+1];
                 memcpy(outqual2,illqual+size+gap1+PEgap,size+gap2);
                 outqual2[size+gap2] = '\0';

                 char* outread1_rc = new char[size+gap1+1];
                 memcpy(outread1_rc,illread_rc,size+gap1);
                 outread1_rc[size+gap1] = '\0';

                 char* outread2_rc = new char[size+gap2+1];
                 memcpy(outread2_rc,illread_rc+size+gap1+PEgap,size+gap2);
                 outread2_rc[size+gap2] = '\0';

                 char* outqual1_rc = new char[size+gap1+1];
                 memcpy(outqual1_rc,illqual_rc,size+gap1);
                 outqual1_rc[size+gap1] = '\0';

                 char* outqual2_rc = new char[size+gap2+1];
                 memcpy(outqual2_rc,illqual_rc+size+gap1+PEgap,size+gap2);
                 outqual2_rc[size+gap2] = '\0';

                 char* read1 = new char[size+1];
                 char* read2 = new char[size+1];

                 char* qual1 = new char[size+1];
                 char* qual2 = new char[size+1];

                 char* read1_rc = new char[size+1];
                 char* read2_rc = new char[size+1];

                 char* qual1_rc = new char[size+1];
                 char* qual2_rc = new char[size+1];
                 
                 int seq_idx1 =0;
                 int seq_idx2 =0;
                 int qual_idx1 =0;
                 int qual_idx2 =0;
                 int seq_idx1_rc =0;
                 int seq_idx2_rc =0;
                 int qual_idx1_rc =0;
                 int qual_idx2_rc =0;

                 for (int readIdx = 0; readIdx < 30 + gap1; readIdx++)
                 {
                        read1[seq_idx1++] = outread1[readIdx];
                        qual1[qual_idx1++] = outqual1[readIdx];
                        read1_rc[seq_idx1_rc++] = outread1_rc[readIdx];
                        qual1_rc[qual_idx1_rc++] = outqual1_rc[readIdx];
                        RD_bases[npos+readIdx]++;
                   // if ((left == false) && readIdx == 9)
                     //   readIdx += gap1;
                   // else if ((left == true) && readIdx == 19)
                        if(readIdx == 19)
                          readIdx += gap1;
                    
                 }

                 for (int readIdx = 0; readIdx < 30 + gap2; readIdx++)
                 {
                        read2[seq_idx2++] = outread2[readIdx];
                        qual2[qual_idx2++] = outqual2[readIdx];
                        read2_rc[seq_idx2_rc++] = outread2_rc[readIdx];
                        qual2_rc[qual_idx2_rc++] = outqual2_rc[readIdx];
                        RD_bases[npos+readIdx+PEgap]++;
                    //if ((left == false) && readIdx == 9)
                      //  readIdx += gap2;
                   // else if ((left == true) && readIdx == 19)
                       // readIdx += gap2;
                      if(readIdx == 9)
                          readIdx += gap2;
                 }

                 read1[size] = '\0';
                 read2[size] = '\0';
                 qual1[size] = '\0';
                 qual2[size] = '\0';

                 read1_rc[size] = '\0';
                 read2_rc[size] = '\0';
                 qual1_rc[size] = '\0';
                 qual2_rc[size] = '\0';
                 
                 if (strand == false)    // rev
                 {
                        rev_random_bases[npos] = true;
                        fprintf(OutFile,"@RD_CHR:%d_POS:%lu_STRAND:REV\n", chr_num, npos + 1);
                        fprintf(OutFile,"%s\t%s\n", read1_rc,read2_rc);
                        fprintf(OutFile,"%s\t%s\n", qual1_rc,qual2_rc);
                        //fprintf(OutFile,"%s\n%s\n", illread_rc,illqual_rc);
                  }
                  else    //fwd
                  {
                        fwd_random_bases[npos] = true;
                        fprintf(OutFile,"@RD_CHR:%d_POS:%lu_STRAND:FWD\n", chr_num,npos+1);
                        //fprintf(OutFile,"%s\n%s\n", illread, illqual);
                        fprintf(OutFile,"%s\t%s\n", read1,read2);
                        fprintf(OutFile,"%s\t%s\n", qual1,qual2);
                   }

               /* if (((strstr(illqual, "#") != NULL) || (strstr(illqual, "+") != NULL) || (strstr(illqual, "-") != NULL)|| mutation == false))
                {
                     if (strand == false)    // rev
                     {
                        rev_random_bases[npos] = true;
                        //makeReverse(illqual_rc, illqual);
                        //makeReverseComplement(illread_rc, illread);
                        //makeReverse(qual1_rc, qual1);
                        //makeReverse(qual2_rc, qual2);
                        //makeReverse(read1_rc, read1);
                        //makeReverse(read2_rc, read2);

                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:REV\n", chr_num, npos + 1);
                        //fprintf(OutFile,"%s\n%s\n", illread, illqual);
                       // fprintf(OutFile,"%s\t%s\n", outread1,outread2);
                        fprintf(OutFile,"%s\t%s\n", read1_rc,read2_rc);
                        fprintf(OutFile,"%s\t%s\n", qual1_rc,qual2_rc);
                        
                        //fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                       
                     }
                     else    //fwd
                     {
                        fwd_random_bases[npos] = true;
                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:FWD\n", chr_num,npos+1);
                        //fprintf(OutFile,"%s\n%s\n", illread, illqual);
                       // fprintf(OutFile,"%s\t%s\n", outread1,outread2);
                        fprintf(OutFile,"%s\t%s\n", read1,read2);
                        fprintf(OutFile,"%s\t%s\n", qual1,qual2);
                      //  fprintf(OutFile,"%s\n", sread.c_str());
                      }
                }
                else
                {
                     if (strand == false)    // rev
                     {
                         // makeReverse(illqual_rc, illqual);
                         // makeReverseComplement(illread_rc, illread);

                         //makeReverse(qual1_rc, qual1);
                         //makeReverse(qual2_rc, qual2);
                         //makeReverse(read1_rc, read1);
                         //makeReverse(read2_rc, read2);
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:REV\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\t%s\n", read1_rc,read2_rc);
                         fprintf(OutFile,"%s\t%s\n", qual1_rc,qual2_rc);
                       
                        // fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                     }
                     else    //fwd
                     {
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:FWD\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\t%s\n", read1,read2);
                         fprintf(OutFile,"%s\t%s\n", qual1,qual2);
                        // fprintf(OutFile,"%s\n+\n%s\n", illread, illqual);
                     }
                 }*/
            }
            read_count++;
            fprintf(stderr, "Read: %lu (%.2f%%)\r", read_count, double(read_count)/double(Totreads)*100);
        }
        numreads--;
     }

     if(dist)        // Implementing distribution
     {
        char SuppReadsFileName[500];
        char PrefixFilename[500];
        char ReadDensityFileName[500];

        strcpy(PrefixFilename,outfileprefix);
        strcat(PrefixFilename,"_CHR");
        strcat(PrefixFilename,chr);

        strcpy(SuppReadsFileName,PrefixFilename);
        strcat(SuppReadsFileName,"Supporting_Reads");

        strcpy(ReadDensityFileName,PrefixFilename);
        strcat(ReadDensityFileName,"ReadDensity");

        queue <int> diffqueue;
        FILE * SRFile = fopen(SuppReadsFileName, "w");
           assert(SRFile != NULL);
        size_t SRcount = 0;
        size_t RDcount = 0;
        size_t SRtotal = 0;
        size_t RDtotal = 0;
        double SRsdtotal = 0;
        double RDsdtotal = 0;

        for (int count = 0; count < used; count++)
        {
            if(SR_bases[count] != 0)
            {
                diffqueue.push(count);
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
                SRcount++;
                SRtotal += SR_bases[count];
                fprintf(SRFile, "%lu\t%c>%c\t%d\n", count + 1, bases[count], snpmap[count],SR_bases[count]);
            }
        }
        fclose(SRFile);
        SRmean = (double)SRtotal/(double)SRcount;
        // calculate sd
        for (int count = 0; count < used; count++)
        {
                SRsdtotal += (SRmean-SR_bases[count])*(SRmean-SR_bases[count]);
        }

        SRstdev = sqrt(SRsdtotal/SRcount);

        while(!diffqueue.empty())
        {
             diffqueue.pop();
        }

        FILE * RDFile = fopen(ReadDensityFileName, "w");
           assert(RDFile != NULL);

        for (int count = 0; count < used; count++)
        {
             RDcount++;
             RDtotal += RD_bases[count];
             fprintf(RDFile, "%lu\t%d\n", count + 1, RD_bases[count]);
        }
        fclose(RDFile);


        RDmean = (double)RDtotal/(double)RDcount;

        for (int count = 0; count < used; count++)
        {
                RDsdtotal += (RDmean-RD_bases[count])*(RDmean-RD_bases[count]);
        }

        RDstdev = sqrt(RDsdtotal/used);
     } 

     used = 0;
     free(bases);
     bases = NULL;

     delete [] fwd_random_bases;
     delete [] rev_random_bases;
     delete [] RD_bases;
     delete [] SR_bases;

     if (mutation_bases != NULL)
     {
        delete [] mutation_bases;
        mutation_bases = NULL;
     }

     if (pmutation_bases != NULL)
     {
        delete [] pmutation_bases;
        pmutation_bases = NULL;
     }
   }


   void genSimIlluminaSNV()
   {
        typedef unsigned int UINT;
        bool stop = false;
        size_t SRCount = 0;
        size_t RDCount = 0;

        bool mutation = true;
        bool strand = true;
        char* mutation_bases;
        char* pmutation_bases;
        map <size_t, string> insmap;
        int *delarray;

        int *SR_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            SR_bases[i] = 0;

        int *RD_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            RD_bases[i] = 0;

        bool * fwd_random_bases = new bool[used];
        for (size_t i = 0; i < used; i++)
            fwd_random_bases[i] = false;

        bool * rev_random_bases = new bool[used];
        for (size_t i = 0; i < used; i++)
            rev_random_bases[i] = false;

        strcpy(outFileName,outfileprefix);
        strcat(outFileName,"CHR_");
        strcat(outFileName,chr);
        if(illumina)
             strcat(outFileName,".fastq");
        else
             strcat(outFileName,".fqx");

         FILE * OutFile = fopen(outFileName, "w");
           assert(OutFile != NULL);

        if(seed)
        {
            seednum =   time(NULL);
            srand ( seednum );
        }
        char nucleotides[4] = {'A','C','G','T'};

        if (mutation)
        {
            size_t snpcount = size_t(used * snppercent)/(size_t)(100);
            size_t inscount = size_t(used * inspercent)/(size_t)(100);
            size_t delcount = size_t(used * delpercent)/(size_t)(100);
            mutation_bases = new char[used];
            pmutation_bases = new char[used];
            string delallele,mnvallele;

            for (size_t i = 0; i < used; i++)
            {
                pmutation_bases[i] = 0;
                mutation_bases[i] = 0;
            }
            fprintf(stderr, "\n");

           while (snpcount > 0)
           {
               bool mut_ok = true;
               fprintf(stderr, "Mutating SNP %lu        \r", snpcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int select = rand()%2;
               size_t randmnvsize = genUniformRandnumber(mnvmin,mnvmax);

               if(select == 0)
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = 'S';
                     pmutation_bases[mut_pos+kk] = 'S';
                     SR_bases[mut_pos] = 1;
                  }
               }
               else
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                    mutation_bases[mut_pos+kk] = 'S';
                    SR_bases[mut_pos] = 1;
                  }
               }
               snpcount--;
          }
         fprintf(stderr, "\n");

         while (inscount > 0)
         {
               string insallele;
               bool mut_ok = true;
               fprintf(stderr, "Mutating Insert %lu        \r", inscount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;
               int select = rand()%2;


               size_t randinssize = genUniformRandnumber(insmin,insmax);
               for (int sz =0; sz < randinssize; sz++)
               {
                    insallele += nucleotides[rand()%4];
               }

               insmap[mut_pos] = insallele;


               if(select == 0)
               {
                  mutation_bases[mut_pos] = '+';
                  pmutation_bases[mut_pos] = '+';
               }
               else
                   mutation_bases[mut_pos] = '+';
               SR_bases[mut_pos] = 1;
               inscount--;
        }
        fprintf(stderr, "\n");

        while (delcount > 0)
        {
               bool mut_ok = true;
               fprintf(stderr, "Mutating Delete %lu        \r", delcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int dsize =5;
               int select = rand()%2;

               size_t randdelsize = genUniformRandnumber(delmin,delmax);

               if(select == 0)
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                     pmutation_bases[mut_pos+kk] = '-';
                  }
               }
               else
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                  }
               }
               SR_bases[mut_pos] = 1;
               delcount--;
        }
        fprintf(stderr, "\n");
      }
        
       ///////////////////////////////////////////////////////////////////////////////////

       size_t pos;
       size_t Totreads = size_t (used*coverage)/(size_t)size;
       size_t numreads = Totreads;

       while(numreads > 0)
       {
          strand = rand() % 2;
          if (random)
          {
             pos = rand() % used;
             if (strand)
             {
                while ((pos >= (used - size - 1)) || (fwd_random_bases[pos] == true))
                    pos = rand() % used;
             }
             else
             {
                while ((pos >= (used - size - 1)) || (rev_random_bases[pos] == true))
                    pos = rand() % used;
             }
          }
         bool ok = true;
         char illqual[1024] = {0}, illread[1024] = {0};
         char illqual_rc[1024] = {0}, illread_rc[1024] = {0};
         char testread[1024] = {0};
         string sread, sillqual;
         int dels = 0;
         char snpallele[200];
         char delallele[200];
         char insallele[200];

         for (long j = 0; j < (size + dels) && ok; j++)
         {
             if (pos + j >= (used - 1))
             {
                ok = false;
                break;
             }
             if ((mutation_bases == NULL) || (mutation_bases[pos + j] == 0))
             {
                switch (bases[pos + j])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                        illread[j] = bases[pos + j];
                        sread += bases[pos + j];
                        illqual[j] = '[';
                        sillqual += '[';
                        break;
                    default:
                        ok = false;
                        break;
                }
             }                                         // nucleotides
             else if(mutation_bases[pos + j] == 'S')
             {
                 int rspos = rand()%4;
                 char mutchar;
                 while(toupper(nucleotides[rspos]) == toupper(bases[pos+j]))
                     rspos = rand()%4;
                 mutchar = nucleotides[rspos];
                 //Key* key = new Key(chr_num,mut_pos,'S');

                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] = mutchar;
                        sread += mutchar;
                        illqual[j] = '#';
                        sillqual += '#';
                        SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == 'S')
                        {
                            illread[j] =  mutchar;
                            sread += mutchar;
                            illqual[j] = '#';
                            sillqual += '#';
                            SR_bases[pos+j]++;
                        }
                        else
                        {
                            illread[j] =  bases[pos+j];
                            sread += bases[pos+j];
                            illqual[j] = '[';
                            sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] = mutchar;
                     sread += mutchar;
                     illqual[j] = '#';
                     sillqual += '#';
                     SR_bases[pos+j]++;
                 }
             }
             else if(mutation_bases[pos + j] == '-')
             {
                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] = mutation_bases[pos + j];
                        //test += mutation_bases[pos + j];
                        illqual[j] = '-';
                        sillqual += '-';
                        dels++;
                        //SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == '-')
                        {
                            illread[j] =  pmutation_bases[pos + j];
                            //test += pmutation_bases[pos + j];
                            illqual[j] = '-';
                            sillqual += '-';
                            dels++;
                        }
                        else
                        {
                           illread[j] =  bases[pos + j];
                           sread += bases[pos + j];
                           illqual[j] = '[';
                           sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] = mutation_bases[pos + j];
                     //test += mutation_bases[pos + j];
                     illqual[j] = '-';
                     sillqual += '-';
                     dels++;
                     //SR_bases[pos+j]++;
                 }
              }
              else if(mutation_bases[pos + j] == '+')
              {
                 if(hetero)
                 {
                     int schr = rand() % 2;
                     if(schr == 0)
                     {
                        illread[j] =  bases[pos + j];
                        sillqual += '[';
                        illread[j] = mutation_bases[pos + j];
                        sread += bases[pos + j];
                        sread += insmap[pos + j];
                        for (int sz =0; sz < insmap[pos + j].size(); sz++)
                        {
                            sillqual += '+';
                        }
                        illqual[j] = '+';
                        //sillqual += '+';
                        //SR_bases[pos+j]++;
                     }
                     else
                     {
                        if(pmutation_bases[pos + j] == '+')
                        {
                            illread[j] =  bases[pos + j];
                            sillqual += '[';
                            illread[j] =  pmutation_bases[pos + j];
                            sread += bases[pos + j];
                            sread += insmap[pos + j];
                           for (int sz =0; sz < insmap[pos + j].size(); sz++)
                            {
                               sillqual += '+';
                            }
                            illqual[j] = '+';

                        }
                        else
                        {
                           illread[j] =  bases[pos + j];
                           sread += bases[pos + j];
                           illqual[j] = '[';
                           sillqual += '[';
                        }
                     }
                 }
                 else
                 {
                     illread[j] =  bases[pos + j];
                     sread += bases[pos + j];
                     sillqual += '[';
                     illread[j] = mutation_bases[pos + j];
                     sread += insmap[pos + j];
                     for (int sz =0; sz < insmap[pos + j].size(); sz++)
                     {
                          sillqual += '+';
                     }
                     illqual[j] = '+';
                     //SR_bases[pos+j]++;
                 }
              }
        }


       ////////////////////////////////////////////////////////////////////////////

        if (ok)
        {
            for (int k =0; k < size; k++)
                RD_bases[pos+k]++;

            int rlen = sread.size();
            strcpy(illqual,sillqual.c_str());

            for (int kk =0; kk< (rlen-size+1); kk++)
            {
                 int dcount = 0;
                 int icount = 0;
                 int id = 0;
                 int inpos =0;
                 while(inpos < kk+1 )
                 {
                     if(illqual[id] == '[')
                        inpos++;
                     else if(illqual[id] == '-')
                        dcount++;
                     else if(illqual[id] == '+')
                        {icount--;inpos++;}
                     id++;
                 }

                 strcpy(illread,substr(kk,size+kk,sread.c_str()));
                 size_t npos =  pos + dcount+kk+icount;


                if (((strstr(illqual, "#") != NULL) || (strstr(illqual, "+") != NULL) || (strstr(illqual, "-") != NULL)|| mutation == false))
                {
                     if (strand == false)    // rev
                     {
                        rev_random_bases[npos] = true;
                        makeReverse(illqual_rc, illqual);
                        makeReverseComplement(illread_rc, illread);
                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:REV\n", chr_num, npos + 1);
                        fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                     }
                     else    //fwd
                     {
                        fwd_random_bases[npos] = true;
                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:FWD\n", chr_num,npos+1);
                        fprintf(OutFile,"%s\n+\n%s\n", illread, illqual);
                        fprintf(OutFile,"%s\n", sread.c_str());
                      }
                }
                else
                {
                     if (strand == false)    // rev
                     {
                         makeReverse(illqual_rc, illqual);
                         makeReverseComplement(illread_rc, illread);
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:REV\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                     }
                     else    //fwd
                     {
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:FWD\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\n+\n%s\n", illread, illqual);
                     }
                 }
            }
            read_count++;
            fprintf(stderr, "Read: %lu (%.2f%%)\r", read_count, double(read_count)/double(Totreads)*100);
        }
        numreads--;
     }

     if(dist)        // Implementing distribution
     {
        char SuppReadsFileName[500];
        char PrefixFilename[500];
        //char ReadDensityFileName[500];

        strcpy(PrefixFilename,outfileprefix);
        strcat(PrefixFilename,"_CHR");
        strcat(PrefixFilename,chr);

        strcpy(SuppReadsFileName,PrefixFilename);
        strcat(SuppReadsFileName,"Supporting_Reads");

        //strcpy(ReadDensityFileName,PrefixFilename);
        //strcat(ReadDensityFileName,"ReadDensity");

        queue <int> diffqueue;
        FILE * SRFile = fopen(SuppReadsFileName, "w");
           assert(SRFile != NULL);
        size_t SRcount = 0;
        size_t RDcount = 0;
        size_t SRtotal = 0;
        size_t RDtotal = 0;
        double SRsdtotal = 0;
        double RDsdtotal = 0;

        for (int count = 0; count < used; count++)
        {
            if(SR_bases[count] != 0)
            {
                diffqueue.push(count);
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
                SRcount++;
                SRtotal += SR_bases[count];
                fprintf(SRFile, "%lu\t%c>%c\t%d\n", count + 1, bases[count], mutation_bases[count],SR_bases[count]);
            }
        }
        fclose(SRFile);
        SRmean = (double)SRtotal/(double)SRcount;
        // calculate sd
        for (int count = 0; count < used; count++)
        {
                SRsdtotal += (SRmean-SR_bases[count])*(SRmean-SR_bases[count]);
        }

        SRstdev = sqrt(SRsdtotal/SRcount);

        while(!diffqueue.empty())
        {
             diffqueue.pop();
        }

        //FILE * RDFile = fopen(ReadDensityFileName, "w");
            //assert(RDFile != NULL);

        for (int count = 0; count < used; count++)
        {
             RDcount++;
             RDtotal += RD_bases[count];
             // fprintf(RDFile, "%lu\t%d\n", count + 1, RD_bases[count]);
        }
        //fclose(RDFile);


        RDmean = (double)RDtotal/(double)RDcount;

        for (int count = 0; count < used; count++)
        {
                RDsdtotal += (RDmean-RD_bases[count])*(RDmean-RD_bases[count]);
        }

        RDstdev = sqrt(RDsdtotal/used);
     }

     used = 0;
     free(bases);
     bases = NULL;

     delete [] fwd_random_bases;
     delete [] rev_random_bases;
     delete [] RD_bases;
     delete [] SR_bases;

     if (mutation_bases != NULL)
     {
        delete [] mutation_bases;
        mutation_bases = NULL;
     }

     if (pmutation_bases != NULL)
     {
        delete [] pmutation_bases;
        pmutation_bases = NULL;
     }
   }


   void genSimCGSNV()
   {
        typedef unsigned int UINT;
        bool stop = false;
        size_t SRCount = 0;
        size_t RDCount = 0;


        bool mutation = true;
        bool strand = true;
        char* mutation_bases;
        char* pmutation_bases;
        map <size_t, string> insmap;
        int *delarray;

        int *SR_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            SR_bases[i] = 0;

        int *RD_bases = new int[used];
        for (size_t i = 0; i < used; i++)
            RD_bases[i] = 0;

      

        strcpy(outFileName,outfileprefix);
        strcat(outFileName,"CHR_");
        strcat(outFileName,chr);
        if(illumina)
             strcat(outFileName,".fastq");
        else
             strcat(outFileName,".fqx");

         FILE * OutFile = fopen(outFileName, "w");
           assert(OutFile != NULL);

        if(seed)
        {
            seednum =   time(NULL);
            srand ( seednum );
        }
           
        char nucleotides[4] = {'A','C','G','T'};
        size_t nchrsize = strlen(bases);
        
        if (mutation)
        {
            size_t snpcount = size_t(used * snppercent)/(size_t)(100);
            size_t inscount = size_t(used * inspercent)/(size_t)(100);
            size_t delcount = size_t(used * delpercent)/(size_t)(100);
            mutation_bases = new char[used];
            pmutation_bases = new char[used];
            string delallele,mnvallele;

            for (size_t i = 0; i < used; i++)
            {
                pmutation_bases[i] = 0;
                mutation_bases[i] = 0;
            }
            
            fprintf(stderr, "\n");

           cout << snpcount << "\t" << inscount << "\t" << delcount << endl;

       
           
           while (snpcount > 0)
           {
               bool mut_ok = true;
               fprintf(stderr, "Mutating SNP %lu        \r", snpcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int select = rand()%2;
               size_t randmnvsize = genUniformRandnumber(mnvmin,mnvmax);
               
               if(select == 0)
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = 'S';
                     pmutation_bases[mut_pos+kk] = 'S';
                     SR_bases[mut_pos] = 1;
                  }
               }
               else
               {
                  for(int kk = 0; kk < randmnvsize; kk++)
                  {
                    mutation_bases[mut_pos+kk] = 'S';
                    SR_bases[mut_pos] = 1;
                  }
               }
               snpcount--;
          }
         fprintf(stderr, "\n");

         while (inscount > 0)
         {
               string insallele;
               bool mut_ok = true;
               fprintf(stderr, "Mutating Insert %lu        \r", inscount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;
               int select = rand()%2;


               size_t randinssize = genUniformRandnumber(insmin,insmax);
               nchrsize += randinssize;
               for (int sz =0; sz < randinssize; sz++)
               {
                    insallele += nucleotides[rand()%4];
               }

               insmap[mut_pos] = insallele;


               if(select == 0)
               {
                  mutation_bases[mut_pos] = '+';
                  pmutation_bases[mut_pos] = '+';
               }
               else
                   mutation_bases[mut_pos] = '+';
               //SR_bases[mut_pos] = 1;
               inscount--;
        }
        fprintf(stderr, "\n");

        while (delcount > 0)
        {
               bool mut_ok = true;
               fprintf(stderr, "Mutating Delete %lu        \r", delcount);
               size_t mut_pos = rand() % used;
               while ((mut_pos <= size_t(size)) || (mut_pos >= (used - size)) || (mutation_bases[mut_pos] != 0))
                 mut_pos = rand() % used;

               for (UINT i = mut_pos - size; i <= mut_pos + size; i++)
               {
                  if (bases[i] == 'N')
                  {
                     mut_ok = false;
                     break;
                  }
               }
               if (mut_ok == false)
                 continue;

               int dsize =5;
               int select = rand()%2;

               size_t randdelsize = genUniformRandnumber(delmin,delmax);

               nchrsize -= randdelsize;
               
               if(select == 0)
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                     pmutation_bases[mut_pos+kk] = '-';
                  }
               }
               else
               {
                  for(int kk = 0; kk < randdelsize; kk++)
                  {
                     mutation_bases[mut_pos+kk] = '-';
                  }
               }
               //SR_bases[mut_pos] = 1;
               delcount--;
        }
        fprintf(stderr, "\n");
      }
       ///////////////////////////////////////////////////////////////////////////////////
       size_t pos;
       size_t Totreads = size_t (used*coverage)/(size_t)size;
       size_t numreads = Totreads;
/////////////////////////////// Updating Chromosome  ///////////////////////////////////////////////////
       string chrombases;
       //char* chrombase = new char(nchrsize +1);
       size_t nrdx =0;
       cout << strlen(bases) << endl;
       for (int i =0; i < used; i++)
       {
           if(mutation_bases[i] == 'S')
           {
               int rspos = rand()%4;
               while(nucleotides[rspos] == bases[i])
                    rspos = rand()%4;
                 chrombases += nucleotides[rspos];
                 //chrombase[nrdx++] = nucleotides[rspos];
           }
           else if (mutation_bases[i] == '+')
           {
                  chrombases += bases[i];
                  //chrombase[nrdx++] = bases[i];
                  chrombases += insmap[i];
                  //strcat(chrombase,insmap[i].c_str());
                  //nrdx += strlen(insmap[i].c_str());
           }
           else if(mutation_bases[i] == 0)
           { 
               chrombases += bases[i];
               //chrombase[nrdx++] = bases[i];
           }
       }
       //chrombase[nchrsize] = '\0';
       cout << chrombases.size() << "\t"<<nchrsize << endl;

       // Mutated Genome sequence

        bool * fwd_random_bases = new bool[nchrsize];
        for (size_t i = 0; i < nchrsize; i++)
            fwd_random_bases[i] = false;

        bool * rev_random_bases = new bool[nchrsize];
        for (size_t i = 0; i < nchrsize; i++)
            rev_random_bases[i] = false;

        
       while(numreads > 0)
       {
          strand = rand() % 2;
          if (random)
          {
             pos = rand() % nchrsize;
             if (strand)
             {
                while ((pos >= (nchrsize - size - 1)) || (fwd_random_bases[pos] == true))
                    pos = rand() % nchrsize;
             }
             else
             {
                while ((pos >= (nchrsize - size - 1)) || (rev_random_bases[pos] == true))
                    pos = rand() % nchrsize;
             }
          }
         bool ok = true;
         char illqual[1024] = {0}, illread[1024] = {0};
         char illqual_rc[1024] = {0}, illread_rc[1024] = {0};
         char testread[1024] = {0};
         string sread, sillqual;
         int dels = 0;
     
         char snpallele[200];
         char delallele[200];
         char insallele[200];

         for (long j = 0; j < (size + dels) && ok; j++)
         {
             if (pos + j >= (nchrsize - 1))
             {
                ok = false;
                break;
             }
             if ((mutation_bases == NULL) || (mutation_bases[pos + j] == 0))
             {
                switch (chrombases[pos + j])
                {
                    case 'A':
                    case 'C':
                    case 'G':
                    case 'T':
                    case 'a':
                    case 'c':
                    case 'g':
                    case 't':
                        illread[j] = chrombases[pos + j];
                        sread += chrombases[pos + j];
                        illqual[j] = '[';
                        sillqual += '[';
                        break;
                    default:
                        ok = false;
                        break;
                }
             }                                         // nucleotides
        }


       ////////////////////////////////////////////////////////////////////////////

        if (ok)
        {
            for (int k =0; k < size; k++)
                RD_bases[pos+k]++;

            int rlen = sread.size();
            strcpy(illqual,sillqual.c_str());

            for (int kk =0; kk< (rlen-size+1); kk++)
            {
                 int dcount = 0;
                 int icount = 0;
                 int id = 0;
                 int inpos =0;
                 while(inpos < kk+1 )
                 {
                     if(illqual[id] == '[')
                        inpos++;
                     else if(illqual[id] == '-')
                        dcount++;
                     else if(illqual[id] == '+')
                        {icount--;inpos++;}
                     id++;
                 }

                 strcpy(illread,substr(kk,size+kk,sread.c_str()));
                 size_t npos =  pos + dcount+kk+icount;


                if (((strstr(illqual, "#") != NULL) || (strstr(illqual, "+") != NULL) || (strstr(illqual, "-") != NULL)|| mutation == false))
                {
                     if (strand == false)    // rev
                     {
                        rev_random_bases[npos] = true;
                        makeReverse(illqual_rc, illqual);
                        makeReverseComplement(illread_rc, illread);
                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:REV\n", chr_num, npos + 1);
                        fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                     }
                     else    //fwd
                     {
                        fwd_random_bases[npos] = true;
                        fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:MISMATCH_STRAND:FWD\n", chr_num,npos+1);
                        fprintf(OutFile,"%s\n+\n%s\n", illread, illqual);
                        fprintf(OutFile,"%s\n", sread.c_str());
                      }
                }
                else
                {
                     if (strand == false)    // rev
                     {
                         makeReverse(illqual_rc, illqual);
                         makeReverseComplement(illread_rc, illread);
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:REV\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\n+\n%s\n", illread_rc,illqual_rc);
                     }
                     else    //fwd
                     {
                         fprintf(OutFile,">RD_CHR:%d_POS:%lu_HIT:PERFECT_STRAND:FWD\n", chr_num, npos + 1);
                         fprintf(OutFile,"%s\n+\n%s\n", illread, illqual);
                     }
                 }
            }
            read_count++;
            fprintf(stderr, "Read: %lu (%.2f%%)\r", read_count, double(read_count)/double(Totreads)*100);
        }
        numreads--;
     }
/*
     if(dist)        // Implementing distribution
     {
        char SuppReadsFileName[500];
        char PrefixFilename[500];
        //char ReadDensityFileName[500];

        strcpy(PrefixFilename,outfileprefix);
        strcat(PrefixFilename,"_CHR");
        strcat(PrefixFilename,chr);

        strcpy(SuppReadsFileName,PrefixFilename);
        strcat(SuppReadsFileName,"Supporting_Reads");

        //strcpy(ReadDensityFileName,PrefixFilename);
        //strcat(ReadDensityFileName,"ReadDensity");

        queue <int> diffqueue;
        FILE * SRFile = fopen(SuppReadsFileName, "w");
           assert(SRFile != NULL);
        size_t SRcount = 0;
        size_t RDcount = 0;
        size_t SRtotal = 0;
        size_t RDtotal = 0;
        double SRsdtotal = 0;
        double RDsdtotal = 0;

        for (int count = 0; count < used; count++)
        {
            if(SR_bases[count] != 0)
            {
                diffqueue.push(count);
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
                SRcount++;
                SRtotal += SR_bases[count];
                fprintf(SRFile, "%lu\t%c>%c\t%d\n", count + 1, bases[count], mutation_bases[count],SR_bases[count]);
            }
        }
        fclose(SRFile);
        SRmean = (double)SRtotal/(double)SRcount;
        // calculate sd
        for (int count = 0; count < used; count++)
        {
                SRsdtotal += (SRmean-SR_bases[count])*(SRmean-SR_bases[count]);
        }

        SRstdev = sqrt(SRsdtotal/SRcount);

        while(!diffqueue.empty())
        {
             diffqueue.pop();
        }

        //FILE * RDFile = fopen(ReadDensityFileName, "w");
            //assert(RDFile != NULL);

        for (int count = 0; count < used; count++)
        {
             RDcount++;
             RDtotal += RD_bases[count];
             // fprintf(RDFile, "%lu\t%d\n", count + 1, RD_bases[count]);
        }
        //fclose(RDFile);


        RDmean = (double)RDtotal/(double)RDcount;

        for (int count = 0; count < used; count++)
        {
                RDsdtotal += (RDmean-RD_bases[count])*(RDmean-RD_bases[count]);
        }

        RDstdev = sqrt(RDsdtotal/used);
     }
*/
     used = 0;
     free(bases);
     bases = NULL;

     delete [] fwd_random_bases;
     delete [] rev_random_bases;
     delete [] RD_bases;
     delete [] SR_bases;

     if (mutation_bases != NULL)
     {
        delete [] mutation_bases;
        mutation_bases = NULL;
     }

     if (pmutation_bases != NULL)
     {
        delete [] pmutation_bases;
        pmutation_bases = NULL;
     }



   }

   void genSummary()
   {
        time_t rawtime;
        time ( &rawtime );

        char SummaryFileName[500];
        strcpy(SummaryFileName,outfileprefix);
        strcat(SummaryFileName,".Summary");
        FILE *SummaryFile = fopen (SummaryFileName,"w");

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
           fprintf (SummaryFile, "All\t%s\n",(all)?"true":"false" );
        }
        fprintf (SummaryFile, "Mutation Data\n");

        fprintf (SummaryFile, "Zygosity\t%s\n",(hetero)?"hetero":"homo" );
        fprintf (SummaryFile, "Seed\t%lu\n",(seed)?(unsigned long)seednum:1 );
        fprintf (SummaryFile, "Read Type\t%s\n",(pair_end)?"Pair End":"Single End" );
        fprintf (SummaryFile, "Sample Type\t%s\n",(illumina)?"Illumina":"CG" );
        fprintf (SummaryFile, "SNP Percentage\t%4.2f\n",snppercent);
        fprintf (SummaryFile, "Delete Precentage\t%4.2f\n",delpercent);
        fprintf (SummaryFile, "Insert Percentage\t%4.2f\n",inspercent);
        fprintf (SummaryFile, "CNV Percentage\t%4.2f\n",cnvpercent);
        fprintf (SummaryFile, "Tranlocation Precentage\t%4.2f\n",transpercent);
        fprintf (SummaryFile, "Inversion Percentage\t%4.2f\n",invpercent);
        fprintf (SummaryFile, "Minimum Delete Size \t%d\n",delmin);
        fprintf (SummaryFile, "Maximum Delete Size \t%d\n",delmax);
        fprintf (SummaryFile, "Minimum Insert Size \t%d\n",insmin);
        fprintf (SummaryFile, "Maximum Insert Size \t%d\n",insmax);
        fprintf (SummaryFile, "Minimum MNV Size \t%d\n",mnvmin);
        fprintf (SummaryFile, "Maximum MNV Size \t%d\n",mnvmax);
        fprintf (SummaryFile, "Minimum CNV Size \t%d\n",mincnv);
        fprintf (SummaryFile, "Maximum CNV Size \t%d\n",maxcnv);
        fprintf (SummaryFile, "Minimum Inversion Size \t%d\n",mininv);
        fprintf (SummaryFile, "Maximum Inversion Size \t%d\n",maxinv);
        fprintf (SummaryFile, "Minimum Translocation Size \t%d\n",mintrans);
        fprintf (SummaryFile, "Maximum Translocation Size \t%d\n",maxtrans);


        fprintf (SummaryFile, "Sequencing Stat\n");

        fprintf (SummaryFile, "Sequencing Error (illumina)\t%4.2f\n",serr);
        fprintf (SummaryFile, "Ambiguous Bases (N %)\t%4.2f\n",N);
        fprintf (SummaryFile, "PolyBase Error(454)\t%4.2f\n",pbaseerr);
        fprintf (SummaryFile, "Read Length\t%d\n",size);
        fprintf (SummaryFile, "Coverage\t%d\n",coverage);
        fprintf (SummaryFile, "Number of Reads\t%lu\n",(unsigned long) read_count);

        fprintf (SummaryFile, "SIM Density\n");

        if(dist)
        {
           fprintf (SummaryFile, "Mean Supporting Reads\t%4.2f\n",SRmean);
           fprintf (SummaryFile, "Mean Read Density\t%4.2f\n",RDmean);
           fprintf (SummaryFile, "Mean Pair End Gap\t%4.2f\n",PEGAPmean);
           fprintf (SummaryFile, "Mean Distance between SNPs\t%4.2f\n",Bases_apartmean);
           fprintf (SummaryFile, "Standard Deviation Supporting Reads\t%4.2f\n",SRstdev);
           fprintf (SummaryFile, "Standard Deviation Read Density\t%4.2f\n",RDstdev);
           fprintf (SummaryFile, "Standard Deviation Pair End Gap\t%4.2f\n",PEGAPstdev);
           fprintf (SummaryFile, "Standard Deviation SNP Distance\t%4.2f\n",Bases_apartstdev);
        }

        fprintf (SummaryFile, "\n");
        fclose (SummaryFile);
   }

   void generateDistribution()
   {
       for (SNVmap::iterator it = snvoffmap.begin(); it != snvoffmap.end(); it++)
       {
          cout << it->first->chr << "\t"<< it->first->offset << "\t" << it->first->snv_type<<"\t" << it->second->mutated_allele<<"\t" << it->second->original_allele << endl;
       }
   }

   void freeMemory()
   {
       pedistmap.clear();
       mutatiogapdist.clear();
       used = 0;
       Filevec.clear();
   }

   void run()
   {
        for (int i =0; i< Filevec.size(); i++)
        {
             strcpy(inputFileName,Filevec[i]);
             loadChrFile();
             // genSimIlluminaSNV();
              genSimCGPairEndSNV();
              // genSimCGSNV();
            // gen_gapbase_distribution();
           /*  if(pair_end == 1)
             {
                if(illumina == 1)
                  // genPairEndillumina();
                else
                   //genPairEndCG();
                if(dist)
                {
                    //gen_pairend_distribution();
                   // gen_gapbase_distribution();
                }
            }
            else
            {
                if(illumina == 1)
                   genSingleEndIllumina();
                else
                  // genSingleEndCG();

                if(dist)
                   // gen_gapbase_distribution();
            }*/
            genSummary();
        }
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

     //gs->genInsertillumina();
     //gs->genDeleteillumina();
     gs->run();
     gs->gen_pairend_distribution();
     gs->generateDistribution();
     //gs->loadChrFile();
     //gs->genSingleEndIllumina();
     gs->freeMemory();

     delete gs;
     return (EXIT_SUCCESS);
}


