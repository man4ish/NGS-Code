#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <SXSynabase.h>
#include <SXOligoSearch.h>
#include <SXNACodes.h>
#include <assert.h>

#include "dist.h"
#include "hit.h"
#include "transalign.h"
#include "align_illumina.h"
using namespace std;
static const char *db = NULL;
static const char *seq_fn = NULL;
static int max_align_threshold = -20;
static int maxns = 6;
static float scoreratio = 0.0;
static const char *dist_fn = NULL;
static const char *hkg_fn = NULL;
static const char *nohit_fn = NULL;
static int gap_open = 12, gap_penalty = 5;
static int margin = 0;
static int minspec = 12;
static float psnp = 0.0;
static int align_report_limit = 10;
static string sampletyp;
static void banner(const char *progname)
{
	printf("SynaBASE Gene-expression mapping tool\n");
	printf("Usage: %s db FASTQ [optional parameters]\n", progname);
	printf("Required parameters:\n");
	printf("\tdb\t\t\tSynaBASE database name\n");
	printf("\tFASTQ\t\t\tInput FASTQ file\n");
	printf("Optional parameters:\n");
	printf("\t-x threshold,margin\tMaximum alignment threshold (default = 20,0)\n");
	printf("\t-dist file name\t\tOutput distribution table file (default is none)\n");
	printf("\t-N Ns\t\t\tMaximum number of Ns (default = 6)\n");
	printf("\t-rr prob\t\tProbability cut-off factor (default = %.3f)\n", scoreratio);
	printf("\t-o gap-open\t\tGap-open penalty (default = %d)\n", gap_open);
	printf("\t-g gap-extent\t\tGap-extent penalty (default = %d)\n", gap_penalty);
	printf("\t-s minspec\t\tMinimum non-Ns value (default = 12)\n");
	printf("\t-p PSNP\t\t\tPrior probability for SNPs (default = %.2f)\n", psnp);
	printf("\t-r report-limit\tLimit alignment reports (default = %d, 0 is unlimited)\n", align_report_limit);
	printf("\t-hkg file name\t\tInput housekeeping file (default is none)\n");
	printf("\t-nohit file name\tOutput file for un-aligned reads (default is none)\n");
        printf("\t-sample_type name\n");
	printf("\nExample:\n\t%s hsRefseqmRnar28sb309 S1.fastq -x 18 -hkg human.rna.fna.2col.len.s.hk2 -sample_type NM> S1.align.txt\n", progname);
}

static void parseCmdLine(int argc, char *argv[])
{
	putchar('#');
	for (int i = 0; i < argc; i++) {
		if (i > 0)
			putchar(' ');
		printf("%s", argv[i]);
	}
	putchar('\n');

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-dist") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -dist requires a filename.\n");
				exit(9);
			}

			dist_fn = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "-x") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -x requires a maximum alignment threshold value.\n");
				exit(9);
			}

			max_align_threshold = -abs(atoi(argv[++i]));
			const char *p = strchr(argv[i], ',');
			if (p)
				margin = atoi(p + 1);

			continue;
		}

		if (strcmp(argv[i], "-N") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -N requires a maximum N value.\n");
				exit(9);
			}

			maxns = abs(atoi(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-p") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -p requires a prior probability for SNPs.\n");
				exit(9);
			}

			psnp = fabs(atof(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-r") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -r requires a alignment-report limit value.\n");
				exit(9);
			}

			align_report_limit  = fabs(atof(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-rr") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -rr requires a probability cut-off factor value.\n");
				exit(9);
			}

			scoreratio  = fabs(atof(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-hkg") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -hkg requires a filename.\n");
				exit(9);
			}

			hkg_fn = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "-o") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -o requires a gap open value.\n");
				exit(9);
			}

			gap_open = abs(atoi(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-g") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -o requires a gap extension value.\n");
				exit(9);
			}

			gap_penalty = abs(atoi(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-s") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -o requires a minspec value.\n");
				exit(9);
			}

			minspec = abs(atoi(argv[++i]));
			continue;
		}

		if (strcmp(argv[i], "-nohit") == 0) {
			if (i >= argc - 1) {
				fprintf(stderr, "Parameter -nohit requires a filename.\n");
				exit(9);
			}

			nohit_fn = argv[++i];
			continue;
		}
              
                if (strcmp(argv[i], "-sample_type") == 0) {
                        if (i >= argc - 1) {
                                fprintf(stderr, "Parameter -sample_type requires a sample_type.\n");
                                exit(9);
                        }

                        sampletyp = argv[++i];
                        continue;
                }


		if (db == NULL)
			db = argv[i];
		else if (seq_fn == NULL)
			seq_fn = argv[i];
		else {
			fprintf(stderr, "Invalid parameter: %s\n", argv[i]);
		}
	}
}

static void TransAlignFASTQ(const char *db, const char *fn)
{
	Align_Illumina p(db);
	p.setMaxNs(maxns);
	p.setDistFilename(dist_fn);
	p.setHkgFilename(hkg_fn);
	p.setGapOpen(gap_open);
	p.setGapPenalty(gap_penalty);
	p.setMaxAlignThreshold(max_align_threshold);
	p.setMargin(margin);
	p.setNohitFilename(nohit_fn);
	p.setPSNP(psnp);
	p.setScoreRatio(scoreratio);
	p.setAlignReportLimit(align_report_limit);
        p.setSampleType(sampletyp);
	p.run(fn);
}

int main(int argc, char *argv[])
{
	parseCmdLine(argc, argv);
	if (db == NULL) {
		banner(argv[0]);
		return 1;
	}

	TransAlignFASTQ(db, seq_fn);
	
	return 0;
}
