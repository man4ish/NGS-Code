#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <vector>
#include <limits>
#include <SXSynabase.h>
#include <NXDefs.h>
#include <fstream>
#include "dist.h"
#include <map>
using namespace std;

extern "C" int compar(const void *a, const void *b)
{
	return strcmp(((DistNode *)a)->gene_name, ((DistNode *)b)->gene_name);
}

extern "C" int compar2(const void *a, const void *b)
{

	if (((DistNode *)a)->align_prob_max != ((DistNode *)b)->align_prob_max) {
		return (((DistNode *)b)->align_prob_max - ((DistNode *)a)->align_prob_max) > 0.0 ? 1 : -1;
	}

	size_t total_count_a = ((DistNode *)a)->unique_count + ((DistNode *)a)->threshold_count + ((DistNode *)a)->balance_count;
	size_t total_count_b = ((DistNode *)b)->unique_count + ((DistNode *)b)->threshold_count + ((DistNode *)b)->balance_count;
	if (total_count_a > 0 && total_count_b > 0) {
		if (((DistNode *)a)->align_prob_sum / total_count_a != ((DistNode *)b)->align_prob_sum / total_count_b) {
			return (((DistNode *)b)->align_prob_sum / total_count_b - ((DistNode *)a)->align_prob_sum / total_count_a) > 0.0 ? 1 : -1;
		}
	}

	if (((DistNode *)a)->align_prob_min != ((DistNode *)b)->align_prob_min) {
		return (((DistNode *)b)->align_prob_min - ((DistNode *)a)->align_prob_min) > 0.0 ? 1 : -1;
	}

	double tfr_a = (double)((DistNode *)a)->numberofsymbols / (double)FRAGMENT_BP_LEN;
	double tfr_b = (double)((DistNode *)b)->numberofsymbols / (double)FRAGMENT_BP_LEN;
	if (tfr_a > 0.0 && tfr_b > 0.0) {
		return ((((DistNode *)b)->align_prob_unique + ((DistNode *)b)->align_prob_threshold + ((DistNode *)b)->align_prob_balance) / tfr_b -
			(((DistNode *)a)->align_prob_unique + ((DistNode *)a)->align_prob_threshold + ((DistNode *)a)->align_prob_balance) / tfr_a)  > 0.0 ? 1 : -1;
	}

	return 0;
}

Dist::Dist(size_t size) :
	allocated(0), used(0)
{
	if (size == 0)
		allocated = DIST_SIZE;
	else
		allocated = size;

	dist_tab = new DistNode[allocated];
	assert(dist_tab != NULL);
}
/////////////////////////////////////////////////////////////////////////////////////////////
 string Dist::getid(string s)
  {
  	
  	int y =s.find("|ref|");
	int z= s.size();
	int len = z-y-6;
	string geneid = s.substr(y+5,len) ;
	return geneid;
  } 
  void Dist::LoadExDB(const char* filename)
   {
    string line;
    ifstream myfile (filename);
    if (myfile.is_open())
     {
       while (! myfile.eof() )
        {
         getline (myfile,line);
       //  cout<<line;
         char *s = &line[0];
         int end = strlen(s) - 1;
         if (end >= 0 && s[end] == '\n')
         s[end] = '\0';
         const char* delimiter=",";
  	     char* pch;
  	     pch = strtok (s,delimiter);
  	     vector<string> a;
         while (pch != NULL)
          {
            a.push_back(pch);
            pch = strtok (NULL, delimiter);  
           } 
          chromosomenum[a[0]]= a[2];
          FirstExpos[a[0]]= a[1];
         }
         myfile.close();
      }
      
      else printf("%s\n", "No file found"); 
   }
   



 string Dist:: getchromosome(string refseqid)
   {
   	string a = chromosomenum[refseqid];
      if(a=="")
     return "XXXX";
  else 
        return a;
       
    }
   
   string Dist::getfirstexonspos(string refseqid)
    {
   	string b = FirstExpos[refseqid];
   if(b!="")
       return b;
    else 
      return "XXXX";
    }


////////////////////////////////////////////////////////////////////////////////////////////
Dist::~Dist()
{
	delete[] dist_tab;
}

DistNode *Dist::SearchAppend(const char *gene_n, bool hkg_fl)
{
	DistNode key, *found;
	strncpy(key.gene_name, gene_n, GENE_NAME_SIZE);
	key.gene_name[GENE_NAME_SIZE - 1] = '\0';
	found = (DistNode *)bsearch(&key, dist_tab, used, sizeof(DistNode), compar);
	if (found)
               {	
		return found;
        	}

	assert(used < allocated);

	// Append a new entry
	key.unique_count = 0;
	key.threshold_count = 0;
	key.balance_count = 0;
	key.align_prob_sum = 0.0;
	key.align_prob_unique = key.align_prob_threshold = key.align_prob_balance = 0.0;
	key.align_prob_max = std::numeric_limits<double>::min();
	key.align_prob_min = std::numeric_limits<double>::max();
	key.gene_len = 0;
	key.numberofsymbols = 0;
	key.hkg_flag = hkg_fl;

	memmove(&dist_tab[used++], &key, sizeof(key));

	// ... and sort it
	qsort(dist_tab, used, sizeof(DistNode), compar);

	// Look for it again
	found = (DistNode *)bsearch(&key, dist_tab, used, sizeof(DistNode), compar);
	assert(found);
	return found;
}

bool Dist::Append(const char *gene_n, t_foffset no_symbols, double align_prob, double align_prob_threshold)
{
	char gn[8192];
	size_t pos_gn = 0;
	while (*gene_n != ' ' && pos_gn < sizeof(gn) - 1)
		gn[pos_gn++] = *gene_n++;
	gn[pos_gn++] = '\0';

	DistNode *found = SearchAppend(gn);
	if (!found)
		return false;

	// Update entry
	found->numberofsymbols = no_symbols;
	if (align_prob == 1.0) {
		found->unique_count++;
		found->align_prob_unique += align_prob;
	}
	else if (align_prob <= align_prob_threshold) {
		found->threshold_count++;
		found->align_prob_threshold += align_prob;
	}
	else {
		found->balance_count++;
		found->align_prob_balance += align_prob;
	}

	found->align_prob_sum += align_prob;
	if (found->align_prob_max < align_prob)
		found->align_prob_max = align_prob;
	if (found->align_prob_min > align_prob)
		found->align_prob_min = align_prob;

	return true;
}

void Dist::Dump(FILE *f)
{
	
	qsort(dist_tab, used, sizeof(DistNode), compar2);
        LoadExDB("Exonlist.txt");
        // LoadExDB("latFirstExonlist.txt");
	int rc;
	rc = fprintf(f, "#target");
	rc = fprintf(f, "\tReadCount(RC1.00)");
	rc = fprintf(f, "\tProbabilitySum(PS1.00)");
	rc = fprintf(f, "\tReadCount(RC1.00..%.2f)", (float)DEFAULT_THRESHOLD);
	rc = fprintf(f, "\tProbabilitySum(PS1.00..%.2f)", (float)DEFAULT_THRESHOLD);
	rc = fprintf(f, "\tReadCount(RC1.00..0.00)");
	rc = fprintf(f, "\tProbabilitySum(PS1.00..0.00)");
	rc = fprintf(f, "\tTargetLength(TL)");
	rc = fprintf(f, "\tTarget:%dbpFragment(TFR)", FRAGMENT_BP_LEN);
	rc = fprintf(f, "\tRC(1.0/TFR)");
	rc = fprintf(f, "\tPS(1.0/TFR)");
	rc = fprintf(f, "\tRC(1.00..%.2f/TFR)", (float)DEFAULT_THRESHOLD);
	rc = fprintf(f, "\tPS(1.00..%.2f/TFR)", (float)DEFAULT_THRESHOLD);
	rc = fprintf(f, "\tRC(1.00..0.00/TFR)");
	rc = fprintf(f, "\tPS(1.00..0.00/TFR)");
	rc = fprintf(f, "\tPS(1.00/TL)");
	rc = fprintf(f, "\tPS(1.00..%.2f/TL)", (float)DEFAULT_THRESHOLD);
	rc = fprintf(f, "\tPS(1.00..0.00/TL)");
	rc = fprintf(f, "\tPmax");
	rc = fprintf(f, "\tPavg");
	rc = fprintf(f, "\tPmin");
	rc = fprintf(f, "\tHKG");
        rc = fprintf(f,"\tchromosome_no");
        rc = fprintf(f,"\tFirst_Exon_Position");
	putc('\n', f);

	for (size_t i = 0; i < used && rc >= 0; i++)
            {
		size_t total_count = dist_tab[i].unique_count + dist_tab[i].threshold_count + dist_tab[i].balance_count;
        	if (total_count == 0)
		continue;	// ignore records with no reads...
		double tfr = (double)dist_tab[i].numberofsymbols / (double)FRAGMENT_BP_LEN;

             string ge = dist_tab[i].gene_name;
             string gey = getid(ge);
             
           rc = fprintf(f, "%s", dist_tab[i].gene_name);
                
                rc = fprintf(f, "\t%lu", dist_tab[i].unique_count);

                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_unique);
                rc = fprintf(f, "\t%lu", dist_tab[i].unique_count + dist_tab[i].threshold_count);
                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold);
                rc = fprintf(f, "\t%lu", dist_tab[i].unique_count + dist_tab[i].threshold_count + dist_tab[i].balance_count);
                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold + dist_tab[i].align_prob_balance);

                rc = fprintf(f, "\t%lu", dist_tab[i].numberofsymbols);
                rc = fprintf(f, "\t%.4lf", tfr);

                rc = fprintf(f, "\t%.4lf", dist_tab[i].unique_count / tfr);
                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_unique / tfr);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].unique_count + dist_tab[i].threshold_count) / tfr);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold) / tfr);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].unique_count + dist_tab[i].threshold_count + dist_tab[i].balance_count) / tfr);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold + dist_tab[i].align_prob_balance) / tfr);

                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_unique / dist_tab[i].numberofsymbols);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold) / dist_tab[i].numberofsymbols);
                rc = fprintf(f, "\t%.4lf", (dist_tab[i].align_prob_unique + dist_tab[i].align_prob_threshold + dist_tab[i].align_prob_balance) / dist_tab[i].numberofsymbols);

                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_max);
                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_sum / total_count);
                rc = fprintf(f, "\t%.4lf", dist_tab[i].align_prob_min);
                rc = fprintf(f, "\t%s", dist_tab[i].hkg_flag ? "HKG" : "-");
                string ab = getchromosome(gey);
                string cd = getfirstexonspos(gey);
                char* aac = &ab[0];
                char* ccd = &cd[0];
                rc = fprintf(f, "\t%s", aac);
                rc = fprintf(f, "\t%s", ccd);
        
putc('\n', f);          
 }
         
	if (rc < 0)
               {
		fprintf(stderr, "Dist::Dump failed: %s\n", strerror(errno));
		return;
        	}
}

bool Dist::LoadHkgFile(const char *hkg_fn)
{
	FILE *f = fopen(hkg_fn, "r");
	if (!f) {
		fprintf(stderr, "Cannot find housekeeping gene file: %s\n", hkg_fn);
		return false;
	}


	DistNode key;
	key.unique_count = 0;
	key.threshold_count = 0;
	key.balance_count = 0;
	key.align_prob_unique = key.align_prob_threshold = key.align_prob_balance = 0.0;
	key.align_prob_sum = 0.0;
	key.align_prob_max = std::numeric_limits<double>::min();
	key.align_prob_min = std::numeric_limits<double>::max();
	key.gene_len = 0;
	key.numberofsymbols = 0;
	char line[1024];
	while (fgets(line, sizeof(line), f)) {
		char *p = key.gene_name, *l = line;
		while (*l && *l != '\t' && *l != '\n') {
			*p++ = *l++;
		}
		*p = '\0';

		assert(strlen(key.gene_name) < sizeof(key.gene_name));

		// skip to 3rd column...
		l++;
		l = strchr(l, '\t');
		l++;

		assert(*l != '\0');
		assert(*l != '\n');
		assert(*l == '-' || (*l == 'H' && *(l + 1) == 'K' && *(l + 2) == 'G'));

		assert(used < allocated);

		key.hkg_flag = *l != '-';

		// Append a new entry
		memmove(&dist_tab[used++], &key, sizeof(key));

	
	}

	fclose(f);

	// ... and sort it
	qsort(dist_tab, used, sizeof(DistNode), compar);



	return false;




}
















































