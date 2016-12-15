#ifndef _DIST_H_
#define _DIST_H_

#define DIST_SIZE 100000
#define DEFAULT_THRESHOLD 0.04
#define GENE_NAME_SIZE 101

#define FRAGMENT_BP_LEN 200
#include <iostream>
#include  <map>
using namespace std;

typedef struct {

size_t unique_count, threshold_count, balance_count, gene_len;
	t_foffset numberofsymbols;
	double align_prob_unique,		// 1.00
		   align_prob_threshold,	// 0.04
		   align_prob_balance;		// everything else
	double align_prob_sum, align_prob_max, align_prob_min;
	char gene_name[GENE_NAME_SIZE];
	bool hkg_flag;
} DistNode;

class Dist {
private:
  
 map<string, string> chromosomenum,FirstExpos;
        
	DistNode *dist_tab;
	size_t allocated, used;
	inline static char *trim(char *s) { size_t len = strlen(s); while (len > 0 && isspace(s[len - 1])) s[--len] = '\0'; return s; }
	DistNode *SearchAppend(const char *gene_n, bool hkg_fl = false);

public:
       
 
        string getchromosome(string);
        string getfirstexonspos(string);
        void LoadExDB(const char*);
	Dist(size_t size = 0);
	~Dist();
	bool Append(const char *gene_n, t_foffset no_symbols, double align_prob, double align_prob_threshold = DEFAULT_THRESHOLD);
	bool LoadHkgFile(const char *hkg_fn);
        string getid(string );  
	void Dump(FILE *f);
};



#endif /* _DIST_H_ */





