#ifndef _TRANSALIGN_H_
#define _TRANSALIGN_H_
#include <string>
using namespace std;
#ifdef MAX
#undef MAX
#endif
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#define ALIGN_THRESHOLD -40

#define MAX_LINE_LENGTH 8192

// see if we can find additional alignments using best alignment scores + BEST_PLUS_SCORE
#define BEST_PLUS_SCORE -3

#define TABULATE_SCORE_DIST
#define TABULATE_IMPROVEMENTS

#if             defined(__GNUC__)
/*
   GCC fix.
   GCC has problems with placement of __attribute__() keywords.
   It is originally meant to be placed following any declarations or
   definitions.
   Specific declarations require them to be placed at the end (eg: class
   constructors & destructors).
   However, they must be placed between return type & function name for
   function definitions.
*/
#if             !defined(__attribute__)
#define         __attribute__(attr)
#endif       /* !defined(__attribute__) */
#endif       /* defined(__GNUC__) */

class bustardtoweight {
	static const int maxb = 30;
	int btow[2*maxb + 1];
	int notbtow[2*maxb + 1];
	float entropy[2*maxb + 1];
	float fastqentropy[2*maxb + 1];
	float normaliser;
	float psnp;

	int bustardtolog(int b) {
		double p;
		p = pow(10., b/10.);
		p = p/(1. + p);
		p += (1 - p) * psnp;
		double logp = log(p)/normaliser;
		return (int) (-0.5 + logp);
	}
	;

	int notbustardtolog(int b) {
		double p;
		p = pow(10., b/10.);
		p = p/(1. + p);
		p = (1.0 - p + psnp) / 3.0;
		return MAX(-10, (int)(-0.5 + log(p)/normaliser));
	}
	;

	float bustardtoentropy(int b) {
		double p;
		p = pow(10., b/10.);
		p = p/(1. + p);
		return p*log(1/p)/log(2);
	}
	;

	float bustardtofastqentropy(int b) {
		double p;
		p = pow(10., b/10.);
		p = p/(1. + p);
		return (p*log(1/p) + 3 * (1 - p)*log(1/(1-p)))/log(2);
	}
	;

public:
	bustardtoweight(float psnp) {
		this->psnp = psnp;
		normaliser = log(psnp + 0.001 - psnp * 0.001)/-10.0; // log(0.001)/10 gives score of -10 for good quality mismatch
		for (int i = -maxb; i <= maxb; i++) {
			btow[maxb + i] = bustardtolog(i);
			notbtow[maxb + i] = notbustardtolog(i);
			entropy[maxb + i] = bustardtoentropy(i);
			fastqentropy[maxb + i] = bustardtofastqentropy(i);
		}
	}

	int btowm(int bustard) {
		if (bustard > maxb) // Limit range to 30  to fit our -10 = good quality mismatch. Also fits with P(SNP) = 1/1000
			bustard = maxb;
		else if (bustard < -maxb)
			bustard = -maxb;
		return btow[maxb + bustard];
	}

	int notbtowm(int bustard) {
		if (bustard > maxb) // Limit range to 30  to fit our -10 = good quality mismatch. Also fits with P(SNP) = 1/1000
			bustard = maxb;
		else if (bustard < -maxb)
			bustard = -maxb;
		return notbtow[maxb + bustard];
	}

	float btoentropy(int bustard) {
		if (bustard > maxb) // Limit range to 30  to fit our -10 = good quality mismatch. Also fits with P(SNP) = 1/1000
			bustard = maxb;
		else if (bustard < -maxb)
			bustard = -maxb;
		return entropy[maxb + bustard];
	}

	float fastqtoentropy(int bustard) {
		if (bustard > maxb) // Limit range to 30  to fit our -10 = good quality mismatch. Also fits with P(SNP) = 1/1000
			bustard = maxb;
		else if (bustard < -maxb)
			bustard = -maxb;
		return fastqentropy[maxb + bustard];
	}

};

class logodds {
   static const int maxlg = 100;
   float normaliser;
   float p[maxlg+1];
   float cvrt(int logP) const {
	  return exp(logP * normaliser);
   }

public:
   void normalise(float psnp) {
	  normaliser = log(psnp + 0.001 - psnp * 0.001 )/-10.0;
	  for(int i = 0; i <= maxlg; i++)
		 p[i] = cvrt(-i);
   }
   float P(int logP) const {
	  // printf("normaliser: %.3f\n", normaliser);
	  if(-logP <= maxlg)
		 return p[-logP];
	  return cvrt(logP);
   }
};

class MyQueryIF : public SXInterface { /// Mechanism Synabase uses to report data.
private:
	OligoList Oligo;
	float sumPAi;
	int topscore;
	const char *queryid;
	const char *query;
	int querylen;

	logodds log2P;

	SXOligoSearch *os;
	SXSynabase *sxsb;

public:
	MyQueryIF(SXOligoSearch *osearch, SXSynabase *sxbase, float psnp);
	APIFUNCINTFATTRH APIFUNCINTFATTRB ~MyQueryIF() APIFUNCINTFATTRT;
	APIFUNCINTFATTRH bool APIFUNCINTFATTRB reportseqref() APIFUNCINTFATTRT;

	inline int setTopScore(int n) { topscore = n; return topscore; }
	inline int getTopScore() const { return topscore; }

	inline const char *setQueryID(const char *s) { queryid = s; return queryid; }
	inline const char *getQueryID() const { return queryid; }
	inline const char *setQuery(const char *s) { query = s; return query; }
	inline const char *getQuery() const { return query; }
	inline int setQueryLen(int n) { querylen = n; return querylen; }
	inline int getQueryLen() const { return querylen; }

	float setSumPAi(float f) { sumPAi = f; return sumPAi; }
	float getSumPAi() const { return sumPAi; }

	inline void ClearResults() { for (unsigned int i = 0; i < Oligo.size(); i++) delete Oligo[i]; Oligo.clear(); } // clear search results
	OligoList& getOligoList() { return Oligo; }

	logodds& getLog2P() { return log2P; }
};

class TransAlign {
public:

private:
	int max_align_threshold, mismatch_penalty, gap_penalty, gap_open, margin;
	double align_prob_threshold;

	SXSynabase *sxsb;
	SXOligoSearch *os;

	MyQueryIF *qryif;

	Dist *dist;
	const char *dist_fn, *hkg_fn, *nohit_fn;

	virtual int read() = 0;
	virtual void write_nohit() = 0;

protected:
	SXNACodes sxnacodes;
	SXOligoSearchWM fwdwm, revwm;
	FILE *fr, *fq;						// files for reading sequences / reads
	FILE *fnohit;

	char *query, *reverse;				// each contains one line of query and reverse read from file
	char *query_id;						// read / probe id
	size_t query_alloc,					// no. of bytes allocated to query
	query_len;						// no. of bytes read into queries
	char *qual, *rev_qual;
        string sampletype;				// each contains one line of quality information
	int threshold;
	float psnp, scoreratio;
	int maxNs;							// The maximum number of Ns letters in a read. Reads with more than this many Ns will not be searched.
	int minspec;
	int align_report_limit;

	int runSearch(MyQueryIF *qryif, SXOligoSearchWM& fwdwm, SXOligoSearchWM& revwm, int threshold);
	void report(const struct hit *h, double align_prob, FILE *out);

	inline static char *chomp(char *s) { char *p = strchr(s, '\r'); if (p) *p = '\0'; p = strchr(s, '\n'); if (p) *p = '\0'; return s; }
	inline static char *trim(char *s) { size_t len = strlen(s); while (len > 0 && isspace(s[len - 1])) s[--len] = '\0'; return s; }
	static inline char *reverse_read(char *str, size_t len) {
		for (size_t i = 0; i < len / 2; i++) {
			char c = str[i];
			str[i] = str[len - i - 1];
			str[len - i - 1] = c;
		}
		
		return str;
	}

public:
       
	TransAlign(const char *db);
	virtual ~TransAlign();
        int getendposition(char*,int);
        int matchscore(char*, char*);
        void concat(char* ,char);
        int getlength(int,char*);
        int getlenquery(char* );
        string Getformat(char* ,char* ,int );  
	const char *setDistFilename(const char *fn) { dist_fn = fn; return dist_fn; }
	int setMaxAlignThreshold(int n) { max_align_threshold = n; return max_align_threshold; }
	int getMaxAlignThreshold() const { return max_align_threshold; }
	int setMaxNs(int n) { maxNs = n; return maxNs; }
	int getMaxNs() const { return maxNs; }
	const char *setHkgFilename(const char *fn) { hkg_fn = fn; return hkg_fn; }
	const char *setNohitFilename(const char *fn) { nohit_fn = fn; return nohit_fn; }
	int setGapOpen(int v) { gap_open = v; return gap_open; }
	int getGapOpen() const { return gap_open; }
	int setGapPenalty(int v) { gap_penalty = v; return gap_penalty; }
	int getGapPenalty() const { return gap_penalty; }
	int setMargin(int v) { margin = v; return margin; }
	int getMargin() const { return margin; }
	int setMinSpec(int v) { minspec = v; return minspec; }
        void setSampleType(string v) { sampletype =v; }
        string  getSampleType() const  {return sampletype;}
	int getMinSpec() const { return minspec; }
	float setPSNP(float v) { psnp = v; return psnp; }
	float getPSNP() const { return psnp; }
	float setScoreRatio(float v) { scoreratio = v; return scoreratio; }
	float getScoreRatio() const { return scoreratio; }
	int setAlignReportLimit(int v) { align_report_limit = v; return align_report_limit; }
	int getAlignReportLimit() const { return align_report_limit; }
	void run(const char *fn, const char *qn = NULL);
};

/* FGETS reads a line from file and optionally terminates if EOF reached */
#define FGETS(buf,fil,fname,EOF_OK) if (fgets(buf, sizeof(buf), fil) == NULL) { if (!EOF_OK) { fprintf(stderr, "Unexpected end-of-file: %s\n", fname); exit(3); } }

#endif /* _TRANSALIGN_H_ */

