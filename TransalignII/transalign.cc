#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <string>
#include <NXDefs.h>
#include <SXSynabase.h>
#include <SXOligoSearch.h>
#include <SXNACodes.h>
#include <assert.h>
#include <map>
#include "dist.h"
#include "hit.h"
#include "transalign.h"
using namespace std;
MyQueryIF::MyQueryIF(SXOligoSearch *osearch, SXSynabase *sxbase, float psnp) : SXInterface()
{
	assert(Pattern != NULL);
	sxsb = sxbase;
	os = osearch;
	topscore = -9999;
	log2P.normalise(psnp);
	queryid = NULL;
}

APIFUNCINTFATTRH APIFUNCINTFATTRB
MyQueryIF::~MyQueryIF()
APIFUNCINTFATTRT
{
}

APIFUNCINTFATTRH bool APIFUNCINTFATTRB
MyQueryIF::reportseqref() // Synabase calls this to report a sequence reference.
APIFUNCINTFATTRT
{
	t_offset offset = SXSeqRefGet()->Offset;
	int score = (int)SXSeqRefGet()->wmscore;
	if (score > topscore)
		topscore = score;
	Oligo.push_back(new hit(sxsb, queryid, query, score, offset, Pattern->PatternLen, os->getalignment()));
	sumPAi += log2P.P(score);
	return true;
};

TransAlign::TransAlign(const char *db) :
	sxsb(NULL), os(NULL), qryif(NULL),
	query(NULL), reverse(NULL),
	query_alloc(0), query_len(0),
	psnp(0.0), scoreratio(0.0),
	qual(NULL),
	rev_qual(NULL),
	maxNs(6), minspec(12),
	max_align_threshold(-20),
	mismatch_penalty(10),
	gap_penalty(5),
	gap_open(12),
	margin(0),
	dist(NULL),
	dist_fn(NULL), hkg_fn(NULL), nohit_fn(NULL),
	fnohit(NULL)
{
	try {
		sxsb = new SXSynabase(db, NXRead);
	}
	catch(std::exception &e)
	{
		fprintf(stderr, "\n%s\n", e.what());
		return;
	}
	catch (...) {
		fprintf(stderr, "\nError: Unable to open database '%s'", db);
		return;
	}
}

TransAlign::~TransAlign()
{
	if (qryif)
		delete qryif;
	if (os)
		delete os;
	if (sxsb)
		delete sxsb;
}

int TransAlign:: matchscore(char* que ,char* tar) {
 int match;
 for(int i=0;i<strlen(que);i++) {
    if(que[i]==tar[i])
   match++;
   }
  return match;
}


int TransAlign::runSearch(MyQueryIF *qryif, SXOligoSearchWM& fwdwm, SXOligoSearchWM& revwm, int threshold)
{
	qryif->ClearResults();

	qryif->setTopScore(-9999);

#if LOGFILE
FILE *logfile;
#endif

	int _max = std::numeric_limits<int>::max();
	for (;;) {
		// fprintf(stderr, "Align threshold: %d, max align threshold: %d\n", threshold, max_align_threshold);
		qryif->setSumPAi(0.0);
		if (threshold <= fwdwm.MinScore()) {
			os->setWeightMatrix(&fwdwm);
			_end = F;
			qryif->setQuery(query);
			qryif->setQueryLen(query_len);

                      #if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d Search #1 threshold = %d, _max = %d [%s:%s]\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec,
    threshold, _max, query, qual);
  fclose(logfile);
}
#endif

			os->Search(sxsb, qryif, threshold, _max);

#if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d Search #1 threshold = %d, _max = %d [%s:%s] completed\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec,
    threshold, _max, query, qual);
  fclose(logfile);
}
#endif

			os->setWeightMatrix(&revwm);
			_end = R;
			qryif->setQuery(reverse);

#if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d Search #2 threshold = %d, _max = %d [%s:%s]\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec,
    threshold, _max, query, qual);
  fclose(logfile);
}
#endif

			os->Search(sxsb, qryif, threshold, _max);

#if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d Search #2 threshold = %d, _max = %d [%s:%s] completed\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec,
    threshold, _max, query, qual);
  fclose(logfile);
}
#endif

		}
		// check SXOligoSearch results...
		const OligoList& Oligo = qryif->getOligoList();
		if (Oligo.size() > 0) {
			// printf("#1 Oligo.size: %lu\n", Oligo.size());
			if (margin == 0 || qryif->getTopScore() - margin >= threshold)
				break;
			threshold = qryif->getTopScore() - margin;

#if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d Before ClearResult()\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec);
  fclose(logfile);
}
#endif

			qryif->ClearResults();

#if LOGFILE
logfile = fopen("transalign.log", "a");
if (logfile) {
  time_t now;
  time(&now);
  struct tm *t = localtime(&now);
  fprintf(logfile, "%2d/%02d/%04d %2d:%02d:%02d After ClearResult()\n",
    t->tm_mday, t->tm_mon + 1, t->tm_year + 1900, t->tm_hour, t->tm_min, t->tm_sec);
  fclose(logfile);
}
#endif

			continue;
		}

		if (threshold <= max_align_threshold)
			break;

		threshold -= mismatch_penalty;
	}

	return threshold;
}
int TransAlign::getlenquery(char * query){
  int lenq = 0;
  while(*query) {
   if(*query != '-') {
     lenq++;
      }
     *query++;
    
   }
return lenq;
}


int TransAlign::getendposition(char* query,int stpos)
 {
  int lenq = stpos -1;
  while(*query) {
   if(*query != '-') {
     lenq++;
      }
     *query++;
   }
return lenq;
 }

int TransAlign::getlength(int h ,char* t)
{  
    int sum =0;
   
     for(int k =1;k<=h;k++)
     {
         if(*(t+k)!='-')
          {         
            sum += 1;                   
          }         
     }
    
     return sum;
};
void TransAlign::concat( char *c, char b) {
                int l = strlen(c);
                c[l] = b;
                c[l+1] = '\0';
        }


void TransAlign::report(const struct hit *h, double align_prob, FILE *out)
{
        string val;
	struct SXSeqHdr sqhdr;
	char q[MAX_LINE_LENGTH], t[MAX_LINE_LENGTH];
	char sequencehdr[MAX_LINE_LENGTH];
	int i;
	const char *qry, *tgt;
	size_t from, to;
	double lowscore = qryif->getLog2P().P(qryif->getTopScore()) * scoreratio / qryif->getSumPAi();
	double Pa = qryif->getLog2P().P(h->score) / qryif->getSumPAi();
	if (Pa < lowscore)
		return;
	sxsb->SXSequenceHdr(h->seqno, &sqhdr);
	qry = h->query;
	tgt = h->seq;
	if (isalpha(sqhdr.sequencehdr[0]))
		strcpy(sequencehdr, sqhdr.sequencehdr);
	else
		strcpy(sequencehdr, sqhdr.sequencehdr + 1);
	char *s = sequencehdr + strlen(sequencehdr);
	while (s != sequencehdr) {
		if (*s == '|') {
			*++s = '\0';
			break;
		}
		s--;
	}

        fprintf(out, ">%s\t", sequencehdr);
 	int gap_open_count = 0, gap_count = 0;
	int mismatches = 0;
	bool gap = false;
	bool gap_open_flag = false;
	t[0] = q[0] = '\0';
	for(i = 0; h->align[i] != '\0'; i++) {
		switch(h->align[i]) {
			case '-':
				if (!gap_open_flag) {
					gap_open_count++;
					gap_open_flag = true;
				}

				t[i] = '-';
				q[i] = *qry++;
				gap = true;
		        	gap_count++;
				break;

			case '+':
				if (!gap_open_flag) {
				gap_open_count++;
					gap_open_flag = true;
				}

				t[i] = *tgt++;
				q[i] = '-';
				gap = true;
				gap_count++;
				break;
			default:
				if (*tgt != *qry)
				mismatches++;

				t[i] = *tgt++;
				q[i] = *qry++;
				gap_open_flag = false;
				break;
		}
	}
	t[i] = q[i]  = '\0';

	from = h->seqoff;
	to = from + h->len - 1;

fprintf(out, "%lu\t%lu\t", from,to);
fprintf(out, "%s\t", h->queryid);	
  int pos;     
  int matchcount =0,qcount =0;
  int start = from;
  if(t[0] == '-') 
  start = start -1;
  bool insertionflag,deletionflag,mismatchflag;
  vector<string > vout;
  int  insum,delsum,missum;
         for(int i =0;i<strlen(q);i++){
                 if(*(q+i)=='-') 
                     {
                       delsum=0; deletionflag=0; mismatchflag=0;insertionflag=1;
                       pos=start+getlength(i,t);   
                       char b1[10],c1[10];
                       sprintf(b1,"%d",pos);
                       string posi,a5,a6;
                       char typ1[200]= "\0";                                              
                    strcat(typ1,b1);
                    concat(typ1,'.');
                    concat(typ1,'-');
                    concat(typ1,*(t+i));
                     vout.push_back(typ1);
                     }
                 else if(*(t+i)=='-')
                 {  
                     if(deletionflag ==1)
                       delsum +=1;
                        else delsum=0;
                        deletionflag=1;mismatchflag=0;insertionflag=0;
                    pos=start+getlength(i,t);
                     char b3[10],c3[10];
                       sprintf(b3,"%d",pos);
                       sprintf(c3,"%d",delsum);
                       char typ3[200]= "\0";
                   if(delsum!=0){
                    strcat(typ3,b3);
                    concat(typ3,'.');
                    strcat(typ3,c3);
                    concat(typ3,*(q+i));
                    concat(typ3,'-');
                     vout.push_back(typ3);
                        }
                     else {
                    strcat(typ3,b3);
                    concat(typ3,'.');
                    concat(typ3,*(q+i));
                    concat(typ3,'-');
                     vout.push_back(typ3);
                     }
                 }
                 else if(*(q+i)!=(*(t+i))) {
                        mismatchflag = 1; delsum=0; insertionflag=0; deletionflag=0;
                        pos=start+getlength(i,t);
                        char b2[10],c2[10];
                        sprintf(b2,"%d",pos);
                        char typ2[200]= "\0"; 
                       strcat(typ2,b2);
                       concat(typ2,'.');
                       concat(typ2,*(q+i));
                       concat(typ2,*(t+i));
                       vout.push_back(typ2);
                 }
                 else {  mismatchflag=0; delsum=0; insertionflag=0; deletionflag=0;matchcount++; }
         }
int lencount;
lencount = getlenquery(q);
int orientation;
if(h->dir=='+')
  {
   orientation = 1;
   fprintf(out, "%lu\t%d\t%d\t", (size_t)1,lencount,lencount);
  }
else {
orientation =0;
fprintf(out, "%d\t%lu\t%d\t", lencount ,(size_t)1,lencount);
}

fprintf(out, "%d\t%d\t%8.6f\t%d\t", h->score,orientation, align_prob, matchcount);
string t5;
if(vout.size()!=0) {
    for( int yy =0; yy< vout.size()-1;yy++) {
       t5 = vout[yy];
       fprintf(out,"%s%s ",t5.c_str(),", ");
     }
  string t6;
t6 = vout[vout.size()-1];
fprintf(out,"%s\t",t6.c_str());
  }
string testrresult = getSampleType();
fprintf(out,"%c\t%s\n",'T',&testrresult[0]);
fflush(out);
}

void TransAlign::run(const char *fn, const char *qn)
{
	try  {
		os = new SXOligoSearch(mismatch_penalty, 3, gap_penalty, gap_open);     // 3 is penalty for ambiguous codes in the reference sequence.
		qryif = new MyQueryIF(os, sxsb, psnp);
		qryif->getLog2P().normalise(psnp);
	}
	catch(std::exception &e)
	{
		fprintf(stderr, "\n%s\n", e.what());
		return;
	}
	catch (...) {
		fprintf(stderr, "\nError: Unable to open database");
		return;
	}

	fr = fopen(fn, "r");
	if (fr == NULL) {
		fprintf(stderr, "Unable for open file %s for reading\n", fn);
		exit(1);
	}
	if (qn) {
		fq = fopen(qn, "r");
		if (fq == NULL) {
			fprintf(stderr, "Unable for open file %s for reading\n", qn);
			exit(1);
		}
	}

	if (dist_fn)
		dist = new Dist;

	if (dist && hkg_fn)
		dist->LoadHkgFile(hkg_fn);

	// Output "no-hit" to file?
	if (nohit_fn)
		fnohit = fopen(nohit_fn, "w");

	size_t unique_count = 0;
	size_t align_count = 0;
	size_t multiple_align_count = 0;
	size_t usable_count = 0;

	time_t t_start;
	time(&t_start);

	for (;;) {
		// call the virtual function read() to read sequences and qualities into query_id, query and qual
		int read_result = read();
		if (read_result < 0)
			break;

		usable_count += read_result;

		try {
			long int from, to;
			t_offset ofs;

			qryif->setQueryID(query_id);

			threshold = runSearch(qryif, fwdwm, revwm, max_align_threshold);
			align_count++;

			OligoList& Oligo = qryif->getOligoList();

			if (Oligo.size() == 0) {
				write_nohit();
				continue;
			}

			if (Oligo.size() == 1)
				unique_count++;
			else if (Oligo.size() > 1)
				multiple_align_count++;
			fflush(stdout);

			int report_count = 0;
			for (OligoList::const_iterator it = Oligo.begin(); it != Oligo.end(); it++) {
				double align_prob = qryif->getLog2P().P((*it)->score) / qryif->getSumPAi();

				if (report_count++ < align_report_limit || align_report_limit == 0)
					report((*it), align_prob, stdout);

				if (dist) {
					struct SXSeqHdr sqhdr;

					sxsb->SXSequenceHdr((*it)->seqno, &sqhdr);
					if (isalpha(sqhdr.sequencehdr[0]))
						dist->Append(sqhdr.sequencehdr, sqhdr.numberofsymbols, align_prob);
					else
						dist->Append(sqhdr.sequencehdr + 1, sqhdr.numberofsymbols, align_prob);
				}
			}
		}
		catch(std::exception &e)
		{
			fprintf(stderr, "\n%s\n", e.what());
			return;
		}
	}

	if (fr) {
		fclose(fr);
		fr = NULL;
	}
	if (fq) {
		fclose(fq);
		fq = NULL;
	}

	if (fnohit) {
		fclose(fnohit);
		fnohit = NULL;
	}

	time_t t_stop;
	time(&t_stop);
	double t_diff = difftime(t_stop, t_start);


       if (dist) {
		if (dist_fn) {
			FILE *fdist = fopen(dist_fn, "w");
			if (fdist) {
				dist->Dump(fdist);
				fclose(fdist);
			}
		}

		delete dist;
		dist = NULL;
	}
}

