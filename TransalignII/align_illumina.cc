#include <math.h>

#include <vector>

#include <SXOligoSearch.h>
#include <SXNACodes.h>
#include <assert.h>

#include "dist.h"
#include "hit.h"
#include "transalign.h"
#include "align_illumina.h"

Align_Illumina::Align_Illumina(const char *db)
	: TransAlign(db)
{
	query_len = 0;
	query_alloc = 1024;
	query_id = new char[query_alloc];
	query = new char[query_alloc];
	qual = new char[query_alloc];
	rev_qual = new char[query_alloc];
	reverse = new char[query_alloc];
	lineno = 0;
}

Align_Illumina::~Align_Illumina()
{
	delete[] reverse;
	delete[] rev_qual;
	delete[] qual;
	delete[] query;
	delete[] query_id;
}

/* Sample data
@SLXA-B7_534_1_1_347_257
CCCCCCCCCCACCACACTTACATACAACATCCTAT
+SLXA-B7_534_1_1_347_257
^^^^^^^^^^^^PUHDD^NUCVPGCJCCGJCQJFJ
@SLXA-B7_534_1_1_791_95
AATTCAGTTGGGATTTTCAATAAAAAAAGTACAAC
+SLXA-B7_534_1_1_791_95
^^^^^^^^^^^^^^^^^^^^^^^[^C^^[XWQSQM
@SLXA-B7_534_1_1_852_118
GTGATCATTTTTCATATGCATCGCAAAAATCACTC
+SLXA-B7_534_1_1_852_118
^^^^^^^^^^^^^[^Y^^^U^^^^VJPQJXWQAQF
@SLXA-B7_534_1_1_615_252
CCCCCCCCCCCCCCCCCCAAAAACAAAAAAACAAA
+SLXA-B7_534_1_1_615_252
^^^^^^^^^^^^^^^^QPHHQQ^CV^JQPSJCOJJ
@SLXA-B7_534_1_1_810_120
TTTTTTATAGAAAAAAAGCAATTTCTCTCACTCAA
+SLXA-B7_534_1_1_810_120
^^^^^^^^^^^^^^^^Q^[^^^^^^^^^[GW[UQQ
*/

/* Pair_Illumina::setSXOligoSearchWM
 * Parameters:
 *   which - 1 = "left", 2 = "right"
 */
bool Align_Illumina::setSXOligoSearchWM(SXOligoSearchWM& fwdwm, SXOligoSearchWM& revwm, char *read_id, char *query, char *reverse, char *quality, size_t query_len)
{
	float readentropy = 0.0;
	int lastgood = -1, trailingNs;
	
	bustardtoweight bw(psnp);
	
	fwdwm.setQuerySequence(query);
	fwdwm.LengthSet(query_len);
	// printf("Pair_Illumina set length: %lu\n", query_len);
    revwm.setQuerySequence(reverse);
	revwm.LengthSet(query_len);

    for(size_t i = 0; i < query_len; i++) {
       int idx, prb, sc, notsc;
       char q = toupper(query[i]);
       if (q == 'N') {
          idx = 99;  // Some invalid value
          sc = notsc = -3;
          readentropy += 2;
       } else {
          prb = quality[i] - '@';                // ASCII coded quality values.
          sc = bw.btowm(prb);
          notsc = bw.notbtowm(prb);
          lastgood = i;
          switch(q) {
             case 'A':
                idx = 0;
                break;
             case 'C':
                idx = 1;
                break;
             case 'G':
                idx = 2;
                break;
             case 'T':
                idx = 3;
                break;
             default:
                idx = 99;
          }
          readentropy += bw.fastqtoentropy(prb);
       }
       for(int b = 0; b < 4; b++) {
          fwdwm.WeightSet(i, b, b == idx ? sc : notsc);
          revwm.WeightSet(query_len - i - 1, 3 - b, b == idx ? sc : notsc);
       }
    }

    trailingNs = query_len - lastgood - 1;

    if (trailingNs != 0) {
		readentropy -= 2 * trailingNs;
		query_len -= trailingNs;
		query[query_len] = '\0';
		fwdwm.LengthSet(query_len);
		revwm.LengthSet(query_len);
	}
	if (readentropy > 2 * (maxNs - trailingNs)) {
		/*
		printf("@%s\n#P%d\n%s\t%s\n# Number of Ns (Read entropy %5.1f bits or %d bases) of the query is above the maximum Ns (%d).\n\n",
				read_id, which, query, quality,
				readentropy + trailingNs * 2, (int)(readentropy / 2.0 + 0.5) + trailingNs, maxNs);
		*/
		return false;
	}
	if (query_len - readentropy / 2 < minspec) {
		/*
		printf("@%s\n#P%d\n# Query length less the number of Ns (Read entropy %5.1f bits or %d bases) of the query is below the minimum length (%d).\n\n",
				read_id, which,
				readentropy + trailingNs * 2, (int)(readentropy / 2.0 + 0.5) + trailingNs, minspec);
		*/
		return false;
	}

    return true;
}

/*
@SNPSTER2_50_3_1_119_588
TGATGCTCTTATTCAGTTCTTTTCCATTGAATAAGT
+SNPSTER2_50_3_1_119_588
hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
@SNPSTER2_50_3_1_121_522
TTCTGCTCTCCCTGTCCGTCTCTTCTCTCTCCTCTT
+SNPSTER2_50_3_1_121_522
hhhhhhhhhhhhhhhhhhFhhhhhhhhhhh]hhhFh
*/

int Align_Illumina::read()
{
	bool rc;

	int read_count = 0;
	for (;;) {
		if (fgets(query_id, query_alloc, fr) == NULL)
			return -1;

		if (fgets(query, query_alloc, fr) == NULL)
			return -1;
		trim(query);

		if (fgets(query_id, query_alloc, fr) == NULL)
			return -1;
		trim(query_id);

		if (fgets(qual, query_alloc, fr) == NULL)
			return -1;
		trim(qual);

		lineno += 4;

		if (query_id[0] != '+') {
			fprintf(stderr, "Expecting '+' as the first character in sequence file but found '%c' instead (line #%lu: \"%s\").\n",
					query_id[0], lineno, query_id);
			return -1;
		}

		memmove(query_id, query_id + 1, strlen(query_id));
		sxnacodes.rc(reverse, query);
		query_len = strlen(query);
		memmove(rev_qual, qual, query_len);
		rev_qual[query_len] = '\0';
		reverse_read(rev_qual, query_len);

		rc = setSXOligoSearchWM(fwdwm, revwm, query_id, query, reverse, qual, query_len);
		read_count++;

		if (rc)
			break;
	}
    
	/* Call LengthSet again as setSXOligoSearchWM may have reduced it (due to Ns at the ends)
	 * This depends on SXOligoSearchWM's behavior of keeping the weight matrix even when the length was reduced
	 */
	fwdwm.LengthSet(query_len);
	revwm.LengthSet(query_len);
	
	return read_count;
}


void Align_Illumina::write_nohit()
{
	if (fnohit) {
		fprintf(fnohit, "@%s\n%s\n+%s\n%s\n",
			query_id, query, query_id, qual);
	}
}

