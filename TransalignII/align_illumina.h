#ifndef _ALIGN_ILLUMINA_H_
#define _ALIGN_ILLUMINA_H_

#include "pair.h"

class Align_Illumina : public TransAlign
{
private:
	size_t lineno;

	int read();
	void write_nohit();

	bool setSXOligoSearchWM(SXOligoSearchWM& fwdwm, SXOligoSearchWM& revwm, char *read_id, char *query, char *reverse, char *quality, size_t query_len);

public:
	Align_Illumina(const char *db);
	virtual ~Align_Illumina();
};

#endif /* _ALIGN_ILLUMINA_H_ */
