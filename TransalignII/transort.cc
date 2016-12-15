#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <map>

class Key {
private:	
	const char *target_id;
	long target_pos;

public:
	inline Key(const char *id) { target_id = id; target_pos = 0; }
	inline const char *getTargetId() const { return target_id; }
	inline long setTargetPos(long v) { return target_pos = v; }
	inline long getTargetPos() const { return target_pos; }
};

struct ltstr
{
	inline bool operator()(const Key* k1, const Key* k2) const
	{
		int rc = strcmp(k1->getTargetId(), k2->getTargetId());
		if (rc < 0)
			return true;

		if (rc > 0)
			return false;

		return k1->getTargetPos() < k2->getTargetPos();
	}
};
typedef std::multimap<Key *, long, ltstr> Results;

static int compare_token(const char *tok, const char *line, int tok_no, char sep = '\t')
{
	int to = 0;

	while (to < tok_no) {
		if (*line == '\0')
			return -1;
		if (*line == sep)
			to++;
		line++;
	}

	while (*tok == *line) {
		tok++;
		line++;
	}

	if (*tok == '\0' && (*line == sep || *line == '\0'))
		return to;

	return -1;
}

static char *extract_token(const char *line, int tok_no, char sep = '\t')
{
	while (tok_no > 0) {
		if (*line == '\0')
			return NULL;
		if (*line == sep)
			tok_no--;
		line++;
	}

	size_t tok_len = 0;
	const char *p = line;
	while (*p && *p != sep) {
		tok_len++;
		p++;
	}

	char *r = new char[tok_len + 1];
	char *r2 = r;
	while (*line && *line != sep) {
		*r2++ = *line++;
	}
	*r2 = '\0';

	return r;
}

static void build_results(Results& r, FILE *f)
{
	int state = 0;
	char line[4096];
	long file_pos;
	Key *key;

	while (!feof(f)) {
		long f_pos = ftell(f);
		fgets(line, sizeof(line), f);

		if (state == 1) {
			if (memcmp(line, "Score =", 7) != 0) {
				fprintf(stderr, "Expecting \"Score =\" instead of %s", line);
				exit(2);
			}
			state++;
			continue;
		}

		if (state == 2) {
			if (memcmp(line, "Query:", 6) != 0) {
				fprintf(stderr, "Expecting \"Query:\" instead of %s", line);
				exit(2);
			}
			state++;
			continue;
		}

		if (state == 3) {
			if (memcmp(line, "      ", 6) != 0) {
				fprintf(stderr, "Expecting \"      \" instead of %s", line);
				exit(2);
			}
			state++;
			continue;
		}

		if (state == 4) {
			if (memcmp(line, "Sbjct:", 6) != 0) {
				fprintf(stderr, "Expecting \"Sbjct:\" instead of %s", line);
				exit(2);
			}

			assert(key != NULL);
			key->setTargetPos(atol(line + 6));

			r.insert(std::pair<Key *, long>(key, file_pos));
			key = NULL;

			state = 0;
			continue;
		}

		// look for "sxoght:" as first token...
		int token_no = compare_token("sxoght:", line, 0);
		if (token_no == 0) {
			if (compare_token("#query", line, 1) < 1) {
				file_pos = f_pos;
				char *token = extract_token(line, 2);
				if (token != NULL) {
					key = new Key(token);
					state = 1;
				}
			}
		}
	}
		
}

static void write_results(const Results& r, FILE *fin, FILE *fout = stdout)
{
	char line[4096];

	for (Results::const_iterator it = r.begin(); it != r.end(); it++) {
		fseek(fin, (*it).second, SEEK_SET);
		do {
			if (fgets(line, sizeof(line), fin) == NULL)
				break;
			fputs(line, fout);
		} while (line[0] != '\0' && line[0] != '\n');
	}
}

static void banner(int argc, char *argv[])
{
	fprintf(stderr, "Usage: %s file\nwhere file is to be sorted\n\n", argv[0]);
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		banner(argc, argv);
		return 1;
	}

	FILE *f = fopen(argv[1], "r");
	if (f != NULL) {
		Results r;
		build_results(r, f);

		write_results(r, f);
		fclose(f);
	}

	return 0;
}

