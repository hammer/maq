#ifndef LH3_READ_H
#define LH3_READ_H

#include <zlib.h>
#include "dword.hh"
#include "const.h"
#include "seq.h"

#define MA_NONAME

typedef struct __longreads_t
{
	int size_l, size_r;
	int n_reads;
	bit8_t **seq;
	char **name;
} longreads_t;

struct __match_data_t;
struct __match_aux_t;

#ifdef __cplusplus
extern "C" {
#endif
	longreads_t *ma_load_reads(gzFile fp_l, int size_l, gzFile fp_r, int size_r);
	void delete_longreads(longreads_t *lr);
	void ma_init_match_data(struct __match_data_t *d, gzFile fp_l, int *size_l, gzFile fp_r, int *size_r, int is_trim,
							const char *adapter, gzFile hits_fp);
	void maq_methy_modify(struct __match_data_t *d, struct __match_aux_t *o);
	static inline int ma_load_1read(gzFile fp, seq_t *seq, char *name)
	{
		int len;
		if (gzread(fp, &len, sizeof(int)) == 0) return -1;
		gzread(fp, name, sizeof(char) * len);
		gzread(fp, &len, sizeof(int));
		if (seq->m < len + 1) {
			seq->m = len + 1;
			seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * seq->m);
		}
		seq->l = len;
		gzread(fp, seq->s, sizeof(char) * len);
		return len;
	}
#ifdef __cplusplus
}
#endif

#endif
