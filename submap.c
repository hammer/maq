#include <zlib.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include "main.h"
#include "maqmap.h"

typedef struct {
	int min_mapQ, max_mm, is_pair, max_sum_err;
} submap_opt_t;

static submap_opt_t *maq_new_submap_opt()
{
	submap_opt_t *p;
	p = (submap_opt_t*)calloc(1, sizeof(submap_opt_t));
	p->min_mapQ = 10;
	p->max_mm = 3;
	p->is_pair = 0;
	p->max_sum_err = 60;
	return p;
}

static bit64_t maq_submap_core(gzFile fpin, gzFile fpout, const submap_opt_t *so)
{
	maqmap_t *mm;
	maqmap1_t m1;
	bit64_t n = 0;
	mm = maqmap_read_header(fpin);
	mm->n_mapped_reads = 0;
	maqmap_write_header(fpout, mm);
	while (gzread(fpin, &m1, sizeof(maqmap1_t))) {
		int is_keep = 1;
		if (m1.map_qual < so->min_mapQ || (m1.info1&0xf) > so->max_mm) is_keep = 0;
		if (so->is_pair && (m1.flag&PAIRFLAG_PAIRED) == 0) is_keep = 0;
		if (m1.info2 > so->max_sum_err) is_keep = 0;
		if (is_keep) {
			gzwrite(fpout, &m1, sizeof(maqmap1_t));
			++n;
		}
	}
	maq_delete_maqmap(mm);
	return n;
}

int maq_submap(int argc, char *argv[])
{
	int c;
	gzFile fpin, fpout;
	submap_opt_t *so = maq_new_submap_opt();
	while ((c = getopt(argc, argv, "q:Q:m:p")) >= 0) {
		switch (c) {
		case 'q': so->min_mapQ = atoi(optarg); break;
		case 'Q': so->max_sum_err = atoi(optarg); break;
		case 'm': so->max_mm = atoi(optarg); break;
		case 'p': so->is_pair = 1; break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   maq submap [options] <out.map> <in.map>\n\n");
		fprintf(stderr, "Options: -q INT      minimum mapping quality [%d]\n", so->min_mapQ);
		fprintf(stderr, "         -Q INT      maximum sum of errors [%d]\n", so->max_sum_err);
		fprintf(stderr, "         -m INT      maximum number of mismatches [%d]\n", so->max_mm);
		fprintf(stderr, "         -p          correctly paired reads only\n\n");
		free(so);
		return 1;
	}
	fpout = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdout), "w") : gzopen(argv[optind], "w");
	fpin = (strcmp(argv[optind+1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind+1], "r");
	maq_submap_core(fpin, fpout, so);
	gzclose(fpin); gzclose(fpout);
	free(so);
	return 0;
}
