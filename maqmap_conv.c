#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <zlib.h>
#include "maqmap.h"
#include "main.h"

typedef struct
{
	bit8_t seq[MAX_READLEN];
	bit8_t size, map_qual, i1, i2, c[2], flag, alt_qual;
	bit32_t seqid, pos;
	int dist;
} maqmap1_oldaux_t;

void maqmap_conv_core(gzFile fpold, gzFile fpnew)
{
	maqmap_t *mm = maq_new_maqmap();
	bit32_t i, n_reads;
	maqmap1_oldaux_t *mo1, mmo1;
	maqmap1_t *m1, mm1;
	int k, len;
	mo1 = &mmo1; m1 = &mm1;
	memset(m1, 0, sizeof(maqmap1_t));
	/* read the header */
	gzread(fpold, &mm->n_ref, sizeof(int));
	if (mm->n_ref == MAQMAP_FORMAT_NEW) {
		fprintf(stderr, "** New map format is detected. No need to convert the format.\n");
		exit(3);
	}
	mm->ref_name = (char**)calloc(mm->n_ref, sizeof(char*));
	for (k = 0; k != mm->n_ref; ++k) {
		gzread(fpold, &len, sizeof(int));
		mm->ref_name[k] = (char*)malloc(len * sizeof(char));
		gzread(fpold, mm->ref_name[k], len);
	}
	gzread(fpold, &n_reads, sizeof(bit32_t));
	mm->n_mapped_reads = n_reads;
	maqmap_write_header(fpnew, mm);
	maq_delete_maqmap(mm);
	/* read and convert records */
	for (i = 0; i != n_reads; ++i) {
		gzread(fpold, mo1, sizeof(maqmap1_oldaux_t));
		memcpy(m1, mo1, sizeof(maqmap1_oldaux_t));
		sprintf(m1->name, "%u", i);
		gzwrite(fpnew, m1, sizeof(maqmap1_t));
	}
}

int ma_mapass2maq(int argc, char *argv[])
{
	gzFile fpold, fpnew;
	if (argc < 3) {
		fprintf(stderr, "maq mapass2maq <mapass2.map> <maq.map>\n");
		return 1;
	}
	fpold = gzopen(argv[1], "r");
	fpnew = gzopen(argv[2], "w");
	assert(fpold && fpnew);
	maqmap_conv_core(fpold, fpnew);
	gzclose(fpold); gzclose(fpnew);
	return 0;
}
