#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "algo.hh"
#include "maqmap.h"
#include "main.h"

#define HEAP_EMPTY 0xffffffffffffffffull

typedef struct
{
	int i;
	bit64_t pos;
	maqmap1_t *m1;
} mapping_heap_t;

inline bool operator < (const mapping_heap_t &a, const mapping_heap_t &b)
{
	return (a.pos > b.pos); // note that this is ">", not "<"
}
// This function will open "n" files at the same time. On most OS, there is a limit.
// This is a O(N log n) algorithm, where N is the total number of reads and n is
// the number of files.
void mapping_merge_core(char *out, int n, char **fn)
{
	gzFile *fp, fpout;
	mapping_heap_t *heap;
	maqmap_t **mm, *mm_out;
	int n_ref;
	
	fpout = (strcmp(out, "-") == 0)? gzdopen(fileno(stdout), "w") : gzopen(out, "w");
	assert(fpout);
	fp = (gzFile*)calloc(n, sizeof(gzFile));
	heap = (mapping_heap_t*)calloc(n, sizeof(mapping_heap_t));
	mm = (maqmap_t**)calloc(n, sizeof(maqmap_t*));
	bit64_t c = 0;
	for (int i = 0; i != n; ++i) {
		mapping_heap_t *h;
		fp[i] = gzopen(fn[i], "r");
		assert(fp[i]);
		// It would be much better if this program can check whether reads are
		// aligned to the same reference. However, I am lazy now. I trust
		// endusers to do the right things.
		mm[i] = maqmap_read_header(fp[i]);
		c += mm[i]->n_mapped_reads;
		h = heap + i;
		h->i = i;
		h->m1 = (maqmap1_t*)malloc(sizeof(maqmap1_t));
		if (maqmap_read1(fp[i], h->m1))
			h->pos = ((bit64_t)h->m1->seqid<<32) | h->m1->pos;
		else h->pos = HEAP_EMPTY;
	}
	// fill mm_out, write to file and then delete it.
	mm_out = maq_new_maqmap();
	n_ref = mm_out->n_ref = mm[0]->n_ref;
	mm_out->n_mapped_reads = c;
	mm_out->ref_name = mm[0]->ref_name;
	maqmap_write_header(fpout, mm_out);
	mm_out->ref_name = 0; mm_out->n_ref = 0;
	maq_delete_maqmap(mm_out);
	// initialize the heap
	int l_record;
	algo_heap_make(heap, n);
	while (heap->pos != HEAP_EMPTY) {
		gzwrite(fpout, heap->m1, sizeof(maqmap1_t));
		if ((l_record = maqmap_read1(fp[heap->i], heap->m1)) != 0) {
			if (l_record != sizeof(maqmap1_t)) {
				fprintf(stderr, "[mapping_mapmerge_core] apparently truncated .map file. Abort!\n");
				exit(1);
			} else if ((int)heap->m1->seqid >= n_ref) {
				fprintf(stderr, "[mapping_mapmerge_core] the %d-th .map file seems to corrupt (%d != %d). Abort!\n",
						heap->i + 1, heap->m1->seqid, mm_out->n_ref);
				exit(1);
			}
			heap->pos = ((bit64_t)heap->m1->seqid<<32) | heap->m1->pos;
		} else heap->pos = HEAP_EMPTY;
		algo_heap_adjust(heap, 0, n);
	}
	// free
	for (int i = 0; i != n; ++i) {
		gzclose(fp[i]);
		free(heap[i].m1);
		maq_delete_maqmap(mm[i]);
	}
	free(mm);
	gzclose(fpout);
	free(fp); free(heap);
}
int ma_mapmerge(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maq mapmerge <out.map> <in1.map> <in2.map> [...]\n");
		return 1;
	}
	mapping_merge_core(argv[1], argc - 2, argv + 2);
	return 0;
}
