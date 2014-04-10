// State-of-art indel detector

#include <assert.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "stdhash.hh"
#include "bfa.h"
#include "main.h"
#include "match.hh"

static int MM_PENALTY = 2;
static int SEARCH_WIN = 5;
static int SKIP_WIN = 10;

typedef struct
{
	int i5, i3;
	bit8_t c3, c5, min;
	float score;
} indel_candidate_t;

// counter - empty:8, 5'-clipping:8, 3'-clipping:8, cross:8
static void fill_counter(bit32_t *counter, int seqid, nst_bfa1_t *l, gzFile fp_map)
{
	maqmap1_t *m1, mm1;
	m1 = &mm1;
	while (maqmap_read1(fp_map, m1)) {
		if (int(m1->seqid) < seqid) continue;
		if (int(m1->seqid) > seqid) break; // WARNING: the last read will be lost
		int begin, end, mid = (m1->pos>>1) + (m1->size>>1);
		int i, j, score, max, max_k, k;
		bit64_t s;
		// right clipping
		i = mid >> 5; j = mid & 0x1f;
		for (k = m1->size>>1, max_k = -1, score = max = 0, s = l->seq[i]; k != m1->size; ++k, ++j) {
			if (j == 32) { // jump to the next word
				s = l->seq[++i];
				j = 0;
			}
			if (int(m1->seq[k]>>6) == int(s>>(31-j<<1)&3)) ++score;
			else score -= MM_PENALTY;
			if (score > max) {
				max = score; max_k = k;
			}
		}
		end = max_k; // > end should be clipped out
		// left clipping
		i = mid >> 5; j = mid & 0x1f;
		for (k = m1->size>>1, max_k = -1, score = max = 0, s = l->seq[i]; k >= 0; --k, --j) {
			if (j == -1) { // jump to the previous word
				s = l->seq[--i];
				j = 31;
			}
			if (int(m1->seq[k]>>6) == int(s>>(31-j<<1)&3)) ++score;
			else score -= MM_PENALTY;
			if (score > max) {
				max = score; max_k = k;
			}
		}
	    begin = max_k; // < begin should be clipped out
		// count
		int pos = m1->pos>>1;
		k = pos + begin - 1;
		if (begin != 0 && k >= 0 && k < l->ori_len && (counter[k]>>16&0xff) < 0xff) counter[k] += 0x10000;
		k = pos + end;
		if (end != m1->size - 1 && k >= 0 && k < l->ori_len && (counter[k]>>8&0xff) < 0xff) counter[k] += 0x100;
		for (k = pos + begin + SEARCH_WIN - 1; k < pos + end - SEARCH_WIN + 1; ++k)
			if (k < l->ori_len && k >= 0 && (counter[k]&0xff) < 0xff)
				++counter[k];
	}
}
static void pickup_indel(bit32_t *counter, nst_bfa1_t *l)
{
	int i, j, k, n;
	long double mean = 0.0;
	indel_candidate_t best, curr;
	memset(&best, 0, sizeof(indel_candidate_t));
	memset(&curr, 0, sizeof(indel_candidate_t));
	for (i = n = 0; i != l->len; ++i) {
		bit64_t mask = l->mask[i];
		for (j = 31; j >= 0; --j) {
			int coor = (i<<5) | (31-j);
			if (coor >= l->ori_len) break;
			if (mask>>(j<<1)&3) {
				++n; mean += counter[coor]&0xff;
			}
		}
	}
	mean /= n; // mean of covered depth
	int depth_thres = int(mean * 0.5 + 0.5);
	best.score = -100.0; best.i3 = best.i5 = 0x7fffffff;
//	for (i = 0; i != l->ori_len; ++i)
//		fprintf(stderr, "%d\t%d\t%d\t%d\n", i+1, counter[i]&0xff, counter[i]>>8&0xff, counter[i]>>16&0xff);
	for (i = 0; i != l->ori_len; ++i) {
		bit32_t c = counter[i];
		if ((c>>8&0xff) || (c>>16&0xff)) { // initiate the search
			int c3 = -1, c5 = -1;
			int i3 = -1, i5 = -1;
			for (k = i; k != l->ori_len && k <= i + SEARCH_WIN; ++k) {
				if (int(counter[k]>>16&0xff) > c5) { c5 = counter[k]>>16&0xff; i5 = k; }
				if (int(counter[k]>>8&0xff) > c3) { c3 = counter[k]>>8&0xff; i3 = k; }
			}
			if (c3 <= 0 || c5 <= 0) continue;
			curr.i3 = i3; curr.i5 = i5;
			curr.c3 = c3; curr.c5 = c5;
			int begin = (i3 < i5)? i3 : i5;
			int end = (i5 > i3)? i5 : i3;
			int min = 0xffff;
			for (k = begin; k <= end; ++k)
				if (int(counter[k]&0xff) < min) min = counter[k]&0xff;
			curr.min = min;
			curr.score = 2.0 * (curr.c3 + curr.c5 - min) / mean;
			if (begin - best.i3 <= SKIP_WIN || begin - best.i5 <= SKIP_WIN) {
				if (curr.score > best.score) best = curr;
			} else {
				if (best.score >= -0.5 && best.min < depth_thres)
					printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\n", l->name, best.i3+1, best.i5-best.i3, best.min,
						   best.c3, best.c5, best.c3+best.c5-best.min);
				best = curr;
			}
		}
	}
}
static void indel_soa_core(FILE *fp_bfa, gzFile fp_map)
{
	nst_bfa1_t *l;
	int seqid = 0, k;
	hash_map_char<int> *hash_map = new hash_map_char<int>;
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names
	for (k = 0; k != mm->n_ref; ++k) hash_map->insert(mm->ref_name[k], k);
	k = mm->n_mapped_reads;
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		if (!hash_map->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		bit32_t *counter = (bit32_t*)calloc(l->ori_len, 4);
		fill_counter(counter, seqid, l, fp_map);
		pickup_indel(counter, l);
		free(counter);
		nst_delete_bfa1(l);
	}
	maq_delete_maqmap(mm);
	delete hash_map;
}

static int usage()
{
	fprintf(stderr, "Usage: maq indelsoa <ref.bfa> <align.map>\n");
	return 1;
}
int ma_indelsoa(int argc, char *argv[])
{
	int c;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc - optind < 2) return usage();
	FILE *fp_bfa = fopen(argv[optind], "r");
	gzFile fp_map = gzopen(argv[optind+1], "r");
	assert(fp_bfa); assert(fp_map);
	indel_soa_core(fp_bfa, fp_map);
	fclose(fp_bfa);
	gzclose(fp_map);
	return 0;
}
