#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <zlib.h>
#include <unistd.h>
#include "match.hh"
#include "stdhash.hh"
#include "main.h"
#include "bfa.h"

#define HIGHQ 30
#define WINSIZE 1000

template<class TYPE>
static long algo_extend(int size, const TYPE *array, int pos, int *begin, int *end)
{
	if (array[pos] < 0) return -1;
	long score, ret, max;
	int max_pos;
	max = score = 0; max_pos = -1;
	for (int i = pos; i != size; ++i) {
		score += array[i];
		if (score > max) {
			max = score;
			max_pos = i;
		}
		if (score < 0) break;
	}
	*end = max_pos + 1; ret = max;
	max = score = 0; max_pos = -1;
	for (int i = pos; i >= 0; --i) {
		score += array[i];
		if (score > max) {
			max = score;
			max_pos = i;
		}
		if (score < 0) break;
	}
	*begin = max_pos; ret += max;
	return ret - array[pos];
}

static void paircov_mask(short *counter, bit8_t *depth, int seqid, nst_bfa1_t *l, gzFile fp_map)
{
	maqmap1_t *m1, mm1;
	m1 = &mm1;
	while (maqmap_read1(fp_map, m1)) {
		if (int(m1->seqid) != seqid) continue;
		if (m1->flag == (PAIRFLAG_PAIRED|PAIRFLAG_FR) && m1->dist > 0) { // pair depth
			short *begin = counter + (m1->pos>>1);
			short *end = begin + m1->dist;
			for (short *p = begin; p != end; ++p)
				if (*p < 0x7fff) ++(*p);
		}
		{ // read depth
			bit8_t *begin = depth + (m1->pos>>1) - 1;
			bit8_t *end = depth + m1->size + (m1->pos>>1);
			for (bit8_t *p = begin; p != end; ++p)
				if (*p < 0xff) ++(*p);
		}
	}
}

static int MIN_PAIRCOV = 1;
static int MIN_N_LEN = 250;

static void paircov_core(FILE *fp_bfa, gzFile fp_map)
{
	nst_bfa1_t *l;
	int seqid = 0, k, i;
	hash_map_char<int> *hash_map = new hash_map_char<int>;
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names
	for (k = 0; k != mm->n_ref; ++k) hash_map->insert(mm->ref_name[k], k);
	k = mm->n_mapped_reads;
	l = nst_load_bfa1(fp_bfa);
	assert(l);
	int ret = hash_map->find(l->name, &seqid);
	assert(ret);
	short *counter = (short*)calloc(l->ori_len, 2);
	bit8_t *depth = (bit8_t*)calloc(l->ori_len, 1);
	paircov_mask(counter, depth, seqid, l, fp_map);
	// calculate mean
	long double mean = 0.0, depth_mean = 0.0;
	for (i = 0; i != l->ori_len; ++i) {
		mean += counter[i];
		depth_mean += depth[i];
	}
	mean /= l->ori_len; depth_mean /= l->ori_len;
	fprintf(stderr, "mean=%Lf,depth_mean=%Lf\n", mean, depth_mean);
	for (i = 0, k = -1; i < l->ori_len; ++i) {
		int is_N = ((l->mask[i>>5]>>((31-(i&0x1f))<<1))&3)? 0 : 1;
		if (is_N) {
			int j;
			for (j = i; j < l->ori_len; ++j)
				if ((l->mask[j>>5]>>((31-(j&0x1f))<<1))&3) break;
			if (j - i > MIN_N_LEN) {
				if (k >= 0) printf("NN\t%s\t%d\t%d\t%d\n", l->name, k+1, i, i - k);
				k = j;
			}
			i = j - 1;
		}
	}
	// detect absolute break point
	for (i = 0; i < l->ori_len; ++i) {
		if (counter[i] > MIN_PAIRCOV) {
			int j;
			for (j = i; j < l->ori_len; ++j)
				if (counter[j] <= MIN_PAIRCOV) break;
			printf("AB\t%s\t%d\t%d\t%d\n", l->name, i+1, j, j - i);
			i = j - 1;
		}
	}
	free(counter); free(depth);
	nst_delete_bfa1(l);
	maq_delete_maqmap(mm);
	delete hash_map;
}
int ma_paircov(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maq paircov [-d %d] [-n %d] <ref.bfa> <align.map>\n",
				MIN_PAIRCOV, MIN_N_LEN);
		return 1;
	}
	FILE *fp_bfa = fopen(argv[1], "r");
	gzFile fp_map = gzopen(argv[2], "r");
	assert(fp_bfa); assert(fp_map);
	paircov_core(fp_bfa, fp_map);
	fclose(fp_bfa);
	gzclose(fp_map);
	return 0;
}

static int abpair_usage()
{
	fprintf(stderr, "maq abpair <input.map>\n");
	return 0;
}
int ma_abpair(int argc, char *argv[])
{
	if (argc == 1) return abpair_usage();
	gzFile fpin = gzopen(argv[1], "r");
	maqmap_t *mm = maqmap_read_header(fpin);
	maqmap1_t *m1, mm1;
	bit64_t n_diffchr, n_nomatch, n_good, n_badori, n_single, n_bad_dist;
	bit64_t n_diffchr2, n_badori2, n_bad_dist2, n_good_alt, n_sw;
	n_diffchr = n_nomatch = n_good = n_badori = n_single = n_bad_dist = n_sw = 0;
	n_diffchr2 = n_badori2 = n_bad_dist2 = n_good_alt = 0;
	fprintf(stdout, "CC\tDC readName chr coor altMapQ\n");
	fprintf(stdout, "CC\tSC readName chr flag altMapQ coor1 coor2\n");
	m1 = &mm1;
	while (maqmap_read1(fpin, m1)) {
		if (m1->flag && m1->alt_qual >= HIGHQ) ++n_good_alt;
		if (m1->flag == 0) ++n_single;
		else if (m1->flag == (PAIRFLAG_PAIRED|PAIRFLAG_FR)) ++n_good;
		else if (m1->flag == (PAIRFLAG_FR|PAIRFLAG_SW)) ++n_sw;
		else if (m1->flag == PAIRFLAG_NOMATCH) ++n_nomatch;
		else if (m1->flag == PAIRFLAG_DIFFCHR) {
			++n_diffchr; if (m1->alt_qual >= HIGHQ) ++n_diffchr2;
			m1->name[MAX_NAMELEN-1] = 0;
			if (m1->alt_qual > 0)
				printf("DC\t%s\t%s\t%d\t%d\n", m1->name, mm->ref_name[m1->seqid],
					   (m1->pos>>1) + 1, m1->alt_qual);
		} else {
			if (m1->flag & PAIRFLAG_FR) { ++n_bad_dist; if (m1->alt_qual >= HIGHQ) ++n_bad_dist2; }
			else { ++n_badori; if (m1->alt_qual >= HIGHQ) ++n_badori2; }
			if (m1->dist > 0 && m1->alt_qual > 0) { // print out bad pairs
				m1->name[MAX_NAMELEN-1] = 0;
				printf("SC\t%s\t%s\t%d\t%d\t%d\t%d\n", m1->name, mm->ref_name[m1->seqid], m1->flag,
					   m1->alt_qual, (m1->pos>>1) + 1, (m1->pos>>1) + m1->size);
			}
		}
	}
	// print out stat in the stderr
	bit64_t n_paired = mm->n_mapped_reads - n_single;
	fprintf(stdout, "CC all the numbers mean 'number of reads', not 'number of pairs'.\n");
	fprintf(stdout, "MM total mapped reads: %llu\n", mm->n_mapped_reads);
	fprintf(stdout, "MM total mapped paired reads: %llu (/%llu = %.3f%%)\n",
			n_paired, mm->n_mapped_reads, 100.0 * n_paired / mm->n_mapped_reads);
	if (n_paired) {
		fprintf(stdout, "MM correctly paired reads: %llu (/%llu = %.3f%%)\n",
				n_good, n_paired , 100.0 * n_good / n_paired);
		fprintf(stdout, "MM reads added by Smith-Waterman alignment: %llu\n", n_sw);
		fprintf(stdout, "MM one end unmapped: %llu (/%llu = %.3f%%)\n",
				n_nomatch, n_paired, 100.0 * n_nomatch / n_paired);
		fprintf(stdout, "MM ends to difference chr: %llu (/%llu = %.3f%%)\n",
				n_diffchr, n_paired, 100.0 * n_diffchr / n_paired);
		fprintf(stdout, "MM bad orientation: %llu (/%llu = %.3f%%)\n",
			n_badori, n_paired, 100.0 * n_badori / n_paired);
		fprintf(stdout, "MM bad distance: %llu (/%llu = %.3f%%)\n",
				n_bad_dist, n_paired, 100.0 * n_bad_dist / n_paired);
		fprintf(stdout, "MM Q%d alt_qual: %llu\n", HIGHQ, n_good_alt);
		fprintf(stdout, "MM Q%d ends to difference chr: %llu\n", HIGHQ, n_diffchr2);
		fprintf(stdout, "MM Q%d bad orientation: %llu\n", HIGHQ, n_badori2);
		fprintf(stdout, "MM Q%d bad distance: %llu\n", HIGHQ, n_bad_dist2);
	}
	maq_delete_maqmap(mm);
	gzclose(fpin);
	return 0;
}
