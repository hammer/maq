#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "bfa.h"
#include "stdaln.h"
#include "match.hh"
#include "read.h"

#define MATCH_RATIO 0.8
#define MATCH_LENGTH 20
#define MAPQ_THRES 20
#define MIN_ALT_QUAL 31

extern "C" {
	int aln_local_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
					   path_t *path, int *path_len, int do_align);
};

// this is the scoring systems:

typedef struct
{
	int n, m;
	int *list;
} simple_list_t;

typedef struct
{
	int n_lname, tot;
	simple_list_t *lists;
} aln_candidate_t;

static void destroy(aln_candidate_t *can)
{
	int i;
	for (i = 0; i != can->n_lname; ++i)
		free(can->lists[i].list);
	free(can->lists); free(can);
}

static inline void add_candidate(aln_candidate_t *can, int seqid, int which, int strand)
{
	simple_list_t *l = can->lists + seqid;
	++can->tot;
	if (l->n == l->m) {
		l->m += 0x10000;
		l->list = (int*)realloc(l->list, sizeof(int) * l->m);
	}
	l->list[l->n++] = which * strand;
}

static int cal_insert_size(const match_aux_t *o, const match_data_t *d, double *avg, double *dev)
{
	int i;
	bit64_t mean, var, n;
	if (o->size_l == 0) return 1;
	mean = var = n = 0;
	for (i = 0; i != d->n_reads>>1; ++i) {
		match_info_t *match = d->match + i;
		int alt_qual = bit64_t((match+1)->last_mm1 < match->last_mm1? (match+1)->last_mm1 : match->last_mm1);
		if ((match->last2 & PAIRFLAG_PAIRED) && match->last_seqid > 0 && alt_qual >= MIN_ALT_QUAL) {
			bit64_t tmp = var;
			int dist = match->last_seqid;
			++n;
			mean += dist;
			var += dist * dist;
			if (var < tmp) {
				mean -= dist; var = tmp; --n;
				fprintf(stderr, "[cal_insert_size] overflow happened.\n");
				break;
			}
		}
	}
	if (n < 20) {
		fprintf(stderr, "[cal_insert_size] too few reads mapped as pairs. Skip Smith-Waterman alignment.\n");
		return 1;
	}
	*avg = (double)mean/n;
	*dev = sqrt((double)var/n - (*avg) * (*avg));
	fprintf(stderr, "[cal_insert_size] %lld read pairs counted. insert size: %lf +/- %lf\n", n, *avg, *dev);
	return 0;
}

// match_info_t will be modified:
// last1 = start position (0 based)
// last2 - stop position
// seqid1 = mate's seqid1
static aln_candidate_t *maq_make_aln_candidate(const match_aux_t *o, match_data_t *d, double avg, double dev)
{
	aln_candidate_t *can = (aln_candidate_t*)calloc(1, sizeof(aln_candidate_t));
	int size[2] = { o->size_l, o->size_r };
	can->n_lname = d->n_lname;
	can->lists = (simple_list_t*)calloc(d->n_lname, sizeof(simple_list_t));
	int k, max_dist, min_dist;
	max_dist = int(avg + 2.0 * dev + 0.5);
	min_dist = int(avg - 2.0 * dev + 0.5);
	if (min_dist < o->size_l+1) min_dist = o->size_l+1;
	if (min_dist < o->size_r+1) min_dist = o->size_r+1;
	if (max_dist > o->max_dist) max_dist = o->max_dist;
	for (k = 0; k != d->n_reads>>1; ++k) {
		match_info_t *match = d->match + (k<<1);
		match_info_t *nomatch = match + 1;
		int tmp_dist;
		int kk = k<<1|1; // kk is the index of the nomatch read
		if (match->i1 == MA_NO_MATCH && nomatch->i1 == MA_NO_MATCH) continue;
		if (match->i1 != MA_NO_MATCH && nomatch->i1 != MA_NO_MATCH) continue;
		if (match->i1 == MA_NO_MATCH) { // swap
			match = nomatch;
			nomatch = match - 1;
			--kk;
		}
		if (bit64_t(match->last_mm1) < MAPQ_THRES) continue;
		if (!o->is_color) { // nucleotide space
			if (match->p1&1) { // matched on the reverse strand
				if ((match->p1>>1) - max_dist > 0) {
					tmp_dist = int(match->p1>>1) - max_dist;
					nomatch->last1 = (tmp_dist < 0)? 0 : tmp_dist;
					tmp_dist = int(match->p1>>1) + size[kk&1] - min_dist;
					nomatch->last2 = (tmp_dist < 0)? 0 : tmp_dist;
					nomatch->seqid1 = match->seqid1;
					add_candidate(can, nomatch->seqid1, kk, 1);
				}
			} else { // the matched end on the forward strand
				tmp_dist = int(match->p1>>1) - size[0] - size[1] + min_dist;
				nomatch->last1 = (tmp_dist < 0)? 0 : tmp_dist;
				tmp_dist = int(match->p1>>1) - size[1-(kk&1)] + max_dist;
				nomatch->last2 = (tmp_dist < 0)? 0 : tmp_dist;
				nomatch->seqid1 = match->seqid1;
				add_candidate(can, nomatch->seqid1, kk, -1);
			}
		} else { // color space
			if (((match->p1^kk)&1) == 1) { // FF-read2-read1 or RR-read1-read2
				if ((match->p1>>1) - max_dist > 0) {
					tmp_dist = int(match->p1>>1) - max_dist;
					nomatch->last1 = (tmp_dist < 0)? 0 : tmp_dist;
					tmp_dist = int(match->p1>>1) + size[kk&1] - min_dist;
					nomatch->last2 = (tmp_dist < 0)? 0 : tmp_dist;
					nomatch->seqid1 = match->seqid1;
					add_candidate(can, nomatch->seqid1, kk, (match->p1&1)? -1 : 1);
				}
			} else {
				tmp_dist = int(match->p1>>1) - size[0] - size[1] + min_dist;
				nomatch->last1 = (tmp_dist < 0)? 0 : tmp_dist;
				tmp_dist = int(match->p1>>1) - size[1-(kk&1)] + max_dist;
				nomatch->last2 = (tmp_dist < 0)? 0 : tmp_dist;
				nomatch->seqid1 = match->seqid1;
				add_candidate(can, nomatch->seqid1, kk, (match->p1&1)? -1 : 1);
			}
		}
	}
	fprintf(stderr, "[maq_make_aln_candidate] %d candidates added\n", can->tot);
	return can;
}

static void align_candidate(const nst_bfa1_t *l, const match_aux_t *o, const simple_list_t *list, const longreads_t *lr,
							match_data_t *d)
{
	bit8_t *ref_seq, *read_seq, *read_qual;
	int i, j, size[2] = { o->size_l, o->size_r };
	int n_ins, n_del, n_sub;
	n_ins = n_del = n_sub = 0;
	ref_seq = (bit8_t*)malloc(l->ori_len);
	read_seq = (bit8_t*)malloc(MAX_READLEN);
	read_qual = (bit8_t*)malloc(MAX_READLEN);
	for (i = 0; i != l->len; ++i) {
		bit64_t seqi = l->seq[i], maski = l->mask[i];
		for (j = 31; j >= 0; --j) {
			int coor = (i << 5) | (31 - j);
			if (coor >= l->ori_len) break;
			ref_seq[coor] = ((maski>>(j<<1)&3) != 3)? 4 : seqi>>(j<<1)&3;
		}
	}
	for (i = 0; i != list->n; ++i) {
		int x = list->list[i];
		int strand = 1;
		if (x < 0) { // reverse strand
			strand = -1; x = -x;
		}
		int s = size[x&1];
		match_info_t *match = d->match + x;
		match_info_t *mate = (x&1)? match - 1: match + 1;
		bit8_t *seq = lr->seq[x];
		// generate read sequence
		if (strand > 0) { // forward strand
			for (j = 0; j < s; ++j) {
				read_seq[j] = (seq[j] == 0)? 4 : seq[j]>>6&3;
				read_qual[j] = seq[j]&0x3f;
			}
		} else { // reverse strand
			if (!o->is_color) { // nucleotide space
				for (j = 0; j < s; ++j) {
					read_seq[s-1-j] = (seq[j] == 0)? 4 : 3 - (seq[j]>>6&3);
					read_qual[s-1-j] = seq[j]&0x3f;
				}
			} else { // color space
				for (j = 0; j < s; ++j) {
					read_seq[s-1-j] = (seq[j] == 0)? 4 : seq[j]>>6&3;
					read_qual[s-1-j] = seq[j]&0x3f;
				}
			}
		}
		path_t *path = (path_t*)malloc((match->last2 - match->last1 + 10 + s) * sizeof(path_t));
		int path_len = 0, begin, end, do_skip = 0;
		begin = match->last1; end = match->last2;
		if (begin < l->ori_len) {
			if (end > l->ori_len) end = l->ori_len;
			if (end - begin > s) aln_local_core(ref_seq + begin, end - begin, read_seq, s, &aln_param_bwa, path, &path_len, 1);
			else do_skip = 1;
		} else do_skip = 1;
		if (path_len < MATCH_LENGTH || (int)match->last1 + path[path_len-1].i - path[path_len-1].j < 0) do_skip = 1;
		if (do_skip) { free(path); continue; }
		int match_length = MATCH_LENGTH;
		if (match_length < int(MATCH_RATIO * s + 0.5)) match_length = int(MATCH_RATIO * s + 0.5);
		if (path->i - path[path_len-1].i >= match_length && path->j - path[path_len-1].j >= match_length) {
			int tmpI = 0, tmpD = 0;
			for (j = path_len-2; j >= 0; --j) {
				path_t *p = path + j;
				if (p->ctype == FROM_I && (p+1)->ctype != FROM_I) ++tmpI;
				if (p->ctype == FROM_D && (p+1)->ctype != FROM_D) ++tmpD;
			}
			if (tmpI + tmpD <= 1) { // just one indel site
				int indel_j = 0, indel_len;
				int mm, err, n_match;
				indel_len = mm = err = 0;
				if (tmpI == 1) ++n_ins;
				else if (tmpD == 1) ++n_del;
				else ++n_sub;
				// get the exact indel coordinate and sum of errors, etc.
				for (j = path_len - 2, n_match = 1; j >= 0; --j) {
					path_t *p = path + j;
					if (p->ctype != FROM_M) {
						if ((p+1)->ctype == FROM_M) indel_j = j + 1;
						++indel_len;
					} else { // this is a match/mismatch, not indel
						++n_match;
						if (ref_seq[begin + p->i - 1] != read_seq[p->j - 1]) { // mismatch
							++mm; err += read_qual[p->j - 1];
						}
					}
				}
				// calculate coordinates, etc.
				path_t *pathb = path + path_len - 1;
				if (tmpD == 1) indel_len = -indel_len;
				match->i1 = 1;
				match->p1 = (bit32_t)(begin + (pathb->i - pathb->j + size[x&1] - 1)) << 1 | (strand>0? 0:1);
				int dist = (mate->p1>>1 > match->p1>>1)? int(mate->p1>>1) - int(match->p1>>1) + size[x&1]
					: -(int(match->p1>>1) - int(mate->p1>>1) + size[1-(x&1)]);
				match->last_seqid = dist;
				match->last1 = (indel_len == 0)? 0 : path[indel_j].j; // indel position
				match->last_mm1 = matches_t(bit8_t(indel_len)<<8) | ((x&1)? (match-1)->last_mm1 : (match+1)->last_mm1);
				match->last_mm2 = ((err <= 0xff)? err:0xff)<<8 | ((mm <= 15)? mm:15); // n_mismatch and sum of err
				// update mate
				mate->last_seqid = -dist;
				// set/update flag
				if (!o->is_color) { // nucleotide space
					match->last2 = PAIRFLAG_SW | PAIRFLAG_FR; // flag
					mate->last2 = PAIRFLAG_FR | PAIRFLAG_PAIRED;
				} else {
					match->last2 = PAIRFLAG_SW | (strand>0? PAIRFLAG_FF : PAIRFLAG_RR);
					mate->last2 = PAIRFLAG_PAIRED | (strand>0? PAIRFLAG_FF : PAIRFLAG_RR);
				}
			}
		}
		free(path);
	}
	free(ref_seq); free(read_seq); free(read_qual);
	//fprintf(stderr, "[align_candidate] %s: (%d, %d, %d) out of %d\n", l->name, n_ins, n_del, n_sub, list->n);
}

void maq_indel_pe(FILE *fp_bfa, const match_aux_t *o, const longreads_t *lr, match_data_t *d)
{
	aln_candidate_t *can;
	clock_t begin = clock();
	double avg, dev;
	nst_bfa1_t *l;
	int ref_id = 0;
	avg = dev = 0.0;
	if (o->is_sw == 0) {
		fprintf(stderr, "[maq_indel_pe] Smith-Waterman alignment is disabled.\n");
		return;
	}
	if (o->size_r == 0 || o->RF_max_dist) {
		fprintf(stderr, "[maq_indel_pe] the indel detector only works with short-insert mate-pair reads.\n");
		return;
	}
	if (cal_insert_size(o, d, &avg, &dev) != 0) return;
	can = maq_make_aln_candidate(o, d, avg, dev);
	fprintf(stderr, "[maq_indel_pe] %d candidates for alignment\n", can->tot);
	rewind(fp_bfa);
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		align_candidate(l, o, can->lists + (ref_id++), lr, d);
		nst_delete_bfa1(l);
	}
	fprintf(stderr, "[maq_indel_pe] CPU time: %.3lf sec\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
	destroy(can);
}
