#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <zlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <sys/resource.h>
#include "algo.hh"
#include "match.hh"
#include "bfa.h"
#include "read.h"
#include "main.h"
#include "dword.hh"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "dummy"
#endif

#define MATCH_INDEX_SHIFT32 8
#define MATCH_INDEX_SHIFT64 (MATCH_INDEX_SHIFT32+32)
#define MATCH_INDEX_SIZE (1ul<<(32-MATCH_INDEX_SHIFT32))

#define MAX_CANDIDATE 32768

static int lsh_global_seed;
static int g_is_first_round, g_processed_percent;
static bit64_t g_n_candidates, g_tot_ref_len, g_processed_len, g_tot_hits, g_last_tot_hits, g_last_n_candidates;

typedef struct { // in attempt to achieve better cache efficiency
	int is_color, size_l, size_r;
	int max_err10, max_dist, min_dist;
	int RF_max_dist;
	char methy_mode;
	gzFile hits_fp;
	int max_hits;
} required_option_t;

typedef struct {
	read_t s, m;
	bit32_t coor, key, k;
	int f;
} ref_candid_t;

typedef struct {
	int n;
	ref_candid_t list[MAX_CANDIDATE];
} ref_canlist_t;

inline int alt_cal_mm(bit64_t x) {
	x = ((x & 0xAAAAAAAAAAAAAAAAull) >> 1)  | (x & 0x5555555555555555ull);
	x = ((x & 0xAAAAAAAAAAAAAAAAull) >> 1)  + (x & 0x5555555555555555ull);
	x = ((x & 0xCCCCCCCCCCCCCCCCull) >> 2)  + (x & 0x3333333333333333ull);
	x = x + (x >> 4) & 0x0F0F0F0F0F0F0F0FLLU;
	x = x + (x >> 8);
	x = x + (x >> 16);
	x = x + (x >> 32) & 0x0000007F;
	return x;
}
inline int alt_cal_err(bit64_t x, bit64_t y) {
	x = (x | ((x & 0xAAAAAAAAAAAAAAAAull)>>1) | ((x&0x5555555555555555ull)<<1)) & y;
	x = ((x & 0xCCCCCCCCCCCCCCCCull) >> 2)  + (x & 0x3333333333333333ull);
	x = x + (x >> 4) & 0x0F0F0F0F0F0F0F0FLLU;
	x = x + (x >> 8);
	x = x + (x >> 16);
	x = x + (x >> 32) & 0x0000007F;
	return x;
}
inline bit64_t alt_shrink(bit64_t x) {
	x = ((x & 0xAAAAAAAAAAAAAAAAull) >> 1)  | (x & 0x5555555555555555ull);
	x = ((x & 0xCCCCCCCCCCCCCCCCull) >> 1)  | (x & 0x3333333333333333ull);
	x = ((x & 0xF0F0F0F0F0F0F0F0ull) >> 2)  | (x & 0x0F0F0F0F0F0F0F0Full);
	x = ((x & 0xFF00FF00FF00FF00ull) >> 4)  | (x & 0x00FF00FF00FF00FFull);
	x = ((x & 0xFFFF0000FFFF0000ull) >> 8)  | (x & 0x0000FFFF0000FFFFull);
	return ((x&0xFFFFFFFF00000000ull)>> 16) | (x & 0x00000000FFFFFFFFull);
}
template<class TYPE>
inline int alt_cal_mm(const dword_t<TYPE> &x) {
	return alt_cal_mm(x.w1) + alt_cal_mm(x.w2);
}
template<class TYPE>
inline int alt_cal_err(const dword_t<TYPE> &x, const dword_t<TYPE> &y) {
	return alt_cal_err(x.w1, y.w1) + alt_cal_err(x.w2, y.w2);
}
template<class TYPE>
inline dword_t<TYPE> alt_shrink(const dword_t<TYPE> &x) {
	return dword_t<TYPE>(alt_shrink(x.w2) << (dword_t<TYPE>::n_bits>>1) | alt_shrink(x.w1));
}

// hash function. 05JAN2006, http://www.cris.com/~Ttwang/tech/inthash.htm
inline bit32_t lsh_hash_64(bit64_t key)
{
	key = (~key) + (key << 18);
	key ^= key >> 31;
	key = (key + (key << 2)) + (key << 4);
	key ^= key >> 11;
	key += key << 6;
	key ^= key >> 22;
	return (bit32_t)key;
}
inline bit32_t lsh_hash_32(bit32_t key)
{
	key = ~key + (key << 15); // key = (key << 15) - key - 1;
	key ^= key >> 12;
	key += key << 2;
	key ^= key >> 4;
	key = (key + (key << 3)) + (key << 11);
	key ^= key >> 16;
	return key;
}
inline bool match_tell_dist(bit32_t last, bit32_t curr, int min, int max, int size, int is_color, int RF_max)
{
	if (RF_max) {
		if (last) {
			if ((last&1) == 0 && (curr&1) == 1) { // FR pair
				int dist = int(curr>>1) - int(last>>1) + size;
				return (dist >= min && dist <= max);
			} else if ((last&1) == 1 && (curr&1) == 0) { // RF pair
				int dist = int(curr>>1) - int(last>>1) + size;
				return (dist >= min && dist <= RF_max);
			}
		}
	} else if (last && (!is_color || ((last^curr)&1) == 0)) {
		int dist = int(curr>>1) - int(last>>1) + size;
		return (dist >= min && dist <= max);
	}
	return false;
}
inline void match_add_pair(pair_info_t *p, bit64_t r1, bit64_t r2, int seqid, const matches_t &mm1, const matches_t &mm2)
{
	bit32_t x = (((r1>>56)+(r2>>56))<<24) | ((r1>>32^r2>>32)&0xffffff);
	if (x <= p->i2) {
		if (x > p->i1) { // then i2 > x > i1
			p->i2 = x;
		} else { // then x <= p->i1
			if (p->p11 != (bit32_t)r1 || p->p12 != (bit32_t)r2 || p->seqid1 != seqid) { // not an identical pair
				p->i2 = p->i1;
				p->i1 = x; p->seqid1 = seqid; p->p11 = (bit32_t)r1; p->p12 = (bit32_t)r2;
				p->m11 = mm1; p->m12 = mm2;
			}
		}
	}
}
static inline void methy_ref(char methy_mode, read_t query[4], int size_l, int size_r, bit64_t seqi, int j)
{
	if (methy_mode == 'c') {
		if ((query[0]&3) == 1) {
			query[0] = (query[1] |= 3);
		} else if ((query[0]&3) == 2) {
			query[2] |= read_t(3)<<((size_l-1)<<1);
			query[3] |= read_t(3)<<((size_r-1)<<1);
		}
	} else if (methy_mode == 'g') {
		if ((query[0]&3) == 2) {
			query[0] = (query[1] ^= 2);
		} else if ((query[0]&3) == 1) {
			query[2] ^= read_t(2)<<((size_l-1)<<1);
			query[3] ^= read_t(2)<<((size_r-1)<<1);
		}
	}
}
static inline void process_candidate(const required_option_t *o, int longseqid, bit64_t tmp_xor,
									 match_data_t *d, ref_canlist_t *list)
{
	int i;
	ref_candid_t *p;
	bit32_t n_reads = d->n_reads;
	g_n_candidates += list->n;
	for (i = 0, p = list->list; i != list->n; ++i, ++p) {
		bit32_t k = p->k, key = p->key, coor = p->coor;
		int f = p->f, m = coor&1;
		match_index_t *index = d->index + f;
		k = index->index[k]; // k is the start position in index->sorted. Here is another cache miss!
		bit128_t *r = index->sorted;
		if (k & 0x80000000u) { // there is only one index hit
			if (r[k & 0x7fffffffu].w2>>32 ^ key) continue; // the sequence is not the same. Cache miss!
			k &= 0x7fffffffu; // get rid of the leading 0x80000000u
		} else { // there are more than one index hit
			for (; (r[k].w2>>32^key) && (r[k].w2>>32^key)>>MATCH_INDEX_SHIFT32 == 0; ++k); // Cache miss!
			if (k == n_reads || (r[k].w2>>32^key)) continue;  // if there is no match, a search will stop here
		}
		// calculate number of mismatches/errors and store them if good enough
		read_t submask, subseq;
		bit64_t subm_shifted, subs_shifted;
		subseq = p->s; submask = ~p->m & index->mask;
		bit32_t signature = bit64_t(subseq) ^ coor ^ tmp_xor; // this will be part of the pseudo-random bits
		subseq &= index->mask;
		int shift = index->shift_seed;
		subs_shifted = bit64_t(subseq>>shift); subm_shifted = bit64_t(submask>>shift);
		for (; key == r[k].w2>>32; ++k) {
			// calculate the number of mismatches in the seed
			bit32_t x = alt_cal_mm((subs_shifted ^ bit64_t(r[k])) | subm_shifted);
			if (x > MAX_MISMATCH || (x == 0 && (f&2))) continue;
			match_info_t *seq = d->match + (bit32_t)r[k].w2; // get the read
			read_t match = (subseq ^ seq->s) | submask | seq->m;
			// calculate the sum of errors
			bit32_t y = alt_cal_err(match, seq->q);
			if (int(y) > o->max_err10) continue;
			x = alt_cal_mm(bit64_t(match>>shift)); // recalculate the number of mismatches in the seed
			if (x > MAX_MISMATCH || (x == 0 && (f&2))) continue;
			if (seq->c[x] < 0xff) ++(seq->c[x]); // update count
			x = (y << 24) | (lsh_hash_32(signature ^ (bit32_t)r[k].w2) & 0xffffff);
			matches_t mm = alt_shrink(match).w1;
			if (o->size_r) { // ***** the following codes are only executed for paired-end data
				int do_append, do_check;
				if (o->is_color) { // color
					do_check = ((m^f)&1) == 0? 1 : 0;
					do_append = !do_check;
				} else { // Illumina
					if (o->RF_max_dist) { // long-insert library
						do_check = do_append = 1;
					} else { // short-insert library
						do_check = m? 1 : 0;
						do_append = !do_check;
					}
				}
				if (do_check) { // check PE
					match_info_t *seq2 = (f&1)? seq - 1 : seq + 1; // get the other end
					if (seq2->last_seqid == longseqid) {
						pair_info_t *pp = d->pair + (((bit32_t)r[k].w2)>>1);
						if (f&1) { // this is the second end
							if (match_tell_dist((bit32_t)seq2->last1, coor, o->min_dist, o->max_dist, o->size_l,
												o->is_color, o->RF_max_dist))
								match_add_pair(pp, seq2->last1, ((bit64_t)x<<32)|coor, longseqid, seq2->mm1, mm);
							if (match_tell_dist((bit32_t)seq2->last2, coor, o->min_dist, o->max_dist, o->size_l,
												o->is_color, o->RF_max_dist))
								match_add_pair(pp, seq2->last2, ((bit64_t)x<<32)|coor, longseqid, seq2->mm2, mm);
						} else { // this is the first end
							if (match_tell_dist((bit32_t)seq2->last1, coor, o->min_dist, o->max_dist, o->size_r,
												o->is_color, o->RF_max_dist))
								match_add_pair(pp, ((bit64_t)x<<32)|coor, seq2->last1, longseqid, mm, seq2->mm1);
							if (match_tell_dist((bit32_t)seq2->last2, coor, o->min_dist, o->max_dist, o->size_r,
												o->is_color, o->RF_max_dist))
								match_add_pair(pp, ((bit64_t)x<<32)|coor, seq2->last2, longseqid, mm, seq2->mm2);
						}
					}
				}
				if (do_append) { // append the reads to the queue
					if (seq->last_seqid != longseqid) {
						seq->last_seqid = longseqid;
						seq->last2 = 0;
						seq->last_mm2 = 0;
						seq->last1 = ((bit64_t)x << 32) | coor;
						seq->last_mm1 = mm;
					} else if (seq->last1 != (((bit64_t)x << 32) | coor)) {
						seq->last2 = seq->last1;
						seq->last_mm2 = seq->last_mm1;
						seq->last1 = ((bit64_t)x << 32) | coor;
						seq->last_mm1 = mm;
					}
				}
			}
			++g_tot_hits;
			if (o->hits_fp) { // print the hits
				bit32_t t = alt_cal_mm(bit64_t(match>>shift));
				if (t < 2 && g_is_first_round == 1 && (int)seq->c[0] + seq->c[1] <= o->max_hits)
					gzprintf(o->hits_fp, "A\t%d\t%d\t%c\t%d\n", (bit32_t)r[k].w2, (coor>>1)+1, (coor&1)? '-' : '+', x>>24);
			}
			// store the match x if it is good enough
			if (x > seq->i2) continue; // x is no better than the second best, reject
			if (x > seq->i1) { // then i2 > x > i1
				if (seq->p1 != coor || seq->seqid1 != longseqid) { // should always stand??
					seq->i2 = x; seq->mm2 = mm;
				}
			} else if (x != seq->i1) { // then x < i1, set the second best as the previous best hit
				// Here is the source of a random error. If x == seq->i1 and the position
				// is not the same as the best, the hit will be ignored. This is not correct,
				// but should happen rarely. The probability that this may happen is 1e-7,
				// determined by the number of random bits in x. Currently, it is 24-bit.
				//
				// This error can be avoided by saving the position of the second best hit.
				// However, this is at the cost of more memory. I do not like this way.
				if (seq->p1 != coor || seq->seqid1 != longseqid) {
					seq->i2 = seq->i1;
					seq->i1 = x; seq->seqid1 = longseqid; seq->p1 = coor; seq->mm1 = mm;
				}
			}
		} // ~for(k)
	} // ~for(i)
	memset(list, 0, sizeof(ref_canlist_t));
}

// This is the core function for the search.
void match_search(const match_aux_t *o, match_data_t *d, const nst_bfa1_t *l, int longseqid, ref_canlist_t *list)
{
	bit64_t tmp_xor;
	int i, j;
	bit32_t n_reads;
	read_t query[4], mask[4], omask, ofilter;
	required_option_t oo;

	oo.is_color = o->is_color; oo.size_l = o->size_l; oo.size_r = o->size_r;
	oo.max_err10 = o->max_err10; oo.max_dist = o->max_dist; oo.min_dist = o->min_dist; oo.RF_max_dist = o->RF_max_dist;
	oo.methy_mode = o->methy_mode; oo.hits_fp = o->hits_fp; oo.max_hits = o->max_hits;
	n_reads = d->n_reads;
	mask[0] = mask[1] = mask[2] = mask[3] = query[0] = query[1] = query[2] = query[3] = 0;
	tmp_xor = longseqid ^ lsh_global_seed;
	for (i = 0; i != l->len; ++i) {
		bit64_t maski, seqi;
		maski = l->mask[i]; seqi = l->seq[i];
		for (j = 31; j >= 0; --j) {
			query[1] = query[0] = (query[0]<<2) | (seqi>>(j<<1) & 0x3); // forward word
			if (!oo.is_color) { // nucleotide space
				query[2] = (query[2]>>2) | (read_t(3 - (seqi>>(j<<1) & 0x3)) << ((oo.size_l-1)<<1)); // reverse word
				query[3] = (query[3]>>2) | (read_t(3 - (seqi>>(j<<1) & 0x3)) << ((oo.size_r-1)<<1)); // reverse word
				if (oo.methy_mode) methy_ref(oo.methy_mode, query, oo.size_l, oo.size_r, seqi, j);
			} else { // color space
				query[2] = (query[2]>>2) | (read_t(seqi>>(j<<1) & 0x3) << ((oo.size_l-1)<<1)); // reverse word
				query[3] = (query[3]>>2) | (read_t(seqi>>(j<<1) & 0x3) << ((oo.size_r-1)<<1)); // reverse word
			}
			mask[1] = mask[0] = (mask[0]<<2) | (maski>>(j<<1) & 0x3);
			mask[2] = (mask[2]>>2) | (read_t(maski>>(j<<1) & 0x3) << ((oo.size_l-1)<<1));
			mask[3] = (mask[3]>>2) | (read_t(maski>>(j<<1) & 0x3) << ((oo.size_r-1)<<1));
			if ((++g_processed_len) * 100 > g_tot_ref_len * g_processed_percent) {
				struct rusage r;
				double sec;
				getrusage(0, &r);
				sec = r.ru_utime.tv_sec + r.ru_stime.tv_sec + (r.ru_utime.tv_usec + r.ru_stime.tv_usec) / 1000000.0;
				fprintf(stderr, "[match_search] %3llu%% processed in %.3lf sec: %llu / %llu = %.3lf\n",
						g_processed_len * 100 / g_tot_ref_len, sec,
						g_tot_hits - g_last_tot_hits, g_n_candidates - g_last_n_candidates,
						(double)(g_tot_hits - g_last_tot_hits) / (g_n_candidates - g_last_n_candidates + 1));
				g_last_tot_hits = g_tot_hits; g_last_n_candidates = g_n_candidates;
				++g_processed_percent;
			}
			if ((mask[0]&0x3) == 0 || (mask[2]&0x3) == 0 || (mask[3]&0x3) == 0) continue; // the first or the last base is N
			if (list->n+8 >= MAX_CANDIDATE) process_candidate(&oo, longseqid, tmp_xor, d, list);
			// loop through all filters
			for (int m = 0; m != 2; ++m) { // 0 for forward strand, and 1 for reverse
				bit32_t coor = (((i << 5) | (31 - j)) << 1) | m;
				for (int f = 0; f != 4; ++f) { // loop through the four filters: read1F, read2F, read1R, read2R
					if (oo.size_r == 0 && (f&1)) continue; // this is not mate-pair read, skip read2F and read2R
					match_index_t *index = d->index + f;
					bit32_t k, key = lsh_hash_64(query[(m<<1)|(f&1)]>>index->shift_seed & index->filter); // the hash key
					k = key>>MATCH_INDEX_SHIFT32; // k is the indexed part of the key
					if ((d->index4[k>>1]&(1<<f<<((k&1)<<2))) == 0) continue; // no hit, skip. Here is the first cache miss!
					ref_candid_t *p = list->list + list->n;
					p->f = f; p->k = k; p->key = key;
					p->s = query[m<<1|(f&1)]; p->m = mask[m<<1|(f&1)];
					p->coor = coor;
					++(list->n);
				}
			}
		}
	}
	process_candidate(&oo, longseqid, tmp_xor, d, list);
}
// index lsh_data_t::sorted to accelerate search
static void match_index_sorted(int m, const bit128_t *a, bit32_t *b)
{
	int i;
	bit32_t curr, begin;
	if (m == 0) {
		fprintf(stderr, "[match_index_sorted] no reasonable reads are available. Exit!\n");
		exit(1);
	}
	for (i = 0; i != MATCH_INDEX_SIZE; ++i) b[i] = MA_NO_INDEX;
	curr = (bit32_t)(a[0].w2>>MATCH_INDEX_SHIFT64);
	begin = 0;
	for (i = 1; i != m; ++i) {
		if (curr != (bit32_t)(a[i].w2>>MATCH_INDEX_SHIFT64)) {
			b[curr] = begin | ((i - begin == 1)? 0x80000000u : 0);
			curr = (bit32_t)(a[i].w2>>MATCH_INDEX_SHIFT64);
			begin = i;
		}
	}
	// now come to the last element
	b[curr] = begin | ((i - begin == 1)? 0x80000000u : 0);
}
// generate "sorted" lists
static void match_index(match_aux_t *o, match_data_t *d, int start)
{
	int i;
	// allocate memory
	d->index4 = (bit8_t*)calloc(MATCH_INDEX_SIZE>>1, 1);
	for (i = 0; i != 4; ++i) {
		d->index[i].filter = o->filters[i+start];
		d->index[i].mask = o->masks[i+start];
		d->index[i].shift_seed = o->shift_seed[i+start];
		if (o->size_r || (i&1) == 0) {
			d->index[i].sorted = (bit128_t*)calloc(d->n_reads + 1, sizeof(bit128_t));
			d->index[i].index = (bit32_t*)calloc(MATCH_INDEX_SIZE, sizeof(bit32_t));
		}
	}
	// fill the sorted array. d->sorted must be in line with o->filters
	for (i = 0; i != 4; ++i) {
		int n = 0;
		match_index_t *index = d->index + i;
		if (o->size_r == 0 && (i&1)) continue;
		for (int j = i & 1; j < d->n_reads; j += 2) { // note "j = i&1" and "j += 2"
			 // discard poor quality reads
			if (alt_cal_mm(d->match[j].m) > MAX_MISMATCH) continue;
			bit64_t tmp = bit64_t(d->match[j].s>>index->shift_seed);
			index->sorted[n++] = bit128_t(tmp, (bit64_t)lsh_hash_64(tmp&index->filter)<<32 | j);
		}
		index->sorted[n] = 0xffffffffull;
		algo_sort(n, index->sorted); // the last element has not been sorted
		match_index_sorted(n, index->sorted, index->index);
		// fill d->index4
		int flag = 1<<i;
		for (int k = 0; k != MATCH_INDEX_SIZE; ++k)
			if (index->index[k] != MA_NO_INDEX)
				d->index4[k>>1] |= flag<<((k&1)<<2);
	}
}

void match_core(const char *fn_out, match_aux_t *o, match_data_t *d, FILE *fp, gzFile fp_bfq_l, gzFile fp_bfq_r)
{
	int i, k;
	nst_bfa1_t *l;
	ref_canlist_t *list;
	
	list = (ref_canlist_t*)calloc(1, sizeof(ref_canlist_t));
	g_tot_ref_len = nst_bfa_len(fp);
	rewind(fp);
	fprintf(stderr, "[match_core] Total length of the reference: %llu\n", g_tot_ref_len);
	gzrewind(fp_bfq_l); gzrewind(fp_bfq_r);
#ifdef _WIN32
	lsh_global_seed = rand();
#else
	lsh_global_seed = lrand48();
#endif
	if (o->methy_mode) maq_methy_modify(d, o);
	d->sum_len = 0;
	for (i = 0; i < o->n_filters; i += 4) {
		g_is_first_round = (i == 0)? 1 : 0;
		fprintf(stderr, "[match_core] round %d/%d...\n", i/4+1, o->n_filters/4);
		fprintf(stderr, "[match_core] making index...\n");
		match_index(o, d, i);
		rewind(fp);
		k = 0; g_processed_len = 0; g_processed_percent = 0;
		while ((l = nst_load_bfa1(fp)) != 0) { // read longseq
			if (o->hits_fp && i == 0) gzprintf(o->hits_fp, "B\t%s\n", l->name); // print reference names
			//if (!o->is_quiet)
			//	fprintf(stderr, "[match_core] processing sequence %s (%d bp)...\n", l->name, l->ori_len);
			if (i == 0) {
				d->sum_len += l->ori_len;
				if ((d->n_lname & 0xff) == 0) // then enlarge the memory block
					d->lname = (char**)realloc(d->lname, sizeof(char*) * (d->n_lname + 0x100));
				d->lname[d->n_lname++] = strdup(l->name);
			}
			match_search(o, d, l, k++, list);
			if (o->hits_fp && i == 0) gzprintf(o->hits_fp, "E\t%s\n", l->name); // print reference names
			nst_delete_bfa1(l);
		}
		// free index
		for (int j = 0; j != 4; ++j) {
			free(d->index[j].sorted); free(d->index[j].index);
		}
		free(d->index4);
	}
	free(list);
	fprintf(stderr, "[match_core] sorting the hits and dumping the results...\n");
	gzFile fpout = (strcmp(fn_out, "-") == 0)? gzdopen(STDOUT_FILENO, "w") : gzopen(fn_out, "w");
	match_data2map(fpout, fp, fp_bfq_l, fp_bfq_r, o, d);
	gzclose(fpout);
}
static int match_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   maq map [options] <out.map> <chr.bfa> <reads_1.bfq> [reads_2.bfq]\n\n");
	fprintf(stderr, "Options: -1 INT      length of the first read (<=%d) [0]\n", MAX_READLEN-1);
	fprintf(stderr, "         -2 INT      length of the second read (<=%d) [0]\n", MAX_READLEN-1);
	fprintf(stderr, "         -m FLOAT    rate of difference between reads and references [0.001]\n");
	fprintf(stderr, "         -e INT      maximum allowed sum of qualities of mismatches [70]\n");
	fprintf(stderr, "         -d FILE     adapter sequence file [null]\n");
	fprintf(stderr, "         -a INT      max distance between two paired reads [250]\n");
	fprintf(stderr, "         -A INT      max distance between two RF paired reads [0]\n");
//	fprintf(stderr, "         -i INT      min distance between two paired reads [0]\n"); // I do not want users to change this.
	fprintf(stderr, "         -n INT      number of mismatches in the first 24bp [2]\n");
	fprintf(stderr, "         -M c|g      methylation alignment mode [null]\n");
	fprintf(stderr, "         -u FILE     dump unmapped and poorly aligned reads to FILE [null]\n");
	fprintf(stderr, "         -H FILE     dump multiple/all 01-mismatch hits to FILE [null]\n");
	fprintf(stderr, "         -C INT      max number of hits to output. >512 for all 01 hits. [250]\n");
	fprintf(stderr, "         -s INT      seed for random number generator [random]\n");
//	fprintf(stderr, "         -P          set alternative mapping quality as pair-2-diff flag\n");
	fprintf(stderr, "         -W          disable Smith-Waterman alignment\n");
#ifndef MAQ_LONGREADS
	fprintf(stderr, "         -N          record mismatch positions (max read length<=55)\n");
#endif
	fprintf(stderr, "         -t          trim all reads (usually not recommended)\n");
	fprintf(stderr, "         -c          match in the colorspace\n");
	fprintf(stderr, "\n");
//	fprintf(stderr, "%lu,%lu,%lu,%lu\n", sizeof(read_t), sizeof(match_info_t), sizeof(pair_info_t), sizeof(mapping1_t));
	return 1;
}
// this is the true main function
int ma_match(int argc, char *argv[])
{
	FILE *fp_bfa, *tmpfp;
	gzFile fp_bfq_l, fp_bfq_r, hits_fp;
	int c, size_l, size_r, n_mismatch, is_sw = 1;
	float mut_rate = 0.001;
	int min_dist, max_dist, RF_max_dist, max_hits, is_color, is_trim, is_mm, is_p2diff, max_err10 = 7;
	match_aux_t *o;
	match_data_t *d;
	char adaptor[1024], dump_file[1024], methy_mode, *fn_out;
	long seed;
#ifdef _WIN32
	seed = time(0);
#else
	seed = lsh_hash_32(time(0)) ^ lsh_hash_32(getpid());
#endif

	min_dist = 0; max_dist = 250;
	RF_max_dist = 0; // for Illumina long insert-size library
	max_hits = 250;
	size_l = size_r = 0; n_mismatch = 2;
	is_color = is_trim = is_mm = is_p2diff = 0;
	methy_mode = 0;
	dump_file[0] = adaptor[0] = '\0';
	hits_fp = 0;
	while ((c = getopt(argc, argv, "1:2:n:m:a:i:d:ctM:u:e:s:NH:C:A:PW")) >= 0) {
		switch (c) {
		case '1': size_l = atoi(optarg); break;
		case '2': size_r = atoi(optarg); break;
		case 'e': max_err10 = int(atoi(optarg) / 10.0 + 0.5); break;
		case 'd':
			tmpfp = fopen(optarg, "r");
			assert(tmpfp);
			fscanf(tmpfp, "%s", adaptor);
			fclose(tmpfp);
			break;
		case 'n': n_mismatch = atoi(optarg); break;
		case 'm': mut_rate = atof(optarg); break;
		case 'a': max_dist = atoi(optarg); break;
		case 'A': RF_max_dist = atoi(optarg); break;
		case 'i': min_dist = atoi(optarg); break;
		case 'u': strcpy(dump_file, optarg); break;
		case 'C': max_hits = atoi(optarg); break;
		case 'H':
			hits_fp = gzopen(optarg, "w");
			assert(hits_fp);
			break;
		case 'c': is_color = 1; break;
		case 't': is_trim = 1; break;
		case 'P': is_p2diff = 1; break;
		case 's': seed = atoi(optarg); break;
		case 'N': is_mm = 1; break;
		case 'M': methy_mode = tolower(optarg[0]); break;
		case 'W': is_sw = 0; break;
		default: return match_usage();
		}
	}
	if (argc - optind < 3) return match_usage();
	if (n_mismatch < 1 || n_mismatch > 3) {
		fprintf(stderr, "[ma_match] please use 1, 2 or 3 with option -n. Abort!\n");
		return 1;
	}
	fp_bfa = fopen(argv[optind+1], "r");
	fp_bfq_l = gzopen(argv[optind+2], "r");
	fp_bfq_r = (optind < argc)? gzopen(argv[optind+3], "r") : 0;
	assert(fp_bfa); assert(fp_bfq_l);
	fn_out = strdup(argv[optind]);
	// core
#ifdef _WIN32
	srand(seed);
#else
	srand48(seed);
#endif
	// load reads into memory
	fprintf(stderr, "-- maq-%s\n", PACKAGE_VERSION);
	d = new_match_data();
	ma_init_match_data(d, fp_bfq_l, &size_l, fp_bfq_r, &size_r, is_trim, adaptor, hits_fp);
	assert(size_l >= 12 && (size_r == 0 || size_r >= 12));
	// initialize parameters
	o = new_match_aux(size_l, size_r, n_mismatch);
	o->q_rate = int(-4.343 * log(mut_rate) + 0.5); // set q_rate
	o->min_dist = min_dist;
	if (o->min_dist < size_l+1 || o->min_dist < size_r+1) {
		o->min_dist = size_l > size_r? size_l+1 : size_r+1;
		fprintf(stderr, "[ma_match] set the minimum insert size as %d.\n", o->min_dist);
	}
	o->max_dist = max_dist;
	o->RF_max_dist = RF_max_dist;
	o->is_color = is_color;
	o->is_mm = is_mm;
	o->is_sw = is_sw;
	o->is_p2diff = is_p2diff;
	o->max_err10 = max_err10;
	o->methy_mode = methy_mode;
	if (dump_file[0]) o->dump_file = strdup(dump_file);
	if (adaptor[0]) o->adapter = strdup(adaptor);
	o->hits_fp = hits_fp;
	o->max_hits = max_hits;
	// core function
	match_core(fn_out, o, d, fp_bfa, fp_bfq_l, fp_bfq_r);
	free(fn_out);
	gzclose(fp_bfq_l);
	if (fp_bfq_r) gzclose(fp_bfq_r);
	fclose(fp_bfa);
	delete_match_data(d);
	delete_match_aux(o);
	return 0;
}
