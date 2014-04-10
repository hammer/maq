#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "read.h"
#include "match.hh"
#include "algo.hh"

// Why "10"? I forget. Probably a wild estimation...
#define END_DIFF_QUAL 10

static float match_norm_c[8] = { 0.5, 0.5, 0.25, 0.11, 0.05, 0.021, 0.009, 0.004 };

static inline int ma_is_paired(const match_aux_t *o, bit32_t pos1, bit32_t pos2)
{
	if (!o->is_color) { // nucleotide space
		if ((pos1&1) == 0) { // the first one is on the forward strand
			if (pos2&1) { // the other on the reverse strand
				int dist = int(pos2>>1) - int(pos1>>1) + o->size_l;
				if (dist >= o->min_dist && dist <= o->max_dist) return 1;
			}
		} else { // the first one in on the reverse strand
			if ((pos2&1) == 0) { // the other on the forward strand
				int dist = int(pos1>>1) - int(pos2>>1) + o->size_r;
				if (dist >= o->min_dist && dist <= o->max_dist) return 1;
			}
		}
		if (o->RF_max_dist) { // long-insert library
			if ((pos1&1) == 1) { // the first is on the reverse
				if ((pos2&1) == 0) { // the second on the forward
					int dist = int(pos2>>1) - int(pos1>>1) + o->size_l;
					if (dist >= o->min_dist && dist <= o->RF_max_dist) return 1;
				}
			} else { // the first is on the forward
				if ((pos2&1) == 1) { // the second on the reverse
					int dist = int(pos1>>1) - int(pos2>>1) + o->size_r;
					if (dist >= o->min_dist && dist <= o->RF_max_dist) return 1;
				}
			}
		}
	} else { // color space
		if ((pos1&1) == 0 && (pos2&1) == 0) { // both on the forward strand
			int dist = int(pos1>>1) - int(pos2>>1) + o->size_r;
			if (dist >= o->min_dist && dist <= o->max_dist) return 1;
		} else if ((pos1&1) == 1 && (pos2&1) == 1) { // both on the reverse strand
			int dist = int(pos2>>1) - int(pos1>>1) + o->size_l;
			if (dist >= o->min_dist && dist <= o->max_dist) return 1;
		}
	}
	return 0;
}
inline int ma_make_pair(const match_aux_t *o, const match_info_t *m1, const match_info_t *m2, pair_info_t *p)
{
	if (p->i1 != MA_NO_MATCH || m1->i1 == MA_NO_MATCH || m2->i1 == MA_NO_MATCH) return 0;
	if (m1->seqid1 == m2->seqid1 && ma_is_paired(o, m1->p1, m2->p1)) {
		p->seqid1 = m1->seqid1;
		p->m11 = m1->mm1; p->m12 = m2->mm1;
		p->p11 = m1->p1;  p->p12 = m2->p1;
		p->i1 = (((m1->i1>>24)+(m2->i1>>24))<<24) | (m1->i1^m2->i1&0xffffff);
		return 1;
	}
	return 0;
}
// match_info_t::last1 MUST be set as the quality before using this function!
static inline int ma_cal_pair_qual(const pair_info_t *pair, const match_info_t *m1, const match_info_t *m2,
								   const bit8_t sl[], const bit8_t sr[], const match_aux_t *o)
{
	if (pair->i1 == MA_NO_MATCH) return -1;
	int q, q_alt = 100;
	if (pair->i2 != MA_NO_MATCH) q_alt = 0;
	q = ((pair->p11 == m1->p1)? m1->last1 : 0) + ((pair->p12 == m2->p1)? m2->last1 : 0);
	if (q_alt < q) q = q_alt;
	return (q <= 99)? q : 99;
}
static inline void ma_dim_end_diff(int size, bit8_t seq[], const matches_t &mm)
{
	int k;
	for (k = 0; k != size; ++k) {
		if ((mm>>(size-1-k)&1) == 0) break;
		if ((seq[k]&0x3f) > END_DIFF_QUAL)
			seq[k] = (seq[k]&0xc0) | END_DIFF_QUAL;
	}
	for (k = size - 1; k >= 0; --k) {
		if ((mm>>(size-1-k)&1) == 0) break;
		if ((seq[k]&0x3f) > END_DIFF_QUAL)
			seq[k] = (seq[k]&0xc0) | END_DIFF_QUAL;
	}
}
static inline int ma_cal_map_qual(int size, bit8_t seq[], const match_aux_t *o, match_info_t *match, const int qs[8])
{
	static int q_miss[] = { 0, 2, 4, 8 };
	int avg, k, tmp;
	int q1, m1, q2, m2, aln_qual, aln_alt_qual;
	int full_m1, full_m2;

	if (match->i1 == MA_NO_MATCH) return -1;
	// calculate sum of errors
	avg = q1 = q2 = full_m1 = full_m2 = 0;
	for (k = size - 1; k >= 0; --k) {
		int q = seq[k]&0x3f;
		if (q > o->q_rate) q = o->q_rate;
		avg += q;
		if (match->mm1>>(size-1-k)&1) { q1 += q; ++full_m1; }
		if (match->mm2>>(size-1-k)&1) { q2 += q; ++full_m2; }
	}
	avg = int((float)avg / size + 0.5);
	// calculate alignment quality
	// part 1: quality difference
	m1 = match->q>>8&0xff; m2 = match->q&0xff;
	if (match->i2 == MA_NO_MATCH) aln_qual = 99;
	else {
		aln_qual = q2 - q1 - o->log_n[match->c[m2]];
		if (m2 - m1 <= 1 && aln_qual > o->q_rate) aln_qual = o->q_rate;
	}
	// part 2: quality for missing hits
	// In principle, a better fomula should be: 4 + avg_diff * C^{24}_k * (k + 1)
	// C^{24}_k is merged to avg. This is about Q13-Q15 for each mismatch.
	if ((avg -= 13) < 0) avg = 0;
	tmp = o->n_mismatch + 1 - m1;
	aln_alt_qual = (tmp >= 0)? q_miss[o->n_mismatch] + avg * tmp + qs[tmp] : 0;
	// calculate final aln_qual
	if (aln_alt_qual < aln_qual) aln_qual = aln_alt_qual;
	if (aln_qual < 0) aln_qual = 0;
	if (aln_qual > 99) aln_qual = 99;
	// write last_mm2
	if (full_m1 > 15) full_m1 = 15;
	if (q1 > 255) q1 = 255;
	match->last_mm2 = q1<<8 | m1<<4 | full_m1;
	return aln_qual;
}
static void mapping_count_single(const match_data_t *d, int qs[], int is_paired)
{
	int i, n;
	for (i = 0; i != 8; ++i) qs[i] = 1; // pseudo-count
	for (i = n = 0; i != d->n_reads; ++i) {
		if (is_paired == 0 && (i&1)) continue;
		match_info_t *m = d->match + i;
		if (m->i1 != MA_NO_MATCH && (m->q>>8&0xff) <= 1) {
			++n;
			if (m->i2 != MA_NO_MATCH) {
				int t = (m->q&0xff) - (m->q>>8&0xff);
				if (t < 0) t = 0;
				++qs[t];
			}
		}
	}
	for (i = 0; i != 8; ++i) {
		qs[i] = int(-4.343 * log((double)qs[i]/n) - 0.5);
		if (qs[i] < 0) qs[i] = 0;
	}
	fprintf(stderr, "[mapping_count_single] %d, %d, %d, %d\n", qs[0], qs[1], qs[2], qs[3]);
}
static inline void ma_set_dist_flag(const pair_info_t *p, match_info_t *m1, match_info_t *m2, int dist)
{
	if (p->i1 != MA_NO_MATCH) {
		m1->last_seqid = dist;
		m2->last_seqid = -dist;
		m1->last2 = m2->last2 = ((p->p11>>1 < p->p12>>1)? 1<<((p->p11&1)<<1|(p->p12&1)) : 1<<((p->p12&1)<<1|(p->p11&1)))
			| PAIRFLAG_PAIRED;
	} else {
		m1->last_seqid = m2->last_seqid = 0; // set dist as zero, by default
		if (m1->i1 != MA_NO_MATCH || m2->i1 != MA_NO_MATCH) {
			if (m1->i1 == MA_NO_MATCH) {
				m2->last2 = PAIRFLAG_NOMATCH;
			} else if (m2->i1 == MA_NO_MATCH) {
				m1->last2 = PAIRFLAG_NOMATCH;
			} else { // then both have a match
				if (m1->seqid1 != m2->seqid1) {
					m1->last2 = m2->last2 = PAIRFLAG_DIFFCHR;
				} else { // both matches are on the same chr
					m1->last_seqid = dist;
					m2->last_seqid = -dist;
					if (dist > 0) { // m1 has smaller coordinate
						m1->last2 = m2->last2 = 1<<((m1->p1&1)<<1 | (m2->p1&1));
					} else { // m2 has smaller coordinate
						m1->last2 = m2->last2 = 1<<((m2->p1&1)<<1 | (m1->p1&1));
					}
				}
			}
		}
	}
}

static inline int cal_mm(bit64_t x) {
	x = ((x & 0xAAAAAAAAAAAAAAAAull) >> 1)  + (x & 0x5555555555555555ull);
	x = ((x & 0xCCCCCCCCCCCCCCCCull) >> 2)  + (x & 0x3333333333333333ull);
	x = ((x & 0xF0F0F0F0F0F0F0F0ull) >> 4)  + (x & 0x0F0F0F0F0F0F0F0Full);
	x = ((x & 0xFF00FF00FF00FF00ull) >> 8)  + (x & 0x00FF00FF00FF00FFull);
	x = ((x & 0xFFFF0000FFFF0000ull) >> 16) + (x & 0x0000FFFF0000FFFFull);
	return ((x&0xFFFFFFFF00000000ull)>> 32) + (x & 0x00000000FFFFFFFFull);
}
template<class TYPE>
inline int cal_mm(const dword_t<TYPE> &x) {
	return cal_mm(x.w1) + cal_mm(x.w2);
}

/*
  In this function, match_info_t will be modified:
    q          -> number of seed mismatches
    m          -> whether the pair is 2-diff (a flag, 0 or 1)
    last1      -> mapping quality
	last_mm1   -> single end mapping quality
	last2      -> pair flag
	last_mm2   -> number of mismatch and sum of mismatch qualities
	last_seqid -> offset of the mate
 */
void match_data2map(gzFile fpout, FILE *fp_bfa, gzFile fp_bfq_l, gzFile fp_bfq_r, const match_aux_t *o, match_data_t *d)
{
	int i, size[2] = { o->size_l, o->size_r };
	longreads_t *lr;
	dword_t<bit64_t> *sorted;
	int qs[8], n, n_mapped;
	extern void maq_indel_pe(FILE *fp_bfa, const match_aux_t *o, const longreads_t *lr, match_data_t *d);
	// *** STEP 1: load reads again
	lr = ma_load_reads(fp_bfq_l, o->size_l, fp_bfq_r, o->size_r);
	// *** STEP 2: calculate read quality
	for (int k = 0; k != d->n_reads; ++k) { // calculate number of mismatches in the seed
		match_info_t *match = d->match + k;
		int m1, m2, s;
		s = size[k&1]; m1 = m2 = 0;
		for (int j = 0; j < SEED_LENGTH && j < s; ++j) {
			if (match->mm1>>(s-1-j)&1) ++m1;
			if (match->mm2>>(s-1-j)&1) ++m2;
		}
		match->q = m1<<8 | m2; // match->q is reused as the number of mismatch in the seed
		match->m = 0;
	}
	mapping_count_single(d, qs, size[1]);
	for (int k = 0; k != d->n_reads; ++k) {
		match_info_t *match = d->match + k;
		// normalize match->c[]
		for (int j = 0; j <= MAX_MISMATCH; ++j) {
			int tmp = int(match->c[j] / (0.5 * o->n_filters) / match_norm_c[j] + 0.5);
			match->c[j] = (tmp < 256)? tmp : 255;
		}
		match->last1 = ma_cal_map_qual(size[k&1], lr->seq[k], o, match, qs); // ->last1 is the SE mapping quality
		match->last_mm1 = match->last1; // ->last_mm1 is the single-end mapping quality
		match->last_seqid = 0; // ->last_seqid should be the offset of the other pair
	}
	// *** STEP 3: calculate pair quality
	if (size[1]) {
		int n1, n2, ql1, ql2, qh1, qh2; // for statistics
		n = n1 = n2 = ql1 = ql2 = qh1 = qh2 = 0;
		for (int k = 0; k != d->n_reads>>1; ++k) {
			pair_info_t *pair = d->pair + k;
			match_info_t *match = d->match + (k<<1);
			int dist;
			{ // calculate the distance between the pair
				int p1, p2;
				if (pair->i1 != MA_NO_MATCH) {
					p1 = (pair->p11&1)? pair->p11>>1 : (pair->p11>>1) - (size[0] - 1);
					p2 = (pair->p12&1)? pair->p12>>1 : (pair->p12>>1) - (size[1] - 1);
				} else { // match or match+1 may be unmapped, but this will be solved later
					p1 = (match->p1&1)? match->p1>>1 : (match->p1>>1) - (size[0] - 1);
					p2 = (match[1].p1&1)? match[1].p1>>1 : (match[1].p1>>1) - (size[1] - 1);
				}
				dist = p2 - p1 + (p2 > p1? 1 : -1);
			}
			if (ma_make_pair(o, match, match+1, pair)) ++n; // find missed pairs
			ma_set_dist_flag(pair, match, match+1, dist); // actually set ->last2 as the flag and ->last_seqid as the dist
			if (pair->i1 == MA_NO_MATCH) continue; // no paired match
			int pair_qual = ma_cal_pair_qual(pair, match, match+1, lr->seq[k<<1], lr->seq[k<<1|1], o);
			// tell whether the pair is at least 2-diff unique
			if (pair->i2 == MA_NO_MATCH) {
				int snm, s;
				s = (o->size_l <= SEED_LENGTH)? 0 : o->size_l - SEED_LENGTH;
				snm = cal_mm(pair->m11>>s);
				s = (o->size_r <= SEED_LENGTH)? 0 : o->size_r - SEED_LENGTH;
				snm += cal_mm(pair->m12>>s);
				match->m = match[1].m = (snm < 2)? 1 : 0;
			}
			// update single ends if the pair can be aligned
			if (pair->seqid1 != match->seqid1 || pair->p11 != match->p1) { // the read is moved!
				if (match->last1 < 30) ++ql1;
				else ++qh1;
				++n1; match->p1 = pair->p11; match->seqid1 = pair->seqid1; match->mm1 = pair->m11;
				match->last_mm1 = match->last1 = 0; // single end Q should be 0, then.
			}
			if ((int)match->last1 < pair_qual) match->last1 = pair_qual; // update the mapping quality
			++match; // match now is the other end
			if (pair->seqid1 != match->seqid1 || pair->p12 != match->p1) { // the read is moved!
				if (match->last1 < 30) ++ql2;
				else ++qh2;
				++n2; match->p1 = pair->p12; match->seqid1 = pair->seqid1; match->mm1 = pair->m12;
				match->last_mm1 = match->last1 = 0; // 0!
			}
			if ((int)match->last1 < pair_qual) match->last1 = pair_qual; // update the mapping quality
		}
		fprintf(stderr, "[match_data2mapping] %d pairs are added afterwards.\n", n);
		fprintf(stderr, "[match_data2mapping] (%d, %d) reads are moved to meet paired-end requirement.\n", n1, n2);
		fprintf(stderr, "[match_data2mapping] quality counts of the first reads: (%d, %d); second reads (%d, %d)\n",
				qh1, ql1, qh2, ql2);
	}
	// we can delete match->pair now
	free(d->pair); d->pair = 0;
	maq_indel_pe(fp_bfa, o, lr, d);
	// *** STEP 4: trim adapter and sort the positions
	if (o->adapter) ma_trim_adapter(o->adapter, d, lr);
	// *** set coorninate of unmapped read in a pair
	if (o->size_r) {
		for (i = 0; i != d->n_reads>>1; ++i) {
			match_info_t *m1 = d->match + (i<<1);
			match_info_t *tmp, *m2 = m1 + 1;
			int do_modify = 0;
			if (m1->i1 == MA_NO_MATCH && m2->i1 != MA_NO_MATCH) do_modify = 1;
			else if (m2->i1 == MA_NO_MATCH && m1->i1 != MA_NO_MATCH) {
				do_modify = 1;
				tmp = m1; m1 = m2; m2 = tmp;
			}
			if (do_modify) {
				m1->i1 = m2->i1; m1->q = 0; m1->m = 0;
				m1->last1 = 0; m1->last_mm1 = 0; m1->last_seqid = 0;
				m1->last2 = PAIRFLAG_SW | PAIRFLAG_NOMATCH; m1->last_mm2 = 0;
				m1->p1 = m2->p1; m1->seqid1 = m2->seqid1;
			}
		}
	}
	// *** sort the positions
	sorted = (dword_t<bit64_t>*)malloc(sizeof(dword_t<bit64_t>) * d->n_reads);
	for (i = n_mapped = 0; i != d->n_reads; ++i) {
		match_info_t *match = d->match + i;
		if (match->i1 != MA_NO_MATCH) {
			ma_dim_end_diff(size[i&1], lr->seq[i], match->mm1); // dim the quality for mismatches at the end of a read
			bit64_t pos = (((match->p1>>1) + 1 - size[i&1]) << 1) | (match->p1 & 1); // this is the start position
			sorted[n_mapped++] = dword_t<bit64_t>((pos<<32) | i, match->seqid1);
		}
	}
	algo_sort(n_mapped, sorted);
	// *** dump no-match and poor-match reads
	if (o->dump_file) {
		fprintf(stderr, "-- Dumping unmapped or poorly mapped reads\n");
		FILE *fp_dump = fopen(o->dump_file, "w");
		for (i = 0; i != d->n_reads; ++i) {
			if (!size[1] && (i&1)) continue;
			match_info_t *match = d->match + i;
			if (match->i1 == MA_NO_MATCH) {
				int s = size[i&1];
				bit8_t *r = lr->seq[i];
				fprintf(fp_dump, "%s\t99\t", lr->name[i>>1]);
				for (int j = 0; j != s; ++j) {
					if (r[j] == 0) fputc('N', fp_dump);
					else fputc("ACGT"[r[j]>>6], fp_dump);
				}
				fputc('\t', fp_dump);
				for (int j = 0; j != s; ++j)
					fputc((r[j]&0x3f) + 33, fp_dump);
				fputc('\n', fp_dump);
			}
		}
		fclose(fp_dump);
	}
	// *** STEP 5: dumping results
	maqmap_t *mm = maq_new_maqmap();
	char new_name[1024];
	mm->n_ref = d->n_lname; mm->ref_name = d->lname; mm->n_mapped_reads = n_mapped;
	maqmap_write_header(fpout, mm);
	free(mm);
	for (i = n = 0; i != n_mapped; ++i) {
		maqmap1_t *m1 = (maqmap1_t*)calloc(1, sizeof(maqmap1_t));
		int k = sorted[i].w1 & 0xfffffffful;
		match_info_t *match = d->match + k; // get the read
		m1->size = size[k&1]; // size of the sequence
		memcpy(m1->seq, lr->seq[k], m1->size); // copy sequence
		memcpy(m1->c, match->c, 2); // only 0-1 are stored
		m1->seqid = match->seqid1; // chromosome ID
		m1->flag = (bit8_t)match->last2; // bit-wise paired information
		m1->pos = (((match->p1>>1) + 1 - m1->size) << 1) | (match->p1 & 1); // position and strand
		if (m1->flag & PAIRFLAG_SW) { // aligned be Smith-Waterman alignment
			m1->map_qual = match->last1; // indel position (0 for no indel)
			m1->alt_qual = match->last_mm1&0xff; // mapping quality of its mate
			m1->seq[MAX_READLEN-1] = match->last_mm1>>8&0xff; // length of indels (>0 for ins, <0 for del, ==0 for no indel)
		} else {
			if (o->is_color) m1->map_qual = ((match->q>>8&0xff) <= 3)? match->last1 : bit64_t(match->last_mm1);
			else m1->map_qual = ((match->q>>8&0xff) <= 2)? match->last1 : bit64_t(match->last_mm1); // mapping quality
			m1->alt_qual = bit64_t(((k&1) == 0)? (((match+1)->last_mm1 < match->last_mm1)?
												  (match+1)->last_mm1 : match->last_mm1)
								   : (((match-1)->last_mm1 < match->last_mm1)? (match-1)->last_mm1 : match->last_mm1)); // altq
			if (o->is_p2diff) m1->alt_qual = match->m&0x1;
			m1->seq[MAX_READLEN-1] = (bit8_t)bit64_t(match->last_mm1); // single end mapping quality
		}
		m1->dist = match->last_seqid; // distance between pairs
		m1->info1 = m1->info2 = 0xff;
		if (match->i1 != MA_NO_MATCH) {
			m1->info1 = match->last_mm2&0xff;
			m1->info2 = match->last_mm2>>8&0xff;
		}
		{ // make read name
			int l = strlen(lr->name[k>>1]);
			strcpy(new_name, lr->name[k>>1]);
			if (new_name[l-1] == '1' && new_name[l-2] == '/') new_name[l-1] += k&1;
			if (l - MAX_NAMELEN + 1 > 0) strcpy(m1->name, new_name + l - MAX_NAMELEN + 1);
			else strcpy(m1->name, new_name);
		}
		if (m1->flag&(PAIRFLAG_PAIRED|PAIRFLAG_SW)) ++n; // count if paired
		if (m1->pos & 1) { // then reverse the read
			if (!o->is_color) { // nucleotide space
				for (int k = 0; k < m1->size>>1; ++k) {
					bit8_t tmp = m1->seq[k];
					m1->seq[k] = ((3 - (m1->seq[m1->size-k-1]>>6)) << 6) | (m1->seq[m1->size-k-1]&0x3f);
					m1->seq[m1->size-k-1] = ((3 - (tmp>>6)) << 6) | (tmp&0x3f);
				}
				if (m1->size&1) m1->seq[m1->size>>1] = ((3 - (m1->seq[m1->size>>1]>>6)) << 6) | (m1->seq[m1->size>>1]&0x3f);
			} else { // color space, no need to do complementary
				for (int k = 0; k < m1->size>>1; ++k) {
					bit8_t tmp = m1->seq[k];
					m1->seq[k] = m1->seq[m1->size-k-1];
					m1->seq[m1->size-k-1] = tmp;
				}
			}
		}
#ifndef MAQ_LONGREADS
		if (o->is_mm) { // record mismatch positions
			bit64_t *p = (bit64_t*)(m1->seq + 55);
			if (m1->pos&1) { // reverse strand, reverse match->mm1
				bit64_t x = match->mm1;
				*p = 0;
				for (int k = 0; k < m1->size; ++k, x >>= 1)
					if (x&1) *p |= 1ull<<(m1->size-1-k);
			} else *p = match->mm1;
		}
#endif
		gzwrite(fpout, m1, sizeof(maqmap1_t));
		free(m1);
	}
	// *** STEP 7: deallocate the memory
	delete_longreads(lr);
	free(sorted);
	fprintf(stderr, "[match_data2mapping] %d out of %d raw reads are mapped with %d in pairs.\n", n_mapped, d->n_reads, n);
	if (o->size_r) fprintf(stderr, "-- (total, isPE, mapped, paired) = (%d, 1, %d, %d)\n", d->n_reads, n_mapped, n);
	else fprintf(stderr, "-- (total, isPE, mapped, paired) = (%d, 0, %d, 0)\n", d->n_reads>>1, n_mapped);
}
