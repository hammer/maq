#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "match.hh"
#include "bfa.h"
#include "main.h"
#include "algo.hh"
#include "assemble.h"
#include "stdhash.hh"

#define WIN_SIZE 7

typedef struct
{
	int refN, baseN, ref_len;
	int called[256], hom[256], het[256];
	int called10[26], hom10[26], het10[26];
} assemble_stat_t;

typedef struct
{
	float esum[4], fsum[4];
	int b[3], bar_e[2];
	bit32_t c[4];
	bit32_t depth, avg01match, unique_ratio;
	bit32_t mapQ_max;
} base_call_aux_t;

typedef struct {
	int is_alt, max_err, min_q;
	int is_single, is_pair_only;
	int max_mm;
	int n_hap;
	float hetero_rate;
	float theta;
	float eta;
} assemble_opt_t;

inline bool operator < (const assemble_posinfo_t &a, const assemble_posinfo_t &b)
{
	return (a.info < b.info);
}

/** helper structure related to doing statistics */
static inline void add_to_stat(bit64_t info, assemble_stat_t *as)
{
	++as->ref_len;
	if (info>>60 == 0xf) ++as->refN;
	if ((info>>56&0xf) == 0xf) ++as->baseN;
	else {
		int j, k;
		j = info>>48&0xff;
		if (j > 99) j = 99;
		k = j/10;
		++as->called[j]; ++as->called10[k];
		if (nst_nt16_count_table[info>>56&0xf] == 2) { ++as->het[j]; ++as->het10[k]; }
		else if (nst_nt16_count_table[info>>56&0xf] == 1 && info>>60 != (info>>56&0xf)) {
			++as->hom[j]; ++as->hom10[k];
		}
	}
}
static void print_stat(assemble_stat_t *as)
{
	int sum_called, sum_hom, sum_het;
	fprintf(stderr, "S0 reference length: %d\n", as->ref_len);
	fprintf(stderr, "S0 number of gaps in the reference: %d\n", as->refN);
	fprintf(stderr, "S0 number of uncalled bases: %d (%.2f)\n", as->baseN, (float)(as->baseN-as->refN)/(as->ref_len-as->refN));
	fprintf(stderr, "CC\n");
	fprintf(stderr, "CC %4s  %10s %10s %10s %10s %10s %10s\nCC\n", "qual", "cover", "cover_sum",
			"hom", "hom_sum", "het", "het_sum");
	sum_called = sum_hom = sum_het = 0;
	for (int i = 99; i >= 0; --i) {
		sum_called += as->called[i];
		sum_hom += as->hom[i];
		sum_het += as->het[i];
		fprintf(stderr, "S1 %.2d    %10d %10d %10d %10d %10d %10d\n", i, as->called[i], sum_called,
				as->hom[i], sum_hom, as->het[i], sum_het);
	}
	fprintf(stderr, "CC\nCC Summary\nCC\n");
	sum_called = sum_hom = sum_het = 0;
	for (int i = 9; i >= 0; --i) {
		sum_called += as->called10[i];
		sum_hom += as->hom10[i];
		sum_het += as->het10[i];
		fprintf(stderr, "S2 %.2d-%.2d %10d %10d %10d %10d %10d %10d\n", i*10, i*10+9, as->called10[i], sum_called,
				as->hom10[i], sum_hom, as->het10[i], sum_het);
	}
}

/** collect information for base calling */
base_call_aux_t *assemble_cns_collect(assemble_pos_t *ap, const assemble_aux_t *aa)
{
	base_call_aux_t *b = (base_call_aux_t*)calloc(1, sizeof(base_call_aux_t));
	b->b[0] = b->b[1] = b->b[2] = -1; // uncalled
	if (ap->n_bases == 0) return b;
	int j, k, depth01, c01, w[8];
	algo_sort(ap->n_bases, ap->bases);
	depth01 = c01 = 0;
	for (k = 0; k != 8; ++k) w[k] = 0;
	b->mapQ_max = 0;
	for (j = ap->n_bases - 1; j >= 0; --j) {
		bit32_t info = ap->bases[j].info;
		if (info>>24 < 4 && (info>>8&0x3f) != 0) info = 4<<24 | info&0xffffff;
		k = info>>16&7;
		if (info>>24 > 0) {
			b->esum[k&3] += aa->fk[w[k]] * (info>>24);
			b->fsum[k&3] += aa->fk[w[k]];
			if (w[k] < 0xff) ++w[k];
			++b->c[k&3];
		}
		if ((info>>19&0x3) <= 1) {
			++depth01; c01 += ap->bases[j].c[info>>19&0x3];
		}
		if (b->mapQ_max < (info&0x7f)) b->mapQ_max = info&0x7f;
	}
	if (b->mapQ_max > 63) b->mapQ_max = 63;
	// find the best three bases;
	float max[3];
	max[0] = max[1] = max[2] = 0.0;
	for (j = 0; j != 4; ++j) {
		if (b->esum[j] > max[0]) {
			max[2] = max[1]; b->b[2] = b->b[1];
			max[1] = max[0]; b->b[1] = b->b[0];
			max[0] = b->esum[j]; b->b[0] = j;
		} else if (b->esum[j] > max[1]) {
			max[2] = max[1]; b->b[2] = b->b[1];
			max[1] = b->esum[j]; b->b[1] = j;
		} else if (b->esum[j] > max[2]) {
			max[2] = b->esum[j]; b->b[2] = j;
		}
	}
	b->bar_e[0] = (b->b[0] >= 0)? (int)(b->esum[b->b[0]] / b->fsum[b->b[0]] + 0.5) : 0;
	b->bar_e[1] = (b->b[1] >= 0)? (int)(b->esum[b->b[1]] / b->fsum[b->b[1]] + 0.5) : 0;
	b->depth = (ap->n_bases >= 0xff)? 0xff : ap->n_bases;
	b->avg01match = bit32_t(16.0 * c01 / depth01 + 0.5);
	if (b->avg01match > 0xfffu) b->avg01match = 0xfff;
	return b;
}
/** base calling and depth calling */
inline bit64_t assemle_cns_call(assemble_pos_t *ap, const assemble_aux_t *aa, bit32_t ref_base)
{
	if (ap->n_bases == 0)
		return ((bit64_t)ref_base<<60) | (0xfull<<56) | (0xffull<<40); // not covered
	base_call_aux_t *bb = assemble_cns_collect(ap, aa);
	bit64_t info = 0;
	int w[3], c[3];
	float q[3];
	bit32_t b[3];
	c[0] = bb->b[0] >= 0 ? bb->c[bb->b[0]] : 0;
	c[1] = bb->b[1] >= 0 ? bb->c[bb->b[1]] : 0;
	if ((c[2] = c[0] + c[1]) > 255) {
		c[0] = int(255.0 * c[0] / c[2] + 0.5);
		c[1] = int(255.0 * c[1] / c[2] + 0.5);
		c[2] = 255;
	}
	q[0] = (c[1] ? bb->esum[bb->b[1]] : 0) + aa->coef[bb->bar_e[1]<<16|c[2]<<8|c[1]];
	q[1] = (c[0] ? bb->esum[bb->b[0]] : 0) + aa->coef[bb->bar_e[0]<<16|c[2]<<8|c[0]];
	q[2] = aa->q_r - 4.343 * aa->lhet[c[1]<<8|c[0]];
	//if (c[1]) fprintf(stderr, "* %f\t%f\t%d\n", bb->esum[bb->b[1]], aa->coef[bb->bar_e[1]<<16|c[2]<<8|c[1]], bb->bar_e[1]);
	//if (c[1] >= 3 && bb->esum[bb->b[1]] >= 80.0) fprintf(stderr, "+ %f\t%f\t%f\n", q[0], q[1], q[2]);
	if (q[0] < 0.0) q[0] = 0.0;
	if (q[1] < 0.0) q[1] = 0.0;
	b[0] = 1<<bb->b[0];
	b[1] = (bb->b[1] >= 0)? 1<<bb->b[1] : 0xf;
	b[2] = b[0] | b[1];
	// actually, the following codes calculate the "rank" of q[0..2]
	w[0] = w[1] = w[2] = -1;
	if (q[0] <= q[1]) { // MUST use "<=" here
		if (q[1] < q[2]) { w[0] = 0; w[1] = 1; w[2] = 2; }
		else { // q[2] <= q[1]
			w[2] = 1;
			if (q[0] < q[2]) { w[0] = 0; w[1] = 2; }
			else { w[0] = 2; w[1] = 0; } 
		}
	} else { // q[1] < q[0]
		if (q[0] < q[2]) { w[0] = 1; w[1] = 0; w[2] = 2; }
		else { // q[2] <= q[0]
			w[2] = 0;
			if (q[1] < q[2]) { w[0] = 1; w[1] = 2; }
			else { w[0] = 2; w[1] = 1; }
		}
	}
	// change to integer and do base-calling
	bit32_t iqual[2];
	iqual[0] = (int)(q[w[1]] - q[w[0]] + 0.5);
	iqual[1] = (int)(q[w[2]] - q[w[1]] + 0.5);
	if (iqual[0] >= 0xff) iqual[0] = 0xff;
	if (iqual[1] >= 0xff) iqual[1] = 0xff;
	if (w[0] != 2 && w[0] == 1) { // in case weird things happen, and they do happen...
		w[0] = 0; w[1] = 1; w[2] = 2; iqual[0] = 3; iqual[1] = 3;
	}
	info = (ref_base<<28) | (b[w[0]]<<24) | (iqual[0]<<16) | (b[w[1]]<<12) | (b[w[2]]<<8) | iqual[1]; // base call
	info <<= 32;
	info |= (bb->avg01match << 20) | (bb->mapQ_max<<8) | bb->depth; // depth call
	free(bb);
	return info;
}
void assemble_core(assemble_aux_t *aa)
{
	nst_bfa1_t *l;
	int k, seqid, half_win = WIN_SIZE/2;
	bit64_t rbuf[WIN_SIZE];
	FILE *fp_bfa = aa->fp_bfa;
	gzFile fpout = aa->fpout, fp_map = aa->fp_map;
	rolling_buf_t *buf = (rolling_buf_t*)calloc(1, sizeof(rolling_buf_t));
	assemble_pos_t *pos = (assemble_pos_t*)calloc(1, sizeof(assemble_pos_t));
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names
	assemble_stat_t *as = (assemble_stat_t*)calloc(1, sizeof(assemble_stat_t));

	hash_map_char<int> *hash = new hash_map_char<int>;
	for (k = 0; k != mm->n_ref; ++k) hash->insert(mm->ref_name[k], k);
	buf->buf = (maqmap1_t*)calloc(ROLLING_BUF_SIZE, sizeof(maqmap1_t));
	k = mm->n_mapped_reads;
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		if (!hash->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		fprintf(stderr, "[assemble_core] Processing %s (%d bp)...\n", l->name, l->ori_len);
		memset(rbuf, 0, sizeof(bit64_t) * WIN_SIZE);
		int len = strlen(l->name) + 1;
		gzwrite(fpout, &len, sizeof(int));
		gzwrite(fpout, l->name, len);
		gzwrite(fpout, &l->ori_len, sizeof(int));
		for (int i = 0; i != l->len; ++i) {
			bit64_t word = l->seq[i];
			bit64_t mask = l->mask[i];
			for (int j = 31; j >= 0; --j) {
				int coor = (i<<5) | (31 - j);
				if (coor >= l->ori_len) break;
				bit32_t ref_base = (mask>>(j<<1)&3)? 1<<(word>>(j<<1)&3) : 0xf;
				assemble_get_pos(seqid, coor, fp_map, buf, pos, aa->max_mm, aa->max_err,
								 aa->min_q, aa->is_single, aa->is_pair_only);
				bit64_t info = assemle_cns_call(pos, aa, ref_base);
				add_to_stat(info, as);
				rbuf[coor%WIN_SIZE] = info;
				if (coor >= half_win) {
					int min = 256;
					bit64_t x;
					info = rbuf[(coor-half_win)%WIN_SIZE];
					for (k = coor + 1 - WIN_SIZE; k < coor - half_win; ++k) {
						x = rbuf[(k+WIN_SIZE)%WIN_SIZE];
						if (min > int(x>>48&0xff)) min = x>>48&0xff;
					}
					for (k = coor - half_win + 1; k <= coor; ++k) {
						x = rbuf[k%WIN_SIZE];
						if (min > int(x>>48&0xff)) min = x>>48&0xff;
					}
					if (min >= 62) min = 62;
					min >>= 1;
					info |= min << 15;
					// something like little endianness, but not exactly
					bit32_t high_info = (bit32_t)(info>>32);
					bit32_t low_info = (bit32_t)info;
					gzwrite(fpout, &low_info, sizeof(bit32_t));
					gzwrite(fpout, &high_info, sizeof(bit32_t));
				}
			}
		}
		for (int coor = (l->ori_len > half_win)? l->ori_len - half_win : 0; coor < l->ori_len; ++coor) {
			bit64_t info = rbuf[coor%WIN_SIZE];
			bit32_t high_info = (bit32_t)(info>>32);
			bit32_t low_info = (bit32_t)info;
			gzwrite(fpout, &low_info, sizeof(bit32_t));
			gzwrite(fpout, &high_info, sizeof(bit32_t));
		}
		nst_delete_bfa1(l);
	}
	print_stat(as);
	delete hash;
	maq_delete_maqmap(mm);
	free(buf->buf); free(buf);
	free(pos->bases); free(pos); free(as);
}

int ma_assemble(int argc, char *argv[])
{
	assemble_aux_t *aa = assemble_parse_opt(argc, argv);
	assemble_core(aa);
	delete_assemble_aux(aa);
	return 0;
}
