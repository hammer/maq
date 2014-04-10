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
#include "glf.h"

typedef struct {
	float esum[4], fsum[4];
	int bar_e[4];
	bit32_t c[4];
	bit32_t mapQ_max;
} glf_call_aux_t;

inline bool operator < (const assemble_posinfo_t &a, const assemble_posinfo_t &b)
{
	return (a.info < b.info);
}

/** collect information for base calling */
static glf1_t *glfgen1_core(assemble_pos_t *ap, const assemble_aux_t *aa, bit8_t ref_base)
{
	glf_call_aux_t *b;
	glf1_t *g = (glf1_t*)calloc(1, sizeof(glf1_t));

	g->ref_base = ref_base;
	if (ap->n_bases == 0) return g;

	int i, j, k, w[8];
	b = (glf_call_aux_t*)calloc(1, sizeof(glf_call_aux_t));
	algo_sort(ap->n_bases, ap->bases);
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
		if (b->mapQ_max < (info&0x7f)) b->mapQ_max = info&0x7f;
	}
	for (j = 0; j != 4; ++j) {
		float tmp1, tmp2;
		tmp1 = tmp2 = 0.0;
		for (k = 0; k != 4; ++k) {
			if (k == j) continue;
			tmp1 += b->esum[k]; tmp2 += b->fsum[k];
		}
		b->bar_e[j] = tmp2 > 0.0? (int)(tmp1 / tmp2 + 0.5) : 0;
	}

	float p[16];
	int c;
	// rescale ->c[]
	for (j = c = 0; j != 4; ++j) c += b->c[j];
	if (c > 255) {
		for (j = 0; j != 4; ++j) b->c[j] = (int)(254.0 * b->c[j] / c + 0.5);
		for (j = c = 0; j != 4; ++j) c += b->c[j];
	}
	for (j = 0; j != 4; ++j) {
		// homozygous
		float sum = 0.0;
		for (k = 0; k != 4; ++k) {
			if (j == k) continue;
			sum += b->esum[k] + (b->bar_e[k]? aa->coef[b->bar_e[k]<<16|c<<8|b->c[k]] : 0);
		}
		p[j<<2|j] = sum;
		// heterozygous
		for (k = j + 1; k < 4; ++k) {
			for (i = 0, sum = 0.0; i != 4; ++i) {
				if (i == j || i == k) continue;
				sum += b->esum[i] + (b->bar_e[i]? aa->coef[b->bar_e[i]<<16|c<<8|b->c[i]] : 0);
			}
			p[j<<2|k] = p[k<<2|j] = -4.343 * aa->lhet[b->c[j]<<8|b->c[k]] + sum;
		}
		//
		for (k = 0; k != 4; ++k)
			if (p[j<<2|k] < 0.0) p[j<<2|k] = 0.0;
	}

	// convert necessary information to glf1_t
	g->ref_base = ref_base; g->max_mapQ = b->mapQ_max;
	g->depth = ap->n_bases;
	float min_p = 1e30;
	for (j = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			if (p[j<<2|k] < min_p) min_p = p[j<<2|k];
	g->min_lk = min_p > 255.0? 255 : (int)(min_p + 0.5);
	for (j = c = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			g->lk[c++] = p[j<<2|k]-min_p > 255.0? 255 : (int)(p[j<<2|k]-min_p + 0.5);

	free(b);
	return g;
}
void glfgen_core(assemble_aux_t *aa)
{
	nst_bfa1_t *l;
	int k, seqid;
	FILE *fp_bfa = aa->fp_bfa;
	gzFile fpout = aa->fpout, fp_map = aa->fp_map;
	rolling_buf_t *buf = (rolling_buf_t*)calloc(1, sizeof(rolling_buf_t));
	assemble_pos_t *pos = (assemble_pos_t*)calloc(1, sizeof(assemble_pos_t));
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names

	hash_map_char<int> *hash = new hash_map_char<int>;
	for (k = 0; k != mm->n_ref; ++k) hash->insert(mm->ref_name[k], k);
	buf->buf = (maqmap1_t*)calloc(ROLLING_BUF_SIZE, sizeof(maqmap1_t));
	k = mm->n_mapped_reads;
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		if (!hash->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		fprintf(stderr, "[glfgen_core] Processing %s (%d bp)...\n", l->name, l->ori_len);
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
				glf1_t *g = glfgen1_core(pos, aa, ref_base);
				gzwrite(fpout, g, sizeof(glf1_t));
				free(g);
			}
		}
		nst_delete_bfa1(l);
	}
	delete hash;
	maq_delete_maqmap(mm);
	free(buf->buf); free(buf);
	free(pos->bases); free(pos);
}

int maq_glfgen(int argc, char *argv[])
{
	assemble_aux_t *aa = assemble_parse_opt(argc, argv);
	glfgen_core(aa);
	delete_assemble_aux(aa);
	return 0;
}
