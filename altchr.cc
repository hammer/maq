#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include "main.h"
#include "const.h"
#include "bfa.h"
#include "stdhash.hh"

typedef struct
{
	hash_map_char<bit64_t> *hash;
	bit64_t **snps;
} snp_array_t;

static void maq_delete_snp_array(snp_array_t *s)
{
	for (int i = 0; i != int(s->hash->size()); ++i)
		free(s->snps[i]);
	free(s->snps);
	delete s->hash;
	free(s);
}

static snp_array_t *maq_load_snp_array(FILE *fp)
{
	char name[256], str1[16], str2[16];
	snp_array_t *sa = (snp_array_t*)calloc(1, sizeof(snp_array_t));
	sa->hash = new hash_map_char<bit64_t>;
	int pos, c;
	char refB, cnsB;
	while (fscanf(fp, "%s%d%s%s", name, &pos, str1, str2) == 4) {
		refB = str1[0]; cnsB = str2[0];
		while ((c = fgetc(fp)) != '\n' && c != EOF);
		bit8_t ref = nst_nt16_table[int(refB)];
		bit8_t cns = nst_nt16_table[int(cnsB)];
//		if (nst_nt16_count_table[ref] != 1 || cns == 0xf || ref == cns) continue;
		assert(nst_nt16_count_table[ref] == 1 && cns != 0xf && ref != cns);
		bit8_t snp;
		if (nst_nt16_count_table[cns] == 2) {
			snp = cns ^ ref;
			if (nst_nt16_count_table[snp] != 1) continue;
		} else if (nst_nt16_count_table[cns] == 1) {
			snp = cns;
		} else continue;
		snp = nst_nt16_nt4_table[ref]<<4 | nst_nt16_nt4_table[snp];
		bit64_t tmp;
		int id, n;
		if (sa->hash->find(name, &tmp)) {
			id = tmp>>32;
			n = tmp&0xffffffffu;
		} else {
			if ((sa->hash->size()&0xff) == 0)
				sa->snps = (bit64_t**)realloc(sa->snps, sizeof(bit64_t*) * (sa->hash->size() + 0x100));
			id = sa->hash->size();
			n = 0;
		}
		if ((n&0xffff) == 0)
			sa->snps[id] = (bit64_t*)realloc(sa->snps[id], sizeof(bit64_t) * (n + 0x10000));
		tmp = bit64_t(id)<<32 | (n+1);
		sa->hash->insert(name, tmp);
		sa->snps[id][n++] = bit64_t(snp)<<32 | (pos-1);
	}
	return sa;
}

void maq_altchr_core(nst_bfa1_t *b, snp_array_t *s)
{
	bit64_t tmp;
	if (!s->hash->find(b->name, &tmp)) return;
	int n = tmp&0xffffffffu;
	bit64_t *snps = s->snps[tmp>>32];
	for (int k = 0; k != n; ++k) {
		bit32_t pos = snps[k]&0xffffffffu;
		bit8_t refB = snps[k]>>36 & 3;
		bit8_t cnsB = snps[k]>>32 & 3;
		int i = pos>>5;
		int j = 31 - (pos&0x1f);
//		fprintf(stderr, "%d\t%d\t%d\n", refB, cnsB, i);
		assert(refB == (b->seq[i]>>(j<<1)&3));
		b->seq[i] = (b->seq[i] & ~(bit64_t(3) << (j<<1))) | bit64_t(cnsB) << (j<<1);
	}
}

int ma_altchr(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "Usage: maq altchr <altchr.bfa> <oldchr.bfa> <cns.snp>\n");
		return 1;
	}
	FILE *fp_bfa, *fp_snp, *fp_alt;
	fp_alt = fopen(argv[1], "w");
	fp_bfa = fopen(argv[2], "r");
	fp_snp = fopen(argv[3], "r");
	assert(fp_alt && fp_bfa && fp_snp);
	snp_array_t *s = maq_load_snp_array(fp_snp);
	nst_bfa1_t *l;
	fclose(fp_snp);
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		maq_altchr_core(l, s);
		int len = strlen(l->name) + 1;
		fwrite(&len, sizeof(int), 1, fp_alt);
		fwrite(l->name, sizeof(char), len, fp_alt);
		fwrite(&l->ori_len, sizeof(int), 1, fp_alt);
		fwrite(&l->len, sizeof(int), 1, fp_alt);
		fwrite(l->seq, sizeof(bit64_t) * l->len, 1, fp_alt);
		fwrite(l->mask, sizeof(bit64_t) * l->len, 1, fp_alt);
	}
	fclose(fp_alt); fclose(fp_bfa);
	maq_delete_snp_array(s);
	return 0;
}
