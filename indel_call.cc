#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include "maqmap.h"
#include "assemble.h"
#include "bfa.h"
#include "stdhash.hh"
#include "algo.hh"
#include "main.h"

#define MIN_INDEL_DIST 10

typedef struct
{
	int coor, score;
	assemble_indelpos_t indelpos;
	int type;
} indel_info_t;

static inline void print_indelpos(const char *name, indel_info_t *ii, nst_bfa1_t *l)
{
	if (ii->indelpos.n_types == 0) return;
	bit64_t x = ii->indelpos.indels[ii->indelpos.n_types-1];
	char indel_bases[MAX_READLEN];
	if (ii->type != 'x') {
		// generate reference sequence
		char s0[33], s1[33];
		int i, j;
		int l_indel = (signed char)(x&0xff);
		for (j = 0, i = (ii->coor > 32)? ii->coor - 32 : 0; i < ii->coor; ++i)
			s0[j++] = (l->mask[i>>5]>>((31-(i&0x1f))<<1)&3)?
				"ACGT"[l->seq[i>>5]>>((31-(i&0x1f))<<1)&3] : 'N';
		s0[j] = 0;
		for (i = ii->coor, j = 0; i < l->ori_len && i < ii->coor + 32; ++i)
			s1[j++] = (l->mask[i>>5]>>((31-(i&0x1f))<<1)&3)?
				"ACGT"[l->seq[i>>5]>>((31-(i&0x1f))<<1)&3] : 'N';
		s1[j] = 0;
		// generate inserted bases
		if (l_indel > 0) {
			for (i = 0; i != l_indel; ++i) {
				int max = -1, max_j = -1;
				for (j = 0; j != 4; ++j) {
					if (ii->indelpos.ins_bases[i][j] > max) {
						max = ii->indelpos.ins_bases[i][j];
						max_j = j;
					}
				}
				if (max == 0) break; // stop
				indel_bases[i] = "ACGT"[max_j];
			}
			indel_bases[i] = 0;
		} else {
			for (i = 0; i < -l_indel; ++i) {
				if (s1[i]) indel_bases[i] = s1[i];
				else break;
			}
			indel_bases[i] = 0;
		}
		// print
		printf("%s\t%d\t%c\t%d\t%d:%s\t%d\t%d\t%s\t%s\t%d\t@", name, ii->coor+1, ii->type, ii->indelpos.n_reads,
			   int((signed char)(x&0xff)), indel_bases, int(x>>16&0xffff), int(x>>32&0xffff), s0, s1, ii->indelpos.n_ungap);
		for (i = 0; i != ii->indelpos.n_ungap; ++i)
			printf("%d,", ii->indelpos.beg_pos[i]);
		printf("\t@");
		for (i = 0; i != ii->indelpos.n_ungap; ++i)
			printf("%d,", ii->indelpos.end_pos[i]);
		printf("\t@");
		for (i = 0; i != ii->indelpos.n_ungap; ++i)
			printf("%d,", ii->indelpos.n_mm[i]);
		printf("\n");
	}
}
int maq_indelpe(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maq indelpe <in.ref.bfa> <in.aln.map>\n");
		return 1;
	}
	FILE *fp_bfa = fopen(argv[1], "r");
	gzFile fp_map = gzopen(argv[2], "r");
	rolling_buf_t *buf = (rolling_buf_t*)calloc(1, sizeof(rolling_buf_t));
	maqmap_t *mm = maqmap_read_header(fp_map);
	hash_map_char<int> *hash_map = new hash_map_char<int>;
	nst_bfa1_t *l;
	int k, seqid;
	indel_info_t curr, last;

	buf->buf = (maqmap1_t*)calloc(ROLLING_BUF_SIZE, sizeof(maqmap1_t));
	for (k = 0; k != mm->n_ref; ++k) hash_map->insert(mm->ref_name[k], k);
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		if (!hash_map->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		last.type = 'x'; // not set
		for (int i = 0; i != l->len; ++i) {
			for (int j = 31; j >= 0; --j) {
				int coor = (i<<5) | (31 - j);
				if (coor >= l->ori_len) break;
				assemble_get_indelpos(seqid, coor, fp_map, buf, &curr.indelpos);
				if (curr.indelpos.n_types == 0) continue; // no potential indels
				int n = curr.indelpos.n_types;
				bit64_t *p = curr.indelpos.indels;
				algo_sort(n, p);
				curr.score = p[n-1]>>48;
				curr.type = '-';
				if (p[n-1]>>48 >= 2) curr.type = '+';
				if ((p[n-1]>>16&0xffff) && (p[n-1]>>32&0xffff)) curr.type = '*';
				curr.coor = coor;
				if (last.type != 'x') {
					if (curr.coor - last.coor < MIN_INDEL_DIST) {
						if (last.score < curr.score) last.type = '.';
						else curr.type = '.';
					}
				}
				print_indelpos(l->name, &last, l);
				memcpy(&last, &curr, sizeof(indel_info_t));
			}
		}
		print_indelpos(l->name, &last, l);
		nst_delete_bfa1(l);
	}
	delete hash_map;
	gzclose(fp_map);
	fclose(fp_bfa);
	return 0;
}

