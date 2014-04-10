#ifndef LH3_ASSEMBLE_H
#define LH3_ASSEMBLE_H

#include <zlib.h>
#include <stdio.h>
#include "maqmap.h"
#include "const.h"

#define ROLLING_BUF_SIZE 0x100000

typedef struct
{
	// map_qual should be no bigger than 99
	// empty:2, qual:6; empty:1, is_present:1, mm:3, strand:1, base:2; empty:2, base_qual:6; map_qual:8
	bit32_t info;
	bit8_t c[2];
	bit8_t pos;
} assemble_posinfo_t;

typedef struct
{
	int n_bases, m_bases;
	assemble_posinfo_t *bases;
} assemble_pos_t;

typedef struct
{
	int n_reads, n_types;
	bit64_t indels[256];
	int n_ungap;
	bit8_t beg_pos[256], end_pos[256], n_mm[256];
	int ins_bases[MAX_READLEN][4];
} assemble_indelpos_t;

typedef struct
{
	int head, tail, is_rounded;
	maqmap1_t *buf;
} rolling_buf_t;

typedef struct
{
	float hetero_rate, theta, eta, maf, q_r;
	int max_mm, max_err, min_q, is_alt, is_single, is_pair_only;
	double *fk, *coef;
	FILE *fp_bfa;
	gzFile fpout, fp_map;

	int n_hap;
	double *lhet;
} assemble_aux_t;

// ref_base:4, base:4, qual:8; base2:4, base3:4, qual2:8; avg01:12, qNei:5, het:1, qMax:6, depth:8
//         60      56      48       44       40       32        20      15     14       8        0

#ifdef __cplusplus
extern "C" {
#endif
	void assemble_get_pos(bit32_t seqid, bit32_t pos, gzFile fpin, rolling_buf_t *ab, assemble_pos_t *ap,
						  int max_mm, int max_err, int min_q, int is_single, int is_pair_only);
	assemble_aux_t *assemble_parse_opt(int argc, char *argv[]);
	void assemble_get_indelpos(bit32_t seqid, bit32_t pos, gzFile fpin, rolling_buf_t *ab, assemble_indelpos_t *ai);
	void delete_assemble_aux(assemble_aux_t *aa);
#ifdef __cplusplus
}
#endif

#endif
