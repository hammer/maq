#ifndef LH3_MATCH_HH
#define LH3_MATCH_HH

#include <zlib.h>
#include <stdio.h>
#include "maqmap.h"
#include "dword.hh"

#ifdef _FASTMAP
#  define SEED_LENGTH 28
#  define MAX_MISMATCH 3
#else
#  define SEED_LENGTH 24
#  define MAX_MISMATCH 4
#endif

#ifdef MAQ_LONGREADS
#  define read_t dword_t<dword_t<bit64_t> >
#  define matches_t dword_t<bit64_t>
#else
#  define read_t dword_t<bit64_t>
#  define matches_t bit64_t
#endif

#define MA_NO_MATCH 0xfffffffful
#define MA_NO_INDEX 0xffffffffu

typedef dword_t<bit64_t> bit128_t;

typedef struct __match_info_t
{
	read_t s, q, m;
	bit32_t i1, p1, i2;
	matches_t mm1, mm2, last_mm1, last_mm2;
 	bit64_t last1, last2;
	int seqid1, last_seqid;
	bit8_t c[MAX_MISMATCH + 1];
} match_info_t;

typedef struct
{
	bit32_t i1, p11, p12;
	bit32_t i2;
	matches_t m11, m12;
	int seqid1;
} pair_info_t;

typedef struct __match_aux_t
{
	int n_mismatch;
	int size_l, size_r;
	int n_filters;
	int max_err10;
	int max_dist, min_dist;
	int RF_max_dist; // for Illumina long insert-size library
	int is_quiet, is_color, is_mm, is_p2diff, is_sw;
	char methy_mode;
	char *adapter, *dump_file;
	int max_hits;
	gzFile hits_fp;
	read_t *masks;
	bit64_t *filters;
	int *shift_seed;
	int log_n[256];
	int q_rate; // mutation rate in phred unit
} match_aux_t;

typedef struct
{
	read_t mask;
	int shift_seed;
	bit64_t filter;
	bit128_t *sorted;
	bit32_t *index;
} match_index_t;

struct __longreads_t;

typedef struct __match_data_t
{
	int n_reads, n_lname;
	bit32_t sum_len;
	match_info_t *match;
	pair_info_t *pair;
	char **lname;
	match_index_t index[4];
	bit8_t *index4;
} match_data_t;

match_aux_t *new_match_aux(int size_l, int size_r, int m);
match_data_t *new_match_data();
void delete_match_data(match_data_t *d);
void delete_match_aux(match_aux_t *o);
void match_data2map(gzFile fpout, FILE *fp_bfa, gzFile fp_bfq_l, gzFile fp_bfq_r, const match_aux_t *o, match_data_t *d);
void ma_trim_adapter(const char adaptor[], const match_data_t *d, struct __longreads_t *lr);

#endif
