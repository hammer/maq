#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include "read.h"
#include "match.hh"
#include "main.h"

#define TRIM_LOWQ 10
#define TRIM_FRACTION 0.5

void delete_longreads(longreads_t *lr)
{
	int i;
	if (lr == 0) return;
	if (lr->seq) {
		for (i = 0; i != lr->n_reads; ++i) free(lr->seq[i]);
		free(lr->seq);
	}
	if (lr->name) {
		for (i = 0; i != (lr->n_reads>>1); ++i) free(lr->name[i]);
		free(lr->name);
	}
	free(lr);
}
longreads_t *ma_load_reads(gzFile fp_l, int size_l, gzFile fp_r, int size_r)
{
	seq_t seq;
	int j, l;
	char name[256];
	longreads_t *lr = (longreads_t*)calloc(1, sizeof(longreads_t));

	fprintf(stderr, "[ma_load_reads] loading reads...\n");
	INIT_SEQ(seq);
	j = 0;
	while ((l = ma_load_1read(fp_l, &seq, name)) >= 0) {
		if (size_l == 0) {
			size_l = l;
			fprintf(stderr, "[ma_load_reads] set length of the first read as %d.\n", size_l);
		}
		if (l < size_l) continue;
		if ((j & 0xfffff) == 0) { // need to enlarge the arrays
			lr->name = (char**)realloc(lr->name, sizeof(char*) * (j + 0x100000));
			lr->seq = (bit8_t**)realloc(lr->seq, sizeof(bit8_t*) * ((j<<1) + 0x200000));
		}
		lr->name[j] = strdup(name);
		lr->seq[j<<1] = (bit8_t*)malloc(size_l);
		memcpy(lr->seq[j<<1], seq.s, size_l);
		if (fp_r) { // read the second read
			assert((l = ma_load_1read(fp_r, &seq, name)) >= 0);
			if (size_r == 0) {
				size_r = l;
				fprintf(stderr, "[ma_load_reads] set length of the second read as %d.\n", size_r);
			}
			lr->seq[(j<<1) | 1] = (bit8_t*)malloc(size_r);
			assert(l >= size_r);
			memcpy(lr->seq[(j<<1) | 1], seq.s, size_r);
			int tl = strlen(name);
			assert(strncmp(name, lr->name[j], tl-1) == 0);
		} else lr->seq[(j<<1) | 1] = 0;
		++j;
	}
	free(seq.s);
	lr->n_reads = j<<1;
	lr->size_l = size_l; lr->size_r = size_r;
	fprintf(stderr, "[ma_load_reads] %d*2 reads loaded.\n", j);
	return lr;
}
static void auto_trim(longreads_t *lr)
{
	int size[2], i, k;
	int *cnt[2];
	size[0] = lr->size_l; size[1] = lr->size_r;
	cnt[0] = (int*)calloc(size[0], sizeof(int));
	cnt[1] = (int*)calloc(size[1], sizeof(int));
	for (i = 0; i != lr->n_reads; ++i) {
		int s = size[i&1];
		bit8_t *seq;
		if (s == 0) continue;
		seq = lr->seq[i];
		for (k = 0; k != s; ++k)
			if ((seq[k]&0x3f) >= TRIM_LOWQ) ++cnt[i&1][k];
	}
	for (i = 0; i != 2; ++i) {
		if (size[i] == 0) continue;
		for (k = size[i] - 1; k >= 0; --k)
			if (2.0 * cnt[i][k] / lr->n_reads >= TRIM_FRACTION) break;
		size[i] = k + 1;
	}
	lr->size_l = size[0]; lr->size_r = size[1];
	fprintf(stderr, "[auto_trim] trimed length: (%d, %d)\n", lr->size_l, lr->size_r);
	free(cnt[0]); free(cnt[1]);
}
match_info_t *ma_longread2read(const longreads_t *lr)
{
	match_info_t *matches;
	read_t mask[2];
	int size[2] = { lr->size_l, lr->size_r };
	int i, n_reads = lr->n_reads;
	fprintf(stderr, "[ma_longread2read] encoding reads... ");
	mask[0] = ~(~read_t(0)<<(size[0]<<1));
	mask[1] = ~(~read_t(0)<<(size[1]<<1));
	matches = (match_info_t*)calloc(n_reads, sizeof(match_info_t));
	assert(matches);
	for (i = 0; i != n_reads; ++i) {
		match_info_t *match;
		read_t s, q, m;
		bit8_t *r = lr->seq[i];
		int cur_size = size[i&1];
		match = matches + i;
		s = q = m = 0;
		for (int k = 0; k != cur_size; ++k) {
			int tmp = (int)r[k];
			s <<= 2; m <<= 2; q <<= 2;
			if (tmp&0x3f) { // a valid residue
				s |= tmp>>6; m |= 3;
			}
			tmp = (int)((tmp&0x3f) / 10.0 + 0.5);
			q |= (tmp < 4)? tmp : 0x3;
		}
		match->s = s; match->q = q; match->m = ~m & mask[i&1];
		match->i1 = match->i2 = MA_NO_MATCH;
		match->seqid1 = match->last_seqid = -1;
	}
	fprintf(stderr, "%d sequences processed.\n", n_reads);
	return matches;
}
void ma_init_match_data(match_data_t *d, gzFile fp_l, int *size_l, gzFile fp_r, int *size_r, int is_trim,
						const char *adapter, gzFile hits_fp)
{
	longreads_t *lr;
	int i;
	lr = ma_load_reads(fp_l, *size_l, fp_r, *size_r);
	if (adapter[0]) ma_trim_adapter(adapter, 0, lr); // prealignment trimming
	if (is_trim) auto_trim(lr);
	d->n_reads = lr->n_reads;
	d->match = ma_longread2read(lr);
	*size_l = lr->size_l; *size_r = lr->size_r;
	if (hits_fp) { // dump read names
		gzprintf(hits_fp, "C\tC comments\n");
		gzprintf(hits_fp, "C\tR readID read_name\n");
		gzprintf(hits_fp, "C\tB chr_name\n");
		gzprintf(hits_fp, "C\tA readID position strand score\n");
		gzprintf(hits_fp, "C\tE chr_name\n");
		for (i = 0; i != (lr->n_reads>>1); ++i)
			gzprintf(hits_fp, "R\t%d\t%s\n", i<<1, lr->name[i]);
	}
	delete_longreads(lr);
	d->pair = 0;
	if (*size_r) {
		d->pair = (pair_info_t*)calloc(d->n_reads>>1, sizeof(pair_info_t));
		for (i = 0; i != (d->n_reads>>1); ++i)
			d->pair[i].i1 = d->pair[i].i2 = MA_NO_MATCH;
	}
}
void maq_methy_modify(match_data_t *d, match_aux_t *o)
{
	int i, k, size[2];
	size[0] = o->size_l; size[1] = o->size_r;
	for (i = 0; i != d->n_reads; ++i) {
		int cur_size = size[i&1];
		match_info_t *m = d->match + i;
		for (k = 0; k < cur_size; ++k) {
			int x = (cur_size-1-k)<<1;
			if (o->methy_mode == 'c') {
				if ((m->s>>x&0x3) == 1) // a C here
					m->s |= read_t(0x3) << x; // C->T
			} else if (o->methy_mode == 'g') {
				if ((m->s>>x&0x3) == 2) // a G here
					m->s ^= read_t(0x2) << x; // G->A
			}
		}			
	}
}
static inline void show_read(int size, read_t s)
{
	for (int i = 0; i < size; ++i) {
		fputc("ACGT"[s>>((size-1-i)<<1)&0x3], stderr);
	}
	fputc('\n', stderr);
}

int ma_bfq2fastq(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maq bfq2fastq <in.bfq> <out.fastq>\n");
		return 1;
	}
	gzFile fpin;
	FILE *fpout;
	seq_t seq;
	char name[256];
	int i, l;
	INIT_SEQ(seq);
	fpin = (strcmp(argv[1], "-") == 0)? gzdopen(STDIN_FILENO, "r") : gzopen(argv[1], "r");
	fpout = (strcmp(argv[2], "-") == 0)? stdout : fopen(argv[2], "w");
	while ((l = ma_load_1read(fpin, &seq, name)) >= 0) {
		fprintf(fpout, "@%s\n", name);
		for (i = 0; i != l; ++i)
			fputc((seq.s[i] == 0)? 'N' : "ACGT"[seq.s[i]>>6&3], fpout);
		fprintf(fpout, "\n+\n");
		for (i = 0; i != l; ++i)
			fputc((seq.s[i]&0x3f) + 33, fpout);
		fputc('\n', fpout);
	}
	gzclose(fpin);
	fclose(fpout);
	return 0;
}

// filter full A reads (for Sanger's "cat" format only, NOT single end or other format)

#define MATCH_A 11
#define NONMATCH_A -59



static void filter_full_A(FILE *fq, int is_pair, int l1)
{
	int l, n_filt = 0, n_tot = 0, l2 = 0;
	seq_t seq, qual;
	char name[256];
	INIT_SEQ(seq); INIT_SEQ(qual);
	while ((l = seq_read_fastq(fq, &seq, &qual, name)) >= 0) {
		int i, score, b = l, e, max, n_A1, n_A, n_N1, n_N, c;
		double r1, r2;
		e = n_A1 = n_A = n_N1 = n_N = 0;
		++n_tot;
		if (l1 == 0) l1 = l>>1;
		if (l2 == 0) l2 = l - l1;
		for (i = score = max = 0; i != l; ++i) {
			qual.s[i] += 33;
			c = toupper(seq.s[i]);
			if (c == 'A' || c == 'N') {
				score += MATCH_A; ++n_A;
				if (c == 'N') ++n_N;
				if (i < l1) {
					++n_A1;
					if (c == 'N') ++n_N1;
				}
			} else score += NONMATCH_A;
			if (score < 0) score = 0;
			if (score > max) {
				max = score; e = i;
			}
		}
		for (i = e, score = max = 0; i >= 0; --i) {
			score += (seq.s[i] == 'A' || seq.s[i] == 'N')? MATCH_A : NONMATCH_A;
			if (score > max) {
				max = score; b = i;
			}
			if (score < 0) break;
		}
		r1 = (double)(e-b+1) / l;
		r2 = (double)n_A / l;
		if (is_pair) {
			if ((n_N1 >= 5 || n_N-n_N1 >= 5) || r1 >= 0.5 || r2 >= 0.9 || (r1 >= 0.45 && r2 >= 0.75)) {
				++n_filt;
				// fprintf(stderr, "%s\t%s\t%s\n", name, seq.s, qual.s);
			} else printf("@%s\n%s\n+\n%s\n", name, seq.s, qual.s);
		} else {
			if (n_N >= 5 || r1 >= 0.9) {
				++n_filt;
			} else printf("@%s\n%s\n+\n%s\n", name, seq.s, qual.s);
		}
	}
	free(seq.s); free(qual.s);
	fprintf(stderr, "[filter_full_A] %d out of %d reads have been sorted out\n", n_filt, n_tot);
}

int maq_catfilter(int argc, char *argv[])
{
	FILE *fp;
	int c, is_pair = 1, l1 = 0;
	while ((c = getopt(argc, argv, "s1:")) >= 0) {
		switch (c) {
		case '1': l1 = atoi(optarg); break;
		case 's': is_pair = 0; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: maq catfilter [-s] [-1 0] <in.fastq>\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? stdin : fopen(argv[optind], "r");
	assert(fp);
	filter_full_A(fp, is_pair, l1);
	fclose(fp);
	return 0;
}
