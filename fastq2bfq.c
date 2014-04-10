#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include "const.h"
#include "main.h"
#include "seq.h"

int64_t fastq2bfq(FILE *fp_fq, const char *fn_bfq, int n_reads)
{
	seq_t seq, qual;
	char name[256], str[1024];
	int l, is_new = 0, l_prefix = 0;
	bit64_t n;
	gzFile *fp = 0;
	INIT_SEQ(seq); INIT_SEQ(qual);
	seq_set_block_size(256);
	n = 0;
	if (n_reads <= 0) strcpy(str, fn_bfq);
	else {
		strcpy(str, fn_bfq);
		if (strcmp(str + strlen(fn_bfq) - 4, ".bfq") == 0) // remove ".bfq" suffix if exist
			str[strlen(fn_bfq) - 4] = '\0';
		l_prefix = strlen(str);
		sprintf(str + l_prefix, "@1.bfq");
	}
	fp = gzopen(str, "w");
	while ((l = seq_read_fastq(fp_fq, &seq, &qual, name)) >= 0) {
		int i, nN;
		bit8_t t, tt;
		if (is_new) {
			sprintf(str + l_prefix, "@%lld.bfq", n+1);
			fp = gzopen(str, "w");
			is_new = 0;
		}
		for (i = 0, nN = 0; i != l; ++i) {
			t = nst_nt4_table[seq.s[i]];
			if (t > 3) ++nN;
			tt = (qual.s[i] > 0x3f)? 0x3f : qual.s[i];
			if (tt == 0) tt = 1;
			seq.s[i] = (t > 3)? 0 : ((t<<6) | tt);
		}
		i = strlen(name) + 1;
		gzwrite(fp, &i, sizeof(int));
		gzwrite(fp, name, sizeof(char) * i);
		gzwrite(fp, &l, sizeof(int));
		gzwrite(fp, seq.s, sizeof(char) * l);
		++n;
		if (n_reads > 0 && n % n_reads == 0) {
			gzclose(fp);
			fprintf(stderr, "-- finish writing file '%s'\n", str);
			fp = 0;
			is_new = 1;
		}
	}
	free(seq.s); free(qual.s);
	if (fp) {
		fprintf(stderr, "-- finish writing file '%s'\n", str);
		gzclose(fp);
	}
	fprintf(stderr, "-- %lld sequences were loaded.\n", n);
	return n;
}
void sol2sanger(FILE *fpin, FILE *fpout)
{
	seq_t seq, qual;
	char name[256];
	int table[128];
	int l;
	/* calculate table */
	for (l = 0; l != 128; ++l) {
		table[l] = (int)(33 + 10.0 * log(1.0 + pow(10.0, (l - 64 + 33) / 10.0)) / log(10.0) + .499);
		if (table[l] >= 126) table[l] = 126;
	}
	INIT_SEQ(seq); INIT_SEQ(qual);
	seq_set_block_size(256);
	while ((l = seq_read_fastq(fpin, &seq, &qual, name)) >= 0) {
		int i;
		fprintf(fpout, "@%s\n%s\n+\n", name, seq.s);
		for (i = 0; i != l; ++i)
			qual.s[i] = table[(int)qual.s[i]];
		fprintf(fpout, "%s\n", qual.s);
	}
	free(seq.s); free(qual.s);
}
int ma_sol2sanger(int argc, char *argv[])
{
	FILE *fpin, *fpout;
	fpin = fpout = 0;
	if (argc < 3) {
		fprintf(stderr, "Usage: maq sol2sanger <in.fastq> <out.fastq>\n");
		return 1;
	}
	fpin = (strcmp(argv[1], "-") == 0)? stdin : fopen(argv[1], "r");
	fpout = (strcmp(argv[2], "-") == 0)? stdout : fopen(argv[2], "w");
	sol2sanger(fpin, fpout);
	fclose(fpin); fclose(fpout);
	return 0;
}
int ma_fastq2bfq(int argc, char *argv[])
{
	FILE *fp_fq;
	int c, n_reads = 0;
	while ((c = getopt(argc, argv, "n:")) >= 0) {
		switch (c) {
		case 'n': n_reads = atoi(optarg); break;
		}
	}
	fp_fq = 0;
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: maq fastq2bfq [-n nreads] <in.fastq> <out.prefix>|<out.bfq>\n");
		return 1;
	}
	fp_fq = (strcmp(argv[optind], "-") == 0)? stdin : fopen(argv[optind], "r");
	fastq2bfq(fp_fq, argv[optind+1], n_reads);
	fclose(fp_fq);
	return 0;
}
