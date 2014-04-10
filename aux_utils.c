#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "const.h"
#include "main.h"

#define LINE_LEN 60

void cns_cns2snp(FILE *fpout, gzFile fp, int is_all, int is_alt)
{
	int i, len, n_high = 0;
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			if (is_all || (high>>28 != (high>>24&0xf) && (high>>24&0xf) != 0xf)) {
				// chr, pos, refB, cnsB, cnsQ, depth, avg01, mapQmax, refDiff
				fprintf(fpout, "%s\t%d\t%c\t%c\t%d\t%d\t%.2f\t%d\t%d", name, i + 1, nst_nt16_rev_table[high>>28],
						nst_nt16_rev_table[high>>24&0xf], high>>16&0xff, low&0xff, (float)(low>>20&0xfff)/16.0,
						low>>8&0x3f, (low>>15&0x1f)<<1);
				if (is_alt) fprintf(fpout, "\t%c\t%d\t%c\n", nst_nt16_rev_table[high>>12&0xf], high&0xff,
									nst_nt16_rev_table[high>>8&0xf]);
				else fputc('\n', fpout);
			}
			if ((low>>8&0x3f) >= 40 && (float)(low>>20&0xfff)/16.0 >= 0.25 && (low&0xff) >= 3) ++n_high;
		}
		free(name);
	}
}
void cns_cns2win(FILE *fpout, gzFile fp, int win_size, const char *chr, int begin, int end, int minQ)
{
	int i, len;
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		int n_cvg, n_snp, n_het, depth, match, gc, is_chr, n_unique, uni_depth;
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		n_cvg = n_snp = n_het = depth = match = gc = 0;
		n_unique = uni_depth = 0;
		is_chr = (chr == 0 || strcmp(chr, name) == 0)? 1 : 0;
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			if (high>>28 != 0xf) {
				if ((high>>16&0xff) >= minQ) {
					++n_cvg;
					if ((low>>20&0xfff) >= 8 && (low>>20&0xfff) <= 32) {
						++n_unique;
						uni_depth += low&0xff;
					}
					if (nst_nt16_count_table[high>>24&0xf] == 2) ++n_het;
					if (high>>28 != (high>>24&0xf)) ++n_snp;
					if (high>>28 == 2 || high>>28 == 4) ++gc;
					depth += low&0xff;
					match += low>>20&0xfff;
				}
			}
			if ((i+1)%win_size == 0 || i+1 == len) { // print
				int pos = i+(win_size>>1)+1;
				if (is_chr && (begin >= end || (pos >= begin && pos <= end))) {
					float alt_depth = (n_unique == 0)? 0.0 : (float)uni_depth/n_unique;
					// chr,pos,snpRate,hetRate,oriDepth,altDepth,01Match,GC%
					if (n_cvg > 0) {
						fprintf(fpout, "%s\t%.6f\t%g\t%g\t%.3f\t%.3f\t%.3f\t%.3f\n", name, pos/1000000.0,
								(float)n_snp/n_cvg, (float)n_het/n_cvg, (float)depth/n_cvg, alt_depth,
								(float)match/n_cvg/16.0, 100.0*gc/n_cvg);
					} else {
						fprintf(fpout, "%s\t%.6f\t%g\t%g\t%.3f\t%.3f\t%.3f\t%.3f\n", name, pos/1000000.0,
								0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
					}
					n_cvg = n_snp = n_het = depth = match = gc = 0;
					n_unique = uni_depth = 0;
				}
			}
		}
		free(name);
	}
}

void cns_cns2ref(FILE *fpout, gzFile fp)
{
	int i, len;
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		fprintf(fpout, ">%s\n", name);
		gzread(fp, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			fputc(nst_nt16_rev_table[high>>28], fpout);
			if ((i+1)%LINE_LEN == 0) fputc('\n', fpout);
		}
		if (i%LINE_LEN != 0) fputc('\n', fpout);
		free(name);
	}
}

void cns_cns2cssnp(FILE *fpout, gzFile fp)
{
	int i, len, q[2], b[2], r[2];
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		r[0] = b[0] = 15; q[0] = 0;
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			int is_snp = 0;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			r[1] = high>>28; b[1] = high>>24&0xf; q[1] = high>>16&0xff;
			if (r[0] != 15 && r[1] != 15 && b[0] != 15 && b[1] != 15 && r[0] != b[0] && r[1] != b[1]
				&& (b[0]&r[0]) == (b[1]&r[1]))
			{
				int c = nst_nt16_count_table[b[0]|r[0]|b[1]|r[1]];
				if (r[0] == r[1]) {
					if (b[0] == b[1]) is_snp = 1;
				} else {
					if (c == 4 || c == 2) is_snp = 1;
				}
			}
			if (is_snp) {
				fprintf(fpout, "%s\t%d\t%c%c\t%c%c\t%d\t%d\n", name, i + 1, nst_nt16_rev_table[r[0]],
						nst_nt16_rev_table[r[1]], nst_nt16_rev_table[b[0]], nst_nt16_rev_table[b[1]], q[0], q[1]);
			}
			r[0] = r[1]; b[0] = b[1]; q[0] = q[1];
		}
		free(name);
	}
}

void cns_cns2fq(FILE *fpout, gzFile fp, bit32_t min_mapq, bit32_t min_depth, bit32_t max_depth, bit32_t min_nQ, float hidden)
{
	int i, len;
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		bit8_t *qual;
		gzread(fp, name, len);
		fprintf(fpout, "@%s", name);
		gzread(fp, &len, sizeof(int));
		qual = (bit8_t*)malloc(len);
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			if (i%LINE_LEN == 0) fputc('\n', fpout);
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			if ((low>>8&0x3f) >= min_mapq && (float)(low>>20&0xfff)/16.0 >= hidden && (low&0xff) >= min_depth
				&& (low&0xff) <= max_depth && (low>>15&0x1f)<<1 >= min_nQ)
			{
				fputc(nst_nt16_rev_table[high>>24&0xf], fpout);
			} else fputc(tolower(nst_nt16_rev_table[high>>24&0xf]), fpout);
			qual[i] = ((high>>16&0xff) + 33 > 126)? 126 : (high>>16&0xff) + 33;
		}
		fprintf(fpout, "\n+");
		for (i = 0; i != len; ++i) {
			if (i%LINE_LEN == 0) fputc('\n', fpout);
			fputc(qual[i], fpout);
		}
		fputc('\n', fpout);
		free(name); free(qual);
	}
}

/* main() functions */

int ma_cns2ref(int argc, char *argv[])
{
	gzFile fp;
	if(argc < 2) {
		fprintf(stderr, "Usage: maq cns2ref <in.cns>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	assert(fp);
	cns_cns2ref(stdout, fp);
	gzclose(fp);
	return 0;
}

int ma_cns2snp(int argc, char *argv[])
{
	gzFile fp;
	if(argc < 2) {
		fprintf(stderr, "Usage: maq cns2snp <in.cns>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	assert(fp);
	cns_cns2snp(stdout, fp, 0, 1);
	gzclose(fp);
	return 0;
}
int ma_cns2cssnp(int argc, char *argv[])
{
	gzFile fp;
	if(argc < 2) {
		fprintf(stderr, "Usage: maq cns2cssnp <in.cns>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	assert(fp);
	cns_cns2cssnp(stdout, fp);
	gzclose(fp);
	return 0;
}
int ma_cns2view(int argc, char *argv[])
{
	gzFile fp;
	if (argc < 2) {
		fprintf(stderr, "Usage: maq cns2view <in.cns>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	assert(fp);
	cns_cns2snp(stdout, fp, 1, 1);
	gzclose(fp);
	return 0;
}
int ma_cns2win(int argc, char *argv[])
{
	gzFile fp;
	int win_size = 1000, begin, end, c, minQ;
	char *chr;
	begin = end = minQ = 0; chr = 0;
	while ((c = getopt(argc, argv, "w:b:e:c:q:")) >= 0) {
		switch (c) {
		case 'w': win_size = atoi(optarg); break;
		case 'b': begin = atoi(optarg); break;
		case 'e': end = atoi(optarg); break;
		case 'c': chr = strdup(optarg); break;
		case 'q': minQ = atoi(optarg); break;
		}
	}
	if (optind >= argc) {
		fprintf(stderr, "Usage: maq cns2win [-w 1000] [-b 0] [-e 0] [-c null] [-q 0] <in.cns>\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	assert(fp);
	cns_cns2win(stdout, fp, win_size, chr, begin, end, minQ);
	gzclose(fp);
	free(chr);
	return 0;
}
int ma_cns2fq(int argc, char *argv[])
{
	gzFile fp;
	double hidden = 0.25;
	int c, min_mapq = 40, min_depth = 3, max_depth = 255, min_nQ = 20;
	while ((c = getopt(argc, argv, "Q:d:D:n:c:")) >= 0) {
		switch (c) {
		case 'Q': min_mapq = atoi(optarg); break;
		case 'd': min_depth = atoi(optarg); break;
		case 'D': max_depth = atoi(optarg); break;
		case 'n': min_nQ = atoi(optarg); break;
		case 'c': hidden = atof(optarg); break;
		default: return 1;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   maq cns2fq [options] <in.cns>\n\n");
		fprintf(stderr, "Options: -Q INT    minimum mapping quality [%d]\n", min_mapq);
		fprintf(stderr, "         -n INT    minimum neighbouring quality [%d]\n", min_nQ);
		fprintf(stderr, "         -d INT    minimum read depth [%d]\n", min_depth);
		fprintf(stderr, "         -D INT    maximum read depth [%d]\n\n", max_depth);
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	assert(fp);
	cns_cns2fq(stdout, fp, min_mapq, min_depth, max_depth, min_nQ, hidden);
	gzclose(fp);
	return 0;
}
