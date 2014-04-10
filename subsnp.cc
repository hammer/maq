#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "const.h"
#include "main.h"
#include "stdhash.hh"

#define MAX_LINE_LEN 65536

hash_map_char<char*> *ma_load_snp(FILE *fp)
{
	char name[256], key[256];
	char *buffer = (char*)calloc(MAX_LINE_LEN, 1);
	int pos, c, i;
	hash_map_char<char*> *hash = new hash_map_char<char*>;
	while (fscanf(fp, "%s%d", name, &pos) == 2) {
		sprintf(key, "%s.%d", name, pos);
		i = 0;
		while ((c = fgetc(fp)) != EOF && c != '\n') buffer[i++] = c;
		buffer[i] = 0;
		hash->insert(key, strdup(buffer));
	}
	free(buffer);
	return hash;
}
hash_set_char *ma_load_snp_set(FILE *fp)
{
	char name[256], key[256];
	int pos, c;
	hash_set_char *hash = new hash_set_char;
	while (fscanf(fp, "%s%d", name, &pos) == 2) {
		sprintf(key, "%s.%d", name, pos);
		while ((c = fgetc(fp)) != EOF && c != '\n');
		hash->insert(key);
	}
	return hash;
}
void ma_free_char_hash(hash_map_char<char*> *hash)
{
	if (hash == 0) return;
	hash_map_char<char*>::iterator iter;
	for (iter = hash->begin(); iter != hash->end(); ++iter) {
		if (iter.isfilled()) free(iter.value());
	}
	delete hash;
}
void cns_snpreg(FILE *fpout, gzFile fp, FILE *fp_snp, bit32_t min_mapq, bit32_t min_depth, bit32_t max_depth, bit32_t min_nQ)
{
	int i, len;
	bit64_t n, nonn, nn, n_thres[10];
	char key[256], *str;
	hash_map_char<char*> *hash = 0;
	if (fp_snp) {
		fprintf(stderr, "-- loading .snp file...\n");
		hash = ma_load_snp(fp_snp);
		fprintf(stderr, "-- %u records loaded.\n", hash->size());
	}
	n = nonn = nn = 0;
	memset(n_thres, 0, sizeof(bit64_t) * 10);
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			if (high>>28 == 0xf) ++nn;
			else ++nonn;
			if ((low>>8&0x3f) >= min_mapq && (float)(low>>20&0xfff)/16.0 >= 0.25 && (low&0xff) >= min_depth
				&& (low&0xff) <= max_depth && (low>>15&0x1f)<<1 >= min_nQ)
			{
				if (hash) {
					sprintf(key, "%s.%d", name, i+1);
					if (hash->find(key, &str))
						fprintf(fpout, "%s\t%d%s\n", name, i+1, str);
				}
				++n;
				int q = high>>16&0xff;
				if (q > 99) q = 99;
				++n_thres[q/10];
			}
		}
		free(name);
	}
	fprintf(stderr, "-- total length of reference: %lld\n", nonn + nn);
	fprintf(stderr, "-- length of 'N' regions on the reference: %lld\n", nn);
	fprintf(stderr, "-- size of good regions: %lld\n", n);
	fprintf(stderr, "-- uncalled: %.2f%% (= %lld / %lld)\n", 100.0 * (1.0 - (float)n/nonn), nonn - n, nonn);
	n = 0;
	for (i = 9; i >= 0; --i) {
		n += n_thres[i];
		fprintf(stderr, "-- %d0-%d9: %lld\n", i, i, n);
	}
	ma_free_char_hash(hash);
}
void cns_subpos(FILE *fpout, gzFile fp, FILE *fp_snp)
{
	int i, len;
	char key[256], *str;
	fprintf(stderr, "[cns_cns2pos] loading .snp file...\n");
	hash_map_char<char*> *hash = ma_load_snp(fp_snp);
	fprintf(stderr, "[cns_cns2pos] %u records loaded.\n", hash->size());
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			sprintf(key, "%s.%d", name, i+1);
			if (hash->find(key, &str)) {
				fprintf(fpout, "%s\t%d\t%c\t%c\t%d\t%d\t%.2f\t%d\t%d%s\n", name, i + 1, nst_nt16_rev_table[high>>28],
						nst_nt16_rev_table[high>>24&0xf], high>>16&0xff, low&0xff, (float)(low>>20&0xfff)/16.0,
						low>>8&0x3f, (low>>15&0x1f)<<1, str);
			}
		}
		free(name);
	}
	ma_free_char_hash(hash);
}
int ma_subpos(int argc, char *argv[])
{
	gzFile fp;
	FILE *fp_snp;
	if(argc < 3) {
		fprintf(stderr, "Usage: maq subpos <in.cns> <.snp>\n");
		return 1;
	}
	fp = gzopen(argv[1], "r");
	fp_snp = fopen(argv[2], "r");
	assert(fp && fp_snp);
	cns_subpos(stdout, fp, fp_snp);
	gzclose(fp);
	fclose(fp_snp);
	return 0;
}

int ma_snpreg(int argc, char *argv[])
{
	gzFile fp;
	FILE *fp_snp = 0;
	int c, min_mapq = 40, min_depth = 3, max_depth = 255, min_nQ = 20;
	while ((c = getopt(argc, argv, "Q:d:D:n:")) >= 0) {
		switch (c) {
		case 'Q': min_mapq = atoi(optarg); break;
		case 'd': min_depth = atoi(optarg); break;
		case 'D': max_depth = atoi(optarg); break;
		case 'n': min_nQ = atoi(optarg); break;
		default: return 1;
		}
	}
	if(optind >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   maq snpreg [options] in.cns [in.snp]\n\n");
		fprintf(stderr, "Options: -Q INT    minimum mapping quality [%d]\n", min_mapq);
		fprintf(stderr, "         -d INT    minimum read depth [%d]\n", min_depth);
		fprintf(stderr, "         -n INT    minimum neighbouring quality [%d]\n", min_nQ);
		fprintf(stderr, "         -D INT    maximum read depth [%d]\n\n", max_depth);
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	if (optind + 1 < argc)
		assert(fp_snp = fopen(argv[optind+1], "r"));
	assert(fp);
	cns_snpreg(stdout, fp, fp_snp, min_mapq, min_depth, max_depth, min_nQ);
	gzclose(fp);
	if (fp_snp) fclose(fp_snp);
	return 0;
}

/* CMD: simucns */

static hash_map_char<bit8_t> *load_simu_snp(FILE *fp)
{
	hash_map_char<bit8_t> *hash = new hash_map_char<bit8_t>;
	char name[256], key[256], a[3], b[3], c[3];
	int pos, i;
	while (fscanf(fp, "%s%d%s%s%s", name, &pos, a, b, c) == 5) {
		if (a[0] == '-' || b[0] == '-') {
			for (i = pos - 5; i <= pos + 5; ++i) {
				sprintf(key, "%s.%d", name, i);
				hash->insert(key, 15);
			}
		} else {
			sprintf(key, "%s.%d", name, pos);
			hash->insert(key, nst_nt16_table[(int)b[0]]);
		}
	}
	return hash;
}
void maq_simucns_core(FILE *fpout, gzFile fp, FILE *fp_snp)
{
	int i, len;
	char key[256];
	bit64_t n_err[10], n_all[10];
	fprintf(stderr, "-- loading .snp file...\n");
	hash_map_char<bit8_t> *hash = load_simu_snp(fp_snp);
	fprintf(stderr, "-- %u records loaded.\n", hash->size());
	memset(n_err, 0, sizeof(bit64_t) * 10);
	memset(n_all, 0, sizeof(bit64_t) * 10);
	while (gzread(fp, &len, sizeof(int))) {
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzread(fp, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			gzread(fp, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			sprintf(key, "%s.%d", name, i+1);
			bit8_t b2, b = high>>28;
			int q10 = (high>>16&0xff) / 10;
			if (q10 > 9) q10 = 9;
			if (hash->find(key, &b2)) b = b2;
			if (b < 15) {
				++n_all[q10];
				if (b != (high>>24&0xf)) ++n_err[q10];
			}
		}
		free(name);
	}
	printf("%4s%13s%13s%15s\n", "cnsQ", "#called", "#wrong", "err_rate");
	for (i = 0; i != 10; ++i)
		printf("%1dx  %13llu%13llu%15e\n", i, n_all[i], n_err[i], (float)n_err[i] / n_all[i]);
	delete hash;
}

int maq_simucns(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: maq simucns <in.cns> <in.true.snp>\n");
		return 1;
	}
	gzFile fp_cns = gzopen(argv[1], "r");
	FILE *fp_snp = fopen(argv[2], "r");
	assert(fp_cns && fp_snp);
	maq_simucns_core(stdout, fp_cns, fp_snp);
	return 0;
}

