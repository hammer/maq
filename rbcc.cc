#include <unistd.h>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "main.h"
#include "const.h"
#include "stdhash.hh"

hash_map_char<int> *cns_load_dbSNP(FILE *fp_dbSNP)
{
	hash_map_char<int> *hash;
	char name[256], str1[256], str2[256];
	float het_rate;
	int qual, pos;
	hash = new hash_map_char<int>;
	hash->resize(100000);
	while (fscanf(fp_dbSNP, "%s%d%s%s%d%f", name, &pos, str1, str2, &qual, &het_rate) == 6) {
		sprintf(str1, "%s.%d", name, pos);
		hash->insert(str1, qual);
	}
	return hash;
}
static inline bit32_t rbcc_core(bit32_t high, bit32_t hrate)
{
	bit32_t b[3], q[2], tmp;
	b[0] = high>>24&0xf; b[1] = high>>12&0xf; b[2] = high>>8&0xf;
	q[0] = high>>16&0xff; q[1] = high&0xff;
	if (nst_nt16_count_table[b[0]] == 2) { // the best call is a het
		q[0] += hrate;
	} else if (nst_nt16_count_table[b[1]] == 2) { // the second best call is a het
		if (q[0] - hrate < 0) { // then the het should be the best call
			tmp = b[0]; b[0] = b[1]; b[1] = tmp; // swap b[0] and b[1]
			q[1] += q[0]; q[0] = hrate - q[0];
		} else { // just adjust qualities
			q[0] -= hrate; q[1] += hrate;
		}
	} else if (nst_nt16_count_table[b[2]] == 2) { // the third is a het
		if (q[0] + q[1] < hrate) { // then the het should be the best call
			tmp = b[0]; b[0] = b[2]; b[2] = b[1]; b[1] = tmp; // swap
			tmp = q[0]; q[0] = hrate - q[0] - q[1]; q[1] = tmp;
		} else if (q[1] < hrate) { // then the het should be the second best call
			tmp = b[1]; b[1] = b[2]; b[2] = tmp;
			tmp = q[0] + q[1] - hrate; q[1] = hrate - q[1]; q[0] = tmp;
		} else q[1] -= hrate; // just adjust qualities
	}
	if (q[0] > 255) q[0] = 255;
	if (q[1] > 255) q[1] = 255;
	return (high&0xf0000000u) | b[0]<<24 | b[1]<<12 | b[2]<<8 | q[0]<<16 | q[1];
}
void cns_rbcc(gzFile fpout, gzFile fp, FILE *fp_dbSNP, float r_ori, float r_good, float r_bad)
{
	int i, len;
	char key[256];
	bit32_t q_good, q_bad;
	assert(r_ori < 1.0 && r_good < 1.0 && r_bad < 1.0);
	q_good = int(4.343 * log(r_good / (1.0-r_good) / (r_ori / (1.0-r_ori))) + 0.5);
	q_bad  = int(4.343 * log(r_bad  / (1.0-r_bad ) / (r_ori / (1.0-r_ori))) + 0.5);
	hash_map_char<int> *hash = cns_load_dbSNP(fp_dbSNP);
	fprintf(stderr, "[cns_rbcc] %d SNPs are loaded.\n", hash->size());
	while (gzread(fp, &len, sizeof(int))) {
		gzwrite(fpout, &len, sizeof(int));
		char *name = (char*)malloc(len);
		gzread(fp, name, len);
		gzwrite(fpout, name, len);
		gzread(fp, &len, sizeof(int));
		gzwrite(fpout, &len, sizeof(int));
		for (i = 0; i != len; ++i) {
			bit32_t low, high;
			int qual;
			gzread(fp, &low, sizeof(bit32_t));
			gzwrite(fpout, &low, sizeof(bit32_t));
			gzread(fp, &high, sizeof(bit32_t));
			sprintf(key, "%s.%d", name, i+1);
			if (hash->find(key, &qual))
				high = rbcc_core(high, qual? q_good : q_bad);
			gzwrite(fpout, &high, sizeof(bit32_t));
		}
		free(name);
	}
	delete hash;
}

int ma_rbcc(int argc, char *argv[])
{
	gzFile fp, fpout;
	FILE *fp_dbSNP;
	int c;
	float r_ori, r_good, r_bad;
	r_ori = 0.001; r_good = 0.1; r_bad = 0.01;
	while ((c = getopt(argc, argv, "r:g:b:")) >= 0) {
		switch (c) {
		case 'r': r_ori = atof(optarg); break;
		case 'g': r_good = atof(optarg); break;
		case 'b': r_bad = atof(optarg); break;
		}
	}
	if(argc - optind < 3) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   maq rbcc [options] <out.cns> <in.cns> <chr_rpt.snp>\n\n");
		fprintf(stderr, "Options: -r FLOAT    original prior [0.001]\n");
		fprintf(stderr, "         -g FLOAT    prior for validated SNPs [0.1]\n");
		fprintf(stderr, "         -b FLOAT    prior for the rest of SNPs [0.01]\n\n");
		return 1;
	}
	fpout = gzopen(argv[optind], "w");
	fp = gzopen(argv[optind+1], "r");
	fp_dbSNP = fopen(argv[optind+2], "r");
	assert(fp && fpout && fp_dbSNP);
	cns_rbcc(fpout, fp, fp_dbSNP, r_ori, r_good, r_bad);
	gzclose(fp);
	gzclose(fpout);
	fclose(fp_dbSNP);
	return 0;
}
