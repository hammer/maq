#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include "bfa.h"
#include "assemble.h"
#include "match.hh"
#include "main.h"
#include "stdhash.hh"

#define MIN_BASEQ 25
#define MAX_BASEQ 33
#define MAX_MM 3

static hash_set_misc<bit64_t> *load_snp(FILE *fp, hash_map_char<int> *names)
{
	char name[256];
	int pos, c, seqid;
	hash_set_misc<bit64_t> *hash = new hash_set_misc<bit64_t>;
	while (fscanf(fp, "%s%d", name, &pos) == 2) {
		if (names->find(name, &seqid))
			hash->insert((bit64_t)seqid<<32 | (pos-1));
		while ((c = fgetc(fp)) != EOF && c != '\n');
	}
	fprintf(stderr, "[load_snp] %u SNPs are loaded.\n", hash->size());
	return hash;
}

static int usage()
{
	fprintf(stderr, "Usage:   maq mapcheck [options] <chr.bfa> <in.map>\n");
	fprintf(stderr, "Options: -s         use single-end mapping qualities\n");
	fprintf(stderr, "         -Q INT     maximum sum of errors [60]\n");
	fprintf(stderr, "         -m INT     maximum number of mismatches [7]\n");
	fprintf(stderr, "         -q INT     minimum mapping quality [41]\n");
	fprintf(stderr, "         -S INT     quality scale [10]\n");
	fprintf(stderr, "         -P FILE    polymorphic sites [null]\n");
	fprintf(stderr, "         -c         print count instead of fraction\n");
	return 1;
}
int ma_mapcheck(int argc, char *argv[])
{
	int c, is_single = 0, max_err = -1, min_mapQ = 41, max_mm = 7, qscale = 10, is_count = 0;
	char *poly_file = 0;
	while ((c = getopt(argc, argv, "sQ:q:m:S:cP:")) >= 0) {
		switch (c) {
		case 's': is_single = 1; break;
		case 'Q': max_err = atoi(optarg); break;
		case 'm': max_mm = atoi(optarg); break;
		case 'q': min_mapQ = atoi(optarg); break;
		case 'S': qscale = atoi(optarg); break;
		case 'c': is_count = 1; break;
		case 'P': poly_file = strdup(optarg); break;
		}
	}
	if (argc - optind < 2) return usage();

	FILE *fp_bfa, *fpout;
	gzFile fp_map;
	nst_bfa1_t *l;
	int seqid;
	bit64_t n_ref, n_bases, n_cov, tot_len, n_ref_nuc[4], n_rea_nuc[4];
	bit64_t n_01bases, n_01ref;
	n_ref = n_bases = tot_len = n_01bases = n_01ref = n_cov = 0;
	n_ref_nuc[0] = n_ref_nuc[1] = n_ref_nuc[2] = n_ref_nuc[3] = 0;
	n_rea_nuc[0] = n_rea_nuc[1] = n_rea_nuc[2] = n_rea_nuc[3] = 0;
	fpout = stdout;
	fp_bfa = fopen(argv[optind], "r");
	fp_map = gzopen(argv[optind+1], "r");
	int k, max_q;
	bit64_t trans[MAX_READLEN][16], qmm[MAX_READLEN][64], qmm2[MAX_READLEN][64], n_nuc[MAX_READLEN][4];
	rolling_buf_t *buf = (rolling_buf_t*)calloc(1, sizeof(rolling_buf_t));
	buf->buf = (maqmap1_t*)calloc(ROLLING_BUF_SIZE, sizeof(maqmap1_t));
	assemble_pos_t *pos = (assemble_pos_t*)calloc(1, sizeof(assemble_pos_t));
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names

	hash_map_char<int> *hash = new hash_map_char<int>;
	for (k = 0; k != mm->n_ref; ++k) hash->insert(mm->ref_name[k], k);

	hash_set_misc<bit64_t> *snp = 0;
	if (poly_file) {
		FILE *f;
		assert(f = fopen(poly_file, "r"));
		snp = load_snp(f, hash);
		fclose(f);
		free(poly_file);
	}

	k = mm->n_mapped_reads;
	max_q = 0; // maximum base quality
	memset(trans, 0, MAX_READLEN * 16 * sizeof(bit64_t)); // ref->read transition
	memset(qmm, 0, MAX_READLEN * 64 * sizeof(bit64_t));   // base-qual count for mismatches
	memset(qmm2, 0, MAX_READLEN * 64 * sizeof(bit64_t));  // base-qual count for all bases
	memset(n_nuc, 0, MAX_READLEN * 4 * sizeof(bit64_t));
	while ((l = nst_load_bfa1(fp_bfa)) != 0) { // the main loop for counting
		if (!hash->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		fprintf(stderr, "[ma_mapcheck] processing %s...\n", l->name);
		++n_ref; tot_len += l->ori_len;
		for (int i = 0; i != l->len; ++i) {
			bit64_t word = l->seq[i];
			bit64_t mask = l->mask[i];
			for (int j = 31; j >= 0; --j) {
				int coor = (i<<5) | (31 - j);
				int depth01 = 0, c01 = 0;
				if (coor >= l->ori_len) break;
				if ((mask>>(j<<1)&3) == 0) continue; // the reference base is N, a gap
				assemble_get_pos(seqid, coor, fp_map, buf, pos, max_mm, max_err, 0, is_single, 0);
				n_bases += pos->n_bases;
				++n_ref_nuc[word>>(j<<1)&3];
				if (pos->n_bases) ++n_cov;
				if (pos->n_bases == 0) continue; // no need to proceed
				if (snp && snp->find((bit64_t)seqid<<32 | coor)) continue; // known SNP site
				for (k = 0; k != pos->n_bases; ++k) {
					assemble_posinfo_t *p = pos->bases + k;
					int bread = p->info >> 16 & 3; // read base
					int bref = word >> (j<<1) & 3; // ref base
					if (p->info>>22&1) ++n_rea_nuc[bread]; // read base is not N, count
					if ((mask>>(j<<1)&3) && (p->info>>22&1) && (p->info&0xff) >= (bit32_t)min_mapQ) { // not N, high mapQ
						if (p->info>>18&1) { // on the reverse strand
							bref = 3 - bref;
							bread = 3 - bread;
						}
						++n_nuc[p->pos][bread];
						if ((p->info>>24&0x3f)/qscale > (bit8_t)max_q) max_q = (p->info>>24&0x3f)/qscale;
						++trans[p->pos][bref<<2 | bread];
						if (bref != bread) ++qmm[p->pos][(p->info>>24&0x3f)/qscale];
						++qmm2[p->pos][(p->info>>24&0x3f)/qscale];
					}
					if ((p->info>>19&0x3) <= 1) {
						++depth01; c01 += p->c[p->info>>19&0x3];
					}
				}
				if (pos->n_bases > 0 && depth01 == c01) { // all the reads can be mapped "uniquely" at first 24bp
					++n_01ref;
					n_01bases += pos->n_bases;
				}
			}
		}
		nst_delete_bfa1(l);
	}
	gzclose(fp_map);
	fclose(fp_bfa);
	maq_delete_maqmap(mm);
	delete snp;
	delete hash;
	free(buf->buf); free(buf);
	free(pos->bases); free(pos);
	// now print out the results
	bit64_t sum, tot_non;
	printf("Number of reference sequences: %llu\n", n_ref);
	for (k = sum = 0; k != 4; ++k) sum += n_ref_nuc[k];
	tot_non = sum;
	printf("Length of reference sequences exlcuding gaps: %llu\n", tot_non);
	printf("Length of gaps in the reference sequences: %llu\n", tot_len - tot_non);
	printf("Length of non-gap regions covered by reads: %llu\n", n_cov);
	printf("Length of 24bp unique regions of the reference: %llu\n", n_01ref);
	printf("Reference nucleotide composition: A: %.2f%%, C: %.2f%%, G: %.2f%%, T: %.2f%%\n", 100.0*n_ref_nuc[0]/sum,
		   100.0*n_ref_nuc[1]/sum, 100.0*n_ref_nuc[2]/sum, 100.0*n_ref_nuc[3]/sum);
	for (k = sum = 0; k != 4; ++k) sum += n_rea_nuc[k];
	printf("Reads nucleotide composition:     A: %.2f%%, C: %.2f%%, G: %.2f%%, T: %.2f%%\n", 100.0*n_rea_nuc[0]/sum,
		   100.0*n_rea_nuc[1]/sum, 100.0*n_rea_nuc[2]/sum, 100.0*n_rea_nuc[3]/sum);
	printf("Average depth across all non-gap regions: %.3lf\n", (double)n_bases / tot_non);
	printf("Average depth across 24bp unique regions: %.3lf\n\n", n_01bases / (tot_non * ((double)n_01ref/n_cov)));
	printf("      A    C    G    T :  AC  AG  AT  CA  CG  CT  GA  GC  GT  TA  TC  TG :");
	for (k = 0; k <= max_q; ++k) printf("%3d?", k);
	printf(" :");
	for (k = 0; k <= max_q; ++k) printf("%3d?", k);
	putchar('\n');
	for (k = 0; k != MAX_READLEN; ++k) {
		int i, j;
		bit64_t sum, sum2, *q, *p = trans[k];
		// scale to 999
		sum2 = 0;
		for (i = 0; i != 4; ++i) {
			sum = 0;
			q = p + (i<<2);
			for (j = 0; j != 4; ++j) sum += q[j];
			sum2 += sum;
			for (j = 0; j != 4; ++j)
				q[j] = sum? int(999.0 * q[j] / sum + 0.5) : 0;
		}
		if (sum2 == 0) break;
		printf("%2d", k+1);
		// composition
		for (i = sum = 0; i != 4; ++i) sum += n_nuc[k][i];
		for (i = 0; i != 4; ++i)
			printf("%5.1f", sum ? 100.0 * n_nuc[k][i] / sum : -1);
		printf(" :");
		// print trans[][]
		for (i = 0; i != 4; ++i)
			for (j = 0; j != 4; ++j)
				if (i != j) printf("%4llu", p[i<<2 | j]);
		printf(" :");
		p = qmm[k]; q = qmm2[k];
		if (!is_count) { // scale to 999 and print
			for (i = 0; i <= max_q; ++i) p[i] = int(999.0 * (p[i] + 0.1/max_q) / (q[i] + 0.1) + 0.5);
			for (i = 0, sum = 0; i <= max_q; ++i) sum += q[i];
			for (i = 0; i <= max_q; ++i) q[i] = int(999.0 * (q[i] + 0.1/max_q) / (sum + 0.1) + 0.5);
			for (i = 0; i <= max_q; ++i) printf("%4llu", q[i]);
			printf(" :");
			for (i = 0; i <= max_q; ++i) printf("%4llu", p[i]);
		} else { // print count
			for (i = 0; i <= max_q; ++i) printf("\t%llu", q[i]);
			printf(" :");
			for (i = 0; i <= max_q; ++i) printf("\t%llu", p[i]);
		}
		putchar('\n');
	}
	return 0;
}
