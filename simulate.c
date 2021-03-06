#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>
#include "main.h"
#include "seq.h"
#include "genran.h"
#include "const.h"
#include "maqmap.h"

#ifndef MAX_READLEN
#define MAX_READLEN 128
#endif
#define MAX_QUAL 63

typedef struct
{
	int l_min, l_max, q_max;
	bit64_t n;
	bit64_t comp_count[MAX_READLEN][5];
	bit64_t qual_count[MAX_READLEN][MAX_QUAL+1];
	double fqual_count[MAX_READLEN][MAX_QUAL+1];
	bit64_t markov_count[MAX_READLEN][MAX_QUAL+1][MAX_QUAL+1];
	double fmarkov_count[MAX_READLEN][MAX_QUAL+1][MAX_QUAL+1];
} fqc_t;

fqc_t *fqc_collect(FILE *fp_fq)
{
	seq_t seq, qual;
	char name[256];
	int i, l, q1, q2;
	fqc_t *fqc = (fqc_t*)calloc(1, sizeof(fqc_t));
	fqc->l_min = 1<<30;
	INIT_SEQ(seq); INIT_SEQ(qual);
	seq_set_block_size(256);
	while ((l = seq_read_fastq(fp_fq, &seq, &qual, name)) >= 0) {
		bit8_t t, tt, prev = 0;
		if (l < fqc->l_min) fqc->l_min = l;
		if (l > fqc->l_max) fqc->l_max = l;
		for (i = 0; i != l; ++i) {
			t = nst_nt4_table[seq.s[i]];
			if (t > 4) t = 4;
			tt = (qual.s[i] > MAX_QUAL)? MAX_QUAL : qual.s[i];
			if (tt > fqc->q_max) fqc->q_max = tt;
			++fqc->comp_count[i][t];
			++fqc->qual_count[i][tt];
			if (i > 0) ++fqc->markov_count[i][prev][tt];
			prev = tt;
		}
		++fqc->n;
	}
	free(seq.s); free(qual.s);
	assert(fqc->l_min == fqc->l_max); // this should not happen
	for (i = 0; i < fqc->l_max; ++i) { // normalize to 1
		bit64_t sum = 0;
		for (q1 = 0; q1 <= fqc->q_max; ++q1)
			sum += fqc->qual_count[i][q1];
		for (q1 = 0; q1 <= fqc->q_max; ++q1)
		    fqc->fqual_count[i][q1] = (double)fqc->qual_count[i][q1] / sum;
		for (q1 = 0; q1 <= fqc->q_max; ++q1) {
			sum = 0;
			for (q2 = 0; q2 <= fqc->q_max; ++q2)
				sum += fqc->markov_count[i][q1][q2];
			for (q2 = 0; q2 <= fqc->q_max; ++q2)
				fqc->fmarkov_count[i][q1][q2] = (double)fqc->markov_count[i][q1][q2] / sum;
		}
	}
	return fqc;
}

int maq_simutrain(int argc, char *argv[])
{
	fqc_t *fqc;
	FILE *fp;
	gzFile *fpout;
	if (argc < 3) {
		fprintf(stderr, "Usage: maq simutrain <simupars.dat> <known_reads.fastq>\n");
		return 1;
	}
	fp = (strcmp(argv[2], "-") == 0)? stdin : fopen(argv[2], "r");
	fpout = gzopen(argv[1], "w");
	assert(fp && fpout);
	fqc = fqc_collect(fp);
	gzwrite(fpout, fqc, sizeof(fqc_t));
	gzclose(fpout);
	fclose(fp);
	free(fqc);
	return 0;
}

static void gen_qual(fqc_t *fqc, char str[], int size)
{
	double ran = ran_uniform();
	double *p = fqc->fqual_count[0];
	double sum;
	int i, k, prev;
	for (k = 0, sum = 0.0; k <= fqc->q_max; ++k) {
		sum += p[k]; // precalculate this sum will be more effficient, but this is not critical...
		if (sum >= ran) break;
	}
	str[0] = 33 + k;
	prev = k;
	for (i = 1; i != size; ++i) {
		double ran = ran_uniform();
		p = fqc->fmarkov_count[i][prev];
		for (k = 0, sum = 0.0; k <= fqc->q_max; ++k) {
			sum += p[k];
			if (sum >= ran) break;
		}
		str[i] = 33 + k;
		prev = k;
	}
	str[size] = 0; // end of a string
}

static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.1;

void maq_mut_diref(const seq_t *seq, int is_hap, seq_t *hap1, seq_t *hap2)
{
	int i;
	seq_t *ret[2];
	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = seq->l; ret[1]->l = seq->l;
	ret[0]->m = seq->m; ret[1]->m = seq->m;
	ret[0]->s = (bit8_t*)calloc(seq->m, 1);
	ret[1]->s = (bit8_t*)calloc(seq->m, 1);
	for (i = 0; i != seq->l; ++i) {
		int c;
		c = ret[0]->s[i] = ret[1]->s[i] = 0x50 | nst_nt4_table[(int)seq->s[i]];
		if ((c&0xf) < 4 && ran_uniform() < MUT_RATE) { // mutation
			if (ran_uniform() >= INDEL_FRAC) { // substitution
				double r = ran_uniform();
				c = (c + (int)(r * 3.0 + 1)) & 3;
				if (is_hap || ran_uniform() < 0.333333) { // hom
					ret[0]->s[i] = ret[1]->s[i] = 0x50|c;
				} else { // het
					ret[ran_uniform()<0.5?0:1]->s[i] = 0x50|c;
				}
			} else { // indel
				if (ran_uniform() < 0.5) { // deletion
					if (is_hap || ran_uniform() < 0.333333) { // hom-del
						ret[0]->s[i] = ret[1]->s[i] = 0x55;
					} else { // het-del
						ret[ran_uniform()<0.5?0:1]->s[i] = 0x55;
					}
				} else { // insersion
					int ins = (int)(ran_uniform() * 4.0);
					if (is_hap || ran_uniform() < 0.333333) { // hom-ins
						ret[0]->s[i] = ret[1]->s[i] = ins << 4 | (c&0xf);
					} else { // het-ins
						ret[ran_uniform()<0.5?0:1]->s[i] = ins << 4 | (c&0xf);
					}
				}
			}
		}
	}
}
void maq_print_mutref(const char *name, const seq_t *seq, seq_t *hap1, seq_t *hap2)
{
	int i;
	for (i = 0; i != seq->l; ++i) {
		int c[3];
		c[0] = nst_nt4_table[(int)seq->s[i]];
		c[1] = hap1->s[i]; c[2] = hap2->s[i];
		if (c[0] >= 4) continue;
		if (c[1] != (c[0]|0x50) || c[2] != (c[0]|0x50)) {
			printf("%s\t%d\t", name, i+1);
			if (c[1] == c[2]) { // hom
				if ((c[1]&0xf) < 5 && c[1]>>4 == 5) { // substitution
					printf("%c\t%c\t-\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
				} else if (c[1]>>4 < 5) { // ins
					printf("-\t%c\t-\n", "ACGTN"[c[1]>>4&0xf]);
				} else if ((c[1]&0xf) == 5) { // del
					printf("%c\t-\t-\n", "ACGTN"[c[0]]);
				} else assert(0);
			} else { // het
				if ((c[1]&0xf) < 5 && c[1]>>4 == 5 && (c[2]&0xf) < 5 && c[2]>>4 == 5) { // substitution
					printf("%c\t%c\t+\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0xf)|1<<(c[2]&0xf)]);
				} else if (c[1]>>4 < 5) { // ins1
					printf("-\t%c\t+\n", "ACGTN"[c[1]>>4&0xf]);
				} else if (c[2]>>4 < 5) { // ins2
					printf("-\t%c\t+\n", "ACGTN"[c[2]>>4&0xf]);
				} else if ((c[1]&0xf) == 5) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else if ((c[2]&0xf) == 5) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else assert(0);
			}
		}
	}
}

void maq_simulate_core(FILE *fpout1, FILE *fpout2, FILE *fp_fa, gzFile fp_par, int is_hap, bit64_t N,
					   int dist, int std_dev, int size_l, int size_r)
{
	fqc_t *fqc;
	seq_t seq, rseq[2];
	bit64_t tot_len, ii;
	int i, l, n_ref, k;
	char name[256], qstr[256], str[1024];
	double Q_table[MAX_QUAL+1];
	int size[2];
	bit8_t tmp_seq[2][256], *target;
	fqc = (fqc_t*)malloc(sizeof(fqc_t));
	l = gzread(fp_par, fqc, sizeof(fqc_t));
	for (i = 0; i <= MAX_QUAL; ++i)
		Q_table[i] = exp(-log(10.0) * i / 10.0);
	INIT_SEQ(seq);
	ran_seed();
	seq_set_block_size(0x1000000);
	assert(size_l <= fqc->l_max && size_r <= fqc->l_max);
	if (size_l == 0) size_l = fqc->l_max;
	if (size_r == 0) size_r = fqc->l_max;
	size[0] = size_l; size[1] = size_r;
	tot_len = n_ref = 0;
	while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		tot_len += l;
		++n_ref;
	}
	fprintf(stderr, "-- %d sequences, total length: %llu\n", n_ref, tot_len);
	rewind(fp_fa);
	while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		bit64_t n_pairs = (bit64_t)((long double)l / tot_len * N + 0.5);
		if (l < dist + 2 * std_dev) {
			fprintf(stderr, "[maq_simulate_core] Each reference sequence should be longer than dist+2*stddev. Abort!\n");
			exit(1);
		}
		maq_mut_diref(&seq, is_hap, rseq, rseq+1);
		maq_print_mutref(name, &seq, rseq, rseq+1);
		for (ii = 0; ii != n_pairs; ++ii) {
			double ran;
			int d, pos, s[2], begin, end, is_flip = 0;
			FILE *fpo[2];
			do {
				ran = ran_normal();
				ran = ran * std_dev + dist;
				if (ran < fqc->l_max) ran = fqc->l_max;
				d = (int)(ran + 0.5);
				pos = (int)((l - d + 1) * ran_uniform());
			} while (pos < 0 || pos >= seq.l);
			// flip or not
			if (ran_uniform() < 0.5) {
				fpo[0] = fpout1; fpo[1] = fpout2;
				s[0] = size[0]; s[1] = size[1];
			} else {
				fpo[1] = fpout1; fpo[0] = fpout2;
				s[1] = size[0]; s[0] = size[1];
				is_flip = 1;
			}
			// generate the read sequences
			target = rseq[ran_uniform()<0.5?0:1].s; // haploid from which the reads are generated
			for (i = pos, k = 0, begin = 0; i < seq.l; ++i) {
				int c = target[i];
				if ((c&0xf) == 5) continue; // deletion
				if (begin == 0) {
					begin = i;
					if ((c&0xf0) < 4) continue; // skip ins at the first base, FIXME!
				}
				if ((c>>4) < 4) {
					tmp_seq[0][k++] = c>>4;
					if (k == s[0]) break;
				}
				tmp_seq[0][k++] = c&0xf;
				if (k == s[0]) break;
			}
			for (i = pos + d - 1, k = 0, end = 0; i >= 0; --i) {
				int c = target[i];
				if ((c&0xf) == 5) continue; // deletion
				if (end == 0) end = i;
				tmp_seq[1][k++] = c&0xf;
				if (k == s[1]) break;
				if ((c>>4) < 4) {
					tmp_seq[1][k++] = c>>4;
					if (k == s[1]) break;
				}
			}
			// start to print out
			if (!is_flip) sprintf(str, "%s_%u_%u_%llx", name, begin+1, end+1, ii);
			else sprintf(str, "%s_%u_%u_%llx", name, begin+1, end+1, ii);
			// print forward read
			fprintf(fpo[0], "@%s/%d\n", str, is_flip+1);
			gen_qual(fqc, qstr, s[0]); // generate qual
			for (i = 0; i < s[0]; ++i) {
				int c = tmp_seq[0][i];
				if (c > 4) c = 4;
				if (c < 4) {
					if (ran_uniform() < Q_table[qstr[i]-33])
						c = (c + (int)(ran_uniform() * 3.0 + 1)) & 3;
				}
				fputc("ACGTN"[c], fpo[0]);
			}
			fprintf(fpo[0], "\n+\n%s\n", qstr);
			// print reverse read
			fprintf(fpo[1], "@%s/%d\n", str, 2-is_flip);
			gen_qual(fqc, qstr, s[1]);
			for (i = 0; i < s[1]; ++i) {
				int c = tmp_seq[1][i];
				if (c > 4) c = 4;
				if (c < 4) {
					c = 3 - c; // complement
					if (ran_uniform() < Q_table[qstr[i]-33])
						c = (c + (int)(ran_uniform() * 3.0 + 1)) & 3;
				}
				fputc("ACGTN"[c], fpo[1]);
			}
			fprintf(fpo[1], "\n+\n%s\n", qstr);
		}
	}
	free(seq.s);
	free(fqc);
}

static int simu_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   maq simulate [options] <read1.out> <read2.out> <ref.fasta> <simupar.dat>\n\n");
	fprintf(stderr, "Options: -d INT        outer distance between the two ends [170]\n");
	fprintf(stderr, "         -s INT        standard deviation [20]\n");
	fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
	fprintf(stderr, "         -1 INT        length of the first read\n");
	fprintf(stderr, "         -2 INT        length of the second read\n");
	fprintf(stderr, "         -r FLOAT      rate of mutations [0.001]\n");
	fprintf(stderr, "         -R FLOAT      fraction of 1bp indels [0.1]\n");
	fprintf(stderr, "         -h            haploid mode\n");
	fprintf(stderr, "\n");
	return 1;
}

int maq_simulate(int argc, char *argv[])
{
	bit64_t N;
	int dist, std_dev, c, size_l, size_r, is_hap = 0;
	FILE *fpout1, *fpout2, *fp_fa;
	gzFile fp_par;
	N = 1000000; dist = 170; std_dev = 20;
	size_l = size_r = 0;
	while ((c = getopt(argc, argv, "d:s:N:1:2:r:R:h")) >= 0) {
		switch (c) {
		case 'd': dist = atoi(optarg); break;
		case 's': std_dev = atoi(optarg); break;
		case 'N': N = atoi(optarg); break;
		case '1': size_l = atoi(optarg); break;
		case '2': size_r = atoi(optarg); break;
		case 'r': MUT_RATE = atof(optarg); break;
		case 'R': INDEL_FRAC = atof(optarg); break;
		case 'h': is_hap = 1; break;
		}
	}
	if (argc - optind < 4) return simu_usage();
	fpout1 = fopen(argv[optind+0], "w");
	fpout2 = fopen(argv[optind+1], "w");
	fp_fa = (strcmp(argv[optind+2], "-") == 0)? stdin : fopen(argv[optind+2], "r");
	fp_par = gzopen(argv[optind+3], "r");
	if (fpout1 == 0)
		fprintf(stderr, "[maq_simulate] cannot write to file '%s'", argv[optind+0]);
	if (fpout2 == 0)
		fprintf(stderr, "[maq_simulate] cannot write to file '%s'", argv[optind+1]);
	if (fp_fa == 0)
		fprintf(stderr, "[maq_simulate] cannot open reference sequence file '%s'", argv[optind+2]);
	if (fp_par == 0)
		fprintf(stderr, "[maq_simulate] cannot open parameter file '%s'\n", argv[optind+3]);
	assert(fpout1 && fpout2 && fp_fa && fp_par);
	maq_simulate_core(fpout1, fpout2, fp_fa, fp_par, is_hap, N, dist, std_dev, size_l, size_r);
	fclose(fpout1);fclose(fpout2); fclose(fp_fa); gzclose(fp_par);
	return 0;
}

/* simustat */

static void simustat_core(gzFile fp, int Q_thres)
{
	static int n_valid_abflags = 5;
	static int valid_abflags[5] = {PAIRFLAG_FF, PAIRFLAG_FR, PAIRFLAG_RF, PAIRFLAG_RR, PAIRFLAG_DIFFCHR};
	bit32_t wc_single[10], wc_pair[10], abpair[256][10], tot_single[10], tot_pair[10];
	maqmap_t *m = maqmap_read_header(fp);
	maqmap1_t *m1, mm1;
	char str[256];
	int i, k, tc[2][2];
	m1 = &mm1;
	memset(wc_single, 0, 40); memset(tot_single, 0, 40);
	memset(wc_pair, 0, 40); memset(tot_pair, 0, 40);
	memset(abpair, 0, 4 * 256 * 10);
	memset(tc[2], 0, 4 * sizeof(int));
	while (maqmap_read1(fp, m1)) {
		int is_correct;
		bit32_t p1, p2;
		if (m1->info2 > 60) continue;
		p1 = p2 = 0xffffffffu;
		strcpy(str, m1->name);
		for (k = 0, i = strlen(str) - 1; i >= 0; --i) {
			if (str[i] == '_') {
				if (k == 0) {
					++k;
					str[i] = 0;
				} else if (k == 1) {
					p2 = atoi(str + i + 1);
					++k;
					str[i] = 0;
				} else if (k == 2) {
					p1 = atoi(str + i + 1);
					++k;
					str[i] = 0;
				}
			}
		}
		assert(p1 != 0xffffffffu && p2 != 0xffffffffu);
		is_correct = 1;
		if (strcmp(str, m->ref_name[m1->seqid])) is_correct = 0;
		else {
			if (p1 < p2) {
				--p1; // zero based
				if (m1->pos&1) { // reverse
					if (abs((m1->pos>>1) - (p2 - m1->size)) > 1) is_correct = 0;
				} else { // forward
					if (abs((m1->pos>>1) - p1) > 1) is_correct = 0;
				}
			} else {
				--p2; // zero based
				if (m1->pos&1) { // reverse
					if (abs((m1->pos>>1) - (p1 - m1->size)) > 1) is_correct = 0;
				} else { // forward
					if (abs((m1->pos>>1) - p2) > 1) is_correct = 0;
				}
			}
		}
		k = m1->seq[MAX_READLEN-1] / 10;
		if (k > 9) k = 9;
		++tot_single[k];
		if (!is_correct) ++wc_single[k];
		k = m1->map_qual / 10;
		if (k > 9) k = 0;
		++tot_pair[k];
		// ++tc[is_correct][m1->info1>>4&1]; // debug only
		if (!is_correct) ++wc_pair[k];
		if (!is_correct && m1->map_qual >= Q_thres) {
			fprintf(stderr, "%s\t%u\t%c\t%d\t%d\t%d\t%d\n", m1->name, m1->pos>>1, (m1->pos&1)? '-' : '+',
					m1->map_qual, m1->seq[MAX_READLEN-1], m1->c[0], m1->c[1]);
		}
		// count abnormal pairs
		if (m1->flag != (PAIRFLAG_FR | PAIRFLAG_PAIRED) && m1->flag != PAIRFLAG_NOMATCH) {
			k = m1->alt_qual / 10;
			if (k > 9) k = 9;
			++abpair[m1->flag][k];
		}
	}
	// print out
	// fprintf(stderr, "%d,%d,%d,%d\n", tc[0][0], tc[0][1], tc[1][0], tc[1][1]);
	printf("%2s%12s / %-12s%12s / %-12s", "Q", "SE_wrong", "SE_tot", "PE_wrong", "PE_tot");
	for (i = 0; i != n_valid_abflags; ++i)
		printf("%12d", valid_abflags[i]);
	printf("\n");
	for (k = 0; k <= 9; ++k) {
		printf("%dx%12d / %-12d%12d / %-12d", k, wc_single[k], tot_single[k], wc_pair[k], tot_pair[k]);
		for (i = 0; i != n_valid_abflags; ++i)
			printf("%12d", abpair[valid_abflags[i]][k]);
		printf("\n");
	}
	maq_delete_maqmap(m);
}
int maq_simustat(int argc, char *argv[])
{
	gzFile fp;
	int Q = 100;
	if (argc < 2) {
		fprintf(stderr, "Usage: maq simustat <simu_align.map> [<Q>=100]\n");
		return 1;
	}
	fp = gzopen(argv[1], "r");
	if (argc >= 3) Q = atoi(argv[2]);
	assert(fp);
	simustat_core(fp, Q);
	gzclose(fp);
	return 0;
}
