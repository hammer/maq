#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "assemble.h"

#ifdef __cplusplus
extern "C" {
#endif
	long double expl(long double);
	long double logl(long double);
#ifdef __cplusplus
}
#endif

/*
  P(<b1,b2>) = \theta \sum_{i=1}^{N-1} 1/i
  P(D|<b1,b2>) = \sum_{k=1}^{N-1} p_k 1/2 [(k/N)^n_2(1-k/N)^n_1 + (k/N)^n1(1-k/N)^n_2]
  p_k = i/k / \sum_{i=1}^{N-1} 1/i
 */
static void cal_het(assemble_aux_t *aa)
{
	int k, n1, n2;
	double sum_harmo; // harmonic sum
	double poly_rate;
	double p1 = 0.0, p3 = 0.0; // just for testing
	
	aa->lhet = (double*)calloc(256 * 256, sizeof(double));
	sum_harmo = 0.0;
	for (k = 1; k <= aa->n_hap - 1; ++k)
		sum_harmo += 1.0 / k;
	for (n1 = 0; n1 < 256; ++n1) {
		for (n2 = 0; n2 < 256; ++n2) {
			long double sum = 0.0;
			double lC = lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1); // \binom{n1+n2}{n1}
			for (k = 1; k <= aa->n_hap - 1; ++k) {
				double pk = 1.0 / k / sum_harmo;
				double log1 = log((double)k/aa->n_hap);
				double log2 = log(1.0 - (double)k/aa->n_hap);
				sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
			}
			aa->lhet[n1<<8|n2] = lC + logl(sum);
			if (n1 == 17 && n2 == 3) p3 = lC + logl(expl(logl(0.5) * 20));
			if (n1 == 19 && n2 == 1) p1 = lC + logl(expl(logl(0.5) * 20));
		}
	}
	poly_rate = aa->hetero_rate * sum_harmo;
	aa->q_r = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));
	fprintf(stderr, "[cal_het] harmonic sum: %lf\n", sum_harmo);
	fprintf(stderr, "[cal_het] het penalty: %.2lf vs. %.2lf\n", aa->q_r,
			-4.343 * log(2.0 * aa->hetero_rate / (1.0 - aa->hetero_rate)));
	fprintf(stderr, "[cal_het] 3 differences out of 20 bases: %.2lf vs. %.2lf\n",
			-4.343 * aa->lhet[17<<8|3], -4.343 * p3);
	fprintf(stderr, "[cal_het] 1 differences out of 20 bases: %.2lf vs. %.2lf\n",
			-4.343 * aa->lhet[19<<8|1], -4.343 * p1);
}

static void cal_coef_ind(assemble_aux_t *aa)
{
	int n, q, k;
	double *lC;

	lC = (double*)calloc(256 * 256, sizeof(double));
	aa->fk = (double*)calloc(256, sizeof(double));
	aa->coef = (double*)calloc(256*256*64, sizeof(double));
	for (n = 1; n != 256; ++n)
		for (k = 1; k <= n; ++k)
			lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
	for (n = 0; n != 256; ++n) aa->fk[n] = 1.0;
	for (q = 1; q != 64; ++q) {
		double e = pow(10.0, -q/10.0);
		double le1 = log(1.0-e);
		for (n = 1; n != 256; ++n) {
			double *coef = aa->coef + (q<<16|n<<8);
			for (k = n; k >= 0; --k)
				coef[k] = -4.343 * (lC[n<<8|k] + (n-k) * le1);
		}
	}
	free(lC);
}

/** initialize the helper structure */
static void cal_coef(assemble_aux_t *aa)
{
	int k, n, q;
	long double sum_a[257], b[256], q_c[256], tmp[256], fk2[256];
	double *lC;

	lC = (double*)calloc(256 * 256, sizeof(double));
	// aa->lhet will be allocated and initialized 
	aa->fk = (double*)calloc(256, sizeof(double));
	aa->coef = (double*)calloc(256*256*64, sizeof(double));
	aa->fk[0] = fk2[0] = 1.0;
	for (n = 1; n != 256; ++n) {
		aa->fk[n] = ((!aa->is_alt)? pow(aa->theta, n) : pow(n+1, -aa->theta)) * (1.0 - aa->eta) + aa->eta;
		fk2[n] = aa->fk[n>>1];
	}
	for (n = 1; n != 256; ++n)
		for (k = 1; k <= n; ++k)
			lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
	for (q = 1; q != 64; ++q) {
		double e = pow(10.0, -q/10.0);
		double le = log(e);
		double le1 = log(1.0-e);
		for (n = 1; n != 256; ++n) {
			double *coef = aa->coef + (q<<16|n<<8);
			sum_a[n+1] = 0.0;
			for (k = n; k >= 0; --k) { // a_k = \sum_{i=k}^n C^n_k \epsilon^k (1-\epsilon)^{n-k}
				sum_a[k] = sum_a[k+1] + expl(lC[n<<8|k] + k*le + (n-k)*le1);
				b[k] = sum_a[k+1] / sum_a[k];
				if (b[k] > 0.99) b[k] = 0.99;
			}
			for (k = 0; k != n; ++k) // c_k
				q_c[k] = -4.343 * fk2[k] * logl(b[k] / e);
			for (k = 1; k != n; ++k) q_c[k] += q_c[k-1]; // \prod_{i=0}^k c_i
			for (k = 0; k <= n; ++k) { // powl() in 64-bit mode seems broken on my Mac OS X 10.4.9
				tmp[k] = -4.343 * logl(1.0 - expl(fk2[k] * logl(b[k])));
				coef[k] = (k? q_c[k-1] : 0) + tmp[k]; 
			}
			if (q == 1000 && n == 20) { // This never happens. I use these codes to look at values when debugging.
				double ttt[257], *tt;
				tt = ttt + 1;
				tt[0] = fk2[k];
				for (k = 1; k != n; ++k) tt[k] = tt[k-1] + fk2[k];
				for (k = 0; k != n; ++k) {
					fprintf(stderr, "%d\t%Lf\t%Lg\t%Lg\t%g\t%Lg\t%Lg\n", k, fk2[k], expl(lC[n<<8|k] + k*le + (n-k)*le1),
							expl(-coef[k]/4.343)*expl(tt[k]*le), coef[k], expl(logl(b[k])*fk2[k]), b[k]);
				}
				exit(0);
			}
		}
	}
	free(lC);
}

static int usage(const assemble_aux_t *aa)
{
	fprintf(stderr, "\nUsage:   maq assemble [options] <out.cns> <chr.bfa> <in.map>\n\n");
	fprintf(stderr, "Options: -r FLOAT    expected rate of heterozygotes [%.3f]\n", aa->hetero_rate);
	fprintf(stderr, "         -t FLOAT    dependency coefficient (theta) [%.2f]\n", aa->theta);
	fprintf(stderr, "         -q INT      minimum mapping quality [%d]\n", aa->min_q);
	fprintf(stderr, "         -Q INT      maximum sum of errors [%d]\n", aa->max_err);
	fprintf(stderr, "         -m INT      maximum number of mismatches [%d]\n", aa->max_mm);
	fprintf(stderr, "         -N INT      number of haplotypes (>=2) [%d]\n", aa->n_hap);
	fprintf(stderr, "         -s          use single-end mapping quality\n");
	fprintf(stderr, "         -p          discard abnormal pairs\n");
//	fprintf(stderr, "         -e FLOAT    minimum independency (eta) [0.03] (EXPERIMENTAL)\n"); // only for me...
//	fprintf(stderr, "         -a          alternative model (EXPERIMENTAL)\n"); // only for me...
	fprintf(stderr, "\n");
	return 1;
}

static assemble_aux_t *new_opt()
{
	assemble_aux_t *aa;
	aa = (assemble_aux_t*)calloc(1, sizeof(assemble_aux_t));
	aa->max_err = 60;
	aa->max_mm = 7;
	aa->hetero_rate = 0.001;
	aa->theta = 0.85;
	aa->eta = 0.03;
	aa->n_hap = 2;
	aa->maf = 0.5; // will be removed in future
	return aa;
}
assemble_aux_t *assemble_parse_opt(int argc, char *argv[])
{
	assemble_aux_t *aa;
	int c;
	aa = new_opt();
	while ((c = getopt(argc, argv, "q:Q:r:t:e:aspm:N:")) >= 0) {
		switch (c) {
		case 'p': aa->is_pair_only = 1; break;
		case 'm': aa->max_mm = atoi(optarg); break;
		case 'q': aa->min_q = atoi(optarg); break;
		case 'Q': aa->max_err = atoi(optarg); break;
		case 'r': aa->hetero_rate = atof(optarg); break;
		case 't': aa->theta = atof(optarg); break;
		case 'e': aa->eta = atof(optarg); break;
		case 'a': aa->is_alt = 1; break;
		case 'N': aa->n_hap = atoi(optarg); break;
		case 's': aa->is_single = 1; break;
		default:
			usage(aa);
			free(aa);
			exit(1);
		}
	}

	if (argc - optind < 3) {
		usage(aa);
		exit(1);
	}
	aa->fpout = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdout), "w") : gzopen(argv[optind], "w");
	aa->fp_bfa = fopen(argv[optind+1], "r");
	aa->fp_map = (strcmp(argv[optind+2], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind+2], "r");
	assert(aa->fpout && aa->fp_bfa && aa->fp_map);

	if (aa->theta >= 1.0) aa->theta = 1.0; // capped at 1.0
	if (aa->theta >= 0.999) cal_coef_ind(aa);
	else cal_coef(aa);
	cal_het(aa);
	return aa;
}

void delete_assemble_aux(assemble_aux_t *aa)
{
	if (aa == 0) return;
	free(aa->lhet); free(aa->fk); free(aa->coef);
	gzclose(aa->fpout); gzclose(aa->fp_map); fclose(aa->fp_bfa);
	free(aa);
}

