#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "read.h"
#include "match.hh"
#include "main.h"
#include "maqmap.h"

match_data_t *new_match_data()
{
	return (match_data_t*)calloc(1, sizeof(match_data_t));
}
void delete_match_data(match_data_t *d)
{
	if (d == 0) return;
	int i;
	// d->index and d->index4 should always be freed in other places 
	if (d->lname) {
		for (i = 0; i < d->n_lname; ++i)
			free(d->lname[i]);
	}
	free(d->lname); free(d->match); free(d->pair);
	free(d);
}
// generate a set of filters. word_length=size, max_mismatch=m
bit64_t *ma_gen_filters(int size, int m, int *l)
{
	int n, i, j, k, c, *tmp;
	bit64_t f, *part, *filter;
	n = m<<1; j = k = 1;
	for (i = n; i >= n - m + 1; --i) {
		j *= i; k *= n - i + 1;
	}
	*l = j / k;
	filter = (bit64_t*)malloc(sizeof(bit64_t) * (*l));
	if (size == 0) { // then return a NULL filter.
		for (i = 0; i != n; ++i) filter[i] = 0;
		return filter;
	}
	part = (bit64_t*)malloc(sizeof(bit64_t) * n);
	tmp = (int*)malloc(sizeof(int) * n);
	for (i = 0, j = size, k = 0; i != n; ++i) {
		tmp[i] = j / (n - i);
		j -= tmp[i];
		k += tmp[i];
		part[i] = bit64_t((1ull<<(tmp[i]<<1)) - 1) << ((k - tmp[i])<<1);
	}
	for (i = c = 0; i < (1<<n); ++i) {
		for (j = 0, k = 0, f = 0ull; j != n; ++j) {
			if (i&1<<j) {
				++k; f |= part[j];
			}
		}
		if (k == m) filter[c++] = f;
	}
	free(part);
	free(tmp);
	// re-order the filters such that the 2n-th filter and 2n+1-th filter is complimentary to each other
	for (i = 0; i < *l - 2; i += 2) {
		if ((filter[i]&filter[i+1]) == 0ull) continue;
		for (k = i + 2; k < *l; ++k) {
			if ((filter[k]&filter[i]) == 0ull) {
				f = filter[k]; filter[k] = filter[i+1]; filter[i+1] = f;
				break;
			}
		}
	}
	return filter;
}
match_aux_t *new_match_aux(int size_l, int size_r, int m)
{
	match_aux_t *o;
	int i;
	o = (match_aux_t*)calloc(1, sizeof(match_aux_t));
	o->size_l = size_l; o->size_r = size_r;
	o->n_mismatch = m;
	o->q_rate = 30; // default q_rate (0.001)
	o->max_dist = 250;
	o->max_hits = 250;
	o->min_dist = 0;
	o->dump_file = 0;
	o->hits_fp = 0;
	o->is_color = o->is_mm = o->is_p2diff = 0;
	o->is_sw = 1;
	o->methy_mode = 0;
	o->adapter = 0;
	o->max_err10 = 7;
	{ // generate filters and masks
		int n_filters;
		bit64_t *filters_l = ma_gen_filters((size_l > SEED_LENGTH)? SEED_LENGTH : size_l, m, &n_filters);
		bit64_t *filters_r = ma_gen_filters((size_r > SEED_LENGTH)? SEED_LENGTH : size_r, m, &n_filters);
		o->filters    = (bit64_t*)calloc(n_filters * 2, sizeof(bit64_t));
		o->masks      = (read_t*)calloc(n_filters * 2, sizeof(read_t));
		o->shift_seed = (int*)calloc(n_filters * 2, sizeof(int));
		for (i = 0; i != n_filters; ++i) {
			o->filters[i<<1]      = filters_l[i];
			o->filters[(i<<1)|1]  = filters_r[i];
			o->masks[i<<1]        = ~(~read_t(0)<<(size_l<<1));
			o->masks[(i<<1)|1]    = ~(~read_t(0)<<(size_r<<1));
			o->shift_seed[i<<1]   = (size_l <= SEED_LENGTH)? 0 : (size_l-SEED_LENGTH)<<1;
			o->shift_seed[i<<1|1] = (size_r <= SEED_LENGTH)? 0 : (size_r-SEED_LENGTH)<<1;
		}
		free(filters_l); free(filters_r);
		o->n_filters = n_filters * 2;
	}
	o->log_n[0] = 0;
	for (i = 1; i != 256; ++i) o->log_n[i] = int(4.343 * log(i) + 0.5);
	return o;
}
void delete_match_aux(match_aux_t *o)
{
	free(o->filters); free(o->masks); free(o->shift_seed); free(o->adapter); free(o->dump_file);
	if (o->hits_fp) gzclose(o->hits_fp);
	free(o);
}

static inline void ma_trim_adapter1(int alen, read_t x, int size, bit8_t seq[], matches_t mm, int counter[], int counter2[],
									int trim_min)
{
	read_t y;
	int j, k, m, trim_len = 0;
	for (k = 0, y = 0; k != size; ++k) y = y<<2 | seq[k]>>6;
	// adapter in the sequence
	for (k = 0; k < size - alen; ++k) {
		for (j = m = 0; j != alen; ++j) {
			if ((y>>((size-1-k-j)<<1)&0x3) != (x>>((alen-1-j)<<1)&0x3)) ++m;
			if (m > 5) break;
		}
		if (m <= 5 && m <= (alen - 1) / 7) {
			trim_len = size;
			++counter[alen];
			goto end_trim;
		}
	}
	// adapter at the 5'-end
	for (j = 1; j < alen - 20; ++j) {
		int end = (size < alen - j)? size : alen - j;
		for (k = m = 0; k < end; ++k) {
			if ((y>>((size-1-k)<<1)&0x3) != (x>>((alen-1-j-k)<<1)&0x3)) ++m;
			if (m > 5) break;
		}
		if (m <= 5 && m <= (alen - j - 1) / 7) {
			trim_len = size;
			++counter[alen];
			goto end_trim;
		}
	}
	// adapter at the 3'-end
	for (j = 0; j != alen - 1; ++j) {
		for (m = k = 0; k != alen - j; ++k) {
			if ((y>>(k<<1)&0x3) != (x>>((j+k)<<1)&0x3)) ++m;
			if (m > 5) break;
		}
		int tmp = alen - j;
		if (m <= 5 && m <= (tmp - 1) / 7) {
			trim_len = tmp;
			++counter[j];
			goto end_trim;
		}
	}

end_trim:
	if (trim_len > trim_min && (mm&((1ull<<trim_len)-1)) != 0) {
		if (trim_len < size) ++counter2[alen-trim_len];
		for (k = size - trim_len; k < size; ++k)
			seq[k] &= 0xc0;
	}
}

// match_info_t::mm1 must be correctly set.
void ma_trim_adapter(const char adaptor[], const match_data_t *d, longreads_t *lr)
{
	read_t x;
	int n, i, j, *counter, *counter2, alen, trim_min = 0;
	int size[2] = { lr->size_l, lr->size_r };
	int len = strlen(adaptor);
	clock_t begin = clock();
	fprintf(stderr, "[ma_trim_adapter] trimming the 3'-adapter sequence...\n");
	counter = (int*)calloc(len + 1, sizeof(int));
	counter2 = (int*)calloc(len + 1, sizeof(int));
	alen = strlen(adaptor);
	for (i = 0, x = 0; i != len; ++i)
		x = x<<2 | (nst_nt4_table[(int)adaptor[i]]&3);
	if (d == 0) trim_min = 4;
	// drop reads
	for (i = 0; i != lr->n_reads; ++i) {
		if (size[i&1] == 0) continue;
		matches_t mm = (d && d->match[i].i1 != MA_NO_MATCH)? d->match[i].mm1 : matches_t(0xffff);
		ma_trim_adapter1(alen, x, size[i&1], lr->seq[i], mm, counter, counter2, trim_min);
	}
	for (j = len-2, n = 0; j >= 0; --j) {
		n += counter2[j];
		fprintf(stderr, "[ma_trim_adapter] match length: %d; count: %d:%d\n", len - j, counter[j], counter2[j]);
	}
	fprintf(stderr, "[ma_trim_adapter] match length: %d+; count: %d\n", len, counter[len]);
	fprintf(stderr, "[ma_trim_adapter] %d reads possibly contains adaptor contamination.\n", n += counter[len]);
	fprintf(stderr, "[ma_trim_adapter] CPU time: %.3f\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
	free(counter); free(counter2);
}

void mapping_stat(gzFile fpin)
{
	int i;
	bit64_t tot_len, tot_reads, flagc[256], qualc_se[26], qualc_pe[26], mmc[16], tot_err;
	maqmap_t *mm;
	maqmap1_t *m1, mm1;
	m1 = &mm1;
	tot_len = tot_reads = tot_err = 0;
	memset(flagc, 0, sizeof(bit64_t) * 256);
	memset(mmc, 0, sizeof(bit64_t) * 16);
	memset(qualc_se, 0, sizeof(bit64_t) * 26);
	memset(qualc_pe, 0, sizeof(bit64_t) * 26);
	mm = maqmap_read_header(fpin);
	while (maqmap_read1(fpin, m1)) {
		int q;
		++tot_reads; tot_len += m1->size;
		++flagc[m1->flag];
		q = m1->seq[MAX_READLEN-1];
		if (q > 99) q = 99;
		++qualc_se[q/10];
		q = m1->map_qual;
		if (q > 99) q = 99;
		++qualc_pe[q/10];
		++mmc[m1->info1&0xf];
		tot_err += m1->info1&0xf;
	}
	maq_delete_maqmap(mm);
	// print
	printf("\n");
	printf("-- Total number of reads: %lld\n", tot_reads);
	printf("-- Sum of read length: %lld\n", tot_len);
	printf("-- Error rate: %lf\n", (double)tot_err/tot_len);
	printf("-- Average read length: %.2lf\n\n", (double)tot_len/tot_reads);
	printf("-- Mismatch statistics:\n\n");
	for (i = 0; i != 16; ++i)
		if (mmc[i]) printf("-- MM %2d    %15lld\n", i, mmc[i]);
	printf("\n");
	printf("-- Mapping quality statistics:\n\n");
	for (i = 0; i != 10; ++i)
		if (qualc_se[i] || qualc_pe[i])
			printf("-- MQ %1d0-%1d9 %15lld %15lld\n", i, i, qualc_se[i], qualc_pe[i]);
	printf("\n");
	printf("-- Flag statistics:\n\n");
	for (i = 0; i != 256; ++i)
		if (flagc[i]) printf("-- FG %3d   %15lld\n", i, flagc[i]);
	printf("\n");
}
int ma_mapstat(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "Usage: maq mapstat <in.map>\n");
		return 1;
	}
	gzFile fp = (strcmp(argv[1], "-") == 0)? gzdopen(STDIN_FILENO, "r") : gzopen(argv[1], "r");
	mapping_stat(fp);
	gzclose(fp);
	return 0;
}
