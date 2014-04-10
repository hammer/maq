#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include "bfa.h"
#include "main.h"
#include "seq.h"

const int nst_color_space_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4};

/*
  AA/CC/GG/TT -> 0 (Blue)
  AC/CA/GT/TG -> 1 (Green)
  AG/GA/CT/TC -> 2 (Orange)
  AT/TA/CG/GC -> 3 (Red)
 */
static void ma_fasta2csfa_core(FILE *fpout, FILE *fpin)
{
	seq_t seq;
	int i, c1, c2, c;
	char name[256], comment[4096];
	INIT_SEQ(seq);
	seq_set_block_size(0x800000); // use longer block size (8M)
	while (seq_read_fasta(fpin, &seq, name, comment) >= 0) {
		fprintf(fpout, ">%s %s", name, comment);
		c1 = nst_nt4_table[(int)seq.s[0]];
		for (i = 1; i < seq.l; ++i) {
			if ((i-1)%60 == 0) fputc('\n', fpout); 
			c2 = nst_nt4_table[(int)seq.s[i]];
			c = (c1 == 4 || c2 == 4)? 4 : nst_color_space_table[((1<<c1)|(1<<c2))&0xf];
			fputc("ACGTN"[c], fpout);
			c1 = c2;
		}
		fputc('\n', fpout);
	}
	free(seq.s);
}

/* fakemut */

typedef struct
{
	bit8_t type, base;
	int pos;
} fakemut_t;

static void ma_fakemut_core(FILE *fpout, FILE *fpin, double rate, double irate)
{
	seq_t seq;
	int n_mut, m_mut;
	fakemut_t *mutarray;
	char name[256], comment[4096];
	srand48(time(0));
	INIT_SEQ(seq);
	seq_set_block_size(0x800000); // use longer block size (8M)
	n_mut = m_mut = 0;
	mutarray = 0;
	while (seq_read_fasta(fpin, &seq, name, comment) >= 0) {
		int i, cur_pos, cur_mut;
		n_mut = 0; // rewind
		fprintf(fpout, ">%s %s\n", name, comment);
		// first round: fill mutarray
		for (i = 0; i < seq.l; ++i) {
			double r = drand48();
			if (nst_nt4_table[(int)seq.s[i]] < 4 && r < rate) {
				bit8_t type, base;
				r = drand48();
				base = (nst_nt4_table[(int)seq.s[i]] + (int)(r * 3.0 + 1)) & 3;
				if (r < irate/2.0) type = 'I'; // insert
				else if (r < irate) type = 'D'; // deletion
				else type = 'S'; // substitution
				if (n_mut == m_mut) { // enlarge
					m_mut += 0x10000;
					mutarray = (fakemut_t*)realloc(mutarray, sizeof(fakemut_t) * m_mut);
				}
				mutarray[n_mut].type = type;
				mutarray[n_mut].base = "ACGT"[base];
				mutarray[n_mut].pos = i;
				++n_mut; // push back
			}
		}
		// second round: output the sequence and mutations
		for (i = cur_pos = cur_mut = 0; i < seq.l; ++i) {
			if (cur_pos && cur_pos%60 == 0) fputc('\n', fpout);
			if (i == mutarray[cur_mut].pos) {
				bit8_t type = mutarray[cur_mut].type;
				if (type == 'S') {
					fputc(mutarray[cur_mut].base, fpout);
					fprintf(stderr, "%s\t%d\t%c\t%c\t99\n", name, cur_pos + 1, mutarray[cur_mut].base, seq.s[i]);
					++cur_pos;
				} else if (type == 'I') {
					fputc(mutarray[cur_mut].base, fpout); // insert before the current base
					fputc(seq.s[i], fpout);
					fprintf(stderr, "%s\t%d\t%c\t-\t99\n", name, cur_pos + 1, mutarray[cur_mut].base);
					cur_pos += 2;
				} else { // else, it is a deletion.
					fprintf(stderr, "%s\t%d\t-\t%c\t99\n", name, cur_pos + 1, seq.s[i]);
				}
				++cur_mut;
			} else {
				fputc(seq.s[i], fpout);
				++cur_pos;
			}
		}
		fputc('\n', fpout);
	}
	free(seq.s); free(mutarray);
}

/* fasta2bfa */

static int nst_fasta2bfa1(FILE *fp_fa, FILE *fp_bfa, seq_t *seq)
{
	int i, len;
	char name[256];
	bit64_t s, m;
	nst_bfa1_t *bfa1;

	len = seq_read_fasta(fp_fa, seq, name, 0);
	if (len < 0) return -1; // no sequence
	bfa1 = nst_new_bfa1();
	bfa1->ori_len = len;
	bfa1->len = len>>5;
	if (len&0x1f) ++(bfa1->len);
	bfa1->seq = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	bfa1->mask = (bit64_t*)malloc(sizeof(bit64_t) * bfa1->len);
	bfa1->name = strdup(name);
	m = s = 0ull;
	for (i = 0; i != len; ++i) {
		int tmp = nst_nt4_table[(int)seq->s[i]];
		s <<= 2; m <<= 2;
		if (tmp < 4) {
			s |= tmp; m |= 0x3;
		} else m &= ~0x3;
		if ((i&0x1f) == 0x1f) { // add to l
			bfa1->seq[i>>5] = s;
			bfa1->mask[i>>5] = m;
		}
	}
	if (len&0x1f) {
		s <<= (32 - (i&0x1f)) << 1; m <<= (32 - (i&0x1f)) << 1;
		bfa1->seq[len>>5] = s;
		bfa1->mask[len>>5] = m;
	}
	i = strlen(bfa1->name) + 1;
	fwrite(&i, sizeof(int), 1, fp_bfa);
	fwrite(bfa1->name, sizeof(char), i, fp_bfa);
	fwrite(&bfa1->ori_len, sizeof(int), 1, fp_bfa);
	fwrite(&bfa1->len, sizeof(int), 1, fp_bfa);
	fwrite(bfa1->seq, sizeof(bit64_t) * bfa1->len, 1, fp_bfa);
	fwrite(bfa1->mask, sizeof(bit64_t) * bfa1->len, 1, fp_bfa);
	nst_delete_bfa1(bfa1);
	return len;
}
int nst_fasta2bfa(FILE *fp_fa, FILE *fp_bfa)
{
	seq_t seq;
	int n = 0;
	INIT_SEQ(seq);
	seq_set_block_size(0x800000); // use longer block size (8M)
	while (nst_fasta2bfa1(fp_fa, fp_bfa, &seq) >= 0) ++n;
	free(seq.s);
	return n;
}
int ma_fasta2bfa(int argc, char *argv[])
{
	FILE *fp_fa, *fp_bfa;
	int n;
	fp_fa = fp_bfa = 0;
	if (argc <= 2) {
		fprintf(stderr, "Usage: maq fasta2bfa <in.fasta> <out.bfa>\n");
		return 1;
	}
	fp_fa = (strcmp(argv[1], "-") == 0)? stdin : fopen(argv[1], "r");
	fp_bfa = (strcmp(argv[2], "-") == 0)? stdin : fopen(argv[2], "w");
	if (fp_fa == 0 || fp_bfa == 0) {
		fprintf(stderr, "ERROR: fail to open file(s).\n");
		return 2;
	}
	n = nst_fasta2bfa(fp_fa, fp_bfa);
	fprintf(stderr, "-- %d sequences have been converted.\n", n);
	fclose(fp_fa); fclose(fp_bfa);
	return 0;
}
int ma_fasta2csfa(int argc, char *argv[])
{
	FILE *fp_fa;
	if (argc < 2) {
		fprintf(stderr, "Usage: maq fasta2csfa <in.fasta>\n");
		return 1;
	}
	fp_fa = (strcmp(argv[1], "-") == 0)? stdin : fopen(argv[1], "r");
	assert(fp_fa);
	ma_fasta2csfa_core(stdout, fp_fa);
	fclose(fp_fa);
	return 0;
}
static int fakemut_usage()
{
	fprintf(stderr, "Usage: maq fakemut [-r 0.001] [-R 0.1] <in.fasta>\n");
	return 1;
}
int ma_fakemut(int argc, char *argv[])
{
	FILE *fp_fa;
	double mutrate = 0.001;
	double indelrate = 0.1;
	int c;
	/* mutrate is the overall mutation rate
	   mutrate*indelrate is the indel rate
	   mutrate*(1-indelrate) is the substitution rate
	 */
	while ((c = getopt(argc, argv, "r:R:")) >= 0) {
		switch (c) {
		case 'r': mutrate = atof(optarg); break;
		case 'R': indelrate = atof(optarg); break;
		}
	}
	if (optind == argc) return fakemut_usage();
	fp_fa = (strcmp(argv[optind], "-") == 0)? stdin : fopen(argv[optind], "r");
	assert(fp_fa);
	ma_fakemut_core(stdout, fp_fa, mutrate, indelrate);
	fclose(fp_fa);
	return 0;
}
