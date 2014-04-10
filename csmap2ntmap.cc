#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "maqmap.h"
#include "stdhash.hh"
#include "main.h"
#include "stdhash.hh"
#include "bfa.h"

/*
  Here is a delicate example. ref_nt=ATTAAC(RBRBG), read_cs=RBBOG. If we
  decode as ATTGAC(RBGOG), there are one color change and one nt change;
  if we decode as ATTAAC(RBRBG), there are two color changes.

  In DP, if color quality is smaller than COLOR_MM, we will use COLOR_MM
  as the penalty; otherwise, we will use color quality as the
  penalty. This means we always prefer two consistent color changes over
  a nt change, but if a color has high quality, we may prefer one nt
  change.

  In the above example, the penalties of the two types of decoding are
  q(B)+25 and q(B)+q(O), respectively. If q(O)>25, we prefer the first;
  otherwise the second. Note that no matter what we choose, the fourth
  base will get a low nt quality.
 */
#define COLOR_MM 19
#define NUCL_MM  25

const int nst_ntnt2cs_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4 };
const int nst_ntcs2nt_table[] = { 0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0 }; // useless at the moment

static inline int get_readseq(const nst_bfa1_t *l, int pos, int size, bit8_t *ntread)
{
	int k;
	bit64_t *s, *m;
	if (pos + size >= l->ori_len) return 1;
	s = l->seq + (pos>>5); m = l->mask + (pos>>5);
	for (k = 0; k <= size; ++k) {
		int i = (31 - ((pos+k)&0x1f)) << 1;
		ntread[k] = ((*m>>i&0x3) == 3)? (*s>>i&0x3) : 4;
		if (31 - ((pos+k)&0x1f) == 0) {
			++s; ++m;
		}
	}
	return 0;
}
static void cs2nt_dp(int size, const bit8_t *nt_ref, const bit8_t *csseq, bit8_t *nt_read, bit8_t *btarray)
{
	int h[8], curr, last;
	int x, y, xmin, hmin, k;

	// recursion: initial value
	if (nt_ref[0] >= 4) memset(h, 0, sizeof(int) << 2);
	else {
		for (x = 0; x != 4; ++x) h[x] = NUCL_MM;
		h[nt_ref[0]] = 0;
	}
	// recursion: main loop
	curr = 1; last = 0;
	for (k = 1; k <= size; ++k) {
		for (x = 0; x != 4; ++x) {
			int min = 0x7fffffff, ymin = 0;
			for (y = 0; y != 4; ++y) {
				int s = h[last<<2|y];
				if (csseq[k-1] && csseq[k-1]>>6 != nst_ntnt2cs_table[1<<x|1<<y])
					s += ((csseq[k-1]&0x3f) < COLOR_MM)? COLOR_MM : (csseq[k-1]&0x3f); // color mismatch
				if (nt_ref[k] < 4 && nt_ref[k] != x) s += NUCL_MM; // nt mismatch
				if (s < min) {
					min = s; ymin = y;
				}
			}
			h[curr<<2|x] = min; btarray[k<<2|x] = ymin;
		}
		last = curr; curr = 1 - curr; // swap
	}
	// back trace
	hmin = 0x7fffffff; xmin = 0;
	for (x = 0; x != 4; ++x) {
		if (h[last<<2|x] < hmin) {
			hmin = h[last<<2|x]; xmin = x;
		}
	}
	nt_read[size] = xmin;
	for (k = size - 1; k >= 0; --k)
		nt_read[k] = btarray[(k+1)<<2 | nt_read[k+1]];
}
static void cal_nt_qual(int size, const bit8_t *nt_read, bit8_t *seq, bit8_t *tarray)
{
	int k, c1, c2;
	bit8_t *t2array = tarray + size;
	// get the color sequence of nt_read
	c1 = nt_read[0];
	for (k = 1; k <= size; ++k) {
		c2 = nt_read[k]; // in principle, there is no 'N' in nt_read[]
		tarray[k-1] = (c1 >= 4 || c2 >= 4)? 4 : nst_ntnt2cs_table[1<<c1 | 1<<c2];
		c1 = c2;
	}
	for (k = 1; k != size; ++k) {
		int q = 0;
		if (tarray[k-1] == seq[k-1]>>6 && tarray[k] == seq[k]>>6) {
			q = int(seq[k-1]&0x3f) + int(seq[k]&0x3f) + 10;
		} else if (tarray[k-1] == seq[k-1]>>6) {
			q = int(seq[k-1]&0x3f) - int(seq[k]&0x3f);
		} else if (tarray[k] == seq[k]>>6) {
			q = int(seq[k]&0x3f) - int(seq[k-1]&0x3f);
		} // else, q = 0
		if (q < 1) q = 1;
		if (q > 63) q = 63;
		t2array[k] = nt_read[k]<<6 | q;
		if (seq[k-1] == 0 && seq[k] == 0) t2array[k] = 0;
	}
	memcpy(seq, t2array+1, size-1); // t2array[0] and t2array[size] are not copied
	seq[size-1] = 0;
}

static void csmap2nt_core(gzFile fpout, FILE *fpbfa, gzFile fpmap)
{
	nst_bfa1_t *l = 0;
	bit32_t seqid = 0;
	int k;
	hash_map_char<bit32_t> *hash_map = new hash_map_char<bit32_t>;
	maqmap_t *mm = maqmap_read_header(fpmap); // skip the sequence names
	maqmap1_t mm1, *m1;
	bit8_t nt_ref[MAX_READLEN+1], nt_read[MAX_READLEN+1], tarray[(MAX_READLEN+1) * 4];
	m1 = &mm1;
	for (k = 0; k != mm->n_ref; ++k) hash_map->insert(mm->ref_name[k], k);
	maqmap_write_header(fpout, mm);
	maqmap_read1(fpmap, m1);
	while ((l = nst_load_bfa1(fpbfa)) != 0) {
		if (!hash_map->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
	    do {
			if (m1->seqid != seqid) break;
			if (get_readseq(l, m1->pos>>1, m1->size, nt_ref)) continue;
			cs2nt_dp(m1->size, nt_ref, m1->seq, nt_read, tarray);
			cal_nt_qual(m1->size, nt_read, m1->seq, tarray);
			m1->pos += 2;
			--m1->size;
			gzwrite(fpout, m1, sizeof(maqmap1_t));
		} while (maqmap_read1(fpmap, m1));
		nst_delete_bfa1(l);
		if (gzeof(fpmap)) break; // end of alignment file
	}
	maq_delete_maqmap(mm);
	delete hash_map;
}

int maq_csmap2nt(int argc, char *argv[])
{
	gzFile fpin, fpout;
	FILE *fp_bfa;
	if (argc < 4) {
		fprintf(stderr, "Usage: maq csmap2nt <out.nt.map> <in.ref.nt.bfa> <in.cs.map>\n");
		return 1;
	}
	fpout = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdout), "w") : gzopen(argv[1], "w");
	fp_bfa = fopen(argv[2], "r");
	fpin  = (strcmp(argv[3], "-") == 0)? gzdopen(fileno(stdin), "r")  : gzopen(argv[3], "r");
	assert(fpout && fpin && fp_bfa);
	csmap2nt_core(fpout, fp_bfa, fpin);
	gzclose(fpin); gzclose(fpout); fclose(fp_bfa);
	return 0;
}
