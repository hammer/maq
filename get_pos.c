#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include "maqmap.h"
#include "assemble.h"

/** fill the rolling buffer */
inline void assemble_fill_buffer(gzFile fpin, rolling_buf_t *ab)
{
	if (gzeof(fpin) || (ab->is_rounded == 1 && ab->head == ab->tail)) return; // the buffer is full
	int n_records, n_bytes;
	if (ab->tail >= ab->head) {
		n_records = ROLLING_BUF_SIZE - ab->tail;
		n_bytes = gzread(fpin, ab->buf + ab->tail, sizeof(maqmap1_t) * n_records);
		ab->tail += n_bytes / sizeof(maqmap1_t);
		if ((bit32_t)n_bytes < sizeof(maqmap1_t) * n_records) return;
		ab->is_rounded = 1; ab->tail = 0;
	}
	if (ab->head == 0) return;
	n_records = ab->head - ab->tail;
	n_bytes = gzread(fpin, ab->buf + ab->tail, sizeof(maqmap1_t) * n_records);
	ab->tail += n_bytes / sizeof(maqmap1_t);
}
/** get the piled bases at a specified position */
void assemble_get_pos(bit32_t seqid, bit32_t pos, gzFile fpin, rolling_buf_t *ab, assemble_pos_t *ap,
					  int max_mm, int max_err, int min_q, int is_single, int is_pair_only)
{
	int i, is_rounded;
	if (max_err < 0) max_err = 60; // default value is set here.
	if (max_mm < 0) max_mm = 7; // default value
	ap->n_bases = 0;
	do {
		assemble_fill_buffer(fpin, ab);
		if (gzeof(fpin) && ab->is_rounded == 0 && ab->head == ab->tail) return; // end of file and empty buffer
		// skip reads that are out of range.
		for (i = ab->head, is_rounded = ab->is_rounded; (is_rounded || i != ab->tail)
				 && (seqid > ab->buf[i].seqid || (seqid == ab->buf[i].seqid && pos >= (ab->buf[i].pos>>1) + ab->buf[i].size));)
		{
			if ((++i) == ROLLING_BUF_SIZE) { i = 0; is_rounded = 0; }
		}
		ab->head = i; ab->is_rounded = is_rounded; // update ab
	} while (ab->is_rounded == 0 && ab->head == ab->tail);
	assemble_fill_buffer(fpin, ab);
	is_rounded = ab->is_rounded;
	// fill ap
	while ((is_rounded || i != ab->tail) && seqid == ab->buf[i].seqid && pos >= ab->buf[i].pos>>1) {
		maqmap1_t *m1 = ab->buf + i;
		bit32_t map_qual = is_single? m1->seq[MAX_READLEN-1] : m1->map_qual;
		int is_pair_rm = (is_pair_only && (m1->flag != 0 && m1->flag != (PAIRFLAG_FR|PAIRFLAG_PAIRED)))? 1 : 0;
		if (m1->info2 <= max_err && (m1->info1&0xf) <= max_mm && map_qual >= (bit32_t)min_q && pos < (m1->pos>>1) + m1->size
			&& (m1->flag&PAIRFLAG_SW) == 0 && !is_pair_rm)
		{
			if (ap->n_bases == ap->m_bases) {
				ap->m_bases += 0x1000;
				ap->bases = (assemble_posinfo_t*)realloc(ap->bases, ap->m_bases * sizeof(assemble_posinfo_t));
			}
			assemble_posinfo_t *p = ap->bases + ap->n_bases;
			bit32_t base_qual = m1->seq[pos - (m1->pos>>1)]&0x3f;
			bit32_t base = m1->seq[pos - (m1->pos>>1)]>>6&0x3;
			bit32_t qual = (base_qual < map_qual)? base_qual : map_qual;
			bit32_t is_present = (m1->seq[pos - (m1->pos>>1)] == 0)? 0 : 1;
			bit32_t mm = m1->info1>>4&0xf;
			if (qual < 1) qual = 1;
			bit32_t strand = m1->pos&1;
			p->info = (qual<<24) | (is_present<<22) | (mm<<19) | (strand<<18) | (base<<16) | (base_qual<<8) | map_qual;
			p->c[0] = m1->c[0]; p->c[1] = m1->c[1];
			p->pos = (m1->pos&1)? m1->size - 1 - (pos - (m1->pos>>1)) : pos - (m1->pos>>1);
			++ap->n_bases;
		}
		if ((++i) == ROLLING_BUF_SIZE) { i = 0; is_rounded = 0; }
	}
}

void assemble_get_indelpos(bit32_t seqid, bit32_t pos, gzFile fpin, rolling_buf_t *ab, assemble_indelpos_t *ai)
{
	int i, is_rounded;

	// initialization
	ai->n_reads = ai->n_types = ai->n_ungap = 0;
	memset(ai->ins_bases[0], 0, sizeof(int) * 4 * MAX_READLEN);

	do {
		assemble_fill_buffer(fpin, ab);
		if (gzeof(fpin) && ab->is_rounded == 0 && ab->head == ab->tail) return; // end of file and empty buffer
		// skip reads that are out of range.
		for (i = ab->head, is_rounded = ab->is_rounded; (is_rounded || i != ab->tail)
				 && (seqid > ab->buf[i].seqid || (seqid == ab->buf[i].seqid && pos >= (ab->buf[i].pos>>1) + ab->buf[i].size));)
		{
			if ((++i) == ROLLING_BUF_SIZE) { i = 0; is_rounded = 0; }
		}
		ab->head = i; ab->is_rounded = is_rounded; // update ab
	} while (ab->is_rounded == 0 && ab->head == ab->tail);
	assemble_fill_buffer(fpin, ab);
	is_rounded = ab->is_rounded;

	// fill ai
	while ((is_rounded || i != ab->tail) && seqid == ab->buf[i].seqid && pos >= ab->buf[i].pos>>1) {
		maqmap1_t *m1 = ab->buf + i;
		if (pos < (m1->pos>>1) + m1->size) {
			++ai->n_reads;
			if ((m1->flag & PAIRFLAG_SW) && (signed char)m1->seq[MAX_READLEN-1] != 0) { // indel
				if (pos == (m1->pos>>1) + m1->map_qual) {
					int j;
					bit8_t b = (bit8_t)m1->seq[MAX_READLEN-1]; // force to be positive
					for (j = 0; j != ai->n_types; ++j)
						if (b == (ai->indels[j]&0xff)) break;
					assert(j < 256);
					if (j < ai->n_types) { // exist
						if (m1->pos&1) {
							ai->indels[j] += 1ull<<48;
							ai->indels[j] += 1ull<<16;
						} else {
							ai->indels[j] += 1ull<<48;
							ai->indels[j] += 1ull<<32;
						}
					} else { // new type
						ai->indels[j] = 1ull<<48 | b | ((m1->pos&1)? 1ull<<16 : 1ull<<32);
						++ai->n_types;
					}
					if ((signed char)b > 0) { // insertion, then get inserted bases
						for (j = 0; j < b; ++j)
							++ai->ins_bases[j][m1->seq[j+m1->map_qual]>>6&3];
					}
				}
			} else { // non-indel
				if (pos > (m1->pos>>1) && ai->n_ungap < 255) {
					ai->n_mm[ai->n_ungap] = m1->info1&0xf;
					ai->beg_pos[ai->n_ungap] = pos - (m1->pos>>1);
					ai->end_pos[ai->n_ungap++] = m1->size - (pos - (m1->pos>>1));
				}
			}
		}
		if ((++i) == ROLLING_BUF_SIZE) { i = 0; is_rounded = 0; }
	}
}
