#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include "match.hh"
#include "bfa.h"
#include "assemble.h"
#include "stdhash.hh"
#include "main.h"

/** main function for 'pileup' */
static int pileup_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   maq pileup [options] <chr.bfa> <align.map>\n\n");
	fprintf(stderr, "Options: -Q INT    maximum sum of errors [60]\n");
	fprintf(stderr, "         -m INT    maximum number of mismatches [7]\n");
	fprintf(stderr, "         -q INT    minimum mapping quality [0]\n");
	fprintf(stderr, "         -l FILE   only output required positions [null]\n");
	fprintf(stderr, "         -s        use single-end mapping qualities\n");
	fprintf(stderr, "         -p        discard abnormal pairs\n");
	fprintf(stderr, "         -d        only show depth\n");
	fprintf(stderr, "         -v        verbose mode\n");
	fprintf(stderr, "         -P        print position on the read\n\n");
	return 1;
}
int ma_pileup(int argc, char *argv[])
{
	extern hash_set_char *ma_load_snp_set(FILE *fp_list);
	int k, c, is_verbose, max_err, max_mm, seqid, min_q, is_single = 0, is_pair_only = 0, is_depth_only, is_show_pos;
	FILE *fp_bfa, *fpout, *fp_list;
	gzFile fp_map;
	nst_bfa1_t *l;
	is_verbose = is_show_pos = is_depth_only = 0; max_mm = max_err = -1; fp_list = 0; min_q = 0;
	while ((c = getopt(argc, argv, "Q:l:q:vsdpm:P")) >= 0) {
		switch (c) {
		case 'd': is_depth_only = 1; break;
		case 'p': is_pair_only = 1; break;
		case 'q': min_q = atoi(optarg); break;
		case 'v': is_verbose = 1; break;
		case 'P': is_show_pos = 1; break;
		case 'l': fp_list = fopen(optarg, "r"); assert(fp_list); break;
		case 'Q': max_err = atoi(optarg); break;
		case 's': is_single = 1; break;
		case 'm': max_mm = atoi(optarg); break;
		default: return pileup_usage();
		}
	}
	if (argc - optind < 2) return pileup_usage();
	fpout = stdout;
	fp_bfa = fopen(argv[optind], "r");
	fp_map = (strcmp(argv[optind+1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind+1], "r");
	assert(fp_bfa && fp_map);
	hash_set_char *hash = 0;
	if (fp_list) {
		hash = ma_load_snp_set(fp_list);
		fclose(fp_list);
	}
	char key[256];
	rolling_buf_t *buf = (rolling_buf_t*)calloc(1, sizeof(rolling_buf_t));
	buf->buf = (maqmap1_t*)calloc(ROLLING_BUF_SIZE, sizeof(maqmap1_t));
	assemble_pos_t *pos = (assemble_pos_t*)calloc(1, sizeof(assemble_pos_t));
	maqmap_t *mm = maqmap_read_header(fp_map); // skip the sequence names
	hash_map_char<int> *hash_map = new hash_map_char<int>;
	for (k = 0; k != mm->n_ref; ++k) hash_map->insert(mm->ref_name[k], k);
	k = mm->n_mapped_reads;
	while ((l = nst_load_bfa1(fp_bfa)) != 0) {
		if (!hash_map->find(l->name, &seqid)) {
			nst_delete_bfa1(l);
			continue;
		}
		for (int i = 0; i != l->len; ++i) {
			bit64_t word = l->seq[i];
			bit64_t mask = l->mask[i];
			for (int j = 31; j >= 0; --j) {
				int coor = (i<<5) | (31 - j);
				if (coor >= l->ori_len) break;
				assemble_get_pos(seqid, coor, fp_map, buf, pos, max_mm, max_err, min_q, is_single, is_pair_only);
				if (hash) {
					sprintf(key, "%s.%d", l->name, coor + 1);
					if (!hash->find(key)) continue;
				}
				fprintf(fpout, "%s\t%d\t%c\t%d", l->name, coor + 1, (mask>>(j<<1)&3)? "ACGT"[word>>(j<<1)&3] : 'N',
						pos->n_bases);
				if (is_depth_only) {
					fputc('\n', fpout);
					continue;
				}
				fprintf(fpout, "\t@");
				for (k = 0; k != pos->n_bases; ++k) {
					assemble_posinfo_t *p = pos->bases + k;
					if ((p->info>>22&1) && (word>>(j<<1)&3) == (p->info>>16&3)) c = (p->info>>18&1)? '.' : ',';
					else if ((p->info>>22&1) == 0) c = (p->info>>18&1)? 'n' : 'N';
					else c = ((p->info>>18&1)? "acgt" : "ACGT")[p->info>>16&3];
					fputc(c, fpout);
				}
				if (is_verbose) {
					fprintf(fpout, "\t@");
					for (k = 0; k != pos->n_bases; ++k)
						fputc((pos->bases[k].info>>8&0x3f) + 33, fpout);
					fprintf(fpout, "\t@");
					for (k = 0; k != pos->n_bases; ++k) {
						c = (pos->bases[k].info&0x7f) + 33;
						fputc((c > 126)? 126 : c, fpout);
					}
				}
				if (is_show_pos) {
					fputc('\t', fpout);
					for (k = 0; k != pos->n_bases; ++k)
						fprintf(fpout, "%d,", pos->bases[k].pos+1);
				}
				fputc('\n', fpout);
			}
		}
		nst_delete_bfa1(l);
	}
	gzclose(fp_map); fclose(fp_bfa);
	maq_delete_maqmap(mm);
	delete hash_map;
	if (hash) delete hash;
	free(buf->buf); free(buf);
	free(pos->bases); free(pos);
	return 0;
}
