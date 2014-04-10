/*
 *  19th Aug 2008, Modified by Colin Hercus to include conversion of novoalign and novopaired reports to maq format.
 *  Released under GPL
 */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include "maqmap.h"
#include "algo.hh"
#include "stdhash.hh"
#include "main.h"

static int DEFAULT_QUAL = 25;
static int log_n[256];
static char **name_conv;

static hash_map_char<bit32_t> *read_list(FILE *fp)
{
	hash_map_char<bit32_t> *hash = new hash_map_char<bit32_t>;
	char str[1024];
	int i = 0;
	while (fscanf(fp, "%s", str) != 0) {
		if (feof(fp)) break;
		char *p = (str[0] == '>')? str+1 : str;
		int c;
		bit32_t x;
		assert(!hash->find(p, &x));
		hash->insert(p, i);
		fscanf(fp, "%s", str);
		if ((i&0xff) == 0)
			name_conv = (char**)realloc(name_conv, sizeof(char*) * (i + 0x100));
		name_conv[i++] = strdup(p);
		while ((c = fgetc(fp)) != EOF && c != '\n');
	}
	return hash;
}
static inline int cal_map_qual(int default_qual, bit32_t *count)
{
	if (count[0] == 1) {
		if (count[1] == 0 && count[2] == 0) return 3 * default_qual;
		if (count[1] == 0) return 2 * default_qual - log_n[count[2]];
		return default_qual - log_n[count[1]];
	}
	if (count[1] == 1) {
		if (count[2] == 0) return 2 * default_qual;
		return default_qual - log_n[count[2]];
	}
	if (count[2] == 1) return default_qual - 3;
	return default_qual;
}
static inline int operator < (const maqmap1_t &a, const maqmap1_t &b)
{
	return (a.seqid < b.seqid) || (a.seqid == b.seqid && a.pos < b.pos);
}

static void eland2maq_core(FILE *fp_list, FILE *fp_eland, gzFile fp)
{
	hash_map_char<bit32_t> *hash;
	// initialize maqmap_t
	maqmap_t *mm = maq_new_maqmap();
	int max = 0, i, l;
	hash = read_list(fp_list);
	mm->n_ref = hash->size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	for (i = 0; i != int(hash->size()); ++i)
		mm->ref_name[i] = strdup(name_conv[i]);
	// initialize log_n
	log_n[0] = -1;
	for (i = 1; i != 256; ++i)
		log_n[i] = (int)(3.434 * log((float)i) + 0.5);
	// read the file
	char str[1024], str2[1024];
	bit8_t tmp_seq[MAX_READLEN];
	bit32_t tmp[4];
	while (fscanf(fp_eland, "%s", str) != 0) {
		if (feof(fp_eland)) break;
		if (mm->n_mapped_reads == (bit64_t)max)
			mm->mapped_reads = (maqmap1_t*)realloc(mm->mapped_reads, sizeof(maqmap1_t) * (max += 0x100000));
		maqmap1_t *m1 = mm->mapped_reads + mm->n_mapped_reads;
		// set name
		strncpy(m1->name, str+1, MAX_NAMELEN-1);
		m1->name[MAX_NAMELEN - 1] = 0;
		// set seq
		fscanf(fp_eland, "%s", str);
		memset(m1->seq, 0, MAX_READLEN);
		m1->size = l = strlen(str);
		for (i = 0; i != l; ++i) {
			tmp[0] = nst_nt4_table[(int)str[i]];
			m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | DEFAULT_QUAL);
		}
		// read flag
		fscanf(fp_eland, "%s", str);
		if (str[0] == 'U') { // there is an alignment
			int n_mis = str[1] - '0';
			fscanf(fp_eland, "%u%u%u%s%u%s", tmp, tmp+1, tmp+2, str, tmp+3, str2);
			if (hash->find(str, &m1->seqid)) {
				if (str2[0] == 'R') {
					for (i = m1->size - 1; i >= 0; --i)
						tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
					memcpy(m1->seq, tmp_seq, m1->size);
				}
				m1->c[0] = tmp[0];
				m1->c[1] = tmp[1];
				m1->flag = 0;
				m1->dist = 0;
				m1->pos = (tmp[3]-1)<<1 | (str2[0] == 'F'? 0 : 1);
				m1->info1 = n_mis<<4 | n_mis;
				m1->info2 = n_mis * DEFAULT_QUAL;
				m1->map_qual = m1->seq[MAX_READLEN-1] = m1->alt_qual = cal_map_qual(DEFAULT_QUAL, tmp);
			}
			++mm->n_mapped_reads;
		}
		while ((i = fgetc(fp_eland)) != EOF && i != '\n');
	}
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(fp, mm);
	gzwrite(fp, mm->mapped_reads, sizeof(maqmap1_t) * mm->n_mapped_reads);
	// free
	maq_delete_maqmap(mm);
	for (i = 0; i != int(hash->size()); ++i) free(name_conv[i]);
	free(name_conv);
	delete hash;
}

static int usage()
{
	fprintf(stderr, "Usage: maq eland2maq [-q qual] <out.map> <in.list> <in.eland>\n");
	return 1;
}

int maq_eland2maq(int argc, char *argv[])
{
	FILE *fp_eland, *fp_list;
	gzFile fp_map;
	int c;
	while ((c = getopt(argc, argv, "q:")) >= 0) {
		switch (c) {
		case 'q': DEFAULT_QUAL = atoi(optarg); break;
		}
	}
	if (optind + 2 >= argc) return usage();
	fp_eland = (strcmp(argv[optind+2], "-") == 0)? stdin : fopen(argv[optind+2], "r");
	fp_list = fopen(argv[optind+1], "r");
	fp_map = gzopen(argv[optind], "w");
	assert(fp_eland && fp_list && fp_map);
	eland2maq_core(fp_list, fp_eland, fp_map);
	fclose(fp_eland); fclose(fp_list); gzclose(fp_map);
	return 0;
}

/* novo2maq 
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  info1: mismatches in the 24bp (higher 4 bits) and mismatches (lower 4 bits)
  info2: sum of errors of the best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
// Recover novo indels while making sure to only pick up indels smith-waterman would find.
// Needleman-Wunsch used in novo may report indels too close to the ends of the reads.
static bool gapped(char * mm, char strand, int rdlen, int &readposn, int &indellen) {
    if(*mm == '\0')
        return false;
    int tposn = 1, qposn = 1, posn = 1;
    int match = 15, mismatch = 30, gapopen = 40, gapextend = 10;
    int maxsc = 0, maxposn = 0;
    bool gapped = false;
    int sc[256], mv[256], qp[256], tp[256];
    sc[0] = qp[0] = tp[0] = '\0';
    mv[0] = '|';
    while(*mm != '\0') {
        int i = 0;
        while(isdigit(*mm))
            i = i* 10 + (*mm++ - '0');
        while(tposn < i) {
            sc[posn] = sc[posn-1] + match;
            mv[posn] = '|';
            qp[posn] = qposn;
            tp[posn] = tposn;
            tposn++;
            qposn++;
            posn++;
        }
        if(sc[tposn-1] > maxsc) {
            maxsc = sc[tposn - 1];
            maxposn = tposn - 1;
        }
        switch(*mm) {
            case '+':
                gapped = true;
                if(mv[posn-1] == '+')
                    sc[posn] = sc[posn-1] - gapextend;
                else
                    sc[posn] = sc[posn-1] - gapopen;
                if(sc[posn] < 0)
                    sc[posn] = 0;
                mv[posn] = '+';
                qp[posn] = qposn;
                tp[posn] = tposn;
                qposn++;
                posn++;
                mm +=3;
                break;
            case '-':
                gapped = true;
                if(mv[posn-1] == '-')
                    sc[posn] = sc[posn-1] - gapextend;
                else
                    sc[posn] = sc[posn-1] - gapopen;
                if(sc[posn] < 0)
                    sc[posn] = 0;
                mv[posn] = '-';
                qp[posn] = qposn;
                tp[posn] = tposn;
                tposn++;
                posn++;
                mm += 3;
                break;
            default:
                sc[posn] = sc[posn-1] - mismatch;
                if(sc[posn] < 0)
                    sc[posn] = 0;
                mv[posn] = '|';
                qp[posn] = qposn;
                tposn++;
                qposn++;
                posn++;
                mm+= 4;
        }
    }
    if(!gapped)
        return false;
    while(qposn <= rdlen) {
        sc[posn] = sc[posn-1] + match;
        mv[posn] = '|';
        qp[posn] = qposn;
        tp[posn] = tposn;
        tposn++;
        qposn++;
        posn++;
    }
    if(sc[tposn-1] > maxsc) {
        maxsc = sc[tposn - 1];
        maxposn = tposn - 1;
    }
    if(maxsc <= 0)
        return false;
    int j = maxposn;
    while(sc[j] > 0 && mv[j] == '|')
        j--;
    if(mv[j] == '|')
        return false;
    indellen = 1;
    while(mv[j-1] == mv[j]) {
        indellen++;
        j--;
    }
//    if(strand == 'F')
        readposn = qp[j] - 1;
//    else
//        readposn = rdlen - qp[j + indellen - 1] + 1;
    if(mv[j] == '-')
        indellen = -indellen;
    return true;
}

#define min(x,y) ((x)<(y)? (x):(y))

static void novo2maq_core(FILE *fp_list, FILE *fp_novo, gzFile fp)
{
	hash_map_char<bit32_t> *hash;
	// initialize maqmap_t
	maqmap_t *mm = maq_new_maqmap();
	int max = 0, i, l;
	hash = read_list(fp_list);
	mm->n_ref = hash->size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	for (i = 0; i != int(hash->size()); ++i)
		mm->ref_name[i] = strdup(name_conv[i]);
	// initialize log_n
	log_n[0] = -1;
	for (i = 1; i != 256; ++i)
		log_n[i] = (int)(3.434 * log((float)i) + 0.5);
	// read the file
    enum {n_readid, n_side, n_seq, n_qual, n_status, n_score, n_quality, n_chr, n_offset, n_strand, n_chr2, n_offset2, n_strand2, n_poly};
    char * novo[n_poly + 1];
    char * novo2[n_poly + 1];
	char str[4096], str2[4096];
	bit8_t tmp_seq[MAX_READLEN];
	bit32_t tmp[4];
	while (fgets( str, 4096,fp_novo ) != NULL) {
		if (feof(fp_novo)) break;
        if(str[0] == '#' || str[0] == '\0')
            continue;
		if (mm->n_mapped_reads + 1 >= (bit64_t)max)      // Make sure room for a pair.
			mm->mapped_reads = (maqmap1_t*)realloc(mm->mapped_reads, sizeof(maqmap1_t) * (max += 0x100000));
		maqmap1_t *m1 = mm->mapped_reads + mm->n_mapped_reads;
        char *np = str;
        char mtstr[1] = "";
        for(i = 0; i <= n_poly; i++ ) {
            novo[i] = np;
            while(*np != '\t' && *np != '\n' && *np != '\0')
                np++;
            if(*np == '\n' || *np == '\0') break;
            *np++ = '\0';
        }
        *np = '\0';
        for(i++;i <= n_poly; i++ ) novo[i] = mtstr;
        
        // Read pair
        bool paired =false;
        if(novo[n_side][0] != 'S') {
            paired = true;
            fgets( str2, 4096, fp_novo );
            np = str2;
            for(i = 0; i <= n_poly; i++ ) {
                novo2[i] = np;
                while(*np != '\t' && *np != '\n')
                    np++;
                if(*np == '\n') break;
                *np++ = '\0';
            }
            *np = '\0';
            for(i++;i <= n_poly; i++ ) novo2[i] = mtstr;
        }
        // read flag
		if (novo[n_status][0] == 'U' || (paired && novo2[n_status][0] == 'U')) { // there is an alignment
            // set name
            strncpy(m1->name, novo[n_readid]+1, MAX_NAMELEN-1);
            m1->name[MAX_NAMELEN - 1] = 0;
            // set seq
            m1->size = l = strlen(novo[n_seq]);
            if(novo[n_qual][0] == '\0')
                for (i = 0; i != l; ++i) {
                    tmp[0] = nst_nt4_table[(int)novo[n_seq][i]];
                    m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | DEFAULT_QUAL);
                }
            else
                for (i = 0; i != l; ++i) {
                    tmp[0] = nst_nt4_table[(int)novo[n_seq][i]];
                    m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | ((int)novo[n_qual][i] - 33));
                }
			int n_mis = (atoi(novo[n_score]) + 15)/30;
            if (novo[n_status][0] != 'U')
                n_mis = 0;
            bool aligned = hash->find(novo[n_chr]+1, &m1->seqid);
            if (novo[n_strand][0] == 'R') {
                for (i = m1->size - 1; i >= 0; --i)
                    tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
                memcpy(m1->seq, tmp_seq, m1->size);
            }

            m1->c[0] = 0;   // 0 mismatches alignments
            m1->c[1] = 0;   // 1 mismatch alignments
            if(n_mis == 0)
                m1->c[0] = 1;   // 0 mismatches alignments
            else
                m1->c[1] = 1;   // 1 mismatch alignments

            m1->flag = 0;   //paired flag  FR, correct pair, 16 good pair, 32 diff chr, 64 one read of pair not mapped, 128 this read not mapped and other read is.
            m1->dist = 0;   // Outer coordinates
            m1->pos = (atoi(novo[n_offset])-1)<<1 | (novo[n_strand][0] == 'F'? 0 : 1);
            m1->info1 = n_mis<<4 | n_mis;
            m1->info2 = atoi(novo[n_score])/3;
            m1->map_qual = m1->seq[MAX_READLEN-1] = min(60, atoi(novo[n_quality]));
            m1->alt_qual = m1->map_qual;  // Lesser of two reads
            if(paired) {
                int dist = 0;
                int flag = 64;
                if (novo[n_status][0] != 'U') {
                    m1->flag = 192;
                    m1->info2 = 0;
                    m1->map_qual = m1->seq[MAX_READLEN-1] = 0;
                    m1->alt_qual = min(60, atoi(novo2[n_quality]));
                } else if (novo2[n_status][0] != 'U') {
                    m1->flag = 64;
                    flag = 192;
                    m1->alt_qual = min((int)m1->map_qual, atoi(novo2[n_quality]));
                } else {
                    if(novo[n_chr2][0] == '\0' 
                            && strcmp(novo[n_quality], novo2[n_quality]) == 0
                            && (dist = abs( atoi(novo[n_offset]) - atoi(novo[n_offset2]))) < 2000) {
                        if(novo[n_strand][0] == 'R')
                            dist = -(dist + strlen(novo[n_seq]));
                        else
                            dist += strlen(novo2[n_seq]);
                        m1->flag = flag = 2 + 16;
                        m1->dist = dist;
                    } else if(novo[n_chr2][0] != '\0')
                        m1->flag = flag = 32;
                    else if(novo[n_strand][0] == 'R' && novo2[n_strand][0] == 'R')
                        m1->flag = flag = 8;
                    else if(novo[n_strand][0] == 'F' && novo2[n_strand][0] == 'F')
                        m1->flag = flag = 1;
                    else
                        m1->flag = flag = 8; 
                    int gapposn, gaplen;
                    if(gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen)) {
                        //gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen);
                        m1->flag = 130;
                        m1->map_qual = gapposn;
                        m1->seq[MAX_READLEN-1] = gaplen;
                    }
                }
                if(aligned && m1->flag != 192)
                    ++mm->n_mapped_reads;
                if(flag == 192)
                    continue;
                m1 = mm->mapped_reads + mm->n_mapped_reads;
                            // set name
                strncpy(m1->name, novo2[n_readid]+1, MAX_NAMELEN-1);
                m1->name[MAX_NAMELEN - 1] = 0;
                // set seq
                m1->size = l = strlen(novo2[n_seq]);
                if(novo2[n_qual][0] == '\0')
                    for (i = 0; i != l; ++i) {
                        tmp[0] = nst_nt4_table[(int)novo2[n_seq][i]];
                        m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | DEFAULT_QUAL);
                    }
                else
                    for (i = 0; i != l; ++i) {
                        tmp[0] = nst_nt4_table[(int)novo2[n_seq][i]];
                        m1->seq[i] = (tmp[0] > 3)? 0 : (tmp[0]<<6 | ((int)novo2[n_qual][i] - 33));
                    }
                int n_mis = (atoi(novo2[n_score]) + 15)/30;
                if (novo2[n_status][0] != 'U')
                    n_mis = 0;
                if (!hash->find(novo2[n_chr]+1, &m1->seqid)) {
                    continue;                                        
                }
                aligned = true;
                if (novo2[n_strand][0] == 'R') {
                    for (i = m1->size - 1; i >= 0; --i)
                        tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
                    memcpy(m1->seq, tmp_seq, m1->size);
                }

                m1->c[0] = 0;   // 0 mismatches alignments
                m1->c[1] = 0;   // 1 mismatch alignments
                if(n_mis == 0)
                    m1->c[0] = 1;   // 0 mismatches alignments
                else
                    m1->c[1] = 1;   // 1 mismatch alignments

                m1->flag = flag;   // paired flag  FR, correct pair, 16 good pair, 32 diff chr, 64 one read of pair not mapped, 128 this read not mapped and other read is.
                m1->dist = -dist;   // Outer coordinates
                m1->pos = (atoi(novo2[n_offset])-1)<<1 | (novo2[n_strand][0] == 'F'? 0 : 1);
                m1->info1 = n_mis<<4 | n_mis;
                m1->info2 = atoi(novo2[n_score])/3;
                m1->map_qual = m1->seq[MAX_READLEN-1] = min(60, atoi(novo2[n_quality]));
                if (novo[n_status][0] == 'U')
                    m1->alt_qual = min(60,min((int)m1->map_qual, atoi(novo[n_quality]))); 
                else
                    m1->alt_qual = m1->map_qual;
                int gapposn, gaplen;
                if(gapped(novo2[n_poly], novo2[n_strand][0], m1->size, gapposn, gaplen)) {
                    //gapped(novo2[n_poly], novo2[n_strand][0], m1->size, gapposn, gaplen);
                    m1->flag = 130;
                    m1->map_qual = gapposn;
                    m1->seq[MAX_READLEN-1] = gaplen;
                }
            } else {
                int gapposn, gaplen;
                if(gapped(novo[n_poly], novo[n_strand][0], m1->size, gapposn, gaplen)) {
                    m1->flag = 130;
//                    m1->c[0] = 0;   // 0 mismatches alignments
//                    m1->c[1] = 0;   // 1 mismatch alignments
                    m1->alt_qual = m1->map_qual;
                    m1->map_qual = gapposn;
                    m1->seq[MAX_READLEN-1] = gaplen;
                }
            }
			if(aligned)
                ++mm->n_mapped_reads;
		}
	}
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(fp, mm);
	gzwrite(fp, mm->mapped_reads, sizeof(maqmap1_t) * mm->n_mapped_reads);
	// free
	maq_delete_maqmap(mm);
	for (i = 0; i != int(hash->size()); ++i) free(name_conv[i]);
	free(name_conv);
	delete hash;
}

static int novo2maqusage()
{
	fprintf(stderr, "Usage: maq novo2maq <out.map> <in.list> <in.novo>\n");
	return 1;
}

int maq_novo2maq(int argc, char *argv[])
{
	FILE *fp_novo, *fp_list;
	gzFile fp_map;
//	int c;
//	while ((c = getopt(argc, argv, "")) >= 0) {
//		switch (c) {
//		case 'q': DEFAULT_QUAL = atoi(optarg); break;
//		}
//	}
    optind = 1;
	if (optind + 2 >= argc) return novo2maqusage();
	fp_novo = (strcmp(argv[optind+2], "-") == 0)? stdin : fopen(argv[optind+2], "r");
	fp_list = fopen(argv[optind+1], "r");
	fp_map = gzopen(argv[optind], "w");
	assert(fp_novo && fp_list && fp_map);
	novo2maq_core(fp_list, fp_novo, fp_map);
	fclose(fp_novo); fclose(fp_list); gzclose(fp_map);
	return 0;
}

/* export2maq */

static int read_len[2], max_dist = 250, is_non_filtered = 0;

static inline int tabread(FILE *fp, char buf[])
{
	int c;
	char *p = buf;
	while ((c = fgetc(fp)) != '\n' && c != '\t' && c != EOF)
		*p++ = c;
	*p = '\0';
	return c;
}
static void export2maq_core(FILE *fp_list, FILE *fp_export, gzFile fp)
{
	hash_map_char<bit32_t> *hash;
	// quality conversion table
	int table[128];
	for (int l = 0; l != 128; ++l) {
		table[l] = (int)(10.0 * log(1.0 + pow(10.0, (l - 64) / 10.0)) / log(10.0) + .499);
		if (table[l] >= 63) table[l] = 63;
		if (table[l] == 0) table[l] = 1;
	}
	// initialize maqmap_t
	maqmap_t *mm = maq_new_maqmap();
	int max = 0, i, l;
	hash = read_list(fp_list);
	mm->n_ref = hash->size();
	mm->ref_name = (char**)malloc(sizeof(char*) * mm->n_ref);
	for (i = 0; i != int(hash->size()); ++i)
		mm->ref_name[i] = strdup(name_conv[i]);
	// read the file
	char str[1024], str2[1024];
	bit8_t tmp_seq[MAX_READLEN];
	int ti[5], c;
	while (fscanf(fp_export, "%s%d%d%d%d%d", str, ti, ti+1, ti+2, ti+3, ti+4) != 0) {
		if (feof(fp_export)) break;
		if (mm->n_mapped_reads == (bit64_t)max)
			mm->mapped_reads = (maqmap1_t*)realloc(mm->mapped_reads, sizeof(maqmap1_t) * (max += 0x100000));
		maqmap1_t *m1 = mm->mapped_reads + mm->n_mapped_reads;
		sprintf(str2, "%s_%d:%d:%d:%d:%d", str, ti[0], ti[1], ti[2], ti[3], ti[4]);
		// set name
		strncpy(m1->name, str2, MAX_NAMELEN-1);
		m1->name[MAX_NAMELEN - 1] = 0;
		// set m1->c[2]
		m1->c[0] = m1->c[1] = 0;
		// set seq and qual
		fgetc(fp_export); // skip '\t'
		c = tabread(fp_export, str); // index
		c = tabread(fp_export, str); // 1 or 2
		int cur_read = atoi(str) - 1;
		c = tabread(fp_export, str); // seq
		c = tabread(fp_export, str2); // qual
		memset(m1->seq, 0, MAX_READLEN);
		m1->size = l = strlen(str);
		// update read_len
		read_len[cur_read] = m1->size;
		if (read_len[1-cur_read] == 0) read_len[1-cur_read] = m1->size;
		for (i = 0; i != l; ++i) {
			int t = nst_nt4_table[(int)str[i]];
			m1->seq[i] = (t > 3)? 0 : (t<<6 | table[(int)str2[i]]);
		}
		c = tabread(fp_export, str); // target
		if (strcmp(str, "NM") && hash->find(str, &m1->seqid)) {
			c = tabread(fp_export, str); // contig
			c = tabread(fp_export, str); // position
			c = tabread(fp_export, str2); // strand
			m1->pos = bit32_t(atoi(str) - 1)<<1 | (str2[0] == 'F'? 0 : 1);
			if (m1->pos&1) { // reverse if necessary
				for (i = m1->size - 1; i >= 0; --i)
					tmp_seq[m1->size-i-1] = (m1->seq[i] == 0)? 0 : (0xc0 - (m1->seq[i]&0xc0)) | (m1->seq[i]&0x3f);
				memcpy(m1->seq, tmp_seq, m1->size);
			}
			// set m1->info1 and m1->info2
			c = tabread(fp_export, str); // mismatches
			int n_mismatch = 0;
			l = strlen(str);
			for (int k = 0; k != l; ++k)
				if (isalpha(str[k])) ++n_mismatch;
			m1->info1 = n_mismatch << 4 | n_mismatch;
			m1->info2 = 10 * n_mismatch;
			if (n_mismatch < 2) m1->c[n_mismatch] = 1;
			// set mapping qualities
			c = tabread(fp_export, str); // SE mapping score
			c = atoi(str);
			if (c > 255) c = 255;
			m1->map_qual = m1->seq[MAX_READLEN-1] = c;
			c = tabread(fp_export, str); // PE mapping score
			c = atoi(str);
			if (c > 255) c = 255;
			if (m1->map_qual < c) m1->map_qual = c;
			c = tabread(fp_export, str); // target of mate
			int is_same_chr = (str[0] == '\0')? 1 : 0;
			c = tabread(fp_export, str); // contig of mate?
			// set m1->dist
			c = tabread(fp_export, str); // offset to mate
			int dist = atoi(str);
			m1->dist = (dist == 0 || !is_same_chr)? 0 : (dist > 0? dist + read_len[1-cur_read] : dist - m1->size);
			// set m1->flag
			c = tabread(fp_export, str); // strand of mate
			m1->flag = 0;
			if (str[0] == 'N') m1->flag |= PAIRFLAG_NOMATCH;
			else {
				if (!is_same_chr) m1->flag |= PAIRFLAG_DIFFCHR;
				else {
					if (dist > 0) m1->flag |= 1 << ((m1->pos&1)<<1 | (str[0] == 'F'? 0 : 1));
					else m1->flag |= 1 << ((str[0] == 'F'? 0 : 1)<<1 | (m1->pos&1));
					if (abs(m1->dist) < max_dist && m1->flag == PAIRFLAG_FR)
						m1->flag |= PAIRFLAG_PAIRED;
				}
			}
			c = tabread(fp_export, str); // filtered or not
			if (is_non_filtered || str[0] == 'Y') ++mm->n_mapped_reads;
		} else while ((c = fgetc(fp_export)) != EOF && c != '\n');			
	}
	algo_sort(mm->n_mapped_reads, mm->mapped_reads);
	maqmap_write_header(fp, mm);
	gzwrite(fp, mm->mapped_reads, sizeof(maqmap1_t) * mm->n_mapped_reads);
	// free
	maq_delete_maqmap(mm);
	for (i = 0; i != int(hash->size()); ++i) free(name_conv[i]);
	free(name_conv);
	delete hash;
}

static int export2maq_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   maq export2maq [options] <out.map> <in.list> <export.txt>\n\n");
	fprintf(stderr, "Options: -1 INT   length of read1 [0]\n");
	fprintf(stderr, "         -2 INT   length of read2 [0]\n");
	fprintf(stderr, "         -a INT   maximum insert size [250]\n");
	fprintf(stderr, "         -n       keep filtered reads\n\n");
	return 1;
}

int maq_export2maq(int argc, char *argv[])
{
	FILE *fp_export, *fp_list;
	gzFile fp_map;
	int c;
	while ((c = getopt(argc, argv, "1:2:a:n")) >= 0) {
		switch (c) {
		case '1': read_len[0] = atoi(optarg); break;
		case '2': read_len[1] = atoi(optarg); break;
		case 'a': max_dist = atoi(optarg); break;
		case 'n': is_non_filtered = 1; break;
		}
	}
	if (optind + 2 >= argc) return export2maq_usage();
	fp_export = (strcmp(argv[optind+2], "-") == 0)? stdin : fopen(argv[optind+2], "r");
	fp_list = fopen(argv[optind+1], "r");
	fp_map = gzopen(argv[optind], "w");
	assert(fp_export && fp_list && fp_map);
	export2maq_core(fp_list, fp_export, fp_map);
	fclose(fp_export); fclose(fp_list); gzclose(fp_map);
	return 0;
}
