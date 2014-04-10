#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <zlib.h>
#include "maqmap.h"
#include "stdhash.hh"
#include "main.h"

#define BUFFER_SIZE 0x10000

typedef struct {
	int n, max;
	maqmap1_t **m1;
} tmp_stack_t;

typedef hash_map_misc<bit64_t, maqmap1_t*> best_hash_t;

static int is_auto_trim_name = 1;

static inline void stack_insert(tmp_stack_t *stack, maqmap1_t *m1)
{
	if (stack->n == stack->max) {
		stack->max += 0x10000;
		stack->m1 = (maqmap1_t**)realloc(stack->m1, sizeof(maqmap1_t) * stack->max);
	}
	stack->m1[stack->n++] = m1;
}
static inline void dump_best(tmp_stack_t *stack, best_hash_t *best_hash, gzFile fpout)
{
	for (int i = 0; i != stack->n; ++i) {
		gzwrite(fpout, stack->m1[i], sizeof(maqmap1_t));
		free(stack->m1[i]);
	}
	stack->n = 0;
	if (best_hash->size() > BUFFER_SIZE)
		best_hash->clear();
}
static inline char *trim_name(char *name, char *buf)
{
	if (is_auto_trim_name) {
		int c1, c2, tl = strlen(name);
		if (tl <= 1) return name;
		c1 = name[tl-1]; c2 = name[tl-2];
		if ((c1 == '1' || c1 == '2') && c2 == '/') { // trim the last character
			strncpy(buf, name, tl-1);
			buf[tl-1] = 0;
			return buf;
		} else return name; // keep the original name
	}
	return name;
}
void maq_rmdup_core(gzFile fpin, gzFile fpout)
{
	maqmap_t *mm;
	maqmap1_t m1;
	int last_seqid = -1;
	hash_set_char *del_set = new hash_set_char;
	best_hash_t *best_hash = new best_hash_t;
	bit32_t last_pos = 0xffffffffu;
	bit64_t n_checked = 0, n_removed = 0;
	char buf[MAX_NAMELEN+1];
	tmp_stack_t stack;
	stack.n = stack.max = 0;
	stack.m1 = 0;
	mm = maqmap_read_header(fpin);
	maqmap_write_header(fpout, mm);
	del_set->resize(4 * BUFFER_SIZE);
	best_hash->resize(3 * BUFFER_SIZE);
	while (gzread(fpin, &m1, sizeof(maqmap1_t)) == sizeof(maqmap1_t)) {
		if (int(m1.seqid) != last_seqid || last_pos != m1.pos>>1) {
			dump_best(&stack, best_hash, fpout); // write the result
			if (int(m1.seqid) != last_seqid) {
				best_hash->clear();
				if (del_set->size()) { // check
					fprintf(stderr, "[maq_rmdup_core] %u unmatched pairs\n", unsigned(del_set->size()));
					del_set->clear();
				}
				last_seqid = m1.seqid;
				fprintf(stderr, "[maq_rmdup_core] processing reference %s...\n", mm->ref_name[m1.seqid]);
			}
		}
		if (m1.dist == 0 || abs(m1.dist) <= m1.size) { // SE, singlet, diff chr or very short insert size
			gzwrite(fpout, &m1, sizeof(maqmap1_t));
		} else if (m1.dist > 0) { // paired, head
			++n_checked;
			bit64_t key = (bit64_t(m1.pos)<<32) | m1.dist;
			maqmap1_t *p;
			if (best_hash->find(key, &p)) { // found in best_hash
				++n_removed;
				if (p->map_qual < m1.map_qual) { // m1 is better
					del_set->insert(trim_name(p->name, buf)); // current best will be removed
					memcpy(p, &m1, sizeof(maqmap1_t)); // replaced as m1
				} else del_set->insert(trim_name(m1.name, buf)); // m1 will be removed
			} else { // not found
				p = (maqmap1_t*)malloc(sizeof(maqmap1_t));
				memcpy(p, &m1, sizeof(maqmap1_t));
				best_hash->insert(key, p);
				stack_insert(&stack, p);
			}
		} else { // paired, tail
			char *s = trim_name(m1.name, buf);
			if (del_set->find(s)) del_set->erase(s);
			else gzwrite(fpout, &m1, sizeof(maqmap1_t));
		}
		last_pos = m1.pos>>1;
	}
	dump_best(&stack, best_hash, fpout);
	maq_delete_maqmap(mm);	
	delete del_set;
	delete best_hash;
	free(stack.m1);
	fprintf(stderr, "[maq_rmdup_core] %lld / %lld = %.4lf\n", n_removed, n_checked, (double)n_removed/n_checked);
}
int maq_rmdup(int argc, char *argv[])
{
	gzFile fpin, fpout;
	if (argc < 3) {
		fprintf(stderr, "Usage: maq rmdup <output.map> <input.map>\n");
		return 1;
	}
	fpout = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdout), "w") : gzopen(argv[1], "w");
	fpin = (strcmp(argv[2], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[2], "r");
	assert(fpin && fpout);
	maq_rmdup_core(fpin, fpout);
	gzclose(fpin); gzclose(fpout);
	return 0;
}
