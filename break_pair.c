#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "maqmap.h"

static void breakpair(const char *prefix)
{
	gzFile fp, fp0, fp1, fp2;
	char str[1024];
	maqmap_t *mm;
	maqmap1_t m1;

	strcpy(str, prefix); strcat(str, ".map");
	assert(fp = gzopen(str, "r"));
	strcpy(str, prefix); strcat(str, ".0.map");
	fp0 = gzopen(str, "w");
	strcpy(str, prefix); strcat(str, ".1.map");
	fp1 = gzopen(str, "w");
	strcpy(str, prefix); strcat(str, ".2.map");
	fp2 = gzopen(str, "w");

	mm = maqmap_read_header(fp);
	mm->n_mapped_reads = 0;
	maqmap_write_header(fp0, mm);
	maqmap_write_header(fp1, mm);
	maqmap_write_header(fp2, mm);
	while (gzread(fp, &m1, sizeof(maqmap1_t)) == sizeof(maqmap1_t)) {
		int l;
		if (m1.map_qual <= 40) continue;
		l = strlen(m1.name);
		if (l >= 2 && m1.name[l-2] == '/' && m1.name[l-1] == '1')
			gzwrite(fp1, &m1, sizeof(maqmap1_t));
		else if (l >= 2 && m1.name[l-2] == '/' && m1.name[l-1] == '2')
			gzwrite(fp2, &m1, sizeof(maqmap1_t));
		else gzwrite(fp0, &m1, sizeof(maqmap1_t));
	}
	maq_delete_maqmap(mm);
	gzclose(fp);
	gzclose(fp0); gzclose(fp1); gzclose(fp2);
}

int maq_breakpair(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: maq breakpair <prefix>\n");
		return 1;
	}
	breakpair(argv[1]);
	return 0;
}
