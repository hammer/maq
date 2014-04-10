#include <stdio.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "dummy"
#endif
#include "main.h"


static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: maq (Mapping and Assembly with Qualities)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   maq <command> [options]\n\n");
	fprintf(stderr, "Key commands:\n");
	fprintf(stderr, "         fasta2bfa   convert FASTA to BFA format\n");
	fprintf(stderr, "         fastq2bfq   convert FASTQ to BFQ format\n");
	fprintf(stderr, "         map         map reads to the reference\n"); // valgrind on 10AUG2007
	fprintf(stderr, "         mapmerge    merge several alignments\n"); // valgrind on 11JUL2007
	fprintf(stderr, "         rmdup       remove pairs with identical outer coordinates (PE)\n");
	fprintf(stderr, "         indelpe     indel calling (PAIRED READS ONLY)\n");
	fprintf(stderr, "         indelsoa    state-of-art homozygous indel detection\n"); // valgrind on 11JUL2007
	fprintf(stderr, "         assemble    call the consensus\n"); // valgrind on 11JUL2007
	fprintf(stderr, "         glfgen      generate .glz consensus\n\n");
	fprintf(stderr, "Format converting:\n");
	fprintf(stderr, "         sol2sanger  convert Solexa FASTQ to standard/Sanger FASTQ\n");
	fprintf(stderr, "         mapass2maq  convert mapass2's map format to maq's map format\n");
	fprintf(stderr, "         bfq2fastq   convert BFQ to FASTQ format\n\n");
	fprintf(stderr, "Information extracting:\n");
	fprintf(stderr, "         mapview     view the mapping alignment\n"); // valgrind on 11JUL2007
	fprintf(stderr, "         mapstat     statistics about a .map file\n");
	fprintf(stderr, "         mapcheck    a QC command\n"); // 11JUL2007: valgrind failed on "-g -O2 -m64"; passed on "-g"
	fprintf(stderr, "         pileup      view the alignment in a 'pileup' like format\n"); // valgrind on 11JUL2007
	fprintf(stderr, "         cns2fq      extract the consensus sequences from a CNS file\n");
	fprintf(stderr, "         cns2snp     extract details from a CNS file at the SNP sites\n");
	fprintf(stderr, "         snpreg      calculate the length of regions where SNPs can be called\n");
	fprintf(stderr, "         cns2view    extract details from a CNS file at all sites\n");
	fprintf(stderr, "         cns2ref     extract the reference sequences from a CNS file\n");
	fprintf(stderr, "         cns2win     extract details in a window\n\n");
	fprintf(stderr, "SOLiD related commands:\n");
	fprintf(stderr, "         fasta2csfa  convert FASTA to colour-space FASTA\n");
	fprintf(stderr, "         csmap2nt    convert colour-space .map to nucleotide .map\n\n");
//	fprintf(stderr, "         cns2cssnp   SNP calling for AB SOLiD\n\n");
	fprintf(stderr, "Simulation related commands:\n");
	fprintf(stderr, "         fakemut     simulate references: randomly generate mutations\n");
	fprintf(stderr, "         simutrain   train parameters for simulation\n");
	fprintf(stderr, "         simulate    simulate reads: randomly generate sequencing errors\n");
	fprintf(stderr, "         simucns     evaluate consensus based on simulation\n");
	fprintf(stderr, "         simustat    evaluate alignment based on simulation\n\n");
	fprintf(stderr, "Miscellaneous/advanced utilities:\n");
	fprintf(stderr, "         submap      extract a region from a map file\n");
	fprintf(stderr, "         mapvalidate validate a .map file\n");
//	fprintf(stderr, "         rbcc        reference based consensus calling\n");
	fprintf(stderr, "         subpos      extract a subset of positions\n");
	fprintf(stderr, "         eland2maq   convert eland alignment to maq\n");
	fprintf(stderr, "         export2maq  convert Solexa's export alignment to maq\n");
	fprintf(stderr, "         novo2maq    convert novoalign/novopaired alignment to maq\n");
#ifdef MAQ_SHOW_EXPERIMENTAL
	fprintf(stderr, "         abpair      show the abnormal pairs (EXPERIMENTAL)\n");
	fprintf(stderr, "         paircov     paired coverage (EXPERIMENTAL)\n");
	fprintf(stderr, "         altchr      change reference according to SNPs (EXPERIMENTAL)\n");
#endif
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fasta2bfa") == 0) return ma_fasta2bfa(argc-1, argv+1);
	else if (strcmp(argv[1], "fastq2bfq") == 0) return ma_fastq2bfq(argc-1, argv+1);
	else if (strcmp(argv[1], "map") == 0) return ma_match(argc-1, argv+1);
	else if (strcmp(argv[1], "mapmerge") == 0) return ma_mapmerge(argc-1, argv+1);
	else if (strcmp(argv[1], "indelpe") == 0) return maq_indelpe(argc-1, argv+1);
	else if (strcmp(argv[1], "indelsoa") == 0) return ma_indelsoa(argc-1, argv+1);
	else if (strcmp(argv[1], "assemble") == 0) return ma_assemble(argc-1, argv+1);
	else if (strcmp(argv[1], "glfgen") == 0) return maq_glfgen(argc-1, argv+1);
	else if (strcmp(argv[1], "sol2sanger") == 0) return ma_sol2sanger(argc-1, argv+1);
	else if (strcmp(argv[1], "mapass2maq") == 0) return ma_mapass2maq(argc-1, argv+1);
	else if (strcmp(argv[1], "bfq2fastq") == 0) return ma_bfq2fastq(argc-1, argv+1);
	else if (strcmp(argv[1], "mapview") == 0) return ma_mapview(argc-1, argv+1);
	else if (strcmp(argv[1], "mapcheck") == 0) return ma_mapcheck(argc-1, argv+1);
	else if (strcmp(argv[1], "pileup") == 0) return ma_pileup(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2ref") == 0) return ma_cns2ref(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2fq") == 0) return ma_cns2fq(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2snp") == 0) return ma_cns2snp(argc-1, argv+1);
	else if (strcmp(argv[1], "snpreg") == 0) return ma_snpreg(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2view") == 0) return ma_cns2view(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2win") == 0) return ma_cns2win(argc-1, argv+1);
	else if (strcmp(argv[1], "fasta2csfa") == 0) return ma_fasta2csfa(argc-1, argv+1);
	else if (strcmp(argv[1], "csmap2nt") == 0) return maq_csmap2nt(argc-1, argv+1);
	else if (strcmp(argv[1], "cns2cssnp") == 0) return ma_cns2cssnp(argc-1, argv+1);
	else if (strcmp(argv[1], "fakemut") == 0) return ma_fakemut(argc-1, argv+1);
	else if (strcmp(argv[1], "simutrain") == 0) return maq_simutrain(argc-1, argv+1);
	else if (strcmp(argv[1], "simulate") == 0) return maq_simulate(argc-1, argv+1);
	else if (strcmp(argv[1], "simucns") == 0) return maq_simucns(argc-1, argv+1);
	else if (strcmp(argv[1], "simustat") == 0) return maq_simustat(argc-1, argv+1);
	else if (strcmp(argv[1], "match") == 0) return ma_match(argc-1, argv+1);
	else if (strcmp(argv[1], "rmdup") == 0) return maq_rmdup(argc-1, argv+1);
	else if (strcmp(argv[1], "submap") == 0) return maq_submap(argc-1, argv+1);
	else if (strcmp(argv[1], "mapvalidate") == 0) return ma_mapvalidate(argc-1, argv+1);
	else if (strcmp(argv[1], "rbcc") == 0) return ma_rbcc(argc-1, argv+1);
	else if (strcmp(argv[1], "mapstat") == 0) return ma_mapstat(argc-1, argv+1);
	else if (strcmp(argv[1], "subpos") == 0) return ma_subpos(argc-1, argv+1);
	else if (strcmp(argv[1], "eland2maq") == 0) return maq_eland2maq(argc-1, argv+1);
	else if (strcmp(argv[1], "export2maq") == 0) return maq_export2maq(argc-1, argv+1);
	else if (strcmp(argv[1], "novo2maq") == 0) return maq_novo2maq(argc-1, argv+1);
	else if (strcmp(argv[1], "abpair") == 0) return ma_abpair(argc-1, argv+1);
	else if (strcmp(argv[1], "paircov") == 0) return ma_paircov(argc-1, argv+1);
	else if (strcmp(argv[1], "altchr") == 0) return ma_altchr(argc-1, argv+1);
	else if (strcmp(argv[1], "catfilter") == 0) return maq_catfilter(argc-1, argv+1);
	else if (strcmp(argv[1], "breakpair") == 0) return maq_breakpair(argc-1, argv+1);
	else {
		fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
		return 2;
	}
	return 0;
}
