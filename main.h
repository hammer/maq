#ifndef LH3_MAIN_H
#define LH3_MAIN_H

#ifdef __cplusplus
extern "C" {
#endif
	int ma_fasta2bfa(int argc, char *argv[]);   /* in fasta2bfa.c */
	int ma_fastq2bfq(int argc, char *argv[]);
	int ma_fasta2csfa(int argc, char *argv[]);  /* in fasta2bfa.c */
	int ma_match(int argc, char *argv[]);
	int ma_mapview(int argc, char *argv[]);     /* in maqmap.c */
	int ma_mapstat(int argc, char *argv[]);
	int ma_mapmerge(int argc, char *argv[]);
	int ma_pileup(int argc, char *argv[]);
	int ma_assemble(int argc, char *argv[]);
	int maq_glfgen(int argc, char *argv[]);     /* in glfgen.cc */
	int ma_cns2ref(int argc, char *argv[]);
	int ma_snpreg(int argc, char *argv[]);      /* in subsnp.cc */
	int ma_cns2fq(int argc, char *argv[]);      /* in aux_utils.c */
	int ma_cns2snp(int argc, char *argv[]);
	int ma_cns2cssnp(int argc, char *argv[]);
	int ma_cns2view(int argc, char *argv[]);
	int ma_cns2win(int argc, char *argv[]);
	int ma_mapcheck(int argc, char *argv[]);
	int ma_rbcc(int argc, char *argv[]);
	int ma_subpos(int argc, char *argv[]);
	int ma_sol2sanger(int argc, char *argv[]);
	int ma_bfq2fastq(int argc, char *argv[]);
	int ma_fakemut(int argc, char *argv[]);     /* in fasta2bfa.c */
	int ma_abpair(int argc, char *argv[]);      /* in pair_stat.cc */
	int ma_paircov(int argc, char *argv[]);     /* in pair_stat.cc */
	int ma_indelsoa(int argc, char *argv[]);    /* in indel_soa.cc */
	int ma_mapass2maq(int argc, char *argv[]);  /* in maqmap_conv.c */
	int ma_altchr(int argc, char *argv[]);      /* in altchr.cc */
	int maq_submap(int argc, char *argv[]);     /* in submap.c */
	int ma_mapvalidate(int argc, char *argv[]); /* in maqmap.c */
	int maq_rmdup(int argc, char *argv[]);      /* in rmdup.cc */
	int maq_simulate(int argc, char *argv[]);   /* in simulate.c */
	int maq_simutrain(int argc, char *argv[]);  /* in simulate.c */
	int maq_simustat(int argc, char *argv[]);   /* in simulate.c */
	int maq_indelpe(int argc, char *argv[]);    /* in indel_call.cc */
	int maq_eland2maq(int argc, char *argv[]);  /* in eland2maq.cc */
	int maq_export2maq(int argc, char *argv[]); /* in eland2maq.cc */
	int maq_novo2maq(int argc, char *argv[]); /* in eland2maq.cc */
	int maq_simucns(int argc, char *argv[]);    /* in subsnp.cc */
	int maq_catfilter(int argc, char *argv[]);  /* in read.cc */
	int maq_csmap2nt(int argc, char *argv[]);   /* in csmap2ntmap.cc */
	int maq_breakpair(int argc, char *argv[]);  /* in break_pair.c */
#ifdef __cplusplus
}
#endif

#endif
