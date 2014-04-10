bin_SCRIPTS = scripts/maq.pl scripts/farm-run.pl scripts/maq_plot.pl scripts/maq_eval.pl
bin_PROGRAMS = maq
maq_SOURCES = main.c const.c seq.c bfa.c read.cc fasta2bfa.c fastq2bfq.c merge.cc \
		match_aux.cc match.cc sort_mapping.cc assemble.cc pileup.cc mapcheck.cc get_pos.c \
		assopt.c aux_utils.c rbcc.cc subsnp.cc pair_stat.cc indel_soa.cc maqmap.c maqmap_conv.c \
		altchr.cc submap.c rmdup.cc simulate.c genran.c indel_pe.cc stdaln.c indel_call.cc \
		eland2maq.cc csmap2ntmap.cc break_pair.c glfgen.cc
maq_LDADD = -lm -lz

bfa.o:bfa.h
main.o:main.h
const.o:const.h
read.o:read.h
seq.o:seq.h
fasta2bfa.o:const.h seq.h
fastq2bfq.o:const.h seq.h
match.o:match.hh dword.hh main.h algo.hh bfa.h read.h
match_aux.o:read.h match.hh main.h
sort_mapping.o:match.hh algo.hh read.h
assemble.o:assemble.h algo.hh main.h bfa.h
assopt:assemble.h
genran.o:genran.h
simulate.o:genran.h
