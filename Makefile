CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -m64 # comment out `-m64' for 32-bit compilation
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-D_FASTMAP -DMAQ_LONGREADS
OBJS=		const.o seq.o bfa.o read.o fasta2bfa.o fastq2bfq.o \
			match_aux.o match.o sort_mapping.o merge.o get_pos.o \
			pileup.o mapcheck.o	assopt.o assemble.o maqmap.o \
			aux_utils.o rbcc.o subsnp.o pair_stat.o indel_soa.o \
			maqmap_conv.o altchr.o submap.o rmdup.o simulate.o \
			genran.o indel_pe.o stdaln.o indel_call.o eland2maq.o \
			csmap2ntmap.o break_pair.o glfgen.o
PROG=		maq
MANPAGE=	maq.1
VERSION=	0.7.1
LIBS=		-lz -lm

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $< -o $@

all:$(PROG) $(MANPAGE)

maq:$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

maq.1:maq.pod
		pod2man --center "Bioinformatics Tools" --release "maq-$(VERSION)" maq.pod > $@

main.o:main.c main.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DPACKAGE_VERSION=\"$(VERSION)\" main.c -o main.o

bfa.o:bfa.h
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

clean:
		rm -f *.o a.out *~ *.a $(PROG) $(MANPAGE)
