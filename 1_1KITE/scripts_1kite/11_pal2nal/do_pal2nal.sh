#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -M k.meusemann@zfmk.de
#$ -m be
#$ -N basic_testjob

cd $HOME/1kite_100_TAXA_PAPER/NEW_016_PAL2NAL_OLD/aa_final_alignments/
for i in *.linsi.aa.fas ; do a=`basename $i .linsi.aa.fas`; perl pal2nal.mod.pl $a.linsi.aa.fas $a.nt.fas -output fasta 1> $a.linsi.nt.fas 2> $a.error; done
