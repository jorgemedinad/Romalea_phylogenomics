#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -m be
#$ -M k.meusemann@zfmk.de
#$ -N mafft

cd $HOME/209_gene_and_headers_renamed/
for i in *.fas; do /share/apps/mafft-linsi_7.123we --quiet $i > $i.linsi.fas; done
