#!/bin/bash
set -xe

## A shell for collect selected taxa for phylogeny
## By Xuankun Li: xli8@memphis.edu

store_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage"
python_script_dir="/scratch/group/songlab/phylo/OR_TE/environment/python_scripts" #location of the python scripts

module load GCC/7.3.0-2.30  CUDA/9.2.88  OpenMPI/3.1.1 Python/2.7.15

mkdir test_tree
for taxa in $(cat tree_samples.txt)
do 
##collect genes by taxa names
	sed -n "/${taxa}/,+1p" ${store_dir}/test.aa.fasta >> test_tree/tree_all_gene.aa.fasta # Modify here
	sed -n "/${taxa}/,+1p" ${store_dir}/test.nt.fasta >> test_tree/tree_all_gene.nt.fasta # Modify here
done

mkdir test_tree/aa
mkdir test_tree/nt
for gene in $(cat gene_names.txt)
do
##collect genes by gene names from a single file and make one file for each gene
	sed -n "/${gene}/,+1p" test_tree/tree_all_gene.aa.fasta > test_tree/aa/${gene}.aa.fasta
	sed -n "/${gene}/,+1p" test_tree/tree_all_gene.nt.fasta > test_tree/nt/${gene}.nt.fasta
done

cp gene_names.txt test_tree/aa/
cp gene_names.txt test_tree/nt/

cd test_tree/aa/
for gene in $(cat gene_names.txt)
do
##remove gene names from the sequence headers
	sed -i "s/|${gene}//" ${gene}.aa.fasta
done

cd ../nt/
for gene in $(cat gene_names.txt)
do
##remove gene names from the sequence headers
	sed -i "s/|${gene}//" ${gene}.nt.fasta
done

#AMAS calculate gene coverage by taxa
python ${python_script_dir}/AMAS.py summary -f fasta -d dna -i *.fasta

mv summary.txt ../

cd ..
rm tree_all_gene.*.fasta
