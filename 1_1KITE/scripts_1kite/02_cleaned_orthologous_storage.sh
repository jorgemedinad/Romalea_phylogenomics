#!/bin/bash
set -xe

module purge
module load GCC/10.2.0  Perl/5.32.0

## A shell for collect and store cleaned orthologous loci
## By Xuankun Li: xli8@memphis.edu

work_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_Taeniopoda4/1KITE_example" # Location of the current working directory
store_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage" # Location of cleaned orthologous storage directory

# 11 Pal2Nal (second half, post-manual check)
cd ${work_dir}/11_pal2nal
mkdir aa_aligned_non_interleaved
mv *.linsi.aa.fas aa_aligned_non_interleaved/
mkdir ErrorFiles
mv *.error ErrorFiles/
mkdir nt_aligned_interleaved
mv *.linsi.nt.fas nt_aligned_interleaved/
mkdir nt_not_aligned
mv *.nt.fas nt_not_aligned/
mkdir nt_aligned_non_interleaved

cd nt_aligned_interleaved
for i in *.linsi.nt.fas; do perl ${work_dir}/fastasingleline.pl $i > $i.noninterleaved; done
mv *.noninterleaved ../nt_aligned_non_interleaved/
cd ../nt_aligned_non_interleaved/
rename .linsi.nt.fas.noninterleaved .linsi.nt.fas *.noninterleaved
# generate: cleaned aa files in directory 'aa_aligned_non_interleaved', cleaned nt files in directory 'nt_aligned_non_interleaved'

# 13a Adjust sequence headers to prepare for Alignment Masking
cd ${work_dir}/13a_rewrite_headers_for_Aliscore/aa_headers_renamed
cp ../../11_pal2nal/aa_aligned_non_interleaved/*.fas .
cp ${work_dir}/13a_rewrite_headers_for_Aliscore/rewrite_headers.pl .
perl rewrite_headers.pl
rename .linsi.aa.fas.fas.cleaned .linsi.aa.fas *.cleaned

cd ../nt_headers_renamed
cp ../../11_pal2nal/nt_aligned_non_interleaved/*.fas .
cp ${work_dir}/13a_rewrite_headers_for_Aliscore/rewrite_headers.pl .
perl rewrite_headers.pl
rename .linsi.nt.fas.fas.cleaned .linsi.nt.fas *.cleaned
# generate: aa files in directory 'aa_headers_renamed', cleaned nt files in directory 'nt_headers_renamed'

# Rename sequence headers and store cleaned orthologous
cd ${work_dir}/13a_rewrite_headers_for_Aliscore/aa_headers_renamed
ls *.fas | sed 's/.linsi.aa.fas//g' > ${work_dir}/gene_names.txt
cd ${work_dir}

for gene in $(cat gene_names.txt)
do
#put genes names with taxa names as format:'280_S36|Her_ATP6'
	sed "/>/s/$/|${gene}/" ${work_dir}/13a_rewrite_headers_for_Aliscore/aa_headers_renamed/${gene}.linsi.aa.fas >> ${store_dir}/Taeniopoda4.aa.fasta # Modify here
done

for gene in $(cat gene_names.txt)
do
#put genes names with taxa names as format:'280_S36|Her_ATP6'
	sed "/>/s/$/|${gene}/" ${work_dir}/13a_rewrite_headers_for_Aliscore/nt_headers_renamed/${gene}.linsi.nt.fas >> ${store_dir}/Taeniopoda4.nt.fasta # Modify here
done
# generate: *.aa.fasta and *.nt.fasta in the directory /scratch/user/hsong/OR_TE/1KITE_cleaning/cleaned_orthologous_storage
