#!/bin/sh
#SBATCH --job-name=1k_auto
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jorgemedinad@tamu.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=10:00:00


## A shell for automatic perform 1KITE pipeline steps 2 to 11
## By Xuankun Li: xli8@memphis.edu

module load GCC/10.2.0  OpenMPI/4.0.5 Perl/5.32.0 MAFFT/7.490-with-extensions

raw_dir="/scratch/group/songlab/phylo/OR_TE/assembly/combined_sample_name" # Location of the sample names
ortho_dir="/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_Taeniopoda2" # Location of the orthograph results 
work_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_Taeniopoda4/1KITE_example" # Location of the current working directory, replace test_cleaning with current directory

# 02 RECONVERT Orthograph into HaMStRad format
cd ${work_dir}
cp 01_Orthograph/* ${ortho_dir}/
cd ${ortho_dir}
bash O2H_header_auto.sh

cd ${work_dir}/02b_species_hamstrad_header
rm -r O2H_header_auto.sh/
rm -r log_and_others/
rm -r orthograph2hamstrad.pl/
# generate: header modified files in directory '02b_species_hamstrad_header'

# 03 Summarize Orthograph results
cd ${work_dir}/03_Summarize_Orthograph
perl summarize_orthograph_results.pl -i ../02b_species_hamstrad_header/ -o . -m
# generate: collect loci from all samples into each 'file_by_loci' in directory '03_Summarize_Orthograph'

# 04 Alignment
cd ${work_dir}
cp 03_Summarize_Orthograph/aa_summarized/* 04_Alignment/aa/
cp 03_Summarize_Orthograph/nt_summarized/* 04_Alignment/nt/
cd 04_Alignment/nt/
rename .nt.summarized.fa .nt.fa *.fa
cd ../aa/
rename .aa.summarized.fa .aa.fa *.fa

for i in *.fa; do mafft-linsi --thread 4 --anysymbol $i > ../aa_aligned/$i.linsi.fas;
done

cd ../aa_aligned/
for i in *.aa.fa.linsi.fas; do perl ${work_dir}/fastasingleline.pl $i > $i.noninterleaved; done
rm *.aa.fa.linsi.fas
rename .aa.fa.linsi.fas.noninterleaved .aa.linsi.fas *.noninterleaved
# generate: aligned, noninterleaved files in directory '04_Alignment/aa_aligned'

# 05 Outlier check (part 1)
cd ${work_dir}/04_Alignment/aa_aligned
cp ${raw_dir}/sample_names_Taeniopoda4.txt ./subjects.txt
cp ../../05_CheckOutlier_results/checker_complete* .

perl checker_complete.1.3.1.2.pl
mv outlier.txt ../../05_CheckOutlier_results/outlier_1.txt
mv log.txt ../../05_CheckOutlier_results/log_1.txt
perl checker_complete.1.3.2.2.pl
mv outlier.txt ../../05_CheckOutlier_results/outlier_2.txt
mv log.txt ../../05_CheckOutlier_results/log_2.txt
# generate: outlier_1.txt, olog_1.txt, and outlier_2.txt, olog_2.txt in directory '05_CheckOutlier_results'

# 06 OAR-pipeline-v1.0-1kite
## STEP A: GENERATE SEPARATE TEXT FILES CONTAINING OUTLIER SEQUENCE-HEADERS
cd ${work_dir}/05_CheckOutlier_results
cp ../06_OAR-pipeline-v1.0-1kite_results/Perl/expand-outlier-file-with-example/expand_outlier_file.pl .
perl expand_outlier_file.pl outlier_2.txt

cp *.outlier ../06_OAR-pipeline-v1.0-1kite_results/outlier-by-gene
cp ../04_Alignment/aa_aligned/*.fas ../06_OAR-pipeline-v1.0-1kite_results/fasta_orig/
cd ../06_OAR-pipeline-v1.0-1kite_results/outlier-by-gene
rename .aa.linsi.outlier .aa.outlier *.outlier
sed -i "s/ /_/g" *.outlier
cd ../fasta_orig/
rename .aa.linsi.fas .aa.fas *.fas
sed -i "s/ /_/g" *.fas

## STEP B: GENERATE TEXT FILE WITH GENE NAMES ONLY
ls *.fas > gene-file-names.txt
cp gene-file-names.txt ../

## STEP C: COMPILE C++ SOURCE CODE
cd ${work_dir}/06_OAR-pipeline-v1.0-1kite_results/C++
module load gcc/8.2.0
chmod u+x compile-all
./compile-all
chmod 777 *
cd ..
cp Perl/refine-alignments.pl .

cd Perl/
chmod 777 *
cd ..
chmod 777 *

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/
export PATH=$PATH:${work_dir}/06_OAR-pipeline-v1.0-1kite_results/C++;
echo $PATH

## STEP E: PERFORM OAR: THE ALIGNMENT REFINEMENT OF THE OUTLIER SEQUENCES
module purge
module load GCC/8.2.0-2.31.1
perl refine-alignments.pl
# generate: iter_*.fas in directory 'mafft-refined-alignments-combined-aligned'

# 07 Outlier check (part 2)
cd ${work_dir}/07_outlier_check_II
cp ../06_OAR-pipeline-v1.0-1kite_results/mafft-refined-alignments-combined-aligned/iter* .
cp ../04_Alignment/aa_aligned/checker_complete* .
cp ../04_Alignment/aa_aligned/subjects.txt .

perl checker_complete.1.3.1.2.pl
mv outlier.txt outlier_1_part2.txt
mv log.txt log_1_part2.txt
perl checker_complete.1.3.2.2.pl
mv outlier.txt outlier_2_part2.txt
mv log.txt log_2_part2.txt
# generate: outlier_2_part2.txt contains all sequence headers of the remaining outliers

# 08.1 Remove outliers (screen transcripts)
cd ${work_dir}/08_remove_outliers/aa/
cp ../../06_OAR-pipeline-v1.0-1kite_results/fasta_orig/*.aa.fas .
cd ../../07_outlier_check_II
rename iter_EOG EOG *.fas
rename .fas .aa.fas *.fas
cp *.fas ../08_remove_outliers/aa
cd ../08_remove_outliers/nt
cp ../../04_Alignment/nt/*.fa .
rename .nt.fa .nt.fas *.fa
sed -i "s/ /_/g" *.fas
cd ../../07_outlier_check_II
grep '|' outlier_2_part2.txt | sed -e 's/^/>/' > outliers2remove.txt

# 08.3 Remove all outliers (aa & nt) files
cp outliers2remove.txt ../08_remove_outliers/aa
cp outliers2remove.txt ../08_remove_outliers/nt
cp ${work_dir}/08_remove_outliers/multifastafilter.pl ../08_remove_outliers/aa
cp ${work_dir}/08_remove_outliers/multifastafilter.pl ../08_remove_outliers/nt
cd ../08_remove_outliers/aa
for i in *.fas; do perl multifastafilter.pl outliers2remove.txt $i; done
cd ../nt
for i in *.fas; do perl multifastafilter.pl outliers2remove.txt $i; done
cd ..
mv aa/*.filtered aa_outliers_removed/
mv nt/*.filtered nt_outliers_removed/
cd aa_outliers_removed
rename .aa.fas.filtered .linsi.aa.fas *.filtered
cd ../nt_outliers_removed/
rename nt.fas.filtered nt.fas *.filtered
# generate: *.linsi.aa.fas in directory 'aa_outliers_removed' and *.nt.fas in directory 'nt_outliers_removed'

# 09 Remove sequences from reference species (aa alignments & nt files)
cd ${work_dir}/09_remove_ref_taxa_results
cp ../08_remove_outliers/aa_outliers_removed/*.fas aa_core_taxa_removed/
cp ../08_remove_outliers/nt_outliers_removed/*.fas nt_core_taxa_removed/
cp multifastafilterNEW.pl aa_core_taxa_removed/
cp multifastafilterNEW.pl nt_core_taxa_removed/
cp RemoveCoreTaxa_DB.sh aa_core_taxa_removed/

# 09a Remove reference species sequences from aa alignment files
cd ${work_dir}/09_remove_ref_taxa_results/aa_core_taxa_removed/
bash RemoveCoreTaxa_DB.sh
for file in EOG*.fas; do perl multifastafilterNEW.pl -v header2remove.txt $file; done
rename .linsi.aa.fas.filtered .linsi.aa.fas *.filtered 

# 09b Remove reference species sequences from nt unaligned files
cp header2remove.txt ../nt_core_taxa_removed
cd ../nt_core_taxa_removed
for file in EOG*.fas; do perl multifastafilterNEW.pl -v header2remove.txt $file; done
rename .nt.fas.filtered .nt.fas *.filtered
# generate: .linsi.aa.fas in directory 'aa_core_taxa_removed' and *.nt.fas in directory 'nt_core_taxa_removed'

# 10 Remove gaps only
cd ${work_dir}/10_remove_gaps_only_results
cp ../09_remove_ref_taxa_results/aa_core_taxa_removed/*.fas .
for i in *.fas; do perl selectSites.pl -s '1-' -x 1 $i > $i.gapsremoved; done

rm *.linsi.aa.fas
for i in *.fas.gapsremoved; do perl ${work_dir}/fastasingleline.pl $i > $i.noninterleaved; done
rm *.fas.gapsremoved
rename .linsi.aa.fas.gapsremoved.noninterleaved .linsi.aa.fas *.noninterleaved
# generate: *.linsi.aa.fas

# 11 Pal2Nal (first half, pre-manual check)
cd ${work_dir}/11_pal2nal
cp ../10_remove_gaps_only_results/*.fas .
cp ../09_remove_ref_taxa_results/nt_core_taxa_removed/*.fas .
for i in *.linsi.aa.fas ; do a=`basename $i .linsi.aa.fas`; perl pal2nal.mod.pl $a.linsi.aa.fas $a.nt.fas -output fasta 1>$a.linsi.nt.fas 2>$a.error; done
# generate: *.linsi.nt.fas cleaned in the same way as *.linsi.aa.fas
