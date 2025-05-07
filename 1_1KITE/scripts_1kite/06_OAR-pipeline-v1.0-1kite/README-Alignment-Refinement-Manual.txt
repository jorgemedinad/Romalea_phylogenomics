Author:  Christoph Mayer, Forschungsmuseum Alexander Koenig, Center for molecular biodiversity research
Project: The pipeline descibed in this manual is part of the 1kite analysis pipelien.


Outlier Alignment Refinement (OAR)
==================================

1.) Prerequisites of the OAR pipeline:
a) Alignments must be present in fasta format that contain the sequences that have to be refined.
The fasta files must obey the following name convention: Genename.aa.fas

Examples:
EOG500005.aa.fas
EOG50000C.aa.fas
EOG50000K.aa.fas
EOG50000R.aa.fas
EOG50001C.aa.fas

Genename can be any string. ".aa.fas" is a fixed string. Changing this file convention requires multiple changes
in multiple files of the pipeline and its is not recommended.

b) Lists of sequences that have been identified as outliers.
The OAR pipeline requires separate files which list the outliers for each
gene individually. Each file, which has to have a name equal to "Genename.aa.outlier",
lists the outliers for this gene. In other words. The file name is equal to the fasta file name,
except that ".fas" is replaced by ".outlier". All outlier files need to exist in the outlier-by-gene folder.
The outlier files are required to have the full sequence names of sequences that shall be refined,
each on a separarte line of the file:

Examples:
***********************
EOG502V72.aa.fas:
***********************
EOG502V72|BMORI_2.0|Manduca_sexta|gi343841330-389
***********************

***********************
EOG502V73.aa.fas:
***********************
EOG502V73|ISCAP_1.1|Nannochorista_sp|contig00656-161
***********************

***********************
EOG508KQ4.aa.fas:
***********************
EOG508KQ4|ISCAP_1.1|Cosmioperla_kuna|C589037-279
EOG508KQ4|ISCAP_1.1|Aretaon_asperrimus|C839012-436
EOG508KQ4|ISCAP_1.1|Gynaikothrips_ficorum|s25612_L_265458_0-433
***********************

++++++++++++++++++++++++++++++++++++++++++
+Optional part of the pipeline:
+
+--> Required script: 
+---> expand_outlier_file.pl (Perl)
+
+The pipeline includes a scipt (called "expand_outlier_file.pl") to create individial files from one file 
+which lists the outliner sequences as follows:
+
+***********************
+outliers.txt
+***********************
+EOG502V72.aa.fas:
+EOG502V72|BMORI_2.0|Manduca_sexta|gi343841330-389
+
+EOG502V73.aa.fas:
+EOG502V73|ISCAP_1.1|Nannochorista_sp|contig00656-161
+
+EOG508KQ4.aa.fas:
+EOG508KQ4|ISCAP_1.1|Cosmioperla_kuna|C589037-279
+EOG508KQ4|ISCAP_1.1|Aretaon_asperrimus|C839012-436
+EOG508KQ4|ISCAP_1.1|Gynaikothrips_ficorum|s25612_L_265458_0-433
+***********************
+
++++++++++++++++++++++++++++++++++++++++++

2.) Required data setup:

The complete OAR procedure is conducted in a separate folder,
which we call OAR-folder.

This folder must contain:
i) a file which lists the fasta file names of all genes.
  This must be a superset of the genes that are refined.
  Genes for which no outlies are found in the outlier files will be skipped,
  so it is safe to provide all gene file names of a data set in this file.
  The file must be called: "gene-file-names.txt"

Example:
*************************
gene-file-names.txt
*************************
EOG500005.aa.fas
EOG50000C.aa.fas
EOG50000K.aa.fas
EOG50000R.aa.fas
EOG50001C.aa.fas
*************************

ii) the refine-alignments.pl script, which does most of the refinement procedure
  and which calls external programs to do special tasks.

iii) a folder called "fasta_orig" which contains the fasta sequence files of the genes that need to be refined. Sequence names need to be as in the sequence name file.

iv) a folder called "outlier-by-gene", which must contains the above mentioned
 files which list the outliers for each gene. Each gene requires its own file as described
 above.



3.) Runnung the refine-alignments.pl script:

--> Required software and scipts 
--->   refine-alignments.pl         (Perl)   # must exist in the OAR folder.
--->   divide_all.pl                (Perl)   # must exist in the system path
--->   divide_core_non_core_taxa.pl (Perl)   # must exist in the system path
--->   seq-extractor-v0.4           (C++)    # must exist in the system path

Description of what the refine-alignments.pl script does.

a) Within the "refine-alignments.pl" the first step is to create separate
files for the core and non-core taxa.
The scipt creates the folders: "non-cores-individual-sequences", "cores-individual-sequences".
Within these foulders the script generates a folder for each gene that will contain the cores and non-cores
sequences as individual sequence files.
These steps require external scripts and programs, namely:
divide_all.pl
which in turn calls for each gene fasta file:
"divide_core_non_core_taxa.pl gene-file-name"
as well as
"seq-extractor-v0.4 --split fasta-file-of-cores"


b) Now "refine-alignments.pl" loops over all genes:

--> Required software and scipts
---> seq-extractor-v0.4                    (C++)    # must exist in the system path
---> remove_gap_only_positions-v1.0        (C++)    # must exist in the system path
---> mafft-linsi (Version used: v6.935b)            # must exist in the system path
---> split-fasta-with-sequence-list        (C++)    # must exist in the system path
---> move_core_taxa_to_top_of_alignment.pl (Perl)   # must exist in the system path
---> combine-align                         (C++)    # must exist in the system path


Step 1 in loop over all genes:
     Read the outlier-file for this gene.
     If there is an outlier for this gene, loop over all outlier sequence names.
     Extract all outlier sequneces to indiviual files.
     Note: There is some redundancy to what we have done above. Above we have extracted all sequences individually already.
     This has historic reasons. Basically this script was designed for different refinement approaches.
     This does not do any harm, but it could be changed in a future version of the script.
     This does waste a small amount of disc space and a few CPU seconds, but not more.    

Step 2 in loop over all genes:
     In 3a) we have extracted the alignment of the core sequences.
     Now we remove the gap only positions in all these core-taxon sequence alignments.      

Step 3 in loop over all genes:
     In this step we align each outlier sequence to its corresponding core-taxon alignment.
     This step requires the mafft program: The system call is "mafft-linsi --preservecase --add  <outlier-fasta-file> <core-taxon-alignment> 
     We chose the linsi algorithm of mafft since this is a local and not a global alignment approach. 

Step 4 in loop over all genes:
     i)   In this step we first create a multiple sequence alignment for this gene in which all outliers have been removed.
     System call: split-fasta-with-sequence-list <full-alignment> <list-of-outlier-sequences>
     ii)  Now we remove the gap only positions in the alignment in which outliers have been removed.
     iii) Next we move the core taxa to the top of the alignment. This mereley changes the order of the sequences in the fasta file.
     
     iv) In a loop over all outlier sequences which have been aligned to the core taxa in (Step 3) we now use the core taxon alignment as 
     a backbone to add the newly aligned outlier sequences to the alignment in which the outliers had been removed.
     After this has been done for all outliers, we have refined the alignemnt for this gene.
     System call: combine-align <iterate-msa> <outlier-to-core-taxon-alignment>



**************************************************
Final list of required programs and scripts:
**************************************************
refine-alignments.pl                  (Perl)  # must exist in the OAR folder.
divide_all.pl                         (Perl)  # must exist in the system path
divide_core_non_core_taxa.pl          (Perl)  # must exist in the system path
seq-extractor-v0.4                    (C++)   # must exist in the system path
remove_gap_only_positions-v1.0        (C++)   # must exist in the system path
mafft  mafft-linsi (Version used: v6.935b)    # must exist in the system path
split-fasta-with-sequence-list        (C++)   # must exist in the system path
move_core_taxa_to_top_of_alignment.pl (Perl)  # must exist in the system path
combine-align                         (C++)   # must exist in the system path









