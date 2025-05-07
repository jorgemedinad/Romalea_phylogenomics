# /***************************************************************************************************
# *  This file is part of the 1kite analysis pipeline. It is
# *  distributed under the following license:
# *  
# *  Copyright (c) 2012-20013 Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
# *  All rights reserved.
# *  
# *  Redistribution and use in source and binary forms, with or without
# *  modification, are permitted provided that the following conditions are met:
# *  1. Redistributions of source code (complete or in parts) must retain
# *     the above copyright notice, this list of conditions and the following disclaimer.
# *  2. Redistributions in binary form must reproduce the above copyright
# *     notice, this list of conditions and the following disclaimer in the
# *     documentation and/or other materials provided with the distribution.
# *  3. All advertising materials mentioning features or any use of this software
# *     e.g. in publications must display the following acknowledgement:
# *     This product includes software developed by Christoph Mayer, Forschungsmuseum
# *     Alexander Koenig, Bonn, Germany.
# *  4. Neither the name of the organization nor the
# *     names of its contributors may be used to endorse or promote products
# *     derived from this software without specific prior written permission.
# *  
# *  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
# *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
# *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# *  
# *  IMPORTANT (needs to be included, if code is redistributed):
# *  Please not that this license is not compatible with the GNU Public License (GPL)
# *  due to paragraph 3 in the copyright. It is not allowed under any
# *  circumstances to use the code of this software in projects distributed under the GPL.
# *  Furthermore, it is not allowed to redistribute the code in projects which are
# *  distributed under a license which is incompatible with one of the 4 paragraphs above.
# *  
# *  This project makes use of code coming from other projects. What follows is a complete
# *  list of files which make use of external code. Please refer to the copyright within
# *  these files.
# *  
# *  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
# *                                See copyright in tclap/COPYRIGHT file for details.	
# ***************************************************************************************************/

#!/usr/bin/perl

use strict;

my @genes;
my @outlier_files;

my $filename_genes="gene-file-names.txt";



# The script first wants to read all names of genes.
# This must be a superset of the genes that are refined.
# Genes for which no outlies are found in the outlier file will be skipped,
# so it is safe to give this script a list of all genes available.

open(FILE, $filename_genes) or die("Unable to open file $filename_genes");
@genes = <FILE>;
close(FILE);

## print @genes;

## File names must be of the form: "EOG500005.aa.fas", where EOG500005 is the
## gene name.
## Extract pure gene names by removing ".aa.fas" from the file names.
foreach (@genes)
{
    s/.aa.fas//;
}

# print @genes;

my @outliers;
chomp @genes;
my $gene;


## Check whether all input and result folders exist. 
## Some folders must be created by the use (e.g. input sequences),
## others can be created with this scipt.

if (1) ## Switch on or off. Switch off if this step has been done.
{
    ## The folder "fasta_orig" must exist and it must contain the original
    ## (mafft) alignment files, that need to be refinded.
    ## This folder must be put in place by the user.
    if ( !(-d "fasta_orig") )
    {
	print "Refinement script requires a folder with the original fasta sequences. Its name must be fasta_orig.\n\n";
	exit(-1);
    }

    ## The folder outlier-by-gene must exist and must contain a text file for
    ## each gene that has an outlier. This file must list all full sequence names
    ## of outliers in this gene.
    ## This folder must be put in place by the user.
    
     if ( !(-d "outlier-by-gene") )
    {
	print "Refinement script requires a folder which contains one file per gene that lists the outliers in this gene. Only these outliers will be refined. Genes for which no such file is present will be skipped in the refinement procedure.\n\n";
	exit(-1);
    }

    ## We need the following folders for intermediate results:
    if ( !(-d "outliers-individual-per-gene") )
    {
	mkdir "outliers-individual-per-gene" or die "Could not create directory: outliers-individual-per-gene";
    }

    if (!(-d "alignment_cores_plus_outlier") )
    {
	mkdir "alignment_cores_plus_outlier" or die "Could not create directory: alignment_cores_plus_outlier";
    }

    if (!(-d "alignments-without-outliers"))
    {
	mkdir "alignments-without-outliers"
	    or die "Cannot create directory alignments-without-outliers"; 
    }

    if (!(-d "mafft-refined-alignments-combined-aligned"))
    {
	mkdir "mafft-refined-alignments-combined-aligned"
	    or die "Cannot create directory mafft-refined-alignments-combined-aligned"; 
    }

    if (!(-d "mafft-refined-alignments-realigned-ranges"))
    {
	mkdir "mafft-refined-alignments-realigned-ranges"
	    or die "Cannot create directory mafft-refined-alignments-realigned-ranges"; 
    }

    if (!(-d "fasta-files-also_contains_extracted_core_taxon_alignments") )
    {
	mkdir "fasta-files-also_contains_extracted_core_taxon_alignments"
	    or die "Could not make directory: fasta-files-also_contains_extracted_core_taxon_alignments";
    }

    ## We copy all fata original fasta files to this folder:
    system("cp -pr fasta_orig/* fasta-files-also_contains_extracted_core_taxon_alignments");
}



#########
##
## Run scipt: divide_all.pl
## in fasta-files-also_contains_extracted_core_taxon_alignments directory.

if (1)
{
    ## For some reson, perl does not like "system("divide_all.pl *.aa.fas")" so we will have to do this expansion with a trick:
    chdir ("fasta-files-also_contains_extracted_core_taxon_alignments");

    my $ls_res = `ls *.aa.fas`;
    $ls_res =~ s/\n/ /g;
##     print "--- $ls_res --- \n";

    system("divide_all.pl $ls_res");
    chdir ("..");
}
## This calls the divide_core_non_core_taxa.pl which separates the cores from the non core taxa.
## Furthermore, this script splits all non-core taxa into individual sequences and places them into
## the non-cores-individual-sequences directory, by creating a directory for each gene and putting all
## non-core sequences of one gene in there.
##
#########


## The following loop runs over all genes, not only the once that need be refined.
## Thus, several actions are done for all genes in the gene-file-names.txt file.
## Alignment that do not contain outliers will not be modified.
## This behaviour is not a problem if it does not consume to many resouces.
## For the 1kite data set, which is not very small, it costs less than 10 seconds
## to do the redundant steps.

foreach $gene (@genes)
{
    ## DEBUG: The pipeline can be run step by step if necessary by switching ON/OFF different tasks.
    ## All tasks that create important output are encapsulated with an if (0/1) {} section.
    ## Only commands needed by all sections are not encapsulated. E.g. commands reading the outliers,
    ## setting variable names and debug output.  

    print "Gene: $gene\n";

    if (! (-d "outliers-individual-per-gene/$gene") )
    {
	mkdir "outliers-individual-per-gene/$gene" or die "Could not create directory: outliers-individual-per-gene/$gene";
    }


    my $input_outliers_filename;
    $input_outliers_filename = "outlier-by-gene/$gene.aa.outlier";

    if (-e "$input_outliers_filename")
    {
	print "${gene}:+:$input_outliers_filename\n";
    }
    else
    {
	print "${gene}:-:$input_outliers_filename\n";
    }
    
    @outliers = ();
    ## Is there an outlier file for this gene? If not, there are no outliers.
    if (-e $input_outliers_filename)
    {
	open(OUTL, $input_outliers_filename) or die("Unable to open file $input_outliers_filename");
	@outliers = <OUTL>;
	close OUTL;

	print @outliers;
	chomp(@outliers);

	my $fasta_file_to_extract_from = "fasta_orig/$gene.aa.fas";
	print "Extract outliers from this file: $fasta_file_to_extract_from\n";

	## First step:
	## Extract outliers:

	## TODO: Some things we do are redundant: We have already extracted all sequences above.
        ##       It is more straight forward to skip the extraction above. -> This would be less rewrite here.

	if (1) ## switch on and off - extract outlier sequences to individual files
	{
	    my $outl;
	    my $outl_fname;
	    
	    foreach $outl (@outliers)
	    {
		$outl_fname = $outl;     ## Sequence name
		$outl_fname =~ tr/|/_/;  ## Sequence name for for which "|" has been converted to "_".
		my $single_seq_name = "outliers-individual-per-gene/$gene/$outl_fname".".fas";

		my $extract_command = "seq-extractor-v0.4 -s \"$outl\" $fasta_file_to_extract_from >  $single_seq_name";  ## Create file consisting of the outlier sequence only.
		print "$extract_command\n";
		
	    system ($extract_command);
	    }
	} ## END switch on and off - extract outlier sequences to individual files

        ## Second step:
	if (1) ## Remove gap only positions in core taxon alignment.
	{
	  my $coreAlig        = "fasta-files-also_contains_extracted_core_taxon_alignments/" . $gene . ".aa_cores.fas";
	  my $coreAlig_new    = "fasta-files-also_contains_extracted_core_taxon_alignments/" . $gene . ".aa_cores.gaponly-removed.fas";

	  my $remove_gap_only_command = "remove_gap_only_positions-v1.0 $coreAlig > $coreAlig_new";

#	  print  "$remove_gap_only_command\n";
	  
	  system($remove_gap_only_command);
	}

	## Third step:
	## Align outlier against cores:
	if (1)
	{
	    my $outl_seqname;
	    my $outl_fname;

	    foreach $outl_seqname (@outliers)
	    {
		$outl_fname = $outl_seqname;
		$outl_fname =~ tr/|/_/;
		$outl_fname =~ tr/\'/_/;

		my $single_seq_name    = "outliers-individual-per-gene/$gene/$outl_fname".".fas";

		my $coreAlig_new       = "fasta-files-also_contains_extracted_core_taxon_alignments/" . $gene . ".aa_cores.gaponly-removed.fas";
		my $coreAlig_plus_seq  = "alignment_cores_plus_outlier/" . $gene . ".aa_cores_plus_".$outl_fname.".fas";


		my $alig_command1 = "mafft --maxiterate 1000 --localpair --preservecase --add $single_seq_name $coreAlig_new   > $coreAlig_plus_seq";
		# my $alig_command1 = "mafft-linsi --preservecase --add $single_seq_name $coreAlig_new   > $coreAlig_plus_seq";


		print ("$alig_command1\n\n");
		system ($alig_command1);
	    }
	}

###****************
	## Fourth step:
	## Next steps:
	## - Determine alignment without outliers
        ## - Remove gap only positions
	## - Combine alignments.
	if (1)
	{
	    ############################
            ## Now we can extract the outlier sequences from the original alignment
	    ############################
	    ## File with complete alignment:    $fasta_file_to_extract_from  ## defined in surrounding block
	    ## File with outliers:              $input_outliers_filename     ## defined in surrounding block
	    ## New file name:
	    my $new_alignment_name 
		= "alignments-without-outliers/$gene.without_outlier.fas";
	    my $new_alignment_name2
		= "alignments-without-outliers/$gene.without_outlier_core_at_top.fas";
	    
	    print "1 $fasta_file_to_extract_from\n";
	    print "2 $input_outliers_filename\n";
	    print "3 $new_alignment_name2\n";

	    system("split-fasta-with-sequence-list $fasta_file_to_extract_from $input_outliers_filename");

	    ## Remove gap only positions in the alignment without outliers.
	    system("remove_gap_only_positions-v1.0 not_in.fas > $new_alignment_name");
	    system("rm -f not_in.fas in.fas");

	    print("move_core_taxa_to_top_of_alignment.pl  $new_alignment_name > $new_alignment_name2\n");
            system("move_core_taxa_to_top_of_alignment.pl  $new_alignment_name > $new_alignment_name2");


######################

# Posponed:
#	    ## Now we should realign lower case regions in the alignment
#	    ## without outliers.


############################################
## Combine - rerun from here:
############################################
	    ## Zusammenstricken:

	    ## Name of alignment file we want to use: $new_alignment_name2
	    ## Name of alignment file we want to work with: (initially a copy of $new_alignment_name2)
	    my $iterate_alignment_name = "mafft-refined-alignments-combined-aligned/iter_${gene}.fas";

	    print "Iterate name: $iterate_alignment_name\n";

	    print( "cp $new_alignment_name2 $iterate_alignment_name");
	    system("cp $new_alignment_name2 $iterate_alignment_name");

	    ## For each gene we come here:
	    ## Now for each outlier in this gene we add the outlier to the
            ## $iterate_alignment_name alignmant.

	    my $outl;
	    my $outl_fname;

	    foreach $outl (@outliers)
	    {
		$outl_fname = $outl;
		$outl_fname =~ tr/|/_/;

		my $coreAlig_plus_seq  = "alignment_cores_plus_outlier/" . $gene . ".aa_cores_plus_".$outl_fname.".fas";

		my $comb_command = "combine-align $iterate_alignment_name $coreAlig_plus_seq > tmp.fas"; 

		print "--> $comb_command\n";

		system($comb_command);

		system("mv tmp.fas $iterate_alignment_name");

	    } ## END foreach outlier
	} ## END this step - alignment without outlier and without gap only positions, combine


	## - Refine unaligned regions
	if (0000)
	{
    	    my $iterate_alignment_name = "mafft-refined-alignments-combined-aligned/iter_${gene}.fas";

	    if (-e  $iterate_alignment_name)
	    {
		print "File exists:  $iterate_alignment_name\n";
	    }
	    else
	    {
		print "ERROR: File does not exist: $iterate_alignment_name\n";
		exit(-1);
	    }

## Now refine unaligned regions:

	    my $log_out = "realign_subalignment.log_out";
	    my $log_err = "realign_subalignment.log_err";

	    open (F_OUT, ">>$log_out");
	    print F_OUT "Gene: $gene\n";
	    close(F_OUT);
		
	    open (F_ERR, ">>$log_err");
	    print F_ERR "Gene: $gene\n";
	    close(F_ERR);

	    my $command_refine = "realign_subalignment range $iterate_alignment_name 1>> $log_out 2>> $log_err";
	    print ("$command_refine\n");
	    system($command_refine);

	    my $move_file = $iterate_alignment_name;
	    $move_file =~ s/\.fas/_real_ranges.fas/;
	    system("mv $move_file mafft-refined-alignments-realigned-ranges");
	}
##	exit(0); ## Stop after first gene
	         ## - test funcionality without running through all genes


    } ## END Outliers exist for this gene
    else ## - there are no outliers for this gene
    {
	
    }





} ## END foreach $gene (@genes)


# my $alig_command1;
# my $alig_command2;
# my $alig_command3;
# my $alig_command4;
# my $view_command;

# foreach my $f (@files)
# {
#    $alig_command1    =  "mafft-linsi --add         $f  EOG5X69Q3.fa > ${f}_mafft_linsi_add.fas";
#    $alig_command2    =  "mafft-linsi --addfragment $f  EOG5X69Q3.fa > ${f}_mafft_linsi_addfragment.fas"; 
#    $alig_command3    =  "mafft-linsi --seed        $f  EOG5X69Q3.fa > ${f}_mafft_linsi_seed.fas"; 
#    $alig_command4    =  "mafft-linsi --addprofile  $f  EOG5X69Q3.fa > ${f}_mafft_linsi_addprofile.fas"; 
#    $view_command     =  "seaview ${f}_mafft_linsi_add.fas";
# #   system($alig_command1);
# #   system($alig_command2);
# #   system($alig_command3);
# #   system($alig_command4);
#    system($view_command);

