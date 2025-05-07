#################################
## AHE data analysis - upstream##
#################################
# By Xuankun Li (xli8@memphis.edu)
# Edited and tested by Austin Baker (ajbaker6@memphis.edu)
# Edited for TAMU Grace Cluster by Hojun Song (hsong@tamu.edu)
# From raw reads to store/grab cleaned orthologous
# Remember to modify directories each time
# TAMU edition heavily relies on loading appropriate modules. Module loading commands (module load xxxx) are included in the bash files. Modules require different dependent modules. 
# If in the future the scripts do not work, check module dependencies.
# Last updated (8/17/2024)

0. Getting into the TAMU Cluster

# First way is to get into the cluster is by using a terminal. Open a terminal and type the following command.
ssh netid@grace.tamu.edu # (e.g. hsong@grace.tamu.edu)
## This will connect you to Grace after entering the password and Duo authentification.
## Navigating the files and folders using terminal is quite effective if you become familiar with how to use certain commands.
## Here are some basic commands.

# This command is used for changing directories. For example, when you go into the terminal, you are in your home directory. 
cd directory_name
# The first thing to do would be change the directory to OR_TE. One useful tip instead of typing out the whole file name would be to type the first few letters and press tab button on your keyboard. 
# The terminal will fill in the rest.
cd /scratch/group/songlab/phylo/OR_TE/
# Once you are in the directory, you can start navigating by listing the files and folders. ls command will list the file names only
ls
# This command will show a detailed list of files with permission status, the owner, date of last edit, file names, and folder names.
ls -l
# To go back to the folder above,
cd ..
# Check your job status
squeue -u netid
# To read the file (usually txt file)
cat filename
# To edit a file
nano filename

## Second way is to use the TAMU Grace OnDemand, which is a web-interface. Open a web browser and connect to: https://portal-grace.hprc.tamu.edu/pun/sys/dashboard. 
## After authentification, you can check files and even run the job. This is often slow.
## Grace clusters provide two directories: /home/netid (e.g. /home/hsong) and /scratch/user/netid (e.g. /scratch/user/hsong). For most, the analysis should be done using the scratch directories.
## We have also created a shared disk with a maximum of 50TB. Now all of the bioinformatics pipeline for this analysis has been moved to this new disk. (/scratch/group/songlab/phylo/OR_TE)
## File transfer can be done using terminal, OnDemand, or Globus.
## The basic information about how to use TAMU Grace cluster can be found at: https://hprc.tamu.edu/wiki/Grace 

1. From raw data to assembly
## We only need to do this whole set up just once.
## Prepare directories and copy over necessary software within your main directory
## All of the files are the directory '/scratch/group/songlab/phylo/OR_TE'
## Create directory called 'environment'. This directory contains programs and scripts required for different pipelines.

mkdir environment
# within environment make a directory called 'trimmomatic'
cd environment
mkdir trimmomatic
cd trimmomatic
# place an adaptor file into this folder (TruSeq3-PE.fa)
# The adapter file looks like the following
######################################
# >PrefixPE/1                        #
# TACACTCTTTCCCTACACGACGCTCTTCCGATCT #
# >PrefixPE/2                        # 
# GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT #
######################################
# In TAMU Grace, trimmomatic works as a module. Make sure '02_AHE_raw_to_assembly_paralleled_TAMU.sh' file includes the following line in the beginning.
module purge
module load Trimmomatic/0.39-Java-11 

# Go back to 'environment'
cd /scratch/group/songlab/phylo/OR_TE/environment/
# also make directory 'software'
mkdir software
# in this directory, we will create a folder, 'Orthograph'
cd software
mkdir Orthograph
# We will install Orthograph in the later stage.

# Another program that we will use is SOAPdenovo2, which is loaded up at a module in the bash file '02_AHE_raw_to_assembly_paralleled_TAMU.sh'. 
# This program requires two dependencies, so make sure that the bash file includes the following modules.
module load GCC/8.3.0
module load OpenMPI/3.1.4
module load SOAPdenovo2/r242

# Go back to 'OR_TE'
cd /scratch/group/songlab/phylo/OR_TE
# Create "raw data" folder to keep all raw data (fq.gz files) in one place. Need to do this only once and keep reusing.
mkdir raw_data
# Within this raw data folder, create each plate-specific data folder. Replace "raw_data_plate_xx" with a specific directory name each time we receive a new sequencing plate result
cd raw_data
mkdir raw_data_plate_xx
# As we add more raw data, this raw data folder will include several subfolders from different sequencing plates.
cd raw_data_plate_xx
# Transfer or copy over raw data into the specific raw_data_plate_xx using wget, Globus, or other file transfer programs.

## Inside each raw_data folder, it is absolutely important that the files are organized into folders. Each folder should be named with sample ID (usually TAMUIC_IGC_XXXXXX). 
## Inside each sample folder, you will have two fq.gz files and sometimes MD5.txt file.
## It is important that these files end in fq.gz, not fastq.gz.
## Please spend time to make sure that all the raw files are correctly organized before launching any pipelines.

# Go back to 'OR_TE'
cd /scratch/group/songlab/phylo/OR_TE
# Create a directory called 'assembly' to perform all filtering and assembly. Need to do this only once and keep reusing.
mkdir assembly
cd assembly

## Within this 'assembly' folder, we will create three separate folders, which will be used in various stages of the pipeline. 

# First folder to create is the 'combined_sample_name' folder (which has already been created, so no need to redo). This is a folder where we keep we keep various 'sample_names_xxx.txt' files. 
# These are simple text files that will be used for the downstream analyses (Step 3 below). We can edit files and add as many files. Make sure each file is specific to the analysis.
#  Prior to running downstream analyses, make sure you create a filed called 'sample_names_XXX.txt' which contains a exact subset of data used for the analysis. We can use a text editor to create this file. 
mkdir combined_sample_name


# The second folder to create is the 'read_data_plate_xx' folder. This is where the unzipped fastaq files will be placed and where the results of the trimmomatic filtering and SOAPdenovo2 assembly will be stored.
# Create a sequencing project specific folder within "assembly", replace 'read_data_plate_xx' with a new directory name. This will essentially correspond with the plate number of the 'raw_data_plate_xx'. 
mkdir read_data_plate_xx

# The third folder we will create is the 'script_01_and_02' folder, which contains the templates for 01 and 02 script used for trimmomatic and SOAPdenovo2. 
# Since these files need to be frequently updated for each analysis, we keep copies of unmodified templates in this folder and we will copy these over to appropriate 'read_data_plate_xx' folder whenever 
# we process newly sequenced samples. These files should not be modified.
mkdir script_01_and_02
cd script_01_and_02
# In this folder, place '01_soap_configure.py' and '02_AHE_raw_to_assembly_paralleled_TAMU.sh'.

## OK. Now you are ready to start your pipeline.

# Copy these two useful scripts into the directory 'read_data_plate_xx'. You can also use the OnDemand web browser to do this.
cd /scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx
cp /scratch/group/songlab/phylo/OR_TE/assembly/script_01_and_02/01_soap_configure.py .
cp /scratch/group/songlab/phylo/OR_TE/assembly/script_01_and_02/02_AHE_raw_to_assembly_paralleled_TAMU.sh .

## Before launching the assembly, you need a file that specifies the samples to process. We will call this file 'sample_names.txt'. This file is a simple text file that can be created in a text editor,
## and should contain a list of sample names that match exactly to the sample folder names in your 'raw_data_plate_xx'. It is important that the names match perfectly. Make sure to use underline (_) instead of a dash (-).
## For example, if you are doing an assembly on 100 samples, the "sample_names.txt" should contain 100 lines matching exactly the names of the samples. After creating this file, upload in the 'read_data_plate_xx' folder.
## Now, the newly created 'read_data_plate_xx' folder should include three files: 'sample_names.txt', '01_soap_configure.py', and '02_AHE_raw_to_assembly_paralleled_TAMU.sh'. For each new analysis, we will need to go through the following process of modifying files.

# Modify file '01_soap_configure.py' line 17, change: read_dir = "/scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx" to reflect the new read_data_plate_xx or to your own directory, e.g "/scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx"
# You can use nano command in terminal to edit the file or the Edit function is the Ondemand browser.
# Modify file '02_AHE_raw_to_assembly_paralleled_TAMU.sh' line 7, change the number of simultaneous jobs (numbers after array=) to call: #SBATCH --array=1-96 (if necessary, based on the number of samples in sample_names.txt, this example uses 6 samples so change to 1-6)
# If your the raw data files are extra large, modify fie line 6 to change the time appropriately (maybe to 8 hours or 12 hours) 
# Modify file '02_AHE_raw_to_assembly_paralleled_TAMU.sh' lines 19-23, change according directory path to:
raw_dir="/scratch/group/songlab/phylo/OR_TE/raw_data/raw_data_plate_xx" #location of raw data - make sure to change to your own directory 
read_dir="/scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx" #location of the fastq reads - make sure to change to your own directory
adaptor_file="/scratch/group/songlab/phylo/OR_TE/environment/trimmomatic/TruSeq3-PE.fa" #location of the adaptor file for trimmomatic  - make sure to change to your own directory 

## After making sure all these files are modified, you are ready to launch the assembly pipeline.

## Run 01_soap_configure.py - This python scrip creates configuration files for trimmomatic and SOAPdenovo2 analysis.
# in your terminal, change directory to:
cd /scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx
# Replace the directory to the appropriate one
# Make sure that the file called "sample_names.txt" is placed in this folder, which contains the correct list of samples.
chmod +x 01_soap_configure.py
# This chmod changes the permission to make sure that the python script will work.
./01_soap_configure.py sample_names.txt

# This will generate *.configure files for SOAP analysis in the same folder. To see if it worked, type ls -l in the terminal, which will show the list of newly created configure files.

## Submit paralleled jobs for assembly
sbatch 02_AHE_raw_to_assembly_paralleled_TAMU.sh
# generate: directory 'contig' with 'graph_prefix_*.contig', all assembly files

## This pipeline will probably take several hours per sample to complete. When a job a properly completed, the pipeline will create a folder called 'contig' which contains a list of 'graph_prefix_samplename.contig' files.
## Normally, each contig file should be a few MBs to tens of MBs. This contig file is a fasta file that can be opened in a text editor.
## The pipeline will also create folders that bear the sample names. A properly executed run will include 31 files in each folder.
## the pipeline will also create many slurm-XXXX.out files that are approximate 1.42KB, which are log files. A successful run will have a slurm output file that says 'TrimmomaticPE: Completed successfully.


(1.5). Download and set up Orthograph (only need to do once)
## Download Orthograph directly from github to the cluster into software directory
cd environment/software
module load git/2.31.1
git clone https://github.com/mptrsen/Orthograph.git
## This is already done and the program is now installed properly in the cluster. Our TAMU specific script calls on the program installed in the module.

## The following steps are unnecessary, but we keep as a legacy. START HERE

## Copy across reference genomes and the cluster of orthologous groups (COGs) file
cd Orthograph/
cp -r /home/sshin4/Orthograph/AHE1/OGS . # OGS folder including three genomes used in 1KITE Coleoptera
cp /home/sshin4/Orthograph/AHE1/3AHE.txt ./OGS/ # COGs for AHE dataset, seems a OrthoDB 7 file, a tab-delimited file
# for above files I always get from lab's previous researcher while I traveled in different 1KITE labs
# I will document how to make those files (especially the COGs) if I got a chance to figure that out

## Create the required database structure (generate a *.sqlite file)
# Reformat COGs
cd OGS
cut -f1,3,4 3AHE.txt > Cole_AHE.table

# Create the required database structure
# check and modify 'peptide_set_for_coleoptera_new.sh'
cp /home/xli8/Orthograph/peptide_set_for_coleoptera_new.sh .
bash peptide_set_for_coleoptera_new.sh
# After submit the job, it will print: 
#...
#Any existing Orthograph table structures with that prefix will be erased. Are you sure (y/n)? y # enter y
#...
#...
#Enter the set name (required; ASCII only, no commas!): AHE # choose whatever name you like, will be the name for *.sqlite file, I used AHE here
#Enter a description for the set (optional but recommended). IMPORTANT: Do not use any commas (,): # I did not put anything here
#...
#Enter the OGS ID for Danaus plexippus or 0 if you don't want to use this OGS: 1 # enter 1
#Enter the OGS ID for Nasonia vitripennis or 0 if you don't want to use this OGS: 2 # enter 2
#Enter the OGS ID for Tribolium castaneum or 0 if you don't want to use this OGS: 3 # enter 3
#...
#Database now ready to run Orthograph.
# generate: AHE.sqlite for next step

## The above steps are unnecessary, but we keep as a legacy. END HERE

## Our Orthoptera specific data are in "6K_set_TAM" folder. To see
cd /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/6K_set_TAM
cd OGS
# we have a SQlite database we created for our project, called "TAMU1.sqlite"
cd /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/6K_set_TAM
cd reference_sequences
# in this folder we have four insect reference genomes used for building our database
cd /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/6K_set_TAM
# We also have 6K_Ortho.txt, which is the reference gene dataset for our workflow.

## The following steps are unnecessary because Orthograph runs well in the cluster, but we keep as a legacy. START HERE
## Get one test run working
# Download and modify orthograph.conf in a text editor

#--------------------------------------------------
# # Important stuff first. These settings are mandatory.
#-------------------------------------------------- 
database-backend   = sqlite
sqlite-database    = /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/AHE.sqlite # path to AHE.sqlite generated in last step
sqlite-program     = sqlite3
input-file         = /home/ajbaker6/contig/graph_prefix_DDM0017.contig # path to one test data
species-name       = DDM0017 # I kept all species name as DDM number
ortholog-set       = AHE # name of ortholog set
output-directory   = /home/ajbaker6/contig # path to output directory
#--------------------------------------------------
# Options. These settings are optional, but may be required if the defaults don't work
# 
#-------------------------------------------------- 
alignment-program    = /public/apps/mafft/7.407/bin/mafft --localpair --maxiterate 1000 --anysymbol
hmmbuild-program     = /public/apps/hmmer/3.2.1/bin/hmmbuild
makeblastdb-program  = /public/apps/blast/2.7.1/bin/makeblastdb
translate-program    = /public/apps/exonerate/2.2.0/bin/fastatranslate
hmmsearch-program    = /public/apps/hmmer/3.2.1/bin/hmmsearch
blast-program        = /public/apps/blast/2.7.1/bin/blastp
exonerate-program    = /public/apps/exonerate/2.2.0/bin/exonerate
# This option may be hidden by a #, in which case delete the #
num-threads                = 24
#--------------------------------------------------
# overwrite original orthograph.conf file with new one in Filezilla
cd ..
cp /home/sshin4/Orthograph/AHE1/Test.sh .
# open Test.sh file in text editor and add # to beginning of line 11, change email address to your own on line 6
# overwrite original Test.sh file in Filezilla
sbatch Test.sh
# this job will be killed because of 'Exceeded job memory limit'. 
# But based on my experience, it is ok, it will make files in ./sets/AHE
# May need to run again after first attempt if it fails

## The above steps are unnecessary because Orthograph runs well in the cluster, but we keep as a legacy. END HERE



2. Run Orthograph
## Starting directory: /scratch/user/hsong/OR_TE
## The idea of this step is to submit paralleled jobs, and there will be three steps.
## This is a bioinformatics pipeline called Orthograph, which was developed as part of the 1KITE project.
## Petersen, M. et al. Orthograph: a versatile tool for mapping coding nucleotide sequences to clusters of orthologous genes. BMC Bioinform. 18, 111, doi:10.1186/s12859-017-1529-8 (2017).
## This pipeline takes the assembly results and finds orthologs, and remove paralogs and duplicates.

## 1) Put orthograph.slurm file in directory 'bin'  ##### STEP 1 ALREADY DONE
# create a directory called 'bin' in root directory
mkdir bin ## Path = /scratch/group/songlab/phylo/OR_TE
cp /scratch/group/songlab/phylo/OR_TE/bin/orthograph.slurm
# Modify file 'orthograph.slurm' lines 11 and 12, change path to your Orthograph command files
time /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph-analyzer -c $WORKDIR/$TAXON/orthograph.conf >> ../$TAXON.ortho_analy.log 2>&1 
time /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph-reporter -c $WORKDIR/$TAXON/orthograph.conf >> ../$TAXON.ortho_report.log 2>&1
# In order to run this slurm file in TAMU Grace, we have included module load commands inside the file.
module purge
module load GCC/10.2.0  OpenMPI/4.0.5 Perl/5.32.0 Orthograph/0.7.1
## These steps are already done, so no need to repeat.

## 2) Prepare directories and configure files for submission
cd /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph

## create a working directory and put contig files in there. Each orthograph analysis will run within the Orthograph directory. In order to organize the data effectively, we will create a series of folders. 

## First we create a folder that keeps all the raw analysis results.
mkdir orthograph_results
cd orthograph_results
# Within this folder (which is already created and no need to redo), we will create plate-specific folders that correspond with raw and read data.
mkdir orthograph_plate_xx
cd orthograph_plate_xx

## In this newly created folder, you will need the following three configuration files.
## The template conf files are located in the '/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/configure_template'
# Copy 'orthograph_conf_prep.sh' and 'orthograph.conf' and 'orthograph.sh' in directory '/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx' from the 'configure_template' folder.
cd /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx
cp /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/configure_template/orthograph_conf_prep.sh .
cp /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/configure_template/orthograph.conf .
cp /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/configure_template/orthograph.sh .
# Modify file 'orthograph.conf' line 42 using nano or OnDemand Editor, change to your directory, 
# for example: input-file         = /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx/XYZ/graph_prefix_XYZ.contig 
# Modify file 'orthograph.conf' line 56, change to your directory, 
# for example: output-directory   = /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx/XYZ
## You need to modify 'orthograph.conf' each time you launch a new Orthograph pipeline.
## You do not need to change 'orthograph_conf_prep.sh' file
## If your assembly files ('graph_prefix_sample_name.contig') are large, you may want to increase the walltime in 'orthograph.sh'
# Modify file 'orthograph.sh' line 3 to 6 hours or 8 hours. This is usually not necessary.

## Your newly created "orthograph_plate_xx" folder has now 3 configuration files, and you now need to copy over the results of the assembly.
## At the end of the SOAPdenovo2 analysis, the program will create a folder called 'contig' inside 'read_data_plate_xx' which will contain 'graph_prefix_sample_name.contig' 
## (where sample_name corresponds to the same name in the 'sample_name.txt' file). Now we need to copy over these files to the newly created 'orthograph_plate_xx' folder. 
cp /scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx/contig/* .
# This copy command will copy everything from the contig folder. But if you want to copy over specific contig files, use the OnDemand Copy function from the appropriate source.

## Although it will be great to be able to launch hundreds orthograph runs in parallel, the GRACE cluster cannot actually handle that many at once. You can reliably include up to 50 contig files
## in the folder and run without failure.

## Prepare directories and files for Orthograph
cd /scratch/group/songlab/phylo/OR_TEE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx
bash orthograph_conf_prep.sh
# This bash script will generate one directory for each sample, including the according contig file and the configure file. You can check what the script had done by typing ls -l.

## 3) Submit paralleled jobs for Orthograph
bash orthograph.sh
# it will submit one job for each sample to the cluster

## You may see 'sbatch: error: Batch job submission failed: Socket timed out on send/recv operation' error but that's ok. Just press 'q'. The batch is running.
## use 'control+c' to cancel the command when submission finished
## Once the jobs are submitted, check the status in terminal by using 'squeue -u netid' to see if the job is working ok. Make sure to change netid to your netid (jorgemedinad)
## Each orthograph run will take several hours. Once completed, the pipeline will create directories with sample names, which contains 'aa' and 'nt' directories and 5 additional files. 
## Also, 4 log files are created per run. If everything goes well, 'samplename.ortho_analy.log' file should be about 642KB, and 'samplename.out.log' file should be 0 bytes.
## If these files are of different sizes, your run might have failed.

## Because you have generated many log files, we will clean results
mkdir log_and_others
mv *.log log_and_others/
mv orth* log_and_others/
# These commands will create a new folder called log_and_others, and move all log files and configuration files to this folder. A good practice for file management.

## Assessing the Orthograph results and figure out how many loci have recovered.
## There are actually two ways. One way is by simply going to the resulting directory of the orthograph output and open up 'nt' or 'aa' folder and count the number of files.
# Go to the result folder, in this case, 'aa'
cd /scratch/group/songlab/phylo/OR_TEE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx/samplename/aa
# This nifty code counts the number of files in a folder.
ls -1 | wc -l
## The second way is use a script. Jackson created several useful scripts and they are in 'useful_scripts' folder.
# Go to the orthograph result folder.
cd /scratch/group/songlab/phylo/OR_TEE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx
# Copy 'get_size_genes.sh' file from the 'useful_scripts' folder.
cp /scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/useful_scripts/get_size_genes.sh .
# Modify 'get_size_genes.sh' line 4 and 5 (and line 25 if necessary) accordingly.
# For example, update path_to_contig="/scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate_xx/contig"
# And update path_to_ortho_results="/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx"
# Run the script
bash get_size_genes.sh
## After running the script (very quick), it generates two tsv files (contig_size.tsv, genes_recovered.tsv). You can open these files in text editor or Excel.
## 'genes_recovered.tsv' lists the number of recovered loci per sample.


3. Orthologous cleaning (1KITE pipeline part 1 in an automatic way) ### STEP 3 ALREADY DONE
## Still recommended to read through the 1KITE pipeline
## These pipelines consist of a series of perl scripts and python scripts to process the results from Orthograph. 
cd /scratch/group/songlab/phylo/OR_TE
mkdir 1KITE_cleaning
cd 1KITE_cleaning
# Create a folder '1KITE_example_template'. This folder contains the entire bioinformatics pipeline to handle Orthograph results into a matrix. 
mkdir 1KITE_example_template
# Create a folder 'cleaned_orthologous_storage'. This folder contains the final results from the 1KITE pipeline. The files within are the raw data for phylogenetic analyses.
mkdir cleaned_orthologous_storage
# Each time 03_select_samples_for_tree.sh is run it will create a folder for the dataset specified in the gene_names.txt and tree_samples.txt files. So, always create and update these two text files.
## These three folders have already been created and there is no need to redo.

## OK. Now you have completed your Orthograph pipeline, and you want to create an alignment for downstream phylogenomic analyses.
## Starting directory: /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning
## Preparations
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning
# Create a new folder for the cleaning 'cleaning_plate_xx'. This is usually analysis specific.
mkdir cleaning_plate_xx
# Go to the newly created folder
cd cleaning_plate_xx
# Create 1KITE_example folder
mkdir 1KITE_example
cd 1KITE_example
# Copy the contents of the entire folder '1KITE_example_template' and paste into this folder. 
cp -r /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/1KITE_example_template/* .
cp /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/1KITE_example_template/* .
# Check the files
ls -l

## Change configuration files.
# Go to the folder
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example

# 1) Modify shell script '01_Orthograph/O2H_header_auto.sh' # For 02 RECONVERT Orthograph into HaMStRad format
# Change work_dir location to the current directory path
cd 01_Orthograph
## Modify 'O2H_header_auto.sh' line 4, update 'test_cleaning' to 'cleaning_plate_xx'
## Modify 'O2H_header_auto.sh' line 8, update 'test_cleaning' to 'cleaning_plate_xx'

# 2) Modify perl script '06_OAR-pipeline-v1.0-1kite_results/Perl/refine-alignments.pl' # For 06 OAR-pipeline-v1.0-1kite
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/test_cleaning/1KITE_example/06_OAR-pipeline-v1.0-1kite_results/Perl
## Update test_cleaning to correct folder name.
## line 167:
system("perl /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/test_cleaning/1KITE_example/06_OAR-pipeline-v1.0-1kite_results/Perl/divide_all.pl $ls_res");
## line 321:
system("perl /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/test_cleaning/1KITE_example/06_OAR-pipeline-v1.0-1kite_results/Perl/move_core_taxa_to_top_of_alignment.pl $new_alignment_name > $new_alignment_name2");
 
# 3) Modify perl script '06_OAR-pipeline-v1.0-1kite_results/Perl/divide_all.pl' # For 06 OAR-pipeline-v1.0-1kite
## Update test_cleaning to correct folder name.
## line 74:
system("perl /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/test_cleaning/1KITE_example/06_OAR-pipeline-v1.0-1kite_results/Perl/divide_core_non_core_taxa.pl $fasta");

# 4) Modify shell script '09_remove_ref_taxa_results/RemoveCoreTaxa_DB.sh' # 09 Remove sequences from reference species (aa alignments & nt files)
## change names to APISU, NVITR, PHUMA, RPROL for the Orthoptera probe set
## This step has also been done. Don't change.
awk -F "|" '{if ( $4=="" && $2 ~ /Tribolium_castaneum/ ) { print $0 } } ' $FASTA >>header2remove.txt
awk -F "|" '{if ( $4=="" && $2 ~ /Danaus_plexippus/ ) { print $0 } } ' $FASTA >>header2remove.txt
awk -F "|" '{if ( $4=="" && $2 ~ /Nasonia_vitripennis/ ) { print $0 } } ' $FASTA >>header2remove.txt

# 5) Modify shell script '01_1KITE_pipeline_auto.sh'
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/test_cleaning/1KITE_example/
# line 4: change to your email address
# If you have more than 50 samples to process, change line 8 time to 12 hours. If you have moren than 100 samples, change to 48 hours.
# lines 16-18, change according directory path to:
raw_dir="/scratch/group/songlab/phylo/OR_TE/assembly/combined_sample_name" # Location of the sample names.  Do not modify
ortho_dir="/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_plate_xx" # Location of the orthograph results 
work_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example" # Location of the current working directory
# line 57, change sample_names.txt to your specific sample_name_xxx.txt in '/scratch/group/songlab/phylo/OR_TE/assembly/combined_sample_name'. This text file should contain those samples you want create an alignment for.

## Submit first job
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example
sbatch 01_1KITE_pipeline_auto.sh

## Manually check results
# After the job is completed, check directory '/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example/11_pal2nal', all *.error should be all 0
# If not, maybe because some sequence headers in *.nt.fas are too long, and been cut by mafft in *.linsi.aa.fas
# Compare the sequence headers in *.nt.fas and *.linsi.aa.fas, check the long ones and modify it manually
# Do the command below again:
# only do the following command if there is an error
for i in *.linsi.aa.fas ; do a=`basename $i .linsi.aa.fas`; perl pal2nal.mod.pl $a.linsi.aa.fas $a.nt.fas -output fasta 1>$a.linsi.nt.fas 2>$a.error; done
# generate: *.linsi.nt.fas cleaned in the same way as *.linsi.aa.fas

## Collect and store cleaned orthologous loci
# Modify shell script '02_cleaned_orthologous_storage.sh'
# line 10:
work_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example" # Location of the current working directory
# line 11:
store_dir="/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage" # Location of storage directory ### THIS PATH REMAINS THE SAME
# line 54: (modify the final output file name)
sed "/>/s/$/|${gene}/" ${work_dir}/13a_rewrite_headers_for_Aliscore/aa_headers_renamed/${gene}.linsi.aa.fas >> ${store_dir}/XXX.aa.fasta  ### change to something like Taeniopoda4.aa.fasta. The number refers to the number of the analysis
# line 60: (modify the final output file name)
sed "/>/s/$/|${gene}/" ${work_dir}/13a_rewrite_headers_for_Aliscore/nt_headers_renamed/${gene}.linsi.nt.fas >> ${store_dir}/XXX.nt.fasta

# Submit second job
bash 02_cleaned_orthologous_storage.sh # Take ~2 mins to run, no need to use sbatch
# generate: *.aa.fasta and *.nt.fasta in the directory /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage

4. Post processing of the alignments
## If everything worked, you will see two folders, 'aa_headers_renamed' and 'nt_headers_renamed', in the '13a_rewrite_headers_for_Aliscore' folder.
## These folder contains thousands of alignment files in fasta format. At this point, it is easier to download the entire results and use the program Geneious to further process.
## In Ondemand browser, go to '/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaning_plate_xx/1KITE_example/13a_rewrite_headers_for_Aliscore'.
## Since most of our analyses will use the nucleotide sequences, download the entire folder 'nt_headers_renamed' to your local computer.
## Usually you will download it as a zip file. Unzip the file and open the folder.
## Of the thousands of fasta files, some will have Zero bytes. These are empty files. Sort the files in the folder by size and delete those.
## Some of the remaining files have only a single sequence while others have multiple sequence alignments. 
## Open Geneious Prime and create a project specific folder. Drag all those fasta files into the newly created folder.
## In Geneious, it would be easy to visualize which files have only one taxon. These are not useful for phylogenetic analyses. Delete all those files.
## Sort the alignments and delete any files that have 2 or 3 samples included. You want to keep only those alignments that have at least 4 terminals.
## Now you have a list of alignments that have at least 4 terminals. Here you have some analytical choices to make.
## Which phylogenetic inferences will be used? Concatenation based maximum likelihood or Bayesian? Coalescence-based ASTRAL?
## How much missing data will you allow for your analyses? 
## After you decide on a particular path, you can export the files out of Geneious for downstream analyses.
## Since most phylogenetic programs use phylip format, you can convert the file format and export out as phylip formats.
## In Geneious, select the files you want to export, go to File > Export > To Multiple Files, which will bring up a window. Change file format to Phylip alignment, specify Export to Folder.
## Click OK, and selected relaxed format.
## Now in your designated folder in your computer, you have hundreds or thousands of alignments ready to be analyzed.

5. ASTRAL
## To build individual trees for Astral, we need to run RAxML on each alignment.
## These files should be uploaded to Grace for further data processing. 

# Go to folder
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/alignment
# Create an analysis specific folder
mkdir test_alignment ### e.g., Taeniopoda_alignment
cd test_alignment 
## Upload these files into this folder using OnDemand browser or Globus.
## template_folder has a bash script for RAxML, raxml.sh
# Copy raxml.sh into your alignment folder
cp /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/alignment/template_folder/raxml.sh .
## Create gene_name.txt. This file lists all of the genes in your analysis. The easiest way to create this file is using a text editor. Copy all of the file names from your local
## alignment folder. Each alignment should look something like 'EOG7B0HCX.linsi.nt.phy'. Your text file now should have a long list of these file names. Using find and replace function,
## delete '.linsi.nt.phy' so that you have only the list of gene name. 
## Each gene_name.txt should contain no more than 500 gene names. This is because the Grace cluster does not let us launch more than 500 RAxML runs simultaneously.
## If you have more than 500 genes, simply create multiple gene_name.txt files, named gene_name1.txt, gene_name2.txt, etc.
## Upload those gene_name.txt files into your alignment folder.
## If gene_name.txt files have been created in windows, run the command "dos2unix [file]" to convert txt files from DOS to UNIX system. Otherwise, you will get errors. 
# Change details within raxml.sh
# Modify line 6 to increase time if you have more than 500 terminals.
# Modify line 7 to match the number of genes. If you have 500 genes, change the array to 1-500. This number has to match exactly the number of genes in the gene_name.txt and cannot go over 500.
# Line 7 must be modify every time the number of genes in a given gene_names.txt file changes.
# Modify line 20 to specify the location of alignemnt
# Modify line 24 to include the appropriate gene_name.txt file. This line must be modified every time a new run with different gene_names.txt is done.
# To run the script
sbatch raxml.sh  ### The SLURM configuration request a lot of SU, but the job itself is fast (less than 15 min). Modify time parameter in the script in case you are running low on SUs.
## This script will launch RAxML runs each with 100 bootstrapping with GTRCAT as a model of nucleotide substitution. The script will create one directory per gene, and also copy over the final tree with
## bootstrap support values to a folder called 'tree'. It will also create many slurm files.
## One way to check the status is by counting the number of files in the 'tree' folder. You can launch additional RAxML runs if the number of trees in the 'tree' folder matches with the number of genes input-file
## gene_name.txt. 
cd /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/alignment/test_alignment/tree
ls -1 | wc -l
# Clean slurm files.
mkdir out
mv *.out out/

## OK. Now you have complete the creation of individual gene trees. Download the entire content of the 'tree' folder to your computer.
# Combined individual tree files into a single file, in terminal. Replace "your_combined" with your own tree name. Open terminal in your personal computer and make sure you change your directory to where the tree files are.
awk '{print $0,"\n"}' *.tre > your_combined.tre

# To contract low support branches, in terminal. Replace "your_combined" with your own tree name. This code collapses any node with bootstrap support lower than 10. This requires that you have newick_utils installed input-file
# a conda environment.
module load Anaconda3/2024.02-1
source activate bio
#conda activate newick_utils # I do not need to activate
nw_ed your_combined.tre 'i & b<=10' o > your_combined_BS10.tre
conda deactivate

# Astral. Make sure you have Astral on your computer. The tree file should in in the same folder as Astral program. Replace "your_combined" with your own tree name. 
java -jar -Xmx3000M astral.5.7.8.jar -i your_combined_BS10.tre -o your_combined_BS10_final.tre 2> your_combined_BS10.log

# Name change. Open terminal. Make sure you have your name change file is in the same folder as your tree. You need a tree file (yourtree.tre) which has only the specimen codes. You will a text file format in the following way, named name_change.txt.
# The format of the name_change.txt is as follows. You should have the list matching the species list.
s/specimen_code/whatever_you_want_to_change_to/
#In terminal, type the following. Open Ubuntu in Windows. Go to working directory using format cd /mnt/c/Users/path/to/file
# Do not tip names and do not use parenthesis
sed -f name_change.txt yourtree.tre > yourtree_renamed.tre

6. Other analyses using concatenated data.

## AMAS.py is a python script that can concatenate gene alignments. 
# To concatenate alignments
cd  /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage/test_tree/nt
python /scratch/group/songlab/phylo/OR_TE/environment/python_scripts/AMAS.py concat -f fasta -d dna -i *fasta
cd  /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage/test_tree/aa
python /scratch/group/songlab/phylo/OR_TE/environment/python_scripts/AMAS.py concat -f fasta -d aa -i *fasta

7. Select samples from stored cleaned data
## This step is a bit redundant, and was originally described in this pipeline, but we are not using this method any more.
## Starting directory: /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning
# Create file tree_samples.txt with list of selected sample numbers in the directory '/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage'
# Create file gene_names.txt with list of selected gene names in the directory '/scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage'
# These genes can be selected from the XXX.aa.fasta file using excel
# Modify 03_select_samples_for_tree.sh
# Need to adjust if using data from multiple plates
# Need to modify lines 16 and 17 before running
# Need to modify lines 7 and 8 to your own directory
# Make sure AMAS.py is in python_scripts folder
mkdir /scratch/group/songlab/phylo/OR_TE/environment/python_scripts
cd /scratch/group/songlab/phylo/OR_TE/environment/python_scripts
cp /scratch/group/songlab/phylo/OR_TE/environment/python_scripts/AMAS.py .

# Move 03_select_samples_for_tree.sh to /scratch/group/songlab/phylo/OR_TE/1KITE_cleaning/cleaned_orthologous_storage

# Submit job
bash 03_select_samples_for_tree.sh
# generate: files by gene in directories test_tree/aa/ and test_tree/nt/, and summary.txt generated by AMAS



## End of AHE upstream pipeline ##


#########################################################
#########################################################
## blast based pipeline for rRNAs and mitochondrial genes
## A separate analysis directly after step 1 'From raw data to assembly'
############################################################################################
(1.5'') Setup envirnoment, download and install softwares
## Four directories in directory 'environment':
# Directory 1: exon_capture_pipeline
1) AddOrReplaceReadGroups.jar
2) CreateSequenceDictionary.jar
3) GenomeAnalysisTK.jar
4) picard-tools-1.88.zip
5) TruSeq3-PE-2.fa
6) vcfutils.pl
7) gatk-4.1.3.0.zip
8) sample_names.txt # sample list, including all taxa names need to analysis
9) Sitophilus_oryzae_kx373615_genes.fasta # gene list, including all the reference genes
10) mt_gene_names.txt # gene name list, including all the reference gene names
# Directory 2: python_scripts
1) ambig_counter_with_dummies.py
2) consensus_maker_LT.py
3) makesomethingNotInterleaved.py
4) mask_low_cov2.py
5) pullexons_without_names_LT_2_zDNA_v2.py
# Directory 3: trimmomatic
1) trimmomatic-0.36.jar
2) TruSeq3-PE.fa
# Directory 4: software
1) bbmap
2) tabix (not in the original version, but if the output *_filtered_phylogenes.vcf.gz is empty, need to install this one)
#For load a software, download *.tar.gz, $ tar xvzf file.tar.gz, ls check the files, less README.md (README.txt)
#1) BBMap_36.92.tar.gz, $ cd docs/ $ less readme.txt, after uncompress could directly use the software
#2) For tabix (https://www.jianshu.com/p/07640917764d):
wget https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
tar xjvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6/
make # setup
# Modify exon_capture_pipeline_stack_resource*.sh step 20: add path to tabix-0.2.6/ for bgzip and tabix


############################################################################################
2'' Run blast based pipeline
# 1) Prepare output directory
## Starting directory: /home/xli8/1KITE_cleaning
mkdir blast_1K_15

# 2) Copy sample_names.txt to the directory 'environment/exon_capture_pipeline'
cp /home/xli8/raw_to_orthograph/1K_15_RTA/sample_names.txt /home/xli8/environment/exon_capture_pipeline

# 3) Modify shell script '01_exon_capture_pipeline_stack_resource_xk_MofU_cluster.sh'
line 4: add your email address
line 10: (Based on the number of samples in file 'sample_names.txt')
#SBATCH --array=1-96
line 22: (modify the location of the output directory: you just made in the last step)
out_dir="/home/xli8/1KITE_cleaning/blast_1K_15" #location of the output directory
line 23: (modify the location of the assemblies: the directory of step 1)
assembly_dir="/home/xli8/raw_to_orthograph/1K_15_RTA" #location of assemblies
## Put '01_exon_capture_pipeline_stack_resource_xk_MofU_cluster.sh' into '/home/xli8/environment/exon_capture_pipeline'

# 4) Modify shell script '02_blast_based_orthologous_storage.sh'
line 7: 
store_dir="/home/xli8/1KITE_cleaning/cleaned_orthologous_storage" # Location of cleaned orthologous storage directory
line 14: (modify the final output file name)
sed "s/>/>${taxon}|/g" ${taxon}.phylogenes_masked_ambig.fa >> ${store_dir}/1K_15.blast.nt.fasta
## Put '02_blast_based_orthologous_storage.sh' into /home/xli8/1KITE_cleaning/blast_1K_15

# 5) Submit first job
cd /home/xli8/environment/exon_capture_pipeline
sbath 01_exon_capture_pipeline_stack_resource_xk_MofU_cluster.sh
# generate: files *.phylogenes_masked_ambig.fa in directory '/home/xli8/1KITE_cleaning/blast_1K_15/cov_10_ambig'

# 6) Submit second job when the first one finished
cd /home/xli8/1KITE_cleaning/blast_1K_15
bash 02_blast_based_orthologous_storage.sh
# generate: file '1K_15.blast.nt.fasta' in directory '/home/xli8/1KITE_cleaning/cleaned_orthologous_storage'


