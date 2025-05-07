#!/usr/bin/bash

# uses PAL2NAL to generate nucleotide alignments for all files in aa directory, 
# using files in nt directory and placing results in output directory
# the aa input files must have the suffix ".linsi.aa.fas",
# the nt input files must have the suffix ".nt.fas"
# the nt alignment output files will have the suffix ".linsi.nt.fas"
# the error files will have the suffix ".linsi.aa.fas.error"
usage="Usage: $0 /PATH/TO/pal2nal.pl AA_DIRECTORY NT_DIRECTORY OUTPUT_DIRECTORY"

# exit unless there are exactly 4 arguments
if [[ $# -ne 4 ]]; then echo $usage; exit; fi

# the first argument must be the path to the PAL2NAL script
pal2nal=$1
if [[ ! $pal2nal =~ pal2nal.pl$ ]]; then 
	echo $usage;
	exit 1;
fi
# and the script must exist
if [[ ! -f $pal2nal ]]; then
	echo "${pal2nal}: No such file or directory"
	exit 1;
fi

# exit unless all directories exist and are dirs
if [[ ! -d $2 ]] ; then echo "Not a directory: $2"; exit; fi
if [[ ! -d $3 ]] ; then echo "Not a directory: $3"; exit; fi
if [[ ! -d $4 ]] ; then echo "Not a directory: $4"; exit; fi

# ask for confirmation from the user
read -r -p "Going to PAL2NAL all files in $2 using nucleotide data from $3, placing results in $4. OK? (y/n) "
if [[ ! $REPLY =~ ^y$ ]]; then exit; fi

# for every file in the input dir
for file in $2/*; do
	# generate nt and output filenames
	aaf=$(basename $file)
	ntf=$3/${aaf%linsi.aa.fas}nt.fas
	outf=$4/${aaf%nt.fas}linsi.nt.fas
	errf=$4/${aaf%nt.fas}.error

	# does the nt file exist? if not, skip this aa file
	if [[ ! -f $ntf ]]; then
		echo "nt file for $aaf does not exist in $ntf!";
		continue;
	fi

	# run PAL2NAL on the aa and nt files, placing results into output dir
	echo "Running $pal2nal on $file and $ntf, results in $outf"
	perl $pal2nal $file $ntf -output fasta 1> $outf 2> $errf

	# if there have been errors, say so
	if [[ -s $errf ]]; then 
		echo "Errors running $pal2nal for $aaf and $ntf, look in $errf"
	fi
done

