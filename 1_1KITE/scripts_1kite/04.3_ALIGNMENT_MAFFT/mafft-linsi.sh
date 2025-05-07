#!/bin/bash

# uses mafft-linsi to align all files in input directory,
# placing results in output directory
usage="Usage: $0 INPUT_DIRECTORY OUTPUT_DIRECTORY"

# exit unless there are exactly 2 arguments
if [[ $# -ne 2 ]]; then echo $usage; exit; fi

# exit unless both directories exist and are dirs
if [[ ! -d $1 ]] ; then echo "Not a directory: $1"; exit; fi
if [[ ! -d $2 ]] ; then echo "Not a directory: $2"; exit; fi

# ask for confirmation from the user
read -r -p "Going to align all files in $1, placing results in $2. OK? (y/n) "
if [[ ! $REPLY =~ ^y$ ]]; then exit; fi

# for every file in the input dir
for file in $1/*; do
	# generate a new filename
	new_fn=$(basename $file .fa).linsi.fas
	# and align the input file, placing the result in the output directory
	echo "Aligning $file into $2/$new_fn"
	mafft-linsi --quiet $file > $2/$new_fn
done
