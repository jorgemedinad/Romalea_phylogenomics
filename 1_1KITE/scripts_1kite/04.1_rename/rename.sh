#!/bin/bash

# renames all files in input directory,
# replacing suffix with new suffix

usage="Usage: $0 INPUT_DIRECTORY suffix new_suffix"

# exit unless we have exactly 3 arguments
if [[ $# -ne 3 ]]; then echo $usage; exit; fi

# exit unless the target is a directory
if [[ ! -d $1 ]] ; then echo "Not a directory: $1"; exit; fi

# prompt for y or n
read -r -p "Going to rename all files in $1, changing '$2' to '$3'. OK? (y/n) "
if [[ ! $REPLY =~ ^y$ ]]; then exit; fi

# for every file in the target directory
for file in $1/*$2*; do
	# generate a new file name, replacing the keyword ($2) with the replacement ($3)
	new_fn=$(basename $file | perl -pne "s/\Q$2\E/$3/") 
	# rename the file
	mv -v $file $1/$new_fn
done
