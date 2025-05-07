#!/usr/bin/bash

# uses selectSites.pl to remove gap-only sites from alignment files

usage="Usage: $0 /PATH/TO/SelectSites.pl FASTAFILE FASTAFILE FASTAFILE ..."

# ignore case in string matches
shopt -s nocasematch

# exit if there are no arguments
if [[ $# -eq 0 ]]; then echo $usage; exit; fi

# the first argument must be the path to the SelectSites.pl script
select=$1
if [[ ! $select =~ SelectSites.pl$ ]]; then 
	echo $usage;
	exit 1;
fi
# and the script must exist
if [[ ! -f $select ]]; then
	echo "${select}: no such file or directory"
	exit 1;
fi

# truncate the arguments list so there are only fasta files
shift

# ask for confirmation from the user
read -r -p "Going to use $select to remove all gap-only sites from alignment files. OK? (y/n) "
if [[ ! $REPLY =~ (y|yes) ]]; then exit; fi

# for every file in the argument list
for file in $@; do
	# save the directory name in a variable
	dirn=$(dirname $file)
	# generate a new file name, appending ".gapsremoved"
	new_fn=${file}.gapsremoved;
	# run SelectSites.pl on the file
	perl $select -s '1-' -x 1 $file > $new_fn
done
