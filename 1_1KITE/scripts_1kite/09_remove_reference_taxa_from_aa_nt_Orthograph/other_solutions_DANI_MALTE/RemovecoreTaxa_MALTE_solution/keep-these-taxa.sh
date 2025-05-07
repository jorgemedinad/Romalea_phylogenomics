#!/bin/bash

# Copyright 2015, Malte Petersen <mptrsen@uni-bonn.de>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function usage {
	echo "Usage: $0 REFERENCE_TAXA_FILE OTHER_TAXA_FILE INPUT_DIR OUTPUT_DIR"
	echo ''
	echo 'Filters input Fasta files with Hamstrad-formatted headers (single-line'
	echo 'sequences!) in INPUT_DIR to contain only sequences from the reference'
	echo 'taxa listed in REFERENCE_TAXA_FILE and the analyzed taxa from'
	echo 'OTHER_TAXA_FILE. Places the result in OUTPUT_DIR.'
	echo ''
	echo 'Note: this script is highly specific and is not suitable for anything'
	echo '      else but the intended purpose.'
	echo 'Note: uses a temporary file named /tmp/taxa-include.txt.'
}

if [ $# -ne 4 ]; then
	usage
	exit 1
fi

INCLUDEFILE=/tmp/taxa-include.txt

# remove pattern file
cat /dev/null > $INCLUDEFILE

# create pattern for reference taxa
while read REFTAXON; do
	echo "^>.\+|$REFTAXON|[^|]\+$" >> $INCLUDEFILE
done < $1

# create pattern for actual taxon
while read TAXON; do
	echo "^>.\+|.\+|$TAXON|.\+$" >> $INCLUDEFILE
done < $2

# filter the input files
for file in $3/*.fa; do
	grep --no-group-separator --after-context=1 --file=$INCLUDEFILE $file > $4/$(basename $file)
done

# remove the pattern file
rm $INCLUDEFILE
