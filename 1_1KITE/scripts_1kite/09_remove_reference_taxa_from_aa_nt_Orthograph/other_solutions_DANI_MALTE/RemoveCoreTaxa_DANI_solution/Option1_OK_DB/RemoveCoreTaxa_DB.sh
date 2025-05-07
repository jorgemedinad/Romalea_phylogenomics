#!/bin/bash

#usage:
#sh RemoveCoreTaxa.sh header2remove.txt

#script looks if:
#4th collumn is empty, 
#and if the second column starts with Ixodes,or Tribolium or ... (you have to change here the names for other core taxa)
#ATTENTION: this works only with core taxa (they have 2 pipes instead of 3)
# and writes every header in the text file "header2remove.txt"

#NOTE 
	#the script must be in the directory, where you fasta files are.
	# the extention of you fasta files must have .fas (or change line 18)


DIR="$( cd "$( dirname "$0" )" && pwd )"
FASTA=$DIR/EOG*fas							    #change here the extension of you files

awk -F "|" '{if ( $4=="" && $2 ~ /Tribolium_castaneum/ ) { print $0 } } ' $FASTA >>header2remove.txt
awk -F "|" '{if ( $4=="" && $2 ~ /Acromyrmex_echinatior/ ) { print $0 } } ' $FASTA >>header2remove.txt
awk -F "|" '{if ( $4=="" && $2 ~ /Acyrthosiphon_pisum/ ) { print $0 } } ' $FASTA >>header2remove.txt
awk -F "|" '{if ( $4=="" && $2 ~ /Pediculus_humanus/ ) { print $0 } } ' $FASTA >>header2remove.txt

