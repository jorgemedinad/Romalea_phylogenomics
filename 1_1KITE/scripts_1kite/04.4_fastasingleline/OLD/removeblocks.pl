#!/usr/bin/perl
use strict;
use warnings;

use Seqload::Fasta;

while (my $gen = <*.fas>) {
	my $datei = Seqload::Fasta->open($gen);
		open (my $gen_new, ">", "$gen.2")
			or die "$!\n";
				
					while (my @headseq = $datei->next_seq) 
						{
								print $gen_new ">".$headseq[0]."\n".$headseq[1]."\n";			
									}
										$datei->close; 	
									}
