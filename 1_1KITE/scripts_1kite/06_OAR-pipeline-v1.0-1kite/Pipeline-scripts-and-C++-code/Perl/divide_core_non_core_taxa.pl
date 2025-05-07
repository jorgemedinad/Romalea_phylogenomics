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

#!/usr/bin/perl -w
use strict;

my $argc = @ARGV;

# Input: One fasta file.
# Output: One fasta file, which includes only those sequences with
#         a specific number of '|' in the sequence name. 

if ($argc < 1 || $argc > 1)
{
    printf STDERR "Usage: divide_core_non_core_taxa.pl fasta-file\n";
    exit(0);
}

my $filename_fas        = $ARGV[0];

my $outf_cores    = $filename_fas;
my $outf_noncores = $filename_fas;

if (! ($filename_fas =~ m/\.fas/))
{
    print "Perl skript requires .fas files.\n";
    exit(-3);
}

$outf_cores    =~ s/\.fas/_cores.fas/;
$outf_noncores =~ s/\.fas/_noncores.fas/;

## print "$outf_in\n";
## print "$outf_notin\n";


open(FI, $filename_fas) or die "Can't open file $filename_fas \n\n";

my $N;
my $in_bool;

my @core_taxa_alig;
my @other_taxa_alig;

while(<FI>)
{
    if (/>/) 	# Found header of next data set
    {
	$in_bool = 0;
    	$N = tr/|/|/;
    	
    	if ($N == 2) ## Core taxon
    	{
	    push(@core_taxa_alig, "$_");
	    $in_bool = 1;
    	}
	else ## Non core taxon
	{
	    push(@other_taxa_alig, "$_");
	}
    }
    else
    {
	if ($in_bool == 1)
	{
	    push(@core_taxa_alig, "$_");
	}
	else
	{
	    push(@other_taxa_alig, "$_");
	}
    }

}

open (CORES, "> $outf_cores");
print CORES (@core_taxa_alig);
close CORES;

open (NONCORES, "> $outf_noncores");
print NONCORES (@other_taxa_alig);
close NONCORES;

