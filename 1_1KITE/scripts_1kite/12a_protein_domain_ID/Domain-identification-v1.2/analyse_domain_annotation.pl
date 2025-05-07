#!/usr/bin/perl

## Run this script with a list of all coordinate files that shall be analysed as parameters.
## A comprehensive analysis can only be obtained if all coordinate files are passed at once.
## Example: analyse_domain_annotation.pl  some-base-name_*
## The best way to call this script in the pipeline is to call: "run-all-stats.sh"

#/***************************************************************************************************
#*  This file is part of the 1kite analysis pipeline. It is
#*  distributed under the following license:
#*  
#*  Copyright (c) 2012-20014 Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
#*  All rights reserved.
#*  
#*  Redistribution and use in source and binary forms, with or without
#*  modification, are permitted provided that the following conditions are met:
#*  1. Redistributions of source code (complete or in parts) must retain
#*     the above copyright notice, this list of conditions and the following disclaimer.
#*  2. Redistributions in binary form must reproduce the above copyright
#*     notice, this list of conditions and the following disclaimer in the
#*     documentation and/or other materials provided with the distribution.
#*  3. All advertising materials mentioning features or any use of this software
#*     e.g. in publications must display the following acknowledgement:
#*     This product/publication includes/makes use of software developed by Christoph Mayer,
#*     Forschungsmuseum Alexander Koenig, Bonn, Germany.
#*  4. Neither the name of the organization nor the
#*     names of its contributors may be used to endorse or promote products
#*     derived from this software without specific prior written permission.
#*  5. If this source code is used in full or in part in programs that analyse data that
#*     is being published, the author Christoph Mayer has to be included in the list of
#*     authors of the publication. This applies to all publications that include data that
#*     is analysed with the help of this source code. If you think this restriction is not justified,
#*     please do not use this code or any program created with this code.
#*      
#*  
#*  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
#*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
#*  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#*  
#*  IMPORTANT (needs to be included, if code is redistributed):
#*  Please not that this license is not compatible with the GNU Public License (GPL)
#*  due to paragraph 3 in the copyright. It is not allowed under any
#*  circumstances to use the code of this software in projects distributed under the GPL.
#*  Furthermore, it is not allowed to redistribute the code in projects which are
#*  distributed under a license which is incompatible with one of the 4 paragraphs above.
#*  
#*  This project makes use of code coming from other projects. What follows is a complete
#*  list of files which make use of external code. Please refer to the copyright within
#*  these files.
#*  
#*  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
#*                                See copyright in tclap/COPYRIGHT file for details.	
#***************************************************************************************************/


use strict;

my $usage = "Usage: $0 <domain result file(s)>\n";

my @argv;
while (@ARGV)
{
    my $arg = shift;
    if ($arg =~ /^-/)
    {
	if ($arg eq "-h")
	{
	    print $usage;
	    exit();
	}	
	else { die $usage }
    }
    else
    {
	push @argv, $arg;
    }
}
push @argv, "-" unless @argv;

my %pfamA_acc_hash = ();
my %pfamB_acc_hash = ();

my %clan_hash = ();

my $anz_gene=0;
my $anz_void=0;
my $anz_pfamA=0;
my $anz_pfamB=0;

my @genelengths;
my @void_lengths;
my @pfamA_lengths;
my @pfamB_lengths;

my $cum_len_gene=0;
my $cum_len_void=0;
my $cum_len_pfamA=0;
my $cum_len_pfamB=0;

my @sp;
my $len;

foreach my $file (@argv)
{
    open DAT, "$file" or die "Couldn't open '$file': $!";
    ++$anz_gene;
    while (<DAT>)
    {
	if (/^#/)
	{
	    next;
	}
	@sp  = split/ /;
	$len = $sp[1]-$sp[0]+1;

	if ($sp[5] eq "void")
	{
	    ++$anz_void;
	    push(@void_lengths, $len);
	    $cum_len_void += $len;
	}
	elsif ($sp[5] eq "Pfam-B")
	{
	    $pfamB_acc_hash{$sp[4]} = 1;
	    ++$anz_pfamB;
	    push(@pfamB_lengths, $len);
	    $cum_len_pfamB += $len;
	}
	else ## PfamA Hit
	{
	    $pfamA_acc_hash{$sp[4]} = 1;
	    if (exists $clan_hash{$sp[6]})
	    {
		++$clan_hash{$sp[6]};
	    }
	    else
	    {
		$clan_hash{$sp[6]} = 1;
	    }
	    ++$anz_pfamA;
	    push(@pfamA_lengths, $len);
	    $cum_len_pfamA += $len;
	}

    } ## while (<DAT>) this gene
    $len = $sp[1];
    push(@genelengths, $len);
    $cum_len_gene += $len;

    close DAT;
}

print "Number of genes: $anz_gene\n";
print "Number of voids: $anz_void\n";
print "Number of pfamA: $anz_pfamA\n";
print "Number of pfamB: $anz_pfamB\n";

my $num_Aacc = scalar keys %pfamA_acc_hash;
my $num_Bacc = scalar keys %pfamB_acc_hash;
my $num_clan = scalar keys %clan_hash;


print "\n";

print "Number of different pfamA acc: $num_Aacc\n"; 
print "Number of different pfamB acc: $num_Bacc\n"; 
print "Number of different clans:     $num_clan\n";

print "\n";

print "Cumulative length of genes:   $cum_len_gene\n";
print "Cumulative length of voids:   $cum_len_void\n";
print "Cumulative length of pfamA:   $cum_len_pfamA\n";
print "Cumulative length of pfamB:   $cum_len_pfamB\n";

print "\n";
my $prop;

if ($cum_len_gene != 0)
{
    $prop = $cum_len_void / $cum_len_gene;
    print "Proportion of void regions:   $prop\n";
    $prop = $cum_len_pfamA / $cum_len_gene;
    print "Proportion of pfamA regions:  $prop\n";
    $prop = $cum_len_pfamB / $cum_len_gene;
    print "Proportion of pfamB regions:  $prop\n";
}

print "\n";

if  ($anz_gene != 0)
{
    $prop = ($anz_pfamA+$anz_pfamB) / $anz_gene;
    print "Mean number of pfamAB per gene: $prop\n";
    $prop = ($anz_pfamA) / $anz_gene;
    print "Mean number of pfamA  per gene: $prop\n";
    $prop = ($anz_pfamB) / $anz_gene;
    print "Mean number of pfamB  per gene: $prop\n";
    $prop = ($anz_void) / $anz_gene;
    print "Mean number of voids  per gene: $prop\n";
}
print "\n";

if ($anz_pfamA+$anz_pfamB)
{
    $prop = ($cum_len_pfamA+$cum_len_pfamB) / ($anz_pfamA+$anz_pfamB);
}
else
{
    $prop = 0;
}
print "Mean length of pfamAB regions: $prop\n";

if ($anz_pfamA)
{
    $prop = $cum_len_pfamA/$anz_pfamA;
}
else
{
    $prop = 0;
}
print "Mean length of pfamA  regions: $prop\n";

if ($anz_pfamB)
{
    $prop = $cum_len_pfamB/$anz_pfamB;
}
else
{
    $prop = 0;
}
print "Mean length of pfamB  regions: $prop\n";

if ($anz_void)
{
    $prop = $cum_len_void/$anz_void;
}
else
{
    $prop = 0;
}
print "Mean length of voids  reginos: $prop\n";



$, = "\n";
    
open OUT, ">genelengths.txt";
print OUT @genelengths;
close OUT;

open OUT, ">voidlengths.txt";
print OUT @void_lengths;
close OUT;

open OUT, ">pfamA_lengths.txt";
print OUT @pfamA_lengths;
close OUT;

open OUT, ">pfamB_lengths.txt";
print OUT @pfamB_lengths;
close OUT;

open OUT, ">clans.txt";
while ( my ($k, $v) = each(%clan_hash) )
{
    print OUT "$k: $v\n";
}
close OUT;

