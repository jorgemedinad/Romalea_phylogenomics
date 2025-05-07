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

#!/usr/bin/perl

use strict;


my @filearray = <>;
chomp(@filearray);
## print @filearray;

my $N;
$N = scalar(@filearray);

print $N;

my $i;
my $line;
my $mode; ## 0, no file open, 1 file open
my $filename;

for ($i=1; $i < $N; ++$i)
{
    $line = $filearray[$i];

    while ($line eq "")
    {
	++$i;
	$line = @filearray[$i];
    }

    $filename = $line;
    $filename =~ s/fas:/outlier/;
    open (MYFILE, "> $filename");

    ++$i;
    while ($i < $N && !($filearray[$i] eq "") )
    {
	$line = $filearray[$i];
	print MYFILE "$line\n";
	++$i;
    }
    close(MYFILE);
}


