/***************************************************************************************************
*  This file is part of the 1kite analysis pipeline. It is
*  distributed under the following license:
*  
*  Copyright (c) 2012-20013 Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany
*  All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  1. Redistributions of source code (complete or in parts) must retain
*     the above copyright notice, this list of conditions and the following disclaimer.
*  2. Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer in the
*     documentation and/or other materials provided with the distribution.
*  3. All advertising materials mentioning features or any use of this software
*     e.g. in publications must display the following acknowledgement:
*     This product includes software developed by Christoph Mayer, Forschungsmuseum
*     Alexander Koenig, Bonn, Germany.
*  4. Neither the name of the organization nor the
*     names of its contributors may be used to endorse or promote products
*     derived from this software without specific prior written permission.
*  
*  THIS SOFTWARE IS PROVIDED BY CHRISTOPH MAYER ''AS IS'' AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHTHOLDER OR ITS ORGANISATION BE LIABLE FOR ANY
*  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*  
*  IMPORTANT (needs to be included, if code is redistributed):
*  Please not that this license is not compatible with the GNU Public License (GPL)
*  due to paragraph 3 in the copyright. It is not allowed under any
*  circumstances to use the code of this software in projects distributed under the GPL.
*  Furthermore, it is not allowed to redistribute the code in projects which are
*  distributed under a license which is incompatible with one of the 4 paragraphs above.
*  
*  This project makes use of code coming from other projects. What follows is a complete
*  list of files which make use of external code. Please refer to the copyright within
*  these files.
*  
*  Files in tclap foler:         Copyright (c) 2003 Michael E. Smoot
*                                See copyright in tclap/COPYRIGHT file for details.	
***************************************************************************************************/

#include "CSequences2.h"
#include <iostream>
#include "CFile/CFile2_1.h"
#include "faststring2.h"
#include "fast-realloc-vector.h"

#include "CSeqNameList.h"

using namespace std;

int split_file(const char *fn, const char *tl)
{
  CFile file;

  cout << "Reading fasta input file: " << fn << endl;

  file.ffopen(fn);

  if (file.fail())
  {
    cout << "Could not open specified file: " << fn << endl;
    exit(-1);
  }

  //  cout << "Line:   " <<  file.line() << endl;
  // cout << "Status: " <<  (int)file.status() << endl;

  CSequences2 seqs(CSequence_Mol::protein);
  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(0);

  seqs.read_from_Fasta_File(file, pflag, 0, -1, false);
  file.ffclose();
  cout << "Found this number of taxa:      " << seqs.GetTaxaNum()  << endl;
  cout << "Found this number of positions: " << seqs.GetPosNum()  << endl;

  CSeqNameList snl(tl, 13);

  int i, n = seqs.GetTaxaNum();

  if (true) // Debug output:
  {
    for (i=0; i<n; ++i)
    {
      cout << "Taxon: " << i << ": " << seqs.get_seq(i)->getName() << endl;
      //    cout << seqs.get_seq(i)->getSeqStr() << endl;
    }
    snl.print(cout, 1); 
  }

  faststring alig_with_outlier    = "in.fas";
  faststring alig_without_outlier = "not_in.fas";

  // Separate the sequences into two files.
  snl.use_list_as_filter(fn, alig_with_outlier.c_str(), alig_without_outlier.c_str() , -1u);

  cerr << "Saved alignmentsto files: in.fas. not_in.fas." << endl;
}



int main(int argc, const char ** argv)
{
  if (argc <3 || argc > 3)
  {
    cout << "Usage:" << endl
	 << argv[0] << " fasta-file sequence-list-file." << endl;
    exit(0);
  }
  else
  {
    cout << "Welcome to the split-fasta-with-sequence-list program." << endl;
    cout.flush();
  }
  
  faststring filename = argv[1];
  faststring taxonList= argv[2];
  split_file(filename.c_str(), taxonList.c_str() );

}

