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

#include<iostream>
#include "CSequences2.h"
#include "faststring2.h"
#include "CFile/CFile2_1.h"
#include "CSequence_Mol2_1.h"
#include <vector>

using namespace std;

bool is_gap_pattern(faststring &str)
{
  unsigned N = str.length();
  unsigned i;

  for (i=0; i<N; ++i)
    if (str[i] != '-')
      return false;
  return true;
}


int main(int argc, char **argv)
{
  //  bool inherit_lower_case_in_refine_range = false;

  // Lower case letter in the alignment:
  // In regions we do not align since the backbone taxa have a gap, we write lower case letters to indicate they are not aligned.
  // If the sequences originate form other sources, such as hmmalign they might already contain lower case regions that are not aligned.
  // If they are present in the first file, they are inherited.
  // If they are in the second file we convert them to upper case symbols if they are aligned to the backbone. They remain lower case
  // if there is no backbone.

  if (argc != 3)
  {
//     if (artc == 4 && faststring(argv[3]) == "inherit" )
//     {
//       inherit_lower_case_in_refine_range = true;
//     }
//     else
    {
      cerr << "Usage: " << argv[0] << " fasta-alignment-file fasta-alignment-file2" << endl
	//	   << "OR" << endl
	//	   << argv[0] << " fasta-alignment-file fasta-alignment-file2 inherit  || Does not change lower case to upper case in existing aligned regions." << endl
	   << endl;
    }
  }

  //  faststring fname1 = "EOG5X69Q3_ZNEVA_2.1_Tricholepidion_gertschi_s12686_L_124975_1-517PPC1186713-580.fas_mafft_linsi_add.fas"; // 2286
  //  faststring fname2 = "EOG5X69Q3_ZNEVA_2.1_Peruphasma_schultei_C1181909-134PPs8514_L_124320_0-149.fas_mafft_linsi_add.fas";      // 2278

  faststring fname1 = argv[1];
  faststring fname2 = argv[2];

  CFile is1, is2;

  is1.ffopen(fname1.c_str());
  is2.ffopen(fname2.c_str());
  
  if (is1.fail())
  {
    cerr << "File: " << fname1 << " does not exist. " << endl;
  }

  if (is2.fail())
  {
    cerr << "File: " << fname2 << " does not exist. " << endl;
  }

  CSequences2 seqs1(CSequence_Mol::protein);
  CSequences2 seqs2(CSequence_Mol::protein);

  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(0);
  seqs1.read_from_Fasta_File(is1, pflag, 0, -1, false);
  seqs2.read_from_Fasta_File(is2, pflag, 0, -1, false);

  faststring pat1, pat2;

  faststring s1, s2;

  vector<faststring> new_alig;

  unsigned N1, N2, i1, i2;
  bool gappat1, gappat2;

  unsigned Ntaxa1 = seqs1.GetTaxaNum();
  unsigned Ntaxa2 = seqs2.GetTaxaNum();

  if (Ntaxa2 > Ntaxa1)
  {
    cerr << "The second alignment must be the one with fewer or equal number of taxa comapred to first alignment." << endl;
  }

  if (Ntaxa1 <= 0 || Ntaxa2 <= 0)
  {
    cerr << "Bad number of taxa found in file. One of the files might be empty." << endl;
    exit(-3);
  }

  unsigned num_cores = Ntaxa2-1;

  cerr << "Number of cores used in combine step: " << num_cores << endl;

  N1 = seqs1.GetPosNum();
  N2 = seqs2.GetPosNum();
  //  N = N1; if (N<N2) N=N2;

  if (N1 <= 0 || N2 <= 0)
  {
    cerr << "Bad number of sequence posiitons found in file. One of the files might be empty." << endl;
    exit(-3);
  }

  unsigned itaxon;

  // Space for existing alignment:
  for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
  {
    new_alig.push_back(faststring());
  }
  // Space for new sequence:
  new_alig.push_back(faststring());

  i1=0;
  i2=0;
  while( i1 < N1 && i2 < N2)
  {
    seqs1.get_partial_pattern(i1, num_cores, pat1);
    seqs2.get_partial_pattern(i2, num_cores, pat2);

    gappat1 = is_gap_pattern(pat1);
    gappat2 = is_gap_pattern(pat2);
  
    if (pat1 == pat2)
    {
      // Copy existing alignment column:
      for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
      {
	new_alig[itaxon].push_back(seqs1.GetChar(itaxon, i1));
      }
      // Add symbol to new sequence:
      if (gappat1 && gappat2)
	new_alig[Ntaxa1].push_back(tolower(seqs2.GetChar(num_cores, i2)) );
      else
	//	if (inherit_lower_case_in_refine_range)
	//	  new_alig[Ntaxa1].push_back(seqs2.GetChar(num_cores, i2) );
	//	else
	  new_alig[Ntaxa1].push_back(toupper(seqs2.GetChar(num_cores, i2)) );

      //     s1.push_back(seqs1.GetChar(12, i1));  // Old code to combine just two files with one additional sequnce.
      //     s2.push_back(seqs2.GetChar(12, i2));

      ++i1;
      ++i2;
    }
    else
    {
      if (gappat1)
      {
	// We have gaps in the cores of the first file.
	// This means there was an insertion introduced by one of the remaining sequences of the first file.
	// Furthermore, we have no corresponding base in the second file.
	// So we copy this column from the first file and add a gap in the new sequence.
	//	cerr << "Gap-pat s1: " << i1 << endl;

	// Copy existing alignment column:
	for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
	{
	  // This is an existing alignment, which we do not transform to lower case. //LOWER
	  new_alig[itaxon].push_back(seqs1.GetChar(itaxon, i1));
	}
	// Add gap to new sequence:
	new_alig[Ntaxa1].push_back('-');

	//	s1.push_back(seqs1.GetChar(12, i1)); // Old code to combine just two files with one additional sequnce.
	//	s2.push_back('-');
	++i1;
      }
      else if (gappat2)
      {
	// We have gaps in the cores of the second file.
	// This means there was an insertion introduced by the new sequence in the second file.
	// We add gaps in the existing alignemnt
	// we add the base in lower case to the new alingment.

	//	cerr << "Gap-pat s2: " << i2 << endl;

	for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
	{
	  new_alig[itaxon].push_back('-');
	}
	
	new_alig[Ntaxa1].push_back( tolower( seqs2.GetChar(num_cores, i2) ) );

	//	s2.push_back(seqs2.GetChar(12, i2)); // Old code to combine just two files with one additional sequnce.
	//	s1.push_back('-');
	++i2;
      }
      else
      {
	cerr << "Non matching alignment columns detected: " << endl;
	cerr << "File 1: " << fname1 << endl;
	cerr << "File 2: " << fname2 << endl;
	cerr << "Columns: "  << i1 << " " << i2 << endl;	
	exit(-5);
      }
    }
  } // End while

  // We should either have i1 == N1 or i2 == N2.
  // The other sequence must be completed:
  while (i1 < N1)
  {
    for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
    {
      new_alig[itaxon].push_back(tolower( seqs1.GetChar(itaxon, i1) ) );
    }
    new_alig[Ntaxa1].push_back('-');
    //    s1.push_back(seqs1.GetChar(12, i1));
    ++i1;
  }

  while (i2 < N2)
  {
    for (itaxon = 0; itaxon < Ntaxa1; ++itaxon)
    {
      new_alig[itaxon].push_back('-');
    }
    // This part does not have a backbone, since the seqs1 had no bases here.
    // Therefore we add bases with lower case.
    new_alig[Ntaxa1].push_back(tolower(seqs2.GetChar(num_cores, i2)  ) );

   //     s2.push_back(seqs2.GetChar(12, i2));
    ++i2;
  }


  // Print alignment:

  unsigned i, N = Ntaxa1;

  for (i=0; i<N; ++i)
  {
    cout << ">" << seqs1.get_Seq_Name(i) << endl;
    cout << new_alig[i]       << endl;
  }

  cout << ">"  << seqs2.get_Seq_Name(num_cores) << endl;
  cout << new_alig[Ntaxa1] << endl;

}
