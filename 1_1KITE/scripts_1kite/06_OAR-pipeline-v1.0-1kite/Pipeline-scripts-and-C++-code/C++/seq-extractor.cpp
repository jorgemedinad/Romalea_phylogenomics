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

//    Program: sequence-chopper
//
//    Christoph Mayer
//    Spezielle Zoologie, Department of Special Zoology
//    Gebaeude: ND 05/785
//    Ruhr Universitaet Bochum
//    44780 Bochum
//    Germany
//

#include "seq-extractor.h"


using namespace std;


unsigned count_sequences_in_file(CFile &is)
{
  faststring line;
  unsigned count=0;

  while (!is.fail())
  {
    is.getline(line);
    if (line[0] == '>')
      ++count;
  }
  return count;
}



void welcome(std::ostream &os)
{
  //  os << sizeof(Cssr) << endl;
  os << endl << endl;
  os << "      Welcome to " << PROGNAME << ", version " << VERSION << ","
     << endl
     << "      the sequence-chopper program." << endl << endl << endl;
}


void print_seq(ostream &os, CSequence_Mol &seq, string use_this_seq_name=string() )
{
  
  if (global_end_pos == -1u)
    global_end_pos = seq.length();

  if (use_this_seq_name.empty() )
  {
    //    os << ">" << global_input_filename 
    os << ">" << seq.getFullName()
       << global_delim          <<  global_seq_name
       << global_delim  << ""  <<  global_begin_pos
       << "-"                 << global_end_pos
       << global_delim  << "+strand"
       << endl;
  }
  else
  {
    os << ">" << use_this_seq_name << endl; 
  }

  const char *data     = seq.getSeqBegin() + global_begin_pos - 1; 
  const char *data_end = seq.getSeqBegin() + global_end_pos;

  if (data > data_end)
    return;

//   if (data < seq.getSeqBegin() )
//     data = seq.getSeqBegin();

  if (data_end > seq.getSeqEnd() )
    data_end = seq.getSeqEnd();

  unsigned interleaved = global_interleaved;

  while (data < data_end)
  {
    os << *data;
    --interleaved;
    if (interleaved == 0)
    {
      os << endl;
      interleaved = global_interleaved;
    }
    ++data;
  }
}


void print_revcomp(ostream &os, CSequence_Mol &seq, string use_this_seq_name=string() )
{
  CDnaString s, r;

  if (global_end_pos == -1u)
    global_end_pos = seq.length();

  if (use_this_seq_name.empty() )
  {
    //    os << ">" << global_input_filename 
         os << ">" << seq.getName()
	 << global_delim          <<  global_seq_name
	 << global_delim  << ""  <<  global_begin_pos
	 << "-"                 << global_end_pos
	 << global_delim  << "-strand"
	 << endl;
  }
  else
  {
    os << ">" << use_this_seq_name << endl; 
  }

  const char *data     = seq.getSeqBegin() + global_begin_pos - 1; 
  const char *data_end = seq.getSeqBegin() + global_end_pos;

  if (data > data_end)
    return;

//   if (data < seq.getSeqBegin() )
//     data = seq.getSeqBegin();

  if (data_end > seq.getSeqEnd() )
    data_end = seq.getSeqEnd();

  s.toupper_assign(data, data_end);
  r.setToReverseComplementOf(s);

  data     = r.begin();
  data_end = r.end();

  unsigned interleaved = global_interleaved;

  while (data < data_end)
  {
    os << *data;
    --interleaved;
    if (interleaved == 0)
    {
      os << endl;
      interleaved = global_interleaved;
    }
    ++data;
  }
}


int main(int argc, char **argv)
{
  // Since we do not convert any ambiguity characters, protein sequences can be read as well.
  CSequence_Mol   seq(CSequence_Mol::dna);
  CFile           infile;
  integer        count_seq=0;

  // Welcome message
  if (global_verbose_flag)
  {
    welcome(cerr);
    cerr << endl << endl;
  }

  // Read command line options
  read_and_init_parameters(argc, argv);

  if (global_verbose_flag)
  {
    print_search_parameters(cerr, "");
    cerr << endl << flush;
  }

  // Open file and check errors:
  infile.ffopen(global_input_filename.c_str());
  if ( !infile.exists() )
  {
    cerr << "Sorry! The input file \n   \"" << global_input_filename
	 << "\"\ncould not be opened. I guess it does not exist."
	 << endl;
    good_bye_and_exit(1);
  }

  if (global_sequence_from < 0 || global_sequence_to < 0)
  {
    // In this case we need to determien the number of sequences in the file:
    unsigned N = count_sequences_in_file(infile);
    infile.rewind();

    //    cout << "Number of sequences in file: " << N << endl;

    if (global_sequence_from < 0) global_sequence_from +=(N+1);
    if (global_sequence_to   < 0) global_sequence_to   +=(N+1);
  }



  cout.setf(ios::fixed);
  cout.precision(3);

  bool global_seq_name_specified          = !global_seq_name.empty();
  bool global_seq_name_contains_specified = !global_seq_name_contains.empty(); 

  while (!infile.eof())
  {
    seq.readRawFastaSequence(infile);
    ++count_seq;

    if ( infile.fail() && !infile.eof() )
    {
      cerr << endl << endl;

      cerr << "An error occurred while reading the input file."
	   << " It might not be a valid fasta file.\n"
	   << "File position: line "
	   << faststring(infile.line()).c_str()
	   << endl << flush;
       return -25;
    }

    if (count_seq < global_sequence_from || count_seq > global_sequence_to)
    {
      if (global_verbose_flag)
      {
	cerr << "Skipping sequence: " << seq.getName() << endl;
      }
    }
    else // Sequences that are potentially extracted
    {
      if (global_split_mode) // Save all sequences in individual files
      {
	string fname = seq.getName();  // File name
	string sname = seq.getName();  // Sequence name

	// File names are not allowed to contain all characters that are allowed in the sequence names:
	replace(fname.begin(), fname.end(), '|', '_');
	replace(fname.begin(), fname.end(), ' ', '_');
	replace(fname.begin(), fname.end(), '\t', '_');
	replace(fname.begin(), fname.end(), '/', '_');
	replace(fname.begin(), fname.end(), '\\', '_');
	//      replace(fname.begin(), fname.end(), '(', '_');
	//      replace(fname.begin(), fname.end(), ')', '_');
	
	if (global_revcomp)
	  fname += "-revcomp";
	
	fname +=".fas";
	
	ofstream os(fname.c_str() );
	
	if (!global_revcomp)
	{
	  print_seq(os, seq, sname);
	}
	else
	{
	  sname += global_delim;
	  sname += "revcomp";
	  print_revcomp(os, seq, sname);
	}
	os << endl;
	
	os.close();
      }
      else if (global_seq_name_specified)
      {
	if  (seq.getName() == global_seq_name)  // extract single sequence mode. Write to cout
	{
	  if (!global_revcomp)
	  {
	    if (global_begin_pos == 1 && global_end_pos >= seq.length() ) // in this case we keep the old sequnce name since we extract complete sequence
	      print_seq(cout, seq, seq.getFullName());
	    else
	      print_seq(cout, seq);  // In this case, the coordinates will be added to the sequnce name.
	  }
	  else
	  {
	    if (global_begin_pos == 1 && global_end_pos >= seq.length() ) // in this case we keep the old sequnce name since we extract complete sequence
	      print_revcomp(cout, seq, seq.getFullName());
	    else
	      print_revcomp(cout, seq);  // In this case, the coordinates will be added to the sequnce name.
	  }
	  cout << endl;
	  break;   // We expect only one exact match, so we break out of the "while (!infile.eof())" loop 
	}
      }
      else if (global_seq_name_contains_specified) // search for partial matches -- multiple matches are possible/expected
      {
	unsigned found;
	faststring sname = seq.getFullName();  // Sequence name
	found=sname.find(global_seq_name_contains.c_str());
	if (found != std::string::npos) // search string has been found in sequence name, so we extract this sequence.
	{
	  if (!global_revcomp)
	  {
	    if (global_begin_pos == 1 && global_end_pos  >= seq.length() ) // in this case we keep the old sequnce name since we extract complete sequence
	      print_seq(cout, seq, seq.getFullName());
	    else
	      print_seq(cout, seq);  // In this case, the coordinates will be added to the sequnce name.
	  }
	  else
	  {
	    if (global_begin_pos == 1 && global_end_pos  >= seq.length() ) // in this case we keep the old sequnce name since we extract complete sequence
	      print_revcomp(cout, seq, seq.getFullName());
	    else
	      print_revcomp(cout, seq);  // In this case, the coordinates will be added to the sequnce name.
	  }
	  cout << endl;
	}
      }
      else  // We print all sequences in from  - to range:
      {
	string sname = seq.getFullName();

	if (!global_revcomp)
	{
	  print_seq(cout, seq, sname);
	}
	else
	{
	  sname += global_delim;
	  sname += "revcomp";
	  print_revcomp(cout, seq, sname);
	}
	cout << endl;
      }
    } // else - skip this sequnce
  }

  infile.ffclose();

  if (global_verbose_flag)
  {
    cerr << endl;
    good_bye_and_exit(0);
  }
}
