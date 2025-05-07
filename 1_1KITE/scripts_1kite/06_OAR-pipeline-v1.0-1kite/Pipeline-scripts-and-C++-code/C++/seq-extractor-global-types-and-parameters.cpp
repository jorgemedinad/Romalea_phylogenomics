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

#include "seq-extractor-global-types-and-parameters.h"
#include "tclap/CmdLine.h"
#include "climits"

using namespace TCLAP;
using namespace std;


string                      global_input_filename;

string                      global_delim = "::";

string                      global_seq_name;
string                      global_seq_name_contains;

integer                     global_begin_pos = 1;
integer                     global_end_pos = INT_MAX;

integer                     global_verbose_flag = 0;

unsigned                    global_interleaved = -1u;
bool                        global_revcomp;
bool                        global_split_mode;

integer                     global_sequence_from = 1;          // if negative count from back
integer                     global_sequence_to = INT_MAX;            // if negative count from back

string                      global_seq_numbers_str;


void good_bye_and_exit(int error)
{
  if (error == 0)
    cout << "## Finished successfully. ##" << endl;
  cerr << endl << PROGNAME << " says goodbye." << endl << endl;
  exit(error);
}



void init_param()
{
}



void read_and_init_parameters(int argc, char** argv)
{
  init_param();

  try
  {
    CmdLine cmd("This program can extract parts of a specified sequence.",
		' ', VERSION);

    ValueArg<integer> sequence_from_Arg("", "firstSeq",
       "Number of first sequence to be processed in this run. Negative numbers are counted from the back of the file.",
	false, global_sequence_from, "int");
    cmd.add( sequence_from_Arg );

    ValueArg<integer> sequence_to_Arg("", "lastSeq",
       "Number of last sequence to be processed in this run.  Negative numbers are counted from the back of the file.",
	false, global_sequence_to, "int");
    cmd.add( sequence_to_Arg );


    ValueArg<unsigned> interleaved_Arg("i", "interleaved",
	"The maximum number of sequence characters per line. Default: infinity",
	false, global_interleaved, "int");
    cmd.add(interleaved_Arg);

    SwitchArg verbosity_Arg("V", "verbosity",
	"Provides additional output, e.g. for debugging or test purposes. "
        "Lists for example the parameters that have been specified, so one can check that these have been "
        "read properly. Default: no additional output.",
	global_verbose_flag);
    cmd.add(verbosity_Arg);

    SwitchArg split_mode_Arg("", "split",
	"Split mode: Extracts all sequences in file and saves them to individual fasta files. "
			    "Some special characters in the sequence names are replaced so that they are valid file names. "
			    "Many of the other options are not valid in this mode.",
	false);
    cmd.add(split_mode_Arg);

    SwitchArg revcomp_Arg("R", "revcomp",
        "Print reverse complement.",
        global_revcomp);
    cmd.add(revcomp_Arg);

    ValueArg<string> seq_name_Arg("s","sname",
				 "The name of the sequence to extract from. Exact match to the main sequence name is required. Note that the main sequence name is the name before any spaces in the sequence namen.",
				 false, "", "sequence name"  );
    cmd.add( seq_name_Arg );

    ValueArg<string> seq_name_contains_Arg("c","contains",
					 "If the substring is contained in the full sequence names, the sequence will be extracted. Any number of sequences may be extracted in the specified sequence range.",
					 false, "", "partial sequence name"  );
    cmd.add( seq_name_contains_Arg );


    ValueArg<integer> beg_pos_Arg("","beg",
       "Beginning of extraction region.",
				  false, global_begin_pos, "int");
    cmd.add( beg_pos_Arg );


    ValueArg<integer> end_pos_Arg("","end",
				  "End of extraction region.",
				  false, global_end_pos, "int");
    cmd.add( end_pos_Arg );


    UnlabeledValueArg<string> filename_Arg( "name4",
       "The FASTA file containing the input sequences.", "",
       "fasta file"  );
    cmd.add( filename_Arg );

    cmd.parse( argc, argv );

    // Assigning parameters to variables:
    global_input_filename          = filename_Arg.getValue();

    global_sequence_from           = sequence_from_Arg.getValue();
    global_sequence_to             = sequence_to_Arg.getValue();

    global_seq_name                = seq_name_Arg.getValue();
    global_seq_name_contains         = seq_name_contains_Arg.getValue();

    global_begin_pos               = beg_pos_Arg.getValue();
    global_end_pos                 = end_pos_Arg.getValue();
    global_interleaved             = interleaved_Arg.getValue();
    global_revcomp                 = revcomp_Arg.getValue();
    global_verbose_flag            = verbosity_Arg.getValue();
    global_split_mode              = split_mode_Arg.getValue();
  }
  catch (ArgException &e)
  {
    cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
    good_bye_and_exit(-1);
  }

  // Checking parameters:

  if (global_split_mode && global_seq_name != "" )
  {
    cerr << "Error: In split mode you cannot specify a sequence name." << endl;
    exit(-3);
  }


}



void print_search_parameters(std::ostream &os, const char *s)
{
  os << s << "Parameter settings:"
     << std::endl;

  os << s << "Input file name:                    " << global_input_filename
     << std::endl;

  os << s << "Sequence name:                      " << global_seq_name
     << std::endl;

  os << s << "Sequence range to process:          " << global_sequence_from
     << " to ";
  if (global_sequence_to == INT_MAX)
    os << "last sequence\n";
  else
    os << global_sequence_to << endl;

  os << s << "Begin position:                     " << global_begin_pos
     << std::endl;

  os << s << "End position:                       " << global_end_pos
     << std::endl;

  
  os << s << "Reverse complement:                 " << (global_revcomp ? "yes" : "no")
     << std::endl;

  os << s << "Verbosity:                          " << (global_verbose_flag ? "yes" : "no") 
     << std::endl;
  os << s << "Global split mode:                  " << (global_split_mode ? "yes" : "no")
     << std::endl;
}


void print_output_creation_message(std::ostream &os, char *s)
{
  os << s << "Output created by " << PROGNAME << ", version " << VERSION
     << std::endl;
  os << s << std::endl;
  print_search_parameters(os, s);
  os << s << std::endl;
}


