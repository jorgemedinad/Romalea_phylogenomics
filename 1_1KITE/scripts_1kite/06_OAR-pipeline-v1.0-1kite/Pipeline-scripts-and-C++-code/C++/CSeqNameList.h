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

#ifndef CSEQNAMELIST_H
#define CSEQNAMELIST_H


#include <iostream>
#include <fstream>
#include "faststring2.h"
#include <map>
#include <set>
#include "CFile/CFile2_1.h"
#include "CSequence_Mol2_1.h"
#include <iterator>
#include <algorithm>


inline void set_to_short_name(faststring &str)
{
  unsigned  pos = str.find(' ');
  if (pos != faststring::npos) 
    str.shorten(pos);
}


class CSeqNameList
{
  // A note on memory management.
  // This class collects names and information about sequences.
  // Sequence names are stored in faststrings. The class keeps pointers to these faststrings.
  // The class keeps a copy of the original string it obtained. Therefore, it allocates memory
  // for the strings and it has to release the memory when they are not used any more.


 public:
  typedef std::map<faststring*, unsigned, less_than_pointer_to_faststring_struct> myMap;
  typedef std::set<faststring>            mySet;
  typedef myMap::iterator CSeqNameListIterator;
  typedef mySet::iterator mySetIterator;

 private:
  faststring infile_name;
  bool       only_names;  // Formerly, all vectors had been initialized to ensure that all fields are accessible.
                          // This is not the case any more. If only_names is true, the vectors defined below can be empty.
                          // TODO: Check that this is handled consistently.

  myMap                            m;  // maps sequence names to id's. id's are indices used to refer to the other sequence data.
  std::vector<faststring*>         full_names;
  std::vector<unsigned>            seq_len_vec;
  std::vector<double>              CG_content;
  std::vector<unsigned>            corrected_len_vec;
  std::vector<unsigned>            repeat_density_bp_mbp;

 public:
  void add(const faststring &full_name, unsigned sl, double CG_c, unsigned cl, unsigned rdens)
  {
    // The basic add routine creates a copy of the sequence name that shall be added.
    faststring *tmp = new faststring(full_name);

    full_names.push_back(tmp);
    seq_len_vec.push_back(sl);
    CG_content.push_back(CG_c);
    corrected_len_vec.push_back(cl);
    repeat_density_bp_mbp.push_back(rdens);

    //    std::cerr << "Adding new seq: Name " << *tmp << " len " << sl << " CG " << CG_c << " cl " << cl<< std::endl; 

    unsigned pos = full_names.size()-1;

    m[tmp] = pos;
  }

  void add_only_name(const char *full_name)
  {
    faststring *tmp = new faststring(full_name);
    full_names.push_back(tmp);

    unsigned pos = full_names.size()-1; 
    m[tmp] = pos;
  }

  void add(const CSeqNameList &snl, unsigned id)
  {
    if (snl.only_names)
      add_only_name( snl.full_names[id]->c_str() );
    else
      add(*(snl.full_names[id]), snl.seq_len_vec[id], snl.CG_content[id], snl.corrected_len_vec[id], snl.repeat_density_bp_mbp[id]);
  }

  bool add_non_redundant(const CSeqNameList &snl, unsigned id)
  {
    if (!snl.only_names)
    {
      if (get_id_short_name(*(snl.full_names[id])) == -1u)
      {
	add(*(snl.full_names[id]), snl.seq_len_vec[id], snl.CG_content[id], snl.corrected_len_vec[id], snl.repeat_density_bp_mbp[id]);
	return true;
      }
      else
      {
	std::cerr << "Detected previous entry of " << *snl.full_names[id] << std::endl;
	return false;
      }
    }
    else // only_names
    {
      if (get_id_short_name(*(snl.full_names[id])) == -1u)
      {
	add_only_name( snl.full_names[id]->c_str() );
	return true;
      }
      else
      {
	std::cerr << "Detected previous entry of " << *snl.full_names[id] << std::endl;
	return false;
      }
    }
    
  }


  void reset()
  {
    // Release memory. Otherwise we have a memory leak;
    unsigned i=0, n=full_names.size();
    while (i<n)
    {
      delete full_names[i];
      ++i;
    }

    infile_name.clear();
    only_names = false;
    m.clear();
    full_names.clear();
    seq_len_vec.clear();
    CG_content.clear();
    repeat_density_bp_mbp.clear();
    corrected_len_vec.clear();
  }

  // remove entry from the "list".
  void remove(faststring full_name)
  {
    CSeqNameListIterator it;
    it = m.find(&full_name);    // We have to take the address operator here. Remember that we use less_than_pointer_to_faststring_struct to compare the elements, so we do no search for the pointer but for the string indeed. 

    unsigned pos = it->second;
    *(full_names[pos]) = "";
    if (!only_names)
    {
      seq_len_vec[pos] = 0;
      CG_content[pos]  = 0;
      corrected_len_vec[pos]   = 0;
      repeat_density_bp_mbp[pos]  = 0;
    }
    m.erase(it);  // Erase this element out of the map.
  }

  // Import sequence names from a fasta file
 CSeqNameList(const char *fasta_file_name, bool read_only_names=false):infile_name(fasta_file_name), only_names(read_only_names)
  {
    CFile infile;

    infile.ffopen(fasta_file_name);

    if ( !infile.exists() )
    {
      std::cerr << "Sorry! The input file \n   \"";
      std::cerr << fasta_file_name;
      std::cerr << "\"\ncould not be opened. I guess it does not exist.\n";
      exit (-101);
    }
    
    CSequence_Mol seq(CSequence_Mol::unknown);

    if (!only_names)
    {
      while (!infile.eof())
      {
	seq.readRawFastaSequence(infile);
	if ( infile.fail() && !infile.eof() )
	{
	  std::cerr << "\n\n";
	  std::cerr << "An error occurred while reading the input file. It might not be a valid fasta file.\n";
	  std::cerr << "File position: line ";
	  std::cerr << faststring(infile.line()) << std::endl;
	  std::cerr << std::flush;
	}
	else
	{
	  faststring name_to_add = seq.getFullName();
	  
	  //	 std::cerr << "Procesing sequence: " << name_to_add << std::endl;
	  
	  seq.compute_numbers_of_residues();
	  
	  unsigned Ambigs      = seq.get_number_of_DNARNA_ambigs();
	  unsigned bases       = seq.get_number_of_DNARNA_bases();
	  double   CGc         = seq.get_CG_content();
	  unsigned rdens       = 0;
	 
	  //	 std::cerr << "Adding 2: " << name_to_add << " slen: " << seq.length() << " CG: " << CGc << " bases: " << bases <<  std::endl;
	  
	  add(name_to_add, seq.length(), CGc, bases, rdens);
	}
      } // End while
    }
    else // Read only names:
    {
      while (!infile.eof())
      {
	seq.readSeqName_ignore_Sequence_data(infile);
	if ( infile.fail() && !infile.eof() )
	{
	  std::cerr << "\n\n";
	  std::cerr << "An error occurred while reading the input file \""<< fasta_file_name << "\" It might not be a valid fasta file.\n";
	  std::cerr << "File position: line ";
	  std::cerr << faststring(infile.line()) << std::endl;
	  std::cerr << std::flush;
	  break; // Stop reading this file. 
	}
	else
	{
	  const char * name_to_add = seq.getFullName();
	  //	 std::cerr << "Procesing sequence: " << name_to_add << std::endl;
	  add_only_name(name_to_add);
        }
      } // End while
    } // End else only_names

      infile.ffclose();
  } // END constructor



  // Import sequence names from a text file:
 CSeqNameList(const char *text_file_name, int dummy):infile_name(text_file_name), only_names(true)
  {
    CFile infile;
    faststring name_to_add;

    infile.ffopen(text_file_name);

    if ( !infile.exists() )
    {
      std::cerr << "Sorry! The input file \n   \"";
      std::cerr << text_file_name;
      std::cerr << "\"\ncould not be opened. I guess it does not exist.\n";
      exit (-101);
    }

    infile.getline(name_to_add);
    while (!infile.eof())
    {
      //	 std::cerr << "Procesing sequence: " << name_to_add << std::endl;
      add_only_name(name_to_add.c_str());
      infile.getline(name_to_add);
    } // End while

    infile.ffclose();
  } // END constructor


  // Default constructor
 CSeqNameList(bool read_only_names=false):infile_name(""), only_names(read_only_names)
  {}


  

  // File must be of the form:
  // sequence name\tdenstiy in bp/Mbp
  void add_repeat_densities(const char * filename)
  {
    std::ifstream is(filename);
    faststring line;
    std::vector<faststring> vf;
    faststring sn;
    unsigned rdens;
    unsigned id;

    if (only_names)
    {
      std::cerr << "Internal error: Reading repeat denstities for names only sequence list.\n";
      exit(-22);
    }


    getline(is, line);
    while (is)
    {
      split(vf, line, "\t");
      if (vf.size()!=2)
      {
	std::cerr << "Critical error when addind repeat densities. Line " << line << " not a valid input.\n";
	exit(-21);
      }
      sn    = vf[0];
      rdens = vf[1].ToUnsigned();

      id = get_id_short_name(sn);

      repeat_density_bp_mbp[id] = rdens;

      getline(is, line);
    }
    is.close();
  }

  ~CSeqNameList()
  {
/*     unsigned i=0, n=full_names.size(); */

/*     while (i<n) */
/*     { */
/*       delete full_names[i]; */
/*       ++i; */
/*     } */
    reset();
  } // END destructor





  void use_list_as_filter(const char *fasta_file_name,
			  const char *out_in_list,
			  const char *out_not_in_list,
			  unsigned char_per_line=50)
  {
    CFile infile;

    infile.ffopen(fasta_file_name);

    if ( !infile.exists() )
    {
      std::cerr << "Sorry! The input file \n   \"";
      std::cerr << fasta_file_name;
      std::cerr << "\"\ncould not be opened. I guess it does not exist.\n";
      exit (-101);
    }
    
    FILE *os_in_list; // (out_in_list);
    FILE *os_not_in_list; // (out_not_in_list);

    os_in_list     = fopen(out_in_list, "w");
    os_not_in_list = fopen(out_not_in_list, "w");

    CSequence_Mol seq(CSequence_Mol::unknown);
     while (!infile.eof())
     {
       seq.readRawFastaSequence(infile);
       if ( infile.fail() && !infile.eof() )
       {
	 std::cerr << "\n\n";
	 std::cerr << "An error occurred while reading the input file. It might not be a valid fasta file.\n";
	 std::cerr << "File position: line ";
	 std::cerr << faststring(infile.line()) << std::endl;
	 std::cerr << std::flush;
       }

       faststring fname = seq.getFullName();

       CSeqNameListIterator it;
       it = m.find(&fname);    // We have to take the address operator here. Remember that we use less_than_pointer_to_faststring_struct to compare the elements, so we do no search for the pointer but for the string indeed. 

       if (it != m.end() )
       {
	 seq.writeSequence_fasta(os_in_list, char_per_line);
       }
       else
       {
	 seq.writeSequence_fasta(os_not_in_list, char_per_line);
       }
       
     } // END while

     infile.ffclose();
     fclose(os_in_list);
     fclose(os_not_in_list);

  } // END use_list_as_filter


  // Reads the specified fasta file.
  // If the sequence name is found in this name list write it to output file.
  // XXXXXXXXXXXX TODO: return value: number of sequences with not only sequence name but also sequence info.
  void write_sequences_in_list(const char *fasta_file_name,
			       const char *out_in_list,
			       const char *write_mode)
  {
    CFile infile;

    infile.ffopen(fasta_file_name);

    if ( !infile.exists() )
    {
      std::cerr << "Sorry! The input file \n   \"";
      std::cerr << fasta_file_name;
      std::cerr << "\"\ncould not be opened. I guess it does not exist.\n";
      exit (-101);
    }
    
    FILE *os_in_list; // (out_in_list);

    os_in_list     = fopen(out_in_list, write_mode);

    CSequence_Mol seq(CSequence_Mol::unknown);
    while (!infile.eof())
    {
      seq.readRawFastaSequence(infile);
      if ( infile.fail() && !infile.eof() )
      {
	std::cerr << "\n\n";
	std::cerr << "An error occurred while reading the input file. It might not be a valid fasta file.\n";
	std::cerr << "File position: line ";
	std::cerr << faststring(infile.line()) << std::endl;
	std::cerr << std::flush;
      }

      // Is the sequence we just read in the list of sequence names??
      // This uses the full name

      faststring fname = seq.getFullName();

      CSeqNameListIterator it;
      it = m.find(&fname);    // We have to take the address operator here. Remember that we use less_than_pointer_to_faststring_struct to compare the elements, so we do no search for the pointer but for the string indeed. 

      if (it != m.end() )
      {
	seq.writeSequence_fasta(os_in_list);
      }
    } // END while
    
    infile.ffclose();
    fclose(os_in_list);
    
  } // END use_list_as_filter









  void print(std::ostream &os, short flag=0)
  {
    // With flag == 0 we only print the size of the list.

    os << "This CSeqNameList has " << full_names.size() << " entries.\n";

    if (flag > 0)
      os << "flag: " << flag << std::endl;

    if (flag == 1)
    {
      unsigned i=0, n=full_names.size();
      
      os << "n: " << n << std::endl;
      os << "only_names: " << only_names << std::endl;

      if (only_names)
      {
	while (i<n)
	{
	  os << "Name: " << *(full_names[i]) 
	     << std::endl;
	  ++i;
	}
      }
      else  // !only_names
      {
	while (i<n)
	{
	  os << "Name: " << *(full_names[i]) 
	     << " " 
	     << seq_len_vec[i] 
	     << "  " 
	     << CG_content[i]
	     << "  "
	     << corrected_len_vec[i]
	     << "  "
	     << repeat_density_bp_mbp[i]
	     << std::endl;
	  ++i;
	}
      }
    }
    else if (flag ==2) // fasta format:
    {
      unsigned i=0, n=full_names.size();

      while (i<n)
      {
	os << ">" << *(full_names[i]) 
	   << std::endl;
	++i;
      }
    }


  }

  unsigned size()
  {
    return full_names.size();
  }


  double get_CG_content(unsigned id)
  {
    if (only_names)
      return 0;
    else
      return CG_content[id];
  }

  const char *get_name(unsigned id)
  {
      return full_names[id]->c_str();
  }

  unsigned get_seq_length(unsigned id)
  {
    if (only_names)
      return 0;
    else
      return seq_len_vec[id];
  }

  unsigned get_corrected_len(unsigned id)
  {
    if (only_names)
      return 0;
    else
      return corrected_len_vec[id];
  }

  unsigned get_repeat_content(unsigned id)
  {
    if (only_names)
      return 0;
    else
      return repeat_density_bp_mbp[id];

  }

  //XXXXXXXXXXXXXXXXXXXX
  unsigned get_id(faststring str)
  {
    CSeqNameListIterator  it_find;
    it_find = m.find(&str);   // We have to take the address operator here. Remember that we use less_than_pointer_to_faststring_struct to compare the elements, so we do no search for the pointer but for the string indeed. 
    if (it_find != m.end())
      return it_find->second;
    else
      return -1u;
  }


  // Search for str_in in the list of sequences. str_in can be the full name or the short name.
  // The target list can contain full or short names.
  unsigned get_id_short_name(faststring str) // No refernce since we change the string
  {
    // We shorten str to the short name.
    // If it is the full name and the target only contains
    // the short name this is important.
    set_to_short_name(str);

    CSeqNameListIterator  it_find;
    it_find = m.lower_bound(&str);  // it_find points to an element that is either equal to str or the first element
                                    // in the map that is greater than str. In case, str contains the short name of the
                                    // target, we are looking for this element.
                                    // The question whether the short name is unique is not adressed here. The id of the 
                                    // first suitable short name is returned.
    if (it_find == m.end())
      return -1u;

    //    if (*(it_find->first) == str)      // str is equal to the string that has been found. 
	    //     return it_find->second;

    // str is not equal to the string that has been found, but it might be equal to its short name.
    unsigned pos1 = it_find->first->find(' ');
    if (pos1 == faststring::npos)
      pos1 = it_find->first->size();
    unsigned pos2 = str.size();

    // The short names differ or the complete string is the short name. Since the short name should have been found above,
    // we have not found the string.
    if (pos1 != pos2) 
      return -1u;                
    
    // If str is longer than the target, it cannot be the short name.
    // If the space is found at a position with index larger thant the size of str, it cannot be the short name of str.
    // It could be that this case is not possible.
    //     if (str.size() > it_find->first->size() || pos > str.size() )
    //        return -1u;

    if (strncmp(str.c_str(), it_find->first->c_str(), pos1) != 0)
      return -1u;
    
    return  it_find->second;
  }

  std::vector<faststring*>& get_vec_of_seq_names()
  {
    return full_names;
  }

  void copy_to_this_if_Full_name_has_a_match(CSeqNameList &a, faststring str)
  {
    CSeqNameListIterator it_beg, it_end;
    it_beg = a.m.begin();
    it_end = a.m.end();

    size_t find_pos;
    unsigned id;

    // Copy a to this:
    while (it_beg != it_end)
    {
      find_pos = it_beg->first->find(str);

      if (find_pos != faststring::npos) // we have a match
      {
	id = it_beg->second;
	add(a, id);
      }
	++it_beg;
    }
  }


  void set_to_union_of(CSeqNameList &a, CSeqNameList &b)
  {
    // The use of set_union does not seem to be better than this solution:

    set_to(a, "set union");

/*     CSeqNameListIterator it_beg, it_end; */
/*     unsigned id; */

/*     infile_name = "set union"; */

/*     it_beg = a.m.begin(); */
/*     it_end = a.m.end(); */

/*     // Copy a to this: */
/*     while (it_beg != it_end) */
/*     { */
/*       id = it_beg->second; */
/*       add(a, id); */
/*       ++it_beg; */
/*     } */

    // Copy b to this, for all elements that have not been added yet:

    add_List_non_redundant(b);
  }

  void set_to(CSeqNameList &a, faststring name)
  {
    reset();
    infile_name = name;
    only_names = a.only_names;
    
    CSeqNameListIterator it_beg, it_end;
    unsigned id;

    it_beg = a.m.begin();
    it_end = a.m.end();

    // Copy a to this:
    while (it_beg != it_end)
    {
      id = it_beg->second;
      add(a, id);
      ++it_beg;
    }
   
  }

  void add_List_non_redundant(CSeqNameList &b)
  {
    CSeqNameListIterator it_beg, it_end;
    unsigned id;

    if (only_names != b.only_names)
    {
      std::cerr << "CSeqNameList objects can only be added to another list of the only_names field match." << std::endl;
    }

    it_beg = b.m.begin();
    it_end = b.m.end();

    while (it_beg != it_end)
    {
      id = get_id_short_name(*(it_beg->first));    // Search for sequence in this.
      if (id == -1u)                               // Sequence is unknown in this, so we add it
      {
	add(b, it_beg->second);
      }
      ++it_beg;
    }  
  }

  bool is_only_names()
  {
    return only_names;
  }


  void set_to_intersection_of(CSeqNameList &a, CSeqNameList &b)
  {
    infile_name = "set intersection";

    // use set_intersection
    //    typedef std::set<faststring> faststringSet;

    //   faststringSet            out_fs;

    mySet   out_fs;

    mySetIterator  out_itr( out_fs.begin() );

    mySet a_short;
    mySet b_short;

    CSeqNameListIterator it1, it2;

    faststring tmp;

    it1 = a.m.begin();
    it2 = a.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      a_short.insert(tmp);
      ++it1;
    }

    it1 = b.m.begin();
    it2 = b.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      b_short.insert(tmp);
      ++it1;
    }

    set_intersection(a_short.begin(), a_short.end(),
		     b_short.begin(), b_short.end(),
		     std::inserter( out_fs, out_itr ) );

    mySetIterator   it_beg, it_end;
    unsigned        id;

    it_beg = out_fs.begin();
    it_end = out_fs.end();

    while (it_beg != it_end)
    {
      id = a.get_id_short_name(*it_beg);
      add(a, id);

      ++it_beg;
    }
  }

  void set_to_diff_of(CSeqNameList &a,CSeqNameList &b)
  {
    infile_name = "set diff";

    typedef std::map<faststring*, unsigned> myMap;
    typedef std::set<faststring>            mySet;
    mySet    out_fs;

    mySetIterator  out_itr( out_fs.begin() );

    mySet a_short;
    mySet b_short;

    CSeqNameListIterator it1, it2;

    faststring tmp;

    it1 = a.m.begin();
    it2 = a.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      a_short.insert(tmp);
      ++it1;
    }

    it1 = b.m.begin();
    it2 = b.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      b_short.insert(tmp);
      ++it1;
    }

    set_difference(a_short.begin(), a_short.end(),
		   b_short.begin(), b_short.end(),
		   std::inserter( out_fs, out_itr ) );

    mySetIterator it_beg, it_end;
    unsigned                id;

    it_beg = out_fs.begin();
    it_end = out_fs.end();

    //    std::cerr << "out_fs" << out_fs.size() << std::endl;

    while (it_beg != it_end)
    {
      id = a.get_id_short_name(*it_beg);
      add(a, id);

      ++it_beg;
    }
  }

  void set_to_sym_diff_of(CSeqNameList &a,CSeqNameList &b)
  {
   infile_name = "set sym diff";

    typedef std::map<faststring*, unsigned> myMap;
    typedef std::set<faststring>            mySet;
    mySet                                   out_fs;

    mySetIterator  out_itr( out_fs.begin() );

    mySet a_short;
    mySet b_short;

    CSeqNameListIterator it1, it2;

    faststring tmp;

    it1 = a.m.begin();
    it2 = a.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      a_short.insert(tmp);
      ++it1;
    }

    it1 = b.m.begin();
    it2 = b.m.end();
    while (it1 != it2)
    {
      tmp = *(it1->first);
      set_to_short_name(tmp);
      b_short.insert(tmp);
      ++it1;
    }

    set_symmetric_difference(a_short.begin(), a_short.end(),
			     b_short.begin(), b_short.end(),
			     std::inserter( out_fs, out_itr ) );

    mySetIterator it_beg, it_end;
    unsigned        id;

    it_beg = out_fs.begin();
    it_end = out_fs.end(); 

    while (it_beg != it_end)
    {
      // The seq can be found either in a or in b.
      id = a.get_id_short_name(*it_beg);
      if (id != -1u)
	add(a, id);
      else
      {
	id = b.get_id_short_name(*it_beg);
	add(b, id);
      }

      ++it_beg;
    }

  }


  // Searches for sequence names found in this in list a. If present in list a add it to list c.
  // If not, search for it in list b. If found add it to d. Oherwise produce error.
  void divide_this_uppon_existence_in_first_two_sets(CSeqNameList &a, CSeqNameList &b,
						     CSeqNameList &c, CSeqNameList &d)
  {
    c.reset();
    d.reset();

    CSeqNameListIterator it_beg, it_end;

    it_beg = m.begin();
    it_end = m.end();

    unsigned id;

    while (it_beg != it_end)
    {
      id = a.get_id_short_name(*(it_beg->first));
      
      if (id != -1u)
      {
	c.add(a, id);
      }
      else
      {
	id = b.get_id_short_name(*(it_beg->first));
	if (id != -1u)
	{
	  d.add(b, id);
	}
	else
	{
	  std::cerr << "Internal problem: Sequence found in neither file, even though I think it should be in one of the files: " << it_beg->first << std::endl;
	}
      }
      ++it_beg;
    }
  }

  // All sequence names are split according to delim. The field with 1-based number field is actracted and stored in the field_set.
  void get_set_of_name_field(unsigned field, char delim, std::set<faststring> &field_set)
  {
    unsigned i, n=full_names.size();
    std::vector<faststring> splitter;

    for (i=0; i<n; ++i) // foreach full name
    {
      split(splitter, *full_names[i], "|");
      if (splitter.size() < field )
      {
	std::cerr << "Warning: In get get_set_of_name_field you tried to extract field # "
		  << field << " after splitting " << full_names[i] << " with deliminator " << delim << std::endl;
      }
      else
      {
	field_set.insert(splitter[field-1]);  // Convert 1-based field number to 0-based index.
      }
    }
  }


};  // END CSeqNameList


#endif
