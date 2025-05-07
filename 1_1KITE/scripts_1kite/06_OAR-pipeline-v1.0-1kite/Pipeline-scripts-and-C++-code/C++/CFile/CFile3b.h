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

#ifndef CFILE_H
#define CFILE_H

#include <fstream>
#include <iostream>
#include <string>
#include "../faststring2.h"
#include <vector>

#define BUFFERSIZE 1000000
#define UNDOs      1


// Changes:
// 22.08.2009: Added getline function for faststring. Disadvantage: addidtional dependence on faststring.h
// 07.01.2011: Changed Version to CFile2_1.h
// 07.01.2011: Now uses the faststring2.h

// 30.04.2011: New version 3.0. a/b use faststring.h/faststring2.h
//             New feature: A file with name STDIN reads from stdin
//                          A call to open_STDIN() has the same effect.
//                          This allows programs that use CFILE to read from STDIN
//                          without any further work.
// 29.02.2012: New function: readFileIntoVectorOfStrings 
// Todo: Check that this new version is not significantly slower than preveious versions.


// Idee: ios::   set Buffer size

// TODO: Old mac format not supported.
//       This requires to allow two successive calls to the internal ungetchar command.


class CFile
{
  std::istream *this_is;
  bool          is_ifstream;

 private:
  char      buffer[BUFFERSIZE];
  char*     buffer_end;
  char*     buffer_pos;
  
  char      __status;
  unsigned  __line;

  //  bool      __eof;
  //  bool      __openFailed;

  unsigned fill_buffer(unsigned overlap)
  {
    unsigned          good_overlap = buffer_end - buffer;
    std::streamsize   n;

    if (good_overlap < overlap)
    {
      overlap = good_overlap;
    }

    if (overlap > 0)
      std::memmove(buffer, buffer_end - overlap, overlap);

    if (is_ifstream)
      ((std::ifstream *) this_is)->read(buffer + overlap, BUFFERSIZE - overlap);
    else
      this_is->read(buffer + overlap, BUFFERSIZE - overlap);
    n = this_is->gcount();

    if ( n == 0 )
    {
      __status |=  __eof_flag;     // Set eof flag
      __status &= ~__good_flag;    // Unset good flag
      __status |=  __fail_flag;    // Setting the fail flag is not always correct. Needs to be unset in routines that read more than one char. 
    }

    buffer_pos = buffer+overlap;
    buffer_end = buffer_pos + n;

    return overlap;
  }

  char getchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if ( buffer_pos == buffer_end )
    {
      fill_buffer(UNDOs);
      if (__status & __fail_flag)
	return '\0';
    }
    return *buffer_pos++;
  }

  void ungetchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if (buffer_pos != buffer)
      --buffer_pos;
    else
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }
  }

 public:

  enum {__eof_flag = 1, __good_flag = 2,  __fail_flag = 4, __bad_flag = 8};

  void open(const char *name)
  {
    ffopen(name);
  }

  void open(std::string name)
  {
    ffopen(name);
  }


  void open_STDIN()
  {
     __status = __good_flag;

     this_is = &std::cin;
     is_ifstream = false;

     if ( this->fail() )
     {
       __status |=  __fail_flag;      // Set fail flag.
       __status &= ~__good_flag;     // Unset good flag.
     }

     __line        = 1;
     buffer_end    = buffer;
     buffer_pos    = buffer;   
  }

  void ffopen(const char *name)
  {
    __status = __good_flag;

    if (strcmp(name, "STDIN")==0)
    {
      this_is = &std::cin;
      is_ifstream = false;
    }
    else
    {
      this_is = new std::ifstream(name);
      is_ifstream = true;
    }

    if ( this_is->fail() )
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }

    __line        = 1;
    buffer_end    = buffer;
    buffer_pos    = buffer;
  }

  void ffopen(std::string name)
  {
    ffopen(name.c_str());
  }

  void ffclose()
  {
    if (is_ifstream)
      ((std::ifstream *)this_is)->close();
  }

  void close()
  {
    ffclose();
  }

  bool exists()  // -- depricated -- do not use this any more - check fail() instead !!!!!!!!
  {
    return !fail(); 
  }

  bool fail()
  {
    return (__status & __fail_flag);
  }

  bool good()
  {
    return (__status & __good_flag);
  }

  unsigned line()
  {
    return __line;
  }

  bool eof()
  {
    return (__status & __eof_flag);
  }

  char status()
  {
    return __status;
  }

  void rewind()
  {
    this_is->clear();
    this_is->seekg (0, std::ios::beg);
  }

  void clear(char s = __good_flag)
  {
    __status = s;
    this_is->clear();
  }

  void ungetchar()
  {
    if (*(buffer_pos-1) == '\n' && !(__status & __fail_flag))
      --__line;
    ungetchar_intern();
  }

  void ignore(int delim = -1)
  {
    while (__status == __good_flag && getchar() != delim ){}
  }

  char peekchar()
  {
    char c = getchar();
    ungetchar();
    return c;
  }

  char peek()
  {
    return peekchar();
  }

  char getchar()
  {
    register char c;

    c = getchar_intern();

    if ( c < 14 && !(__status & __fail_flag) )
    {
      if (  c == '\r' )
      {
	// Overwrite the last reading position that contained the \r with \n
	*(buffer_pos-1) = '\n';           // Should always be valid! Works for 1 Undo
	c = getchar_intern();
	if ( c != '\n' && !(__status & __fail_flag) )
	{
	  // 	std::cerr << "Old mac file format currently not supported." << std::endl;
	  ungetchar();    /* old mac format, else dos format     */
	}
	c = '\n';
	++__line;
      }
      else if ( c == '\n')
	++__line;
    }
    return c;  
  }


  char getrawchar()
  {
    return getchar_intern();
  }

  void getline(faststring& str, char delim='\n')
  {
    char c = getchar();
    str.clear();
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && str.size() > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(std::string& str, char delim='\n')
  {
    char c = getchar();
    str="";
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && str.size() > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(char* cstr, unsigned n, char delim='\n')
  {
    char     c = getchar();
    unsigned i = 0;

    while ( !(__status & __fail_flag ) && i < n-1 &&  c != delim)
    {
      cstr[i] = c;
      ++i;
      c = getchar(); 
    }
    if ((__status & __fail_flag) && i > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
    cstr[i] = '\0';
  }

  void readFileIntoString(faststring &str)
  {
    char c;
    
    str.clear();
    
    c = getchar();
    while (!(__status & __fail_flag))
    {
      str.push_back(c);
      c = getchar();
    }
  }

  void readFileIntoVectorOfStrings(std::vector<faststring> &vos)
  {
    faststring tmp;

    getline(tmp);
    while (!(__status & __fail_flag))
    {
      vos.push_back(tmp);
      getline(tmp);
    }
  }
  

  char lastchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }

    /* It could be a '\r' in mac format */
    if (c == '\r')
      c = '\n';

    return c;
  }

  char lastrawchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }
    return c;
  }


};




#endif
