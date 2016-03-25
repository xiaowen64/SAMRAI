//
// File:	ParallelBuffer.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Parallel I/O class buffer to manage parallel ostreams output
//

#include "tbox/ParallelBuffer.h"

#include <string>
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif
#include "tbox/Utilities.h"

#ifndef NULL
#define NULL 0
#endif

#define DEFAULT_BUFFER_SIZE (128)


namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* Construct a parallel buffer object.  The object will require further	*
* initialization to set up I/O streams and the prefix string.		*
*									*
*************************************************************************
*/

ParallelBuffer::ParallelBuffer()
{
   d_active        = true;
   d_prefix        = string();
   d_ostream1      = NULL;
   d_ostream2      = NULL;
   d_buffer        = NULL;
   d_buffer_size   = 0;
   d_buffer_ptr    = 0;
}

/*
*************************************************************************
*									*
* The destructor deallocates internal data buffer.  It does not modify  *
* the output streams.                                                   *
*									*
*************************************************************************
*/

ParallelBuffer::~ParallelBuffer()
{
   if (d_buffer) delete [] d_buffer;
}

/*
*************************************************************************
*									*
* Activate or deactivate the output stream.  If the stream has been	*
* deactivated, then deallocate the internal data buffer.                *
*									*
*************************************************************************
*/

void ParallelBuffer::setActive(bool active)
{
   if (!active && d_buffer) {
      delete [] d_buffer;
      d_buffer = NULL;
      d_buffer_size = 0;
      d_buffer_ptr = 0;
   }
   d_active = active;
}
     
/*
*************************************************************************
*									*
* Set the prefix that begins every new output line.                     *
*									*
*************************************************************************
*/

void ParallelBuffer::setPrefixString(const string &text)
{
   d_prefix = text;
}

/*
*************************************************************************
*									*
* Set the primary output stream.					*
*									*
*************************************************************************
*/

void ParallelBuffer::setOutputStream1(ostream *stream)
{
   d_ostream1 = stream;
}

/*
*************************************************************************
*									*
* Set the secondary output stream.					*
*									*
*************************************************************************
*/

void ParallelBuffer::setOutputStream2(ostream *stream)
{
   d_ostream2 = stream;
}

/*
*************************************************************************
*									*
* Output a string to the output stream by invoking the                  *
* outputString(string,length) method.                                   *
*									*
*************************************************************************
*/

void ParallelBuffer::outputString(const string &text)
{
   outputString(text, text.length());
}

/*
*************************************************************************
*									*
* Write a text string of the specified length to the output stream.	*
* Note that the string data is accumulated into the internal output	*
* buffer until an end-of-line is detected.				*
*									*
*************************************************************************
*/

void ParallelBuffer::outputString(const string &text, const int length)
{
   if ((length > 0) && d_active) {

      /*
       * If we need to allocate the internal buffer, then do so
       */

      if (!d_buffer) {
         d_buffer      = new char[DEFAULT_BUFFER_SIZE];
         d_buffer_size = DEFAULT_BUFFER_SIZE;
         d_buffer_ptr  = 0;
      }

      /*
       * If the buffer pointer is zero, then prepend the prefix if not empty
       */

      if ((d_buffer_ptr == 0) && !d_prefix.empty()) {
         copyToBuffer(d_prefix, d_prefix.length());
      }

      /*
       * Search for an end-of-line in the string 
       */

      int eol_ptr = 0;
      for ( ; (eol_ptr < length) && (text[eol_ptr] != '\n'); eol_ptr++) 
	 NULL_STATEMENT;

      /*
       * If no end-of-line found, copy the entire text string but no output
       */

      if (eol_ptr == length) {
         copyToBuffer(text, length);

      /*
       * If we found end-of-line, copy and output; recurse if more chars
       */

      } else {
         const int ncopy = eol_ptr+1;
         copyToBuffer(text, ncopy);
         outputBuffer();
         if (ncopy < length) {
            outputString(text.substr(ncopy), length-ncopy);
         }
      }
   }
}

/*
*************************************************************************
*									*
* Copy data from the text string into the internal output buffer.	*
* If the internal buffer is not large enough to hold all of the string	*
* data, then allocate a new internal buffer.				*
*									*
*************************************************************************
*/

void ParallelBuffer::copyToBuffer(const string &text, const int length)
{
   /*
    * First check whether we need to increase the size of the buffer
    */

   if (d_buffer_ptr+length > d_buffer_size) {
      const int new_size = Utilities::imax(d_buffer_ptr+length, 2*d_buffer_size);
      char *new_buffer = new char[new_size];

      if (d_buffer_ptr > 0) {
         (void) strncpy(new_buffer, d_buffer, d_buffer_ptr);
      }
      delete [] d_buffer;

      d_buffer      = new_buffer;
      d_buffer_size = new_size;
   }

   /*
    * Copy data from the input into the internal buffer and increment pointer
    */

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_buffer_ptr+length <= d_buffer_size);
#endif

   strncpy(d_buffer+d_buffer_ptr, text.c_str(), length);
   d_buffer_ptr += length;
}

/*
*************************************************************************
*									*
* Output buffered stream data to the active output streams and reset	*
* the buffer pointer to its empty state.				*
*									*
*************************************************************************
*/

void ParallelBuffer::outputBuffer()
{
   if (d_buffer_ptr > 0) {
      if (d_ostream1) {
         d_ostream1->write(d_buffer, d_buffer_ptr);
         d_ostream1->flush();
      }
      if (d_ostream2) {
         d_ostream2->write(d_buffer, d_buffer_ptr);
         d_ostream2->flush();
      }
      d_buffer_ptr = 0;
   }
}

/*
*************************************************************************
*									*
* Synchronize the parallel buffer and write string data.  This routine	*
* is called from streambuf.						*
*									*
*************************************************************************
*/

int ParallelBuffer::sync()
{
   const int n = pptr() - pbase();
   if (n > 0) outputString(pbase(), n);
   return(0);
}

/*
*************************************************************************
*									*
* Write the specified number of characters into the output stream.	*
* This routine is called from streambuf.  If this routine is not	*
* provided, then overflow() is called instead for each character.	*
*									*
* We only define this for g++ and KCC since those libraries define the	*
* streamsize type.  Note that this routine is not required; it only	*
* offers some efficiency over overflow().				*
*									*
*************************************************************************
*/

#if !defined(__INTEL_COMPILER) && (defined(__GNUG__) || defined(__KCC))
streamsize ParallelBuffer::xsputn(const string &text, streamsize n)
{
   sync();
   if (n > 0) outputString(text, n);
   return(n);
}
#endif

/*
*************************************************************************
*									*
* Write a single character into the parallel buffer.  This routine is	*
* called from streambuf.						*
*									*
*************************************************************************
*/

int ParallelBuffer::overflow(int ch)
{
   const int n = pptr() - pbase();
   if (n && sync()) {
      return(EOF);
   }
   if (ch != EOF) {
      char character[2];
      character[0] = ch;
      character[1] = 0;
      outputString(character, 1);
   }
   pbump(-n);
   return(0);
}

#ifdef _MSC_VER
// Should never read from here
int ParallelBuffer::underflow() 
{
   return EOF;
}
#endif

}
}

