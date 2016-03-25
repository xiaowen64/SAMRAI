//
// File:	FileStream.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Simple class to read/write files in XDR format for portability
//

#ifndef included_tbox_FileStream
#define included_tbox_FileStream

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif
#ifndef included_tbox_XDRStream
#include "tbox/XDRStream.h"
#endif


namespace SAMRAI {
   namespace tbox {


/**
 * Class FileStream reads/writes files using XDRStream
 * to convert data into a machine independent format required for
 * portability across machine architectures.
 *
 * There are two constructors to a file stream: (1) provide a FILE pointer
 * to an existing file, or (2) provide a file name and let FileStream
 * open the object.  By default, FileStream calls fclose() on the file
 * only for case (2).  This default behavior can be changed by making
 * a call to closeFileOnExit().
 *
 * @see tbox::XDRStream
 */

class FileStream : public XDRStream
{
public:
   enum StreamMode { Read, Write, Append };

   /**
    * Open a file with the specified filename and the stream mode.  The
    * mode can be either FileStream::Read, FileStream::Write,
    * or FileStream::Append.  If the filename cannot be opened, then
    * the program aborts with an error.  By default, the filename will be
    * closed when the destructor for FileStream is called.
    */
   FileStream(const char *filename, const StreamMode mode);

   /**
    * Use an existing file for the XDR stream.  The stream mode can be
    * either FileStream::Read or FileStream::Write.  The append
    * mode is treated as a write.  By default, the file is NOT closed when
    * the destructor for FileStream is called.
    */
   FileStream(FILE *file, const StreamMode mode);

   /**
    * The virtual destructor for a file stream.  This optionally closes
    * the file depending on the status of the closeFileOnExit() flag.
    */
   virtual ~FileStream();

   /**
    * Flag to close or not close the file when the file streamfinishes.
    * The default value for this flag is true if the filename constructor
    * is used but false if the FILE constructor is used.
    */
   void closeFileOnExit(const bool flag);

   /**
    * Return the FILE data structure corresponding to the file stream.
    */
   FILE *getFILE();

private:
   FileStream(const FileStream&);	// not implemented
   void operator=(const FileStream&);	// not implemented

   XDR d_xdr_stream;
   FILE *d_FILE;
   bool d_close_on_exit;

};


}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/FileStream.I"
#endif
#endif
