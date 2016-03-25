//
// File:	FileStream.C
// Package:	SAMRAI communication and data transfer package
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Simple class to read/write files in XDR format for portability
//

#include "tbox/FileStream.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "tbox/FileStream.I"
#endif

namespace SAMRAI {
   namespace tbox {


/*
*************************************************************************
*									*
* Constructor and destructor for the file stream.  The constructor	*
* probably should not throw an exception if the file cannot be opened,	*
* but this class will disappear later, so this behavior is OK for now.	*
*									*
*************************************************************************
*/

FileStream::FileStream(const char *filename, const StreamMode mode)
{
   d_close_on_exit = true;

   const char *fmode = ((mode == FileStream::Read)  ? "rb" :
                        (mode == FileStream::Write) ? "wb" : "ab");
   if (!(d_FILE = fopen(filename, fmode))) {
      TBOX_ERROR("FileStream: Unable to open file ``" 
                 << filename << "''...\n");
   }

#ifdef HAVE_XDR
   const xdr_op xop = ((mode==FileStream::Read) ? XDR_DECODE : XDR_ENCODE);
   xdrstdio_create(&d_xdr_stream, d_FILE, xop);
   setXDRStream(&d_xdr_stream);
#endif

}

FileStream::FileStream(FILE *file, const StreamMode mode)
{
   d_close_on_exit = false;

   d_FILE = file;

#ifdef HAVE_XDR
   const xdr_op xop = ((mode==FileStream::Read) ? XDR_DECODE : XDR_ENCODE);
   xdrstdio_create(&d_xdr_stream, d_FILE, xop);
   setXDRStream(&d_xdr_stream);
#endif
}

FileStream::~FileStream()
{
#ifdef HAVE_XDR
#ifndef LACKS_PROPER_XDR_HEADER
   xdr_destroy(&d_xdr_stream);
#else
   if (d_xdr_stream.x_ops->x_destroy)
      (*(void(*)(XDR*))(d_xdr_stream.x_ops->x_destroy))(&d_xdr_stream);
#endif
#endif

   if (d_close_on_exit) {
      (void) fclose(d_FILE);
   }
}

}
}
