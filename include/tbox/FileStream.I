//
// File:	FileStream.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Simple class to read/write files in XDR format for portability
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
void FileStream::closeFileOnExit(const bool flag)
{
   d_close_on_exit = flag;
}

inline
FILE *FileStream::getFILE()
{
   return(d_FILE);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

