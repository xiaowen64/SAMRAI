//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/stream/FileStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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

