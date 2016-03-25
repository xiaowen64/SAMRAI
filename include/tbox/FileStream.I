//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/stream/FileStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

