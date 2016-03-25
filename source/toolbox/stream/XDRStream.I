//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/stream/XDRStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description:	Stream class that converts into XDR for portable communication
//

#ifndef NULL
#define NULL (0)
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {

inline
XDRStream::XDRStream()
{
   d_xdr_stream = (XDR *) NULL;
}

inline
void XDRStream::setXDRStream(XDR *xdrs)
{
   d_xdr_stream = xdrs;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
