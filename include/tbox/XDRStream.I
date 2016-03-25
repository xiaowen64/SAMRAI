//
// File:	XDRStream.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
