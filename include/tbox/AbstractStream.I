//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/stream/AbstractStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Abstract base class for data stream abstraction
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
AbstractStream::AbstractStream()
{
}

inline
int AbstractStream::roundXDR(const int n)
{
   const int xdr_bytes = 4;
   return(xdr_bytes*((n+xdr_bytes-1)/xdr_bytes));
}

inline
int AbstractStream::sizeofBool(const int n)
{
   return(roundXDR(n));
}

inline
int AbstractStream::sizeofChar(const int n)
{
   return(roundXDR(n));
}

inline
int AbstractStream::sizeofDoubleComplex(const int n)
{
   return(roundXDR(2*n*sizeof(double)));
}

inline
int AbstractStream::sizeofDouble(const int n)
{
   return(roundXDR(n*sizeof(double)));
}

inline
int AbstractStream::sizeofFloat(const int n)
{
   return(roundXDR(n*sizeof(float)));
}

inline
int AbstractStream::sizeofInt(const int n)
{
   return(roundXDR(n*sizeof(int)));
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
