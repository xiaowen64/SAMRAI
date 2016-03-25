//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/stream/AbstractStream.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1917 $
// Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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
