//
// File:	AbstractStream.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
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
