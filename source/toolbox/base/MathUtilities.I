//
// File:	MathUtilities.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	MathUtilities routines to set up handlers and get signaling NaNs
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

#include <float.h>
#include <math.h>
#include <limits.h>

namespace SAMRAI {
   namespace tbox {

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getZero()
{
   return(s_zero);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getOne()
{
   return(s_one);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getSignalingNaN()
{
   return(s_signaling_nan);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getMax()
{
   return(s_max);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getMin()
{
   return(s_min);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getEpsilon()
{
   return(s_epsilon);
}

template<class TYPE> inline
TYPE MathUtilities<TYPE>::getUndefined()
{
   return(s_undefined);
}

template<class TYPE> inline
bool MathUtilities<TYPE>::isNaN(const TYPE &value)
{
   return(false);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::Min( TYPE a, TYPE b )
{
   return(a < b ? a : b);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::Max( TYPE a, TYPE b )
{
   return(a > b ? a : b);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
