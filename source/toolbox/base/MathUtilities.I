//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/base/MathUtilities.I $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2041 $
// Modified:	$LastChangedDate: 2008-03-11 15:10:41 -0700 (Tue, 11 Mar 2008) $
// Description:	Utilities class to access common POSIX constants and math ops
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

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
bool MathUtilities<TYPE>::isNaN(const TYPE& value)
{
   (void) value;
   return(false);
}

template<class TYPE> inline
bool MathUtilities<TYPE>::equalEps(const TYPE& a, const TYPE& b)
{
   return(a == b);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::Min(TYPE a, TYPE b)
{
   return(a < b ? a : b);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::Max(TYPE a, TYPE b)
{
   return(a > b ? a : b);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::Abs(TYPE value)
{
   return(value);
}

template<class TYPE> inline 
TYPE MathUtilities<TYPE>::round(TYPE x)
{
   return x;
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
