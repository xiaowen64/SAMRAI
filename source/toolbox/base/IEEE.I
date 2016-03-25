//
// File:	IEEE.I
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	IEEE routines to set up handlers and get signaling NaNs
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
float IEEE::getSignalingFloatNaN()
{
   return(s_signaling_nan_float);
}

inline
double IEEE::getSignalingNaN()
{
   return(s_signaling_nan_double);
}

inline
int IEEE::getINT_MAX()
{
   return(s_int_max);
}

inline
int IEEE::getINT_MIN()
{
   return(s_int_min);
}

inline
float IEEE::getFLT_MAX()
{
   return(s_flt_max);
}

inline
float IEEE::getFLT_MIN()
{
   return(s_flt_min);
}

inline
float IEEE::getFLT_EPSILON()
{
   return(s_flt_epsilon);
}

inline
double IEEE::getDBL_MAX()
{
   return(s_dbl_max);
}

inline
double IEEE::getDBL_MIN()
{
   return(s_dbl_min);
}

inline
double IEEE::getDBL_EPSILON()
{
   return(s_dbl_epsilon);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}
