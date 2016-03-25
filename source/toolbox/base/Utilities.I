//
// File:        Utilities.I
// Package:     SAMRAI toolbox
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: A collection of trivial utility functions such as min and max
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

namespace SAMRAI {
   namespace tbox {


inline
int Utilities::imin(const int a, const int b)
{
   return(a < b ? a : b);
}

inline
int Utilities::imax(const int a, const int b)
{
   return(a > b ? a : b);
}

inline
int Utilities::iabs(const int a)
{
   return(a > 0 ? a : -a);
}

inline
float Utilities::fmin(const float a, const float b)
{
   return(a < b ? a : b);
}

inline
float Utilities::fmax(const float a, const float b)
{
   return(a > b ? a : b);
}

inline
float Utilities::fabs(const float a)
{
   return(a > 0.0 ? a : -a);
}

inline
double Utilities::dmin(const double a, const double b)
{
   return(a < b ? a : b);
}

inline
double Utilities::dmax(const double a, const double b)
{
   return(a > b ? a : b);
}

inline
double Utilities::dabs(const double a)
{
   return(a > 0.0 ? a : -a);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

