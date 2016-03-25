//
// File:	MathUtilitiesSpecial.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	MathUtilities routines to set up handlers and get signaling NaNs
//

#include "tbox/MathUtilities.h"

#include <float.h>
#include <math.h>
#include <limits.h>

#include "tbox/Complex.h"

/*
 * Floating point exception handling.  
 * The following lines setup exception handling header files on 
 * systems other than solaris.
 */
#if defined(HAVE_EXCEPTION_HANDLING) 
#include <stdlib.h>
#include <stdio.h>
#include <fpu_control.h>
#include <signal.h>
#endif

/*
 * The following lines setup exception handling headers on the Sun.  If we
 * use Sun's native compiler, just pull in the <sunmath.h> include file.
 * If we are under solaris but use a different compiler (e.g. KCC, g++)
 * we have to explicitly define the functions that <sunmath.h> defines,
 * since we don't have access to this file.
 */
#ifdef __SUNPRO_CC
#include <sunmath.h>
#endif

namespace SAMRAI {
   namespace tbox {

/*
 *  Settings for the various signaling NaNs on different systems
 */
#if !defined(FLT_SNAN_IS_BROKEN)  
#define SAMRAI_FLT_SNAN   FLT_SNAN
#elif !defined(FLT_MAX_IS_BROKEN)
#define SAMRAI_FLT_SNAN   FLT_MAX
#else
#define SAMRAI_FLT_SNAN   NAN
#endif

#if !defined(DBL_SNAN_IS_BROKEN)
#define SAMRAI_DBL_SNAN   DBL_SNAN
#elif !defined(DBL_MAX_IS_BROKEN)
#define SAMRAI_DBL_SNAN   DBL_MAX
#else
#define SAMRAI_DBL_SNAN   NAN
#endif

template<> bool   MathUtilities<bool>::s_zero           = false;
template<> bool   MathUtilities<bool>::s_one            = true;
template<> bool   MathUtilities<bool>::s_signaling_nan  = false;
template<> bool   MathUtilities<bool>::s_max            = true;
template<> bool   MathUtilities<bool>::s_min            = false;
template<> bool   MathUtilities<bool>::s_epsilon        = true;
template<> bool   MathUtilities<bool>::s_undefined      = false;

template<> char   MathUtilities<char>::s_zero           = 0;
template<> char   MathUtilities<char>::s_one            = 1;
template<> char   MathUtilities<char>::s_signaling_nan  = CHAR_MAX;
template<> char   MathUtilities<char>::s_max            = CHAR_MAX;
template<> char   MathUtilities<char>::s_min            = CHAR_MIN;
template<> char   MathUtilities<char>::s_epsilon        = 1;
template<> char   MathUtilities<char>::s_undefined      = (char) 0xff;

template<> int    MathUtilities<int>::s_zero           = 0;
template<> int    MathUtilities<int>::s_one            = 1;
template<> int    MathUtilities<int>::s_signaling_nan   = INT_MAX;
template<> int    MathUtilities<int>::s_max             = INT_MAX;
template<> int    MathUtilities<int>::s_min             = INT_MIN;
template<> int    MathUtilities<int>::s_epsilon         = 1;
template<> int    MathUtilities<int>::s_undefined       = INT_MAX;

template<> float  MathUtilities<float>::s_zero          = 0.0;
template<> float  MathUtilities<float>::s_one           = 1.0;
template<> float  MathUtilities<float>::s_signaling_nan = SAMRAI_FLT_SNAN;
template<> float  MathUtilities<float>::s_max           = FLT_MAX;
template<> float  MathUtilities<float>::s_min           = FLT_MIN;
template<> float  MathUtilities<float>::s_epsilon       = FLT_EPSILON;
template<> float  MathUtilities<float>::s_undefined     = SAMRAI_FLT_SNAN;

template<> double MathUtilities<double>::s_zero          = 0.0;
template<> double MathUtilities<double>::s_one           = 1.0;
template<> double MathUtilities<double>::s_signaling_nan = SAMRAI_DBL_SNAN;
template<> double MathUtilities<double>::s_max           = DBL_MAX;
template<> double MathUtilities<double>::s_min           = DBL_MIN;
template<> double MathUtilities<double>::s_epsilon       = DBL_EPSILON;
template<> double MathUtilities<double>::s_undefined     = SAMRAI_DBL_SNAN;

template<> dcomplex   MathUtilities<dcomplex>::s_zero             = dcomplex(0.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_one              = dcomplex(1.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_signaling_nan  = dcomplex(SAMRAI_DBL_SNAN,SAMRAI_DBL_SNAN);
template<> dcomplex   MathUtilities<dcomplex>::s_max            = dcomplex(DBL_MAX,DBL_MAX);
template<> dcomplex   MathUtilities<dcomplex>::s_min            = dcomplex(DBL_MIN,DBL_MIN);
template<> dcomplex   MathUtilities<dcomplex>::s_epsilon        = dcomplex(DBL_MIN,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_undefined      = dcomplex(SAMRAI_DBL_SNAN,SAMRAI_DBL_SNAN);


template<>
bool MathUtilities<float>::isNaN(const float &value)
{
   int i = isnan(value);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

template<> 
bool MathUtilities<double>::isNaN(const double &value)
{
   int i = isnan(value);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

template<> 
dcomplex MathUtilities<dcomplex>::Min( dcomplex a, dcomplex b )
{
   return(norm(a) < norm(b) ? a : b);
}

template<> 
dcomplex MathUtilities<dcomplex>::Max( dcomplex a, dcomplex b )
{
   return(norm(a) > norm(b) ? a : b);
}

template<> 
bool MathUtilities<bool>::Rand(bool low, bool width)
{
   return mrand48() > 0 ? true : false;
}

template<> 
char MathUtilities<char>::Rand(char low, char width)
{
   return (char)((double)width * drand48()) + low;
}

template<> 
int MathUtilities<int>::Rand(int low, int width)
{
   return (int)((double)width * drand48()) + low;
}

template<> 
float MathUtilities<float>::Rand(float low, float width)
{
   return width * drand48() + low;
}

template<> 
double MathUtilities<double>::Rand(double low, double width)
{
   return width * drand48() + low;
}

template<> 
dcomplex MathUtilities<dcomplex>::Rand(dcomplex low, dcomplex width)
{
   double real_part = real(width) * drand48() + real(low);
   double imag_part = imag(width) * drand48() + imag(low);
   return dcomplex(real_part, imag_part);
}

}
}

