//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/templates/special/MathUtilitiesSpecial.C $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1777 $
// Modified:	$LastChangedDate: 2007-12-13 16:51:06 -0800 (Thu, 13 Dec 2007) $
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

template<> char   MathUtilities<char>::s_zero           = 0;
template<> char   MathUtilities<char>::s_one            = 1;
template<> char   MathUtilities<char>::s_signaling_nan  = CHAR_MAX;
template<> char   MathUtilities<char>::s_max            = CHAR_MAX;
template<> char   MathUtilities<char>::s_min            = CHAR_MIN;
template<> char   MathUtilities<char>::s_epsilon        = 1;

template<> int    MathUtilities<int>::s_zero           = 0;
template<> int    MathUtilities<int>::s_one            = 1;
template<> int    MathUtilities<int>::s_signaling_nan   = INT_MAX;
template<> int    MathUtilities<int>::s_max             = INT_MAX;
template<> int    MathUtilities<int>::s_min             = INT_MIN;
template<> int    MathUtilities<int>::s_epsilon         = 1;

template<> float  MathUtilities<float>::s_zero          = 0.0;
template<> float  MathUtilities<float>::s_one           = 1.0;
template<> float  MathUtilities<float>::s_signaling_nan = SAMRAI_FLT_SNAN;
template<> float  MathUtilities<float>::s_max           = FLT_MAX;
template<> float  MathUtilities<float>::s_min           = FLT_MIN;
template<> float  MathUtilities<float>::s_epsilon       = FLT_EPSILON;

template<> double MathUtilities<double>::s_zero          = 0.0;
template<> double MathUtilities<double>::s_one           = 1.0;
template<> double MathUtilities<double>::s_signaling_nan = SAMRAI_DBL_SNAN;
template<> double MathUtilities<double>::s_max           = DBL_MAX;
template<> double MathUtilities<double>::s_min           = DBL_MIN;
template<> double MathUtilities<double>::s_epsilon       = DBL_EPSILON;

template<> dcomplex   MathUtilities<dcomplex>::s_zero             = dcomplex(0.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_one              = dcomplex(1.0,0.0);
template<> dcomplex   MathUtilities<dcomplex>::s_signaling_nan  = dcomplex(SAMRAI_DBL_SNAN,SAMRAI_DBL_SNAN);
template<> dcomplex   MathUtilities<dcomplex>::s_max            = dcomplex(DBL_MAX,DBL_MAX);
template<> dcomplex   MathUtilities<dcomplex>::s_min            = dcomplex(DBL_MIN,DBL_MIN);
template<> dcomplex   MathUtilities<dcomplex>::s_epsilon        = dcomplex(DBL_MIN,0.0);

template<>
bool MathUtilities<float>::isNaN(const float& value)
{
   int i = isnan(value);
   return( (i != 0) ? true : false );
}

template<> 
bool MathUtilities<double>::isNaN(const double& value)
{
   int i = isnan(value);
   return( (i != 0) ? true : false );
}

template<>
bool MathUtilities<dcomplex>::isNaN(const dcomplex& value)
{
   int i_re = isnan( real(value) );
   int i_im = isnan( imag(value) );
   return( ( (i_re != 0) || (i_im !=0) ) ? true : false );
}

template<>
bool MathUtilities<float>::equalEps(const float& a, const float& b)
{
   float absmax = MathUtilities<float>::Max(
                     MathUtilities<float>::Abs(a),
                     MathUtilities<float>::Abs(b) );
   float numerator = MathUtilities<float>::Abs(a-b);
   float denomenator =
      MathUtilities<float>::Max(absmax,
           MathUtilities<float>::s_epsilon);

   return( numerator/denomenator < sqrt(MathUtilities<float>::s_epsilon) );
}

template<>
bool MathUtilities<double>::equalEps(const double& a, const double& b)
{
   double absmax = MathUtilities<double>::Max(
                      MathUtilities<double>::Abs(a),
                      MathUtilities<double>::Abs(b) );
   double numerator = MathUtilities<double>::Abs(a-b);
   double denomenator =
      MathUtilities<double>::Max(absmax,
           MathUtilities<double>::s_epsilon);

   return( numerator/denomenator < sqrt(MathUtilities<double>::s_epsilon) );
}

template<>
bool MathUtilities<dcomplex>::equalEps(const dcomplex& a, const dcomplex& b)
{
   double a_re = real(a);
   double a_im = imag(a);
   double b_re = real(b);
   double b_im = imag(b);

   return( MathUtilities<double>::equalEps(a_re,b_re) && 
           MathUtilities<double>::equalEps(a_im,b_im) );
}

template<> 
dcomplex MathUtilities<dcomplex>::Min(dcomplex a, dcomplex b)
{
   return(norm(a) < norm(b) ? a : b);
}

template<> 
dcomplex MathUtilities<dcomplex>::Max(dcomplex a, dcomplex b)
{
   return(norm(a) > norm(b) ? a : b);
}

template<>
int MathUtilities<int>::Abs(int a)
{
   return(a > 0 ? a : -a);
}

template<>
float MathUtilities<float>::Abs(float a)
{
   return(a > 0.0 ? a : -a);
}

template<>
double MathUtilities<double>::Abs(double a)
{
   return(a > 0.0 ? a : -a);
}

template<> 
bool MathUtilities<bool>::Rand(const bool& low, const bool& width)
{
   (void) low;
   (void) width;
   return( mrand48() > 0 ? true : false );
}

template<> 
char MathUtilities<char>::Rand(const char& low, const char& width)
{
   return static_cast<char>( static_cast<double>(width) * drand48() ) + low;
}

template<> 
int MathUtilities<int>::Rand(const int& low, const int& width)
{
   return static_cast<int>( static_cast<double>(width) * drand48() ) + low;
}

template<> 
float MathUtilities<float>::Rand(const float& low, const float& width)
{
   return width * drand48() + low;
}

template<> 
double MathUtilities<double>::Rand(const double& low, const double& width)
{
   return width * drand48() + low;
}

template<> 
dcomplex MathUtilities<dcomplex>::Rand(const dcomplex& low, const dcomplex& width)
{
   double real_part = real(width) * drand48() + real(low);
   double imag_part = imag(width) * drand48() + imag(low);
   return dcomplex(real_part, imag_part);
}

}
}

