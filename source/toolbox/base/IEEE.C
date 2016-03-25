//
// File:	IEEE.C
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	IEEE routines to set up handlers and get signaling NaNs
//

#include "tbox/IEEE.h"
#include "tbox/MPI.h"
#include <float.h>
#include <math.h>
#include <limits.h>

/*
 * Floating point exception handling.  
 * 
 * The following lines setup exception handling header files.
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

#ifdef DEBUG_NO_INLINE
#include "tbox/IEEE.I"
#endif

namespace SAMRAI {
   namespace tbox {


/*
 * Create the function invoked when an exception is tripped. 
 */
#if defined(HAVE_EXCEPTION_HANDLING)
static void error_action(int error) 
{
   fprintf(stderr, "floating point exception -- program abort!\n");
   abort();
   MPI::abort();
}
#endif

/*
 *  Settings for the various signaling NaNs on different systems
 */

#if !defined(FLT_SNAN_IS_BROKEN)  
float  IEEE::s_signaling_nan_float  = FLT_SNAN;
#elif !defined(FLT_MAX_IS_BROKEN)
float  IEEE::s_signaling_nan_float  = FLT_MAX;
#else
float  IEEE::s_signaling_nan_float  = NAN;
#endif

#if !defined(DBL_SNAN_IS_BROKEN)
double  IEEE::s_signaling_nan_double  = DBL_SNAN;
#elif !defined(DBL_MAX_IS_BROKEN)
double  IEEE::s_signaling_nan_double  = DBL_MAX;
#else
double  IEEE::s_signaling_nan_double  = NAN;
#endif

int    IEEE::s_int_max = INT_MAX;
int    IEEE::s_int_min = INT_MIN;
float  IEEE::s_flt_max = FLT_MAX;
float  IEEE::s_flt_min = FLT_MIN;
float  IEEE::s_flt_epsilon = FLT_EPSILON;
double IEEE::s_dbl_max = DBL_MAX;
double IEEE::s_dbl_min = DBL_MIN;
double IEEE::s_dbl_epsilon = DBL_EPSILON;
   


/*
*************************************************************************
*									*
* Set up the IEEE exception handlers so that normal IEEE exceptions	*
* will cause a program abort.  How this is done varies wildly from	*
* architecture to architecture. 					*
*************************************************************************
*/

void IEEE::setupExceptionHandlers()
{

#if defined(HAVE_EXCEPTION_HANDLING)
   unsigned short fpu_flags = _FPU_DEFAULT;          
   fpu_flags &= ~_FPU_MASK_IM;  /* Execption on Invalid operation */
   fpu_flags &= ~_FPU_MASK_ZM;  /* Execption on Division by zero  */
   fpu_flags &= ~_FPU_MASK_OM;  /* Execption on Overflow */
   _FPU_SETCW(fpu_flags);
   signal(SIGFPE, error_action);
#endif
}

/*
*************************************************************************
*									*
* Initialize float and double values to the signaling nan.              *
* Initialize int to INT_MAX.                                            *
*									*
*************************************************************************
*/

void IEEE::setNaN(float &f)
{  
   f = s_signaling_nan_float;
}

void IEEE::setNaN(double &d)
{  
   d = s_signaling_nan_double;
}

/*
*************************************************************************
*									*
* Initialize float and double arrays to signaling NaNs.			*
* Initialize int array to INT_MAX.                                      *
*									*
*************************************************************************
*/

void IEEE::initializeArrayToSignalingNaN(Array<float>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_signaling_nan_float;
   }
}

void IEEE::initializeArrayToSignalingNaN(Array<double>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_signaling_nan_double;
   }
}

void IEEE::initializeArrayToINT_MAX(Array<int>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_int_max;
   }
}

void IEEE::initializeArrayToINT_MIN(Array<int>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_int_min;
   }
}

void IEEE::initializeArrayToFLT_MAX(Array<float>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_flt_max;
   }
}

void IEEE::initializeArrayToFLT_MIN(Array<float>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_flt_min;
   }
}

void IEEE::initializeArrayToDBL_MAX(Array<double>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_dbl_max;
   }
}

void IEEE::initializeArrayToDBL_MIN(Array<double>& data)
{
   for (int i = 0; i < data.getSize(); i++) {
      data[i] = s_dbl_min;
   }
}

void IEEE::initializeArrayToSignalingNaN(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_signaling_nan_float;
   }
}

void IEEE::initializeArrayToSignalingNaN(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_signaling_nan_double;
   }
}

void IEEE::initializeArrayToINT_MAX(int *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_int_max;
   }
}

void IEEE::initializeArrayToINT_MIN(int *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_int_min;
   }
}

void IEEE::initializeArrayToFLT_MAX(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_flt_max;
   }
}

void IEEE::initializeArrayToFLT_MIN(float *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_flt_min;
   }
}

void IEEE::initializeArrayToDBL_MAX(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_dbl_max;
   }
}

void IEEE::initializeArrayToDBL_MIN(double *data, const int n)
{
   for (int i = 0; i < n; i++) {
      data[i] = s_dbl_min;
   }
}

/*
*************************************************************************
*									*
* Return whether or not the value is a NaN.     	                *
*									*
*************************************************************************
*/

bool IEEE::isNaN(const float &f) 
{
   int i = isnan(f);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}

bool IEEE::isNaN(const double &d) 
{
   int i = isnan(d);
   if (i != 0) {
     return(true);
   } else {
     return(false);
   }
}


}
}
