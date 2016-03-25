//
// File:	MathUtilities.h
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	IEEE routines to set up handlers and get signaling NaNs
//

#ifndef included_tbox_MathUtilities
#define included_tbox_MathUtilities

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class MathUtilities is a utility for managing basic math functions
 * and the initialization of data to signaling NaNs, and for managing
 * Posix constants like INT_MAX, FLT_MAX, DBL_MAX, etc.  Signaling
 * NaNs force a trap if they are used in a numerical operation, so
 * they are a useful way to track uninitialized floating point data.
 * Setting integer values to INT_MAX is a useful way to track
 * uninitialized integer values.
 *
 * The implementation of this class depends heavily on the particular
 * computer architecture and how it implements floating point arithmetic
 * and hardware traps.  
 * 
 * This class has a similiar purpose to @see tbox::IEEE.  The primary
 * purpose for replication is due to the use of variable type names
 * in the method names in the IEEE class.  This prevents templating
 * based on type.  This class is templated on type so it can be 
 * used by templated classes.
 */

template<class TYPE> class MathUtilities
{

  public:

   /**
    * @brief Get the value 0
    *
    */
   static TYPE getZero();

   /**
    * @brief Get the value 1
    *
    */
   static TYPE getOne();


   /**
    * @brief Get the IEEE signaling NaN on architectures that
    * support it.  
    *
    * Using this value in a numerical expression will
    * cause a program abort.
    *
    * Valid for float/double only.  For non float/double will return 
    * false.
    */
   static TYPE getSignalingNaN();

   /**
    * @brief Indicates whether the supplied value is NaN.
    *
    * Valid for float/double only.  For non float/double will return 
    * false.
    *
    * @param value Value to test
    */
   static bool isNaN(const TYPE &value);

   /**
    * @brief Get max for the templated type.
    */
   static TYPE getMax();

   /**
    * @brief Get min for the templated type.
    */
   static TYPE getMin();

   /**
    * @brief Get epsilon for the templated type.
    */
   static TYPE getEpsilon();

   /**
    * @brief Get value to set for undefined data.  
    *
    * Typicaly used for initialziation during debugging.
    */
   static TYPE getUndefined();

   /**
    * @brief Compute the minimum of a and b.
    * 
    * @param a
    * @param b
    */
   static TYPE Min( TYPE a, TYPE b );

   /**
    * @brief Compute the maximum of a and b.
    * 
    * @param a
    * @param b
    */
   static TYPE Max( TYPE a, TYPE b );

   /**
    * @brief Generate a random value from low to low+width.
    * 
    * @param low   Starting value for range
    * @param width Width of the range.
    */
   static TYPE Rand(TYPE low, TYPE width);

private:
   static TYPE  s_zero;
   static TYPE  s_one;
   static TYPE  s_signaling_nan;
   static TYPE  s_max;
   static TYPE  s_min;
   static TYPE  s_epsilon;
   static TYPE  s_undefined;

};

/*
 * Template specializations.
 */
template<> bool MathUtilities<float>::isNaN(const float &value);
template<> bool MathUtilities<double>::isNaN(const double &value);

template<> bool     MathUtilities<bool>::Rand(bool low, bool width);
template<> char     MathUtilities<char>::Rand(char low, char width);
template<> int      MathUtilities<int>::Rand(int low, int width);
template<> float    MathUtilities<float>::Rand(float low, float width);
template<> double   MathUtilities<double>::Rand(double low, double width);
template<> dcomplex MathUtilities<dcomplex>::Rand(dcomplex low, dcomplex width);

template<> dcomplex MathUtilities<dcomplex>::Min( dcomplex a, dcomplex b );
template<> dcomplex MathUtilities<dcomplex>::Max( dcomplex a, dcomplex b );

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/MathUtilities.I"
#endif

#endif
