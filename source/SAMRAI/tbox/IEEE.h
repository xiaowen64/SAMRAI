/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   IEEE routines to set up handlers and get signaling NaNs 
 *
 ************************************************************************/

#ifndef included_tbox_IEEE
#define included_tbox_IEEE

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Complex.h"

namespace SAMRAI {
namespace tbox {

/*!
 * Class IEEE is a utility providing rotuines for managing IEEE trap
 * handlers and data set to signaling NaNs.  Signaling NaNs force
 * a trap if they are used in a numerical operation, so they are a
 * useful way to track uninitialized floating point data.  Signaling
 * NaN's may only be used for double and float data (and the real
 * ans imaginary parts of dcomplex data) and so operations are
 * provided here for those types only.
 *
 * IMPORTANT: To properly trap operations based on signaling NaN values,
 *            the routine IEEE::setupFloatingPointExceptionHandlers()
 *            must be called.  This is normally done in the
 *            SAMRAIManager::startup() routine.
 *
 * Note that all operations provided by this class (except for setting
 * up exception handling) are implemented in @see tbox::MathUtilities.
 * Operations are provided by this class since it is not templated on
 * data type and so calling the operations provided here may be easier
 * in some cases, such as in codes built based on earlier versions
 * of SAMRAI.  See the tbox::MathUtilities header file for details
 * about the routines.
 *
 * @see tbox::MathUtilities
 */

struct IEEE {
   /*!
    * Set up IEEE exception handlers so that normal IEEE exceptions will
    * cause a program abort.  This is useful for tracking down errors.
    * Note, however, that this may cause problems if your code relies on
    * IEEE underflow or overflow traps.
    */
   static void
   setupFloatingPointExceptionHandlers();

   /*!
    * Get the IEEE float signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static float
   getSignalingFloatNaN();

   /*!
    * Get the IEEE double signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static double
   getSignalingNaN();

   /*!
    * Get the dcomplex value with real and imaginary parts set to the
    * IEEE double signaling NaN on architectures that support it.
    * Using this value in a numerical expression will cause a program abort.
    */
   static dcomplex
   getSignalingComplexNaN();

   /*!
    * Set supplied float value to the signaling NaN.
    */
   static void
   setNaN(
      float& f);

   /*!
    * Set supplied double value to the signaling NaN.
    */
   static void
   setNaN(
      double& d);

   /*!
    * Set real and imaginary parts of supplied dcomplex value to the
    * double signaling NaN.
    */
   static void
   setNaN(
      dcomplex& dc);

   /*!
    * Initialize an array of floats to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      Array<float>& array);

   /*!
    * Initialize an array of doubles to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      Array<double>& array);

   /*!
    * Initialize an array of dcomplex to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      Array<dcomplex>& array);

   /*!
    * Initialize an array of floats to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      float* array,
      int n = 1);

   /*!
    * Initialize an array of doubles to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      double* array,
      int n = 1);

   /*!
    * Initialize an array of dcomplex to signaling NaNs.  Before using this
    * array in any operation, the NaN value should be reset.  Otherwise,
    * an unrecoverable exception will result (as long as floating point
    * exception handling is supported by the compiler).
    */
   static void
   initializeArrayToSignalingNaN(
      dcomplex* array,
      int n = 1);

   /*!
    * Return true if the supplied float value is NaN; else, false.
    */
   static bool
   isNaN(
      const float& f);

   /*!
    * Return true if the supplied double value is NaN; else, false.
    */
   static bool
   isNaN(
      const double& d);

   /*!
    * Return true if if either real and imaginary part of the supplied
    * dcomplex value is NaN; else, false.
    */
   static bool
   isNaN(
      const dcomplex& dc);
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/IEEE.I"
#endif
#endif
