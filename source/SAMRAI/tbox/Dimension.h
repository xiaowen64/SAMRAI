/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Dimension class for abstracting dimension 
 *
 ************************************************************************/

#ifndef included_tbox_Dimension
#define included_tbox_Dimension

#include "SAMRAI/SAMRAI_config.h"

#include <iostream>
#include <limits>

/*
 * These forward declarations are obviously bad and creates a
 * dependency in the packages that violates the general nesting of
 * SAMRAI packages.  Was needed as IntVector needed the default ctor
 * in order to make the library work reasonable close to previous
 * versions and some performance issues would result if we did not
 * allow this.
 *
 * It would be very good to come up with something better than this.
 */
namespace SAMRAI {
namespace hier {

class IntVector;
}

namespace pdat {
template<class TYPE>
class ArrayData;
}

}

namespace SAMRAI {
namespace tbox {

class DatabaseBox;

/**
 * Class Dimension is used to represent the dimension of a SAMRAI
 * object.
 *
 * The maximum dimension is set at compile time using a flag to the
 * configure script.  This is used to allocate arrays in some lower
 * level classes such as IntVector.  If dynamic memory allocation is
 * used the performance impact is significant; a maximum dimension
 * allows for stack based memory allocation in performance critical
 * classes at the expense of wasting storage for objects with
 * dimension less than the maximum dimension.
 *
 * A class is used rather than a simple short or integer to provide
 * enhanced type safety.
 *
 */

class Dimension
{
public:

   /**
    * @brief Default constructor, creating an invalid dimension object
    * that must be set using setValue() before using.
    */
   Dimension();

   /**
    * Constructor for Dimension, object is built using the specified dimension
    *
    * Note that the constructor is "explicit" thus making automatic
    * type conversions from integers impossible.  This is intentionally to avoid
    * unintended conversions.
    *
    * When dimensional assertion checking is active an assert is
    * thrown when dim < 1 or dim > getMaxDimension() value specified when
    * the library is configured (defaults to 3).  dim also cannot be
    * the getInvalidDimension() (the largest unsigned short value).
    *
    */
   explicit Dimension(
      const unsigned short& dim);

   /**
    * Construct a dimension equal to the argument.
    */
   Dimension(
      const Dimension& dimension);

   /**
    * Returns true if Dimension is valid.
    *
    * A valid Dimension != 0; != getInvalidDimension(), 
    * and <= getMaxDimension().
    *
    */
   bool
   isValid() const;

   /**
    * Returns true if Dimension is initialized (not set 
    * to getInvalidDimension()).
    *
    * Uninitialized dimensions may be set using setValue().
    */
   bool
   isInitialized() const;

   /*!
    * @brief Set the value of the dimension.
    *
    * This method can be called at most once per object.
    * It may be called only if the current value is invalid.
    */
   void setValue( const unsigned short &dim );

   /**
    * Equality operator.
    */
   bool
   operator == (
      const Dimension& rhs) const;

   /**
    * Inequality operator.
    */
   bool
   operator != (
      const Dimension& rhs) const;

   /**
    * Greater than operator.
    */
   bool
   operator > (
      const Dimension& rhs) const;

   /**
    * Greater than or equal operator.
    */
   bool
   operator >= (
      const Dimension& rhs) const;

   /**
    * Less than operator.
    */
   bool
   operator < (
      const Dimension& rhs) const;

   /**
    * Less than or equal operator.
    */
   bool
   operator <= (
      const Dimension& rhs) const;

   /**
    * Returns the dimension of the Dimension as an unsigned short.
    *
    * The method is provided to allow sizing of arrays based on the
    * dimension and for iteration.  In general this should not be
    * used for comparisons, the Dimension comparison operations are
    * better suited for that purpose.
    */
   unsigned short getValue() const {
      return d_dim;
   }

   /**
    * Returns the maximum dimension for the currently compiled library
    * as an unsigned short.
    *
    * When the SAMRAI library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    *  double array[SAMRAI::tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
    *
    * The value must be >= 1 and < numeric_limits<unsigned short>::max()
    */
   static const unsigned short MAXIMUM_DIMENSION_VALUE = 
      SAMRAI_MAXIMUM_DIMENSION;
   static unsigned short
   getMaxDimValue();

   /**
    * Returns the maximum dimension for the currently compiled library
    * as a Dimension object.
    *
    * When the SAMRAI library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    */
   static const Dimension&
   getMaxDimension();

   /**
    * An invalid dimension value as a Dimension object.
    */
   static const Dimension& 
   getInvalidDimension();

   /**
    * An invalid dimension value as an unsigned short.
    *
    * Currently this value is numeric_limits<unsigned short>::max() but
    * use this symbol as it is more readable.
    *
    */
   static unsigned short 
   getInvalidDimValue();

   /**
    * Output operator for debugging and error messages.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Dimension& rhs);

private:


   /**
    * Assignment operator is private to prevent dimensions
    * from being assigned.  This was done to improve type
    * safety.
    */
   Dimension&
   operator = (
      const Dimension& rhs);

   unsigned short d_dim;

   static Dimension s_maximum_dimension;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/tbox/Dimension.I"
#endif

#endif
