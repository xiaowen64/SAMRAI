/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI 
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

SAMRAI_INLINE_KEYWORD
Dimension::Dimension()
   : d_dim(Dimension::getInvalidDimValue()) 
{
}

SAMRAI_INLINE_KEYWORD Dimension::Dimension(
   const unsigned short& dim):d_dim(dim) {
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::getMaxDimValue()));
}

SAMRAI_INLINE_KEYWORD Dimension::Dimension(
   const Dimension& dimension):d_dim(dimension.d_dim) {
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::getMaxDimValue()));
}

SAMRAI_INLINE_KEYWORD Dimension& Dimension::operator = (
   const Dimension& rhs) {
   d_dim = rhs.d_dim;
   return *this;
}

SAMRAI_INLINE_KEYWORD bool Dimension::isValid() const {
   /*
    * Note we don't check for != getInvalidDimension() because it is
    * maximum value for a unsigned short so this condition will
    * be caught by the < getMaxDimension() clause.
    *
    * If for some reason getInvalidDimension() implementation is changed
    * may need to be fixed.
    */
   return (d_dim != 0) && (d_dim <= Dimension::getMaxDimValue());
}

SAMRAI_INLINE_KEYWORD bool Dimension::isInitialized() const {
   return d_dim != Dimension::getInvalidDimValue();
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator == (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim == rhs.d_dim;
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator != (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim != rhs.d_dim;
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator > (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim > rhs.d_dim;
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator >= (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim >= rhs.d_dim;
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator < (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim < rhs.d_dim;
}

SAMRAI_INLINE_KEYWORD bool Dimension::operator <= (
   const Dimension& rhs) const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   TBOX_DIM_ASSERT_CHECK_DIM(rhs);
   return d_dim <= rhs.d_dim;
}

#ifdef SGS
SAMRAI_INLINE_KEYWORD Dimension::operator unsigned short () const {
   TBOX_DIM_ASSERT_CHECK_DIM(*this);
   return d_dim;
}
#endif

/*
 * The following methods are so constructed in order to get 
 * passed the static initialization order issue.  
 * This guarantees that initialization will occur on first access.
 */
SAMRAI_INLINE_KEYWORD 
unsigned short 
Dimension::getMaxDimValue() {
   return SAMRAI_MAXIMUM_DIMENSION;
}

SAMRAI_INLINE_KEYWORD 
const Dimension&
Dimension::getMaxDimension() {
   static Dimension dim(SAMRAI_MAXIMUM_DIMENSION);
   return dim;
}

SAMRAI_INLINE_KEYWORD
unsigned short 
Dimension::getInvalidDimValue() {
   static unsigned short invalid = 
      std::numeric_limits<unsigned short>::max();
   return invalid;
}

SAMRAI_INLINE_KEYWORD
const Dimension& 
Dimension::getInvalidDimension()
{
   static Dimension invalidDim(Dimension::getInvalidDimValue());
   return invalidDim;
}

}
}
