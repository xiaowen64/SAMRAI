/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and COPYING.LESSER.
*
* Copyright:     (c) 1997 - 2011 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#ifndef included_pdat_IntegerAttributeId_h
#define included_pdat_IntegerAttributeId_h

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Attribute identifying class for integer-valued attributes
 */
class IntegerAttributeId
{
public:
   IntegerAttributeId(
      int value);
   IntegerAttributeId(
      const IntegerAttributeId& other);
   ~IntegerAttributeId();
   IntegerAttributeId&
   operator = (
      const IntegerAttributeId& rhs);
   bool
   operator == (
      const IntegerAttributeId& other) const;
   bool
   operator != (
      const IntegerAttributeId& other) const;
   int
   operator () () const;

   friend class Attributes;
private:
   IntegerAttributeId();

   int d_val;
}; // end class IntegerAttributeId.

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/IntegerAttributeId.I"
#endif

#endif
