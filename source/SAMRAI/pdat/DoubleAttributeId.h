/********************************************************************** 
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997 - 2011 Lawrence Livermore National Security, LLC
 * Description:   pdat
 **********************************************************************/
#ifndef included_pdat_DoubleAttributeId_h
#define included_pdat_DoubleAttributeId_h

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Attribute identifying class for double-valued attributes
 */
class DoubleAttributeId
{
public:
   explicit DoubleAttributeId(int value);
   DoubleAttributeId(const DoubleAttributeId& other);
   ~DoubleAttributeId();
   DoubleAttributeId& operator= (const DoubleAttributeId& rhs);
   bool operator==(const DoubleAttributeId& other) const; 
   bool operator!=(const DoubleAttributeId& other) const;
   int operator()() const; 

   friend class Attributes;
private:
   DoubleAttributeId();
   int d_val;
}; // end class DoubleAttributeId

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/pdat/DoubleAttributeId.I"
#endif

#endif
