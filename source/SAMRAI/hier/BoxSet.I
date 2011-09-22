/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {
#if 0
SAMRAI_INLINE_KEYWORD
BoxSet::BoxSet():
   d_set()
{
}

SAMRAI_INLINE_KEYWORD
BoxSet::BoxSet(
   const BoxSet& other):
   d_set(other.d_set)
{
}

SAMRAI_INLINE_KEYWORD
BoxSet::~BoxSet()
{
}

SAMRAI_INLINE_KEYWORD
BoxSet& BoxSet::operator = (
   const BoxSet& rhs)
{
   d_set = rhs.d_set;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxSet::operator == (
   const BoxSet& rhs) const
{
   if (size() != rhs.size()) {
      return false;
   }
   return std::equal(d_set.begin(), d_set.end(),
      rhs.d_set.begin(), Box::id_equal());
}

SAMRAI_INLINE_KEYWORD
bool BoxSet::operator != (
   const BoxSet& rhs) const
{
   return !this->operator == (rhs);
}

SAMRAI_INLINE_KEYWORD
void BoxSet::swap(
   BoxSet& other)
{
   d_set.swap(other.d_set);
}

SAMRAI_INLINE_KEYWORD
void BoxSet::separatePeriodicImages(
   std::vector<Box>& real_mapped_box_vector,
   std::vector<Box>& periodic_image_mapped_box_vector) const
{
   if (!empty()) {

      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      real_mapped_box_vector.reserve(real_mapped_box_vector.size() + size());
      for (const_iterator ni = begin(); ni != end(); ++ni) {
         const Box& mapped_box = *ni;
         if (mapped_box.getPeriodicId() == zero_shift_number) {
            real_mapped_box_vector.push_back(*ni);
         } else {
            periodic_image_mapped_box_vector.push_back(*ni);
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
#endif
}
}
