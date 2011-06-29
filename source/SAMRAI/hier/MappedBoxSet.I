/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std 
 *
 ************************************************************************/
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
MappedBoxSet::MappedBoxSet():
   d_set()
{
}

SAMRAI_INLINE_KEYWORD
MappedBoxSet::MappedBoxSet(
   const MappedBoxSet& other):
   d_set(other.d_set)
{
}

SAMRAI_INLINE_KEYWORD
MappedBoxSet::~MappedBoxSet()
{
}

SAMRAI_INLINE_KEYWORD
MappedBoxSet& MappedBoxSet::operator = (
   const MappedBoxSet& rhs)
{
   d_set = rhs.d_set;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool MappedBoxSet::operator == (
   const MappedBoxSet& rhs) const
{
   return d_set == rhs.d_set;
}

SAMRAI_INLINE_KEYWORD
bool MappedBoxSet::operator != (
   const MappedBoxSet& rhs) const
{
   return d_set != rhs.d_set;
}

SAMRAI_INLINE_KEYWORD
void MappedBoxSet::swap(
   MappedBoxSet& other)
{
   d_set.swap(other.d_set);
}

SAMRAI_INLINE_KEYWORD
void MappedBoxSet::separatePeriodicImages(
   std::vector<MappedBox>& real_mapped_box_vector,
   std::vector<MappedBox>& periodic_image_mapped_box_vector) const
{
   if (!empty()) {

      const MappedBox& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                     first_element.getDim())->
                                  getZeroShiftNumber());

      real_mapped_box_vector.reserve(real_mapped_box_vector.size() + size());
      for (const_iterator ni = begin(); ni != end(); ++ni) {
         const MappedBox& mapped_box = *ni;
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

}
}
