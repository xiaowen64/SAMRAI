/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real MappedBoxes in a MappedBoxSet. 
 *
 ************************************************************************/
#ifndef included_hier_RealMappedBoxConstIterator_C
#define included_hier_RealMappedBoxConstIterator_C

#include "SAMRAI/hier/RealMappedBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

RealMappedBoxConstIterator::RealMappedBoxConstIterator(
   const MappedBoxSet& mapped_boxes):
   d_mapped_boxes(&mapped_boxes)
{
   d_ni = d_mapped_boxes->begin();
   while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage()) {
      ++d_ni;
   }
}

RealMappedBoxConstIterator::~RealMappedBoxConstIterator()
{
   d_mapped_boxes = NULL;
}

bool RealMappedBoxConstIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_ni != d_mapped_boxes->end() &&
          !d_ni->isPeriodicImage();
}

RealMappedBoxConstIterator& RealMappedBoxConstIterator::operator = (
   const RealMappedBoxConstIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_ni = r.d_ni;
   return *this;
}

const MappedBox& RealMappedBoxConstIterator::operator * () const
{
   return *d_ni;
}

const MappedBox *RealMappedBoxConstIterator::operator -> () const
{
   return &(*d_ni);
}

bool RealMappedBoxConstIterator::operator == (
   const RealMappedBoxConstIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes && d_ni == r.d_ni;
}

bool RealMappedBoxConstIterator::operator != (
   const RealMappedBoxConstIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes || d_ni != r.d_ni;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

RealMappedBoxConstIterator &RealMappedBoxConstIterator::operator ++ ()
{
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage());
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

RealMappedBoxConstIterator RealMappedBoxConstIterator::operator ++ (
   int)
{
   RealMappedBoxConstIterator saved = *this;
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage());
   return saved;
}

}
}
#endif
