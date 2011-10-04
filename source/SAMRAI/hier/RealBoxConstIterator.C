/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real Boxes in a BoxSet.
 *
 ************************************************************************/
#ifndef included_hier_RealBoxConstIterator_C
#define included_hier_RealBoxConstIterator_C

#include "SAMRAI/hier/RealBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

RealBoxConstIterator::RealBoxConstIterator(
   const BoxSet& mapped_boxes):
   d_mapped_boxes(&mapped_boxes),
   d_ni(d_mapped_boxes->orderedBegin())
{
   while (d_ni != d_mapped_boxes->orderedEnd() && d_ni->isPeriodicImage()) {
      ++d_ni;
   }
}

RealBoxConstIterator::~RealBoxConstIterator()
{
   d_mapped_boxes = NULL;
}

bool RealBoxConstIterator::isValid() const
{
   return d_mapped_boxes != NULL &&
          d_ni != d_mapped_boxes->orderedEnd() &&
          !d_ni->isPeriodicImage();
}

RealBoxConstIterator& RealBoxConstIterator::operator = (
   const RealBoxConstIterator& r)
{
   d_mapped_boxes = r.d_mapped_boxes;
   d_ni = r.d_ni;
   return *this;
}

const Box& RealBoxConstIterator::operator * () const
{
   return *d_ni;
}

const Box *RealBoxConstIterator::operator -> () const
{
   return &(*d_ni);
}

bool RealBoxConstIterator::operator == (
   const RealBoxConstIterator& r) const
{
   return d_mapped_boxes == r.d_mapped_boxes && d_ni == r.d_ni;
}

bool RealBoxConstIterator::operator != (
   const RealBoxConstIterator& r) const
{
   return d_mapped_boxes != r.d_mapped_boxes || d_ni != r.d_ni;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator& RealBoxConstIterator::operator ++ ()
{
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->orderedEnd() && d_ni->isPeriodicImage());
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator RealBoxConstIterator::operator ++ (
   int)
{
   RealBoxConstIterator saved = *this;
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->orderedEnd() && d_ni->isPeriodicImage());
   return saved;
}

}
}
#endif
