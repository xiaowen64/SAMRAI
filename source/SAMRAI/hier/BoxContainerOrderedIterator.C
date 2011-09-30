/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerOrderedIterator_C
#define included_hier_BoxContainerOrderedIterator_C

#include "SAMRAI/hier/BoxContainerOrderedIterator.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerOrderedIterator.I"
#endif

namespace SAMRAI {
namespace hier {

BoxContainerOrderedIterator::BoxContainerOrderedIterator(
   BoxContainer& container,
   bool from_start):
   d_box_container(&container),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end())
{
}

BoxContainerOrderedIterator::BoxContainerOrderedIterator(
   const BoxContainer& container,
   bool from_start):
   d_box_container(&container),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end())
{
}

BoxContainerOrderedIterator& BoxContainerOrderedIterator::operator = (
   const BoxContainerOrderedIterator& rhs)
{
   if (this != &rhs) {
      d_box_container = rhs.d_box_container;
      d_set_iter = rhs.d_set_iter;
   }
   return *this;
}

}
}

#endif
