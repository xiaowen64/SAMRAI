/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerIterator_C
#define included_hier_BoxContainerIterator_C

#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxContainer.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerIterator.I"
#endif

namespace SAMRAI {
namespace hier {

BoxContainerIterator::BoxContainerIterator(
   BoxContainer& container,
   bool from_start):
   d_list_iter(from_start ? container.d_list.begin() :
               container.d_list.end()),
   d_set_iter(from_start ? container.d_set.begin() :
               container.d_set.end()),
   d_ordered(container.d_ordered)
{
}

BoxContainerIterator& BoxContainerIterator::operator = (
   const BoxContainerIterator& rhs)
{
   if (this != &rhs) {
      d_list_iter = rhs.d_list_iter;
      d_set_iter = rhs.d_set_iter;
      d_ordered = rhs.d_ordered;
   }
   return *this;
}

}
}

#endif
