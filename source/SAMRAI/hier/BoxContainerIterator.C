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
   d_box_container(&container),
   d_list_iter(from_start ? container.begin().d_list_iter :
               container.end().d_list_iter)
{
}

BoxContainerIterator& BoxContainerIterator::operator = (
   const BoxContainerIterator& rhs)
{
   if (this != &rhs) {
      d_box_container = rhs.d_box_container;
      d_list_iter = rhs.d_list_iter;
   }
   return *this;
}

}
}

#endif
