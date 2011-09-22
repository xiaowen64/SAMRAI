/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An iterator over containers of boxes.
 *
 ************************************************************************/

#ifndef included_hier_BoxContainerSetConstIterator_C
#define included_hier_BoxContainerSetConstIterator_C

#include "SAMRAI/hier/BoxContainerSetConstIterator.h"
#include "SAMRAI/hier/BoxContainer.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxContainerSetConstIterator.I"
#endif

namespace SAMRAI {
namespace hier {

BoxContainerSetConstIterator::BoxContainerSetConstIterator(
   const BoxContainer& container,
   bool from_start):
   d_box_container(&container),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end())
{
}

BoxContainerSetConstIterator& BoxContainerSetConstIterator::operator = (
   const BoxContainerSetConstIterator& rhs)
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
