/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch classes
 *
 ************************************************************************/

#ifndef included_hier_PatchFactory_C
#define included_hier_PatchFactory_C

#include "SAMRAI/hier/PatchFactory.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/PatchFactory.I"
#endif
namespace SAMRAI {
namespace hier {

PatchFactory::~PatchFactory()
{
}

boost::shared_ptr<Patch> PatchFactory::allocate(
   const Box& mapped_box_level_mapped_box,
   boost::shared_ptr<PatchDescriptor> descriptor) const
{
   return boost::shared_ptr<Patch>(new Patch(mapped_box_level_mapped_box,
                                             descriptor));
}

}
}
#endif
