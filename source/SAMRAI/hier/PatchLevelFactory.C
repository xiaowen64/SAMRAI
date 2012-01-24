/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch level objects
 *
 ************************************************************************/

#ifndef included_hier_PatchLevelFactory_C
#define included_hier_PatchLevelFactory_C

#include "SAMRAI/hier/PatchLevelFactory.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/PatchLevelFactory.I"
#endif

namespace SAMRAI {
namespace hier {

PatchLevelFactory::~PatchLevelFactory()
{
}

boost::shared_ptr<PatchLevel> PatchLevelFactory::allocate(
   const BoxLevel& mapped_box_level,
   const boost::shared_ptr<GridGeometry> grid_geometry,
   const boost::shared_ptr<PatchDescriptor> descriptor,
   boost::shared_ptr<PatchFactory> factory) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(mapped_box_level, *grid_geometry);
   PatchLevel* pl =
      new PatchLevel(mapped_box_level,
         grid_geometry,
         descriptor,
         factory);
   return boost::shared_ptr<PatchLevel>(pl);
}

boost::shared_ptr<PatchLevel> PatchLevelFactory::allocate(
   boost::shared_ptr<tbox::Database> database,
   const boost::shared_ptr<GridGeometry> grid_geometry,
   const boost::shared_ptr<PatchDescriptor> descriptor,
   const ComponentSelector& component_selector,
   boost::shared_ptr<PatchFactory> factory,
   const bool defer_boundary_box_creation) const
{
   PatchLevel* pl =
      new PatchLevel(database,
         grid_geometry,
         descriptor,
         factory,
         component_selector,
         defer_boundary_box_creation);
   return boost::shared_ptr<PatchLevel>(pl);
}

}
}
#endif
