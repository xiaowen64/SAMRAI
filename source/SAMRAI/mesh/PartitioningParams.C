/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Utilitiy for breaking boxes during partitioning.
 *
 ************************************************************************/

#ifndef included_mesh_PartitioningParams_C
#define included_mesh_PartitioningParams_C

#include "SAMRAI/mesh/PartitioningParams.h"

namespace SAMRAI {
namespace mesh {


PartitioningParams::PartitioningParams(
   const hier::BaseGridGeometry &grid_geometry,
   const hier::IntVector &ratio_to_level_zero,
   const hier::IntVector &min_size,
   const hier::IntVector &max_size,
   const hier::IntVector &bad_interval,
   const hier::IntVector &cut_factor ) :
   d_min_size(min_size),
   d_max_size(max_size),
   d_bad_interval(bad_interval),
   d_cut_factor(cut_factor),
   d_load_comparison_tol(0)
{
   for ( int bid(0); bid<grid_geometry.getNumberBlocks(); ++bid ) {
      grid_geometry.computePhysicalDomain(
         d_block_domain_boxes[hier::BlockId(bid)], ratio_to_level_zero, hier::BlockId(bid));
   }
}


}
}


#endif
