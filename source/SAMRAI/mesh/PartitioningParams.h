/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2014 Lawrence Livermore National Security, LLC
 * Description:   Parameters in load balancing.
 *
 ************************************************************************/

#ifndef included_mesh_PartitioningParams
#define included_mesh_PartitioningParams

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BaseGridGeometry.h"

#include <map>

namespace SAMRAI {
namespace mesh {



/*!
 * @brief Light weight class holding parameters generally used
 * in partitioning.
 */

class PartitioningParams
{
public:

   PartitioningParams(
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
      d_load_comparison_tol(1e-8)
      {
         for ( int bid(0); bid<grid_geometry.getNumberBlocks(); ++bid ) {
            grid_geometry.computePhysicalDomain(
               d_block_domain_boxes[hier::BlockId(bid)], ratio_to_level_zero, hier::BlockId(bid));
         }
      }

   double getMinLoad() const {
      return static_cast<double>(d_min_size.getProduct());
   }

   const hier::IntVector &getMinBoxSize() const {
      return d_min_size;
   }

   const hier::IntVector &getMaxBoxSize() const {
      return d_max_size;
   }

   const hier::BoxContainer &getDomainBoxes( const hier::BlockId &bid ) const {
      return d_block_domain_boxes.find(bid)->second;
   }

   const hier::IntVector &getBadInterval() const {
      return d_bad_interval;
   }

   const hier::IntVector &getCutFactor() const {
      return d_cut_factor;
   }

   const tbox::Dimension &getDim() const {
      return d_min_size.getDim();
   }

   const double &getLoadComparisonTol() const {
      return d_load_comparison_tol;
   }

private:

   std::map<hier::BlockId,hier::BoxContainer> d_block_domain_boxes;
   hier::IntVector d_min_size;
   hier::IntVector d_max_size;
   hier::IntVector d_bad_interval;
   hier::IntVector d_cut_factor;
   double d_load_comparison_tol;

};

}
}

#endif
