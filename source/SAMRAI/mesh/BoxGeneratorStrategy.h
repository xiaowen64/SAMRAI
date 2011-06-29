/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box generation routines. 
 *
 ************************************************************************/

#ifndef included_mesh_BoxGeneratorStrategy
#define included_mesh_BoxGeneratorStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/hier/MappedBoxLevel.h"
#include "SAMRAI/hier/PatchLevel.h"

namespace SAMRAI {
namespace mesh {

/**
 * Class BoxGeneratorStrategy is an abstract base class that defines
 * a Strategy pattern interface for operations to build boxes that cover a
 * collection of tagged cells on a single AMR patch hierarchy level.
 *
 * @see hier::PatchLevel
 */

class BoxGeneratorStrategy:public tbox::DescribedClass
{
public:
   /**
    * Default constructor.
    */
   BoxGeneratorStrategy();

   /**
    * Virtual destructor.
    */
   virtual ~BoxGeneratorStrategy();

   /*!
    * @brief Cluster tags into non-overlapping boxes.
    *
    * Create a set of boxes that covers all integer tags on the patch
    * level that match the specified tag value.  Each box should be at
    * least as large as the given minimum size and the tolerances will
    * be met.
    *
    * The efficiency tolerance is a threshold value for the fraction
    * of tagged cells in each box.  If this percentage is below the
    * tolerance, the box will continue to be split into smaller boxes.
    *
    * The combine tolerance is a threshold value for the sum of the
    * volumes of two boxes into which a box may be potentially split.
    * If ratio of that sum and the volume of the original box, the box
    * will not be split.
    */
   virtual void
   findBoxesContainingTags(
      hier::MappedBoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const tbox::Pointer<hier::PatchLevel> tag_level,
      const int tag_data_index,
      const int tag_val,
      const hier::Box& bound_box,
      const hier::IntVector& min_box,
      const double efficiency_tol,
      const double combine_tol,
      const hier::IntVector& max_gcw,
      const hier::BlockId& block_id) const = 0;

private:
   // The following are not implemented:
   BoxGeneratorStrategy(
      const BoxGeneratorStrategy&);
   void
   operator = (
      const BoxGeneratorStrategy&);

};

}
}
#endif
