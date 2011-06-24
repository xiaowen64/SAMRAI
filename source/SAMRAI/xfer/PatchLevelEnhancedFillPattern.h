/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Level fill pattern for enhanced connectivity 
 *
 ************************************************************************/

#ifndef included_xfer_PatchLevelEnhancedFillPattern
#define included_xfer_PatchLevelEnhancedFillPattern

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/tbox/DescribedClass.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Class PatchLevelEnhancedFillPattern is an implementation of the
 * abstract base class PatchLevelFillPattern.
 *
 * This class is used by the MultiblockRefineSchedule to restrict filling to
 * only patches in ghost regions across an ehanced connectivity block boundary.
 * It is intended for users who wish to handle the filling of data at these
 * singularities separately from the filling of all other data.
 *
 * @see xfer::MultiblockRefineSchedule
 */

class PatchLevelEnhancedFillPattern:public PatchLevelFillPattern
{
public:
   /*!
    * @brief Default constructor
    */
   PatchLevelEnhancedFillPattern();

   /*!
    * @brief Destructor
    */
   virtual ~PatchLevelEnhancedFillPattern();

   /*!
    * @brief Compute the mapped boxes representing the region that will
    *        be filled by a MultiblockRefineSchedule.
    *
    * This is currently unimplemented until MappedBoxLevel is
    * multiblock-aware.  An error will occur if this is called.
    *
    * @param[in] dst_mapped_box_level  destination level
    * @param[in] dst_to_dst            Connector of destination to itself
    * @param[in] dst_to_src            Connector of destination to source
    * @param[in] src_to_dst            Connector of source to destination
    * @param[in] fill_ghost_width      Ghost width being filled by refine
    *                                  schedule
    * @param[out] fill_mapped_boxes    Output set of MappedBoxes to be filled
    * @param[out] dst_to_fill_edges    Output NeighborhoodSet between
    *                                  dst_mapped_box_level and
    *                                  and fill_mapped_boxes
    */
   void
   computeFillMappedBoxesAndNeighborhoodSets(
      hier::MappedBoxSet& fill_mapped_boxes,
      hier::NeighborhoodSet& dst_to_fill_edges,
      const hier::MappedBoxLevel& dst_mapped_box_level,
      const hier::Connector& dst_to_dst,
      const hier::Connector& dst_to_src,
      const hier::Connector& src_to_dst,
      const hier::IntVector& fill_ghost_width);

   /*!
    * @brief Return true because source patch owner cannot compute fill
    * boxes across block boundaries.
    *
    * @return  true.
    */
   bool
   needsToCommunicateDestinationFillBoxes() const;

   /*!
    * @brief Tell RefineSchedule not to communicate data directly from source
    * to destination level.
    *
    * With this fill pattern, the RefineSchedule does not attempt to fill
    * as much of the destination level as possible from the source level at
    * the same level of resolution.  By returning 'false', this method
    * tells the RefineSchedule to skip that step.
    *
    * @return false
    */
   bool
   doesSourceLevelCommunicateToDestination() const;

   /*!
    * @brief Return the maximum number of fill boxes.
    *
    * This will not return a valid value until
    * computeFillMappedBoxesAndNeighborhoodSets is fully implemented.  An
    * error will occur if this method is called
    */
   int
   getMaxFillBoxes() const;

   /*!
    * @brief Returns true because this fill pattern fills coarse fine ghost
    * data.
    */
   bool
   fillingCoarseFineGhosts() const;

   /*!
    * @brief Returns true because this fill pattern is specialized for
    * enhanced connectivity only.
    */
   bool
   fillingEnhancedConnectivityOnly() const;

private:
   PatchLevelEnhancedFillPattern(
      const PatchLevelEnhancedFillPattern&);           // not implemented

   void
   operator = (
      const PatchLevelEnhancedFillPattern&);           // not implemented

   /*!
    * @brief Shorthand typedef.
    */
   typedef hier::Connector::NeighborSet NeighborSet;

   /*!
    * @brief Maximum number of fill boxes across all destination patches.
    */
   int d_max_fill_boxes;
};

}
}

#ifdef SAMRAI_INLINE
#include "SAMRAI/xfer/PatchLevelEnhancedFillPattern.I"
#endif
#endif
