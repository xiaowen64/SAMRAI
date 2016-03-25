//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/boxgeometry/SideOverlap.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	hier::Box intersection information for side centered objects
//

#ifndef included_pdat_SideOverlap
#define included_pdat_SideOverlap

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxList.h"
#include "BoxOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class SideOverlap<DIM> represents the intersection between two side 
 * centered geometry boxes.  It is a subclass of hier::BoxOverlap<DIM> and records
 * the portions of index space that needs to be copied between two objects
 * with side centered geometry.
 *
 * @see hier::BoxOverlap
 * @see pdat::SideOverlap
 */

template<int DIM> class SideOverlap : public hier::BoxOverlap<DIM>
{
public:
   /**
    * The constructor takes the list of boxes and the source offset between
    * the source and destination index spaces.  This information is used later
    * in the generation of communication schedules.
    */
   SideOverlap(
      const hier::BoxList<DIM> boxes[DIM], const hier::IntVector<DIM>& src_offset);

   /**
    * The virtual destructor does nothing interesting except deallocate
    * box data.
    */
   virtual ~SideOverlap<DIM>();

   /**
    * Return whether there is an empty intersection between the two
    * side centered boxes.  This method over-rides the virtual function
    * in the hier::BoxOverlap<DIM> base class.
    */
   virtual bool isOverlapEmpty() const;

   /**
    * Return the list of boxes (in side centered index space) that
    * constitute the intersection.  The boxes are given in the
    * destination coordinate space and must be shifted by
    * -(getSourceOffset()) to lie in the source index space.  The axis
    * argument represents which axis is desired: X=0, Y=1, and Z=2.
    * This method over-rides the virtual function in the
    * hier::BoxOverlap<DIM> base class.
    */
   virtual const hier::BoxList<DIM>& getDestinationBoxList(const int axis) const;

   /**
    * Return the offset between the destination and source index spaces.
    * The destination index space is the source index space shifted
    * by this amount.
    */
   virtual const hier::IntVector<DIM>& getSourceOffset() const;

private:
   bool d_is_overlap_empty;
   hier::IntVector<DIM> d_offset;
   hier::BoxList<DIM> d_dst_boxes[DIM];

};


}
}
#ifndef DEBUG_NO_INLINE
#include "SideOverlap.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideOverlap.C"
#endif
