//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/boxgeometry/CellOverlap.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	hier::Box intersection information for cell centered objects
//

#ifndef included_pdat_CellOverlap
#define included_pdat_CellOverlap

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxList.h"
#include "BoxOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/**
 * Class CellOverlap<DIM> represents the intersection between two cell
 * centered geometry boxes.  It is a subclass of hier::BoxOverlap<DIM>
 * and records the portions of index space that needs to be copied
 * between two objects with cell centered geometry.  Note that
 * CellOverlap does NOT compute the overlap of the arguments. It
 * stores the arguments as given and assumes that they already
 * represent an overlap previously computed.
 *
 * @see hier::BoxOverlap
 * @see pdat::CellOverlap
 */

template<int DIM> class CellOverlap : public hier::BoxOverlap<DIM>
{
public:
   /**
    * The constructor takes the list of boxes and the source offset between
    * the source and destination index spaces.  This information is used later
    * in the generation of communication schedules.
    */
   CellOverlap(const hier::BoxList<DIM>& boxes,
                     const hier::IntVector<DIM>& src_offset);

   /**
    * The virtual destructor does nothing interesting except deallocate
    * box data.
    */
   virtual ~CellOverlap<DIM>();

   /**
    * Return whether there is an empty intersection between the two
    * cell centered boxes.  This method over-rides the virtual function
    * in the hier::BoxOverlap<DIM> base class.
    */
   virtual bool isOverlapEmpty() const;

   /**
    * Return the list of boxes (in cell centered index space) that constitute
    * the intersection.  The boxes are given in the destination coordinate
    * space and must be shifted by -(getSourceOffset()) to lie in the source
    * index space.
    */
   virtual const hier::BoxList<DIM>& getDestinationBoxList() const;

   /**
    * Return the offset between the destination and source index spaces.
    * The destination index space is the source index space shifted
    * by this amount.
    */
   virtual const hier::IntVector<DIM>& getSourceOffset() const;

  /**
   * Output the boxes in the overlap region.
   */
   virtual void print(std::ostream& os) const;

private:
   bool d_is_overlap_empty;
   hier::IntVector<DIM> d_offset;
   hier::BoxList<DIM> d_dst_boxes;
   int d_count;

};


}
}
#ifndef DEBUG_NO_INLINE
#include "CellOverlap.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellOverlap.C"
#endif
