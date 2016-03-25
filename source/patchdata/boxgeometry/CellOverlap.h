//
// File:	CellOverlap.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box intersection information for cell centered objects
//

#ifndef included_pdat_CellOverlap
#define included_pdat_CellOverlap

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxList
#include "BoxList.h"
#endif
#ifndef included_hier_BoxOverlap
#include "BoxOverlap.h"
#endif
#ifndef included_hier_IntVector
#include "IntVector.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class CellOverlap<DIM> represents the intersection between two cell 
 * centered geometry boxes.  It is a subclass of hier::BoxOverlap<DIM> and records
 * the portions of index space that needs to be copied between two objects
 * with cell centered geometry.
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
   virtual void print(ostream& os) const;

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
