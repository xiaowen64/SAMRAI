//
// File:	CellGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for cell centered objects
//

#ifndef included_pdat_CellGeometry
#define included_pdat_CellGeometry

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_BoxGeometry
#include "BoxGeometry.h"
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
 * Class CellGeometry<DIM> manages the mapping between the AMR index space
 * and the cell-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between cell-
 * centered box geometries.  That is, cell geometry objects calculate the
 * cell-centered data residing in the intersection of two boxes defining
 * regions of index space on an AMR patch hierarchy.  Since the AMR index 
 * space is also cell-centered, cell-centered data storage maps directly 
 * onto a cell geometry box with no modifications. 
 *
 * @see hier::BoxGeometry
 * @see pdat::CellOverlap
 */

template<int DIM> class CellGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the cell geometry object given the box and ghost cell width.
    */
   CellGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~CellGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source box geometry
    * and the destination box geometry.  Refer to the box geometry class
    * for a detailed description of calculateOverlap().
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /**
    * Return the box extents for this cell centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this cell centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /**
    * Convert an AMR abstract box into a cell geometry box.  The lower
    * index is the same, but the upper index is one greater in each
    * dimension.
    */
   static hier::Box<DIM> toCellBox(const hier::Box<DIM>& box);

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have cell centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const CellGeometry<DIM>& dst_geometry,
      const CellGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   CellGeometry(const CellGeometry<DIM>&);	// not implemented
   void operator=(const CellGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "CellGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellGeometry.C"
#endif
