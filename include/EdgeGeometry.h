//
// File:	EdgeGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for edge centered objects
//

#ifndef included_pdat_EdgeGeometry
#define included_pdat_EdgeGeometry

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
 * Class EdgeGeometry<DIM> manages the mapping between the AMR index space
 * and the edge-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between edge-
 * centered box geometries.  That is, edge geometry objects calculate the 
 * edge-centered data residing in the intersection of two boxes defining 
 * regions of index space on an AMR patch hierarchy.  For example, given 
 * a three-dimensional box [l0:u0,l1:u1,l2:u2], the indices for a 
 * three-dimensional edge data object run as follows:
 * 


 * - \b X edges [l0:u0,l1:u1+1,l2:u2+1]
 * - \b Y edges [l0:u0+1,l1:u1,l2:u2+1]
 * - \b Z edges [l0:u0+1,l1:u1+1,l2:u2]
 * 


 * Recall that edge data is defined so that the edges associated with a
 * given coordinate direction are those whose tangent vector lies in that
 * direction.
 *
 * Note that the intersection between two edge-centered boxes can be 
 * complicated, since edge geometries contain indices on the edges of 
 * boxes.  Thus, there may be overlap between two boxes, even though 
 * the boxes do not intersect in the AMR index space.
 *
 * @see hier::BoxGeometry
 * @see pdat::EdgeOverlap
 */

template<int DIM> class EdgeGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the edge geometry object given the box and ghost cell width.
    */
   EdgeGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~EdgeGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source edge box
    * geometry object and the destination box geometry.  Refer to the
    * box geometry class for a detailed description of calculateOverlap().
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /**
    * Return the box extents for this edge centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this edge centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /**
    * Convert an AMR abstract box into a edge geometry box.  The box indices
    * are cyclically shifted such that the edge direction is first.  The edge
    * direction runs from the corresponding lower index to the upper index
    * plus one.  All other indices run as in the original box.  The axes
    * are given by X=0, Y=1, and Z=2.
    */
   static hier::Box<DIM> toEdgeBox(const hier::Box<DIM>& box, const int axis);

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have edge centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const EdgeGeometry<DIM>& dst_geometry,
      const EdgeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   EdgeGeometry(const EdgeGeometry<DIM>&);	// not implemented
   void operator=(const EdgeGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "EdgeGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "EdgeGeometry.C"
#endif
