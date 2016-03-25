//
// File:	SideGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for side centered objects
//

#ifndef included_pdat_SideGeometry
#define included_pdat_SideGeometry

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
 * Class SideGeometry<DIM> manages the mapping between the AMR index space
 * and the side-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between side-
 * centered box geometries.  That is, side geometry objects calculate the 
 * side-centered data residing in the intersection of two boxes defining
 * regions of index space on an AMR patch hierarchy.  For example, given 
 * a three-dimensional box [l0:u0,l1:u1,l2:u2], the indices for a 
 * three-dimensional side data object run as follows:
 * 


 * - \b X sides [l0:u0+1,l1:u1  ,l2:u2  ]
 * - \b Y sides [l0:u0  ,l1:u1+1,l2:u2  ]
 * - \b Z sides [l0:u0  ,l1:u1  ,l2:u2+1]
 * 


 * Recall that side data is defined so that the sides associated with a
 * given coordinate direction are those whose normal vector lies in that
 * direction.  Also, side data and face data storage are similar.  However,
 * the ordering of the indices is what differentiates the side data classes
 * from the face data classes.
 *
 * Note that the intersection between two side-centered boxes can be
 * complicated since side geometries contain indices on the sides of
 * the boxes.  Thus, there may be overlap between two boxes, even though
 * the boxes do not intersect in the AMR index space.
 *
 * Note that it is possible to manage side data objects where the data
 * is allocated for sides associated with a single coordinate direction only.
 * All side geometry operations involving more than one side geometry object
 * are only defined for the case where all objects have the same data
 * allocation.  If this is not the case, an assertion will result.
 *
 * @see hier::BoxGeometry
 * @see pdat::SideOverlap
 */

template<int DIM> class SideGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct a side geometry object given the box, ghost cell width,
    * and information about which coordinate directions are allocated.
    */
   SideGeometry(const hier::Box<DIM>& box,
                      const hier::IntVector<DIM>& ghosts,
                      const hier::IntVector<DIM>& directions);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~SideGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source side box
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
    * Return the box extents for this side centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this side centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /**
    * Return constant reference to vector describing which coordinate
    * directions managed by this side geometry object.  A vector entry 
    * of zero indicates that this object will not perform operations 
    * involving the corresponding coordinate direction.  A non-zero value 
    * indicates otherwise.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

   /**
    * Convert an AMR abstract box into a side geometry box.  The box indices
    * are cyclically shifted such that the side direction is first.  The side
    * direction runs from the corresponding lower index to the upper index
    * plus one.  All other indices run as in the original box.  The axes
    * are given by X=0, Y=1, and Z=2.
    */
   static hier::Box<DIM> toSideBox(const hier::Box<DIM>& box, const int axis);

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have side centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const SideGeometry<DIM>& dst_geometry,
      const SideGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   SideGeometry(const SideGeometry<DIM>&);	// not implemented
   void operator=(const SideGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;
   hier::IntVector<DIM> d_directions;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "SideGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SideGeometry.C"
#endif
