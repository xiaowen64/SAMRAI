//
// File:	OutersideGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for outerside centered objects
//

#ifndef included_pdat_OutersideGeometry
#define included_pdat_OutersideGeometry

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

template<int DIM> class SideGeometry;

/**
 * Class OutersideGeometry<DIM> manages the mapping between the AMR index 
 * space and the outerside-centered geometry index space.  It is a 
 * subclass of hier::BoxGeometry<DIM> and it computes intersections between two
 * outerside-centered box geometries or between an outerside geometry and 
 * a side box geometry.  Outerside data differs from side data in that, given
 * a box, an outerside data object represents side-centered data living only 
 * on the boundary of the box.  However, for the sides over which outerside
 * data is defined, the outerside geometry index space is the same as the side
 * geometry index space.  For example, given a three-dimensional box
 * [l0:u0,l1:u1,l2:u2], the indices for a three-dimensional outerside data  
 * object run as follows:
 * 


 * - \b X sides (lower and upper) [l1:u1,l2:u2]
 * - \b Y sides (lower and upper) [l0:u0,l2:u2]
 * - \b Z sides (lower and upper) [l0:u0,l1:u1]
 * 


 * Recall that side data is defined so that the sides associated with a
 * given coordinate direction are those whose normal vector lies in that
 * direction.  Outerside data matches this conventions.
 *
 * @see hier::BoxGeometry
 * @see pdat::SideGeometry
 * @see pdat::SideOverlap
 */

template<int DIM> class OutersideGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the outerside geometry object given the box and ghost
    * cell width.
    */
   OutersideGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~OutersideGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source outerside
    * geometry object (or a side geometry object) and the destination
    * object (this).  Refer to the box geometry class for a detailed
    * description of calculateOverlap().
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /**
    * Return the box extents for this outerside box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this outerside box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where the source
    * has outerside geometry and the destination side geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const SideGeometry<DIM>& dst_geometry,
      const OutersideGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   OutersideGeometry(const OutersideGeometry<DIM>&); // not implemented
   void operator=(const OutersideGeometry<DIM>&);	    // not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OutersideGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OutersideGeometry.C"
#endif
