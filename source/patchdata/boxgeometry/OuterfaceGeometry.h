//
// File:	OuterfaceGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for outerface centered objects
//

#ifndef included_pdat_OuterfaceGeometry
#define included_pdat_OuterfaceGeometry

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

template<int DIM> class FaceGeometry;

/**
 * Class OuterfaceGeometry<DIM> manages the mapping between the AMR index 
 * space and the outerface-centered geometry index space.  It is a 
 * subclass of hier::BoxGeometry<DIM> and it computes intersections between two 
 * outerface-centered box geometries or between an outerface geometry and 
 * a face box geometry.  Outerface data differs from face data in that, given 
 * a box, an outerface data object represents face-centered data living only 
 * on the boundary of the box.  However, for the faces over which outerface
 * data is defined, the outerface geometry index space is the same as the face
 * geometry index space.  For example, given a three-dimensional box 
 * [l0:u0,l1:u1,l2:u2], the indices for a three-dimensional outerface data  
 * object run as follows:
 * 


 * - \b X faces (lower and upper) [l1:u1,l2:u2]
 * - \b Y faces (lower and upper) [l2:u2,l0:u0]
 * - \b Z faces (lower and upper) [l0:u0,l1:u1]
 * 


 * Recall that face data is defined so that the faces associated with a
 * given coordinate direction are those whose normal vector lies in that
 * direction.  Also, face data indices are permuted so that the leading
 * dimension of each array corresponds to the direction of the faces. 
 * Outerface data matches these conventions.
 *
 * @see hier::BoxGeometry
 * @see pdat::FaceGeometry
 * @see pdat::FaceOverlap
 */

template<int DIM> class OuterfaceGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the outerface geometry object given the box and ghost
    * cell width.
    */
   OuterfaceGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~OuterfaceGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source outerface
    * geometry object (or a face geometry object) and the destination
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
    * Return the box extents for this outerface box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this outerface box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where the source
    * has outerface geometry and the destination face geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const FaceGeometry<DIM>& dst_geometry,
      const OuterfaceGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   OuterfaceGeometry(const OuterfaceGeometry<DIM>&); // not implemented
   void operator=(const OuterfaceGeometry<DIM>&);	    // not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuterfaceGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuterfaceGeometry.C"
#endif
