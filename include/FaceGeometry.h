//
// File:	FaceGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for face centered objects
//

#ifndef included_pdat_FaceGeometry
#define included_pdat_FaceGeometry

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
 * Class FaceGeometry<DIM> manages the mapping between the AMR index space
 * and the face-centered geometry index space.  It is a subclass of 
 * hier::BoxGeometry<DIM> and it computes intersections between face-
 * centered box geometries.  That is, face geometry objects calculate the 
 * face-centered data residing in the intersection of two boxes defining 
 * regions of index space on an AMR patch hierarchy.  For example, given 
 * a three-dimensional box [l0:u0,l1:u1,l2:u2], the indices for a 
 * three-dimensional face data object run as follows:
 * 


 * - \b X faces [l0:u0+1,l1:u1,l2:u2]
 * - \b Y faces [l1:u1+1,l2:u2,l0:u0]
 * - \b Z faces [l2:u2+1,l0:u0,l1:u1]
 * 


 * Recall that face data is defined so that the faces associated with a
 * given coordinate direction are those whose normal vector lies in that
 * direction.  Also, face data indices are permuted so that the leading
 * dimension of each array corresponds to the direction of the faces.
 * Side data classes provide the same data storage as face classes however
 * the indices are not permuted for side data.
 *
 * Note that the intersection between two face-centered boxes can be 
 * complicated since face geometries contain indices on the faces of 
 * the boxes.  Thus, there may be overlap between two boxes, even though 
 * the boxes do not intersect in the AMR index space.
 *
 * @see hier::BoxGeometry
 * @see pdat::FaceOverlap
 */

template<int DIM> class FaceGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the face geometry object given the box and ghost cell width.
    */
   FaceGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~FaceGeometry<DIM>();

   /**
    * Compute the overlap in index space between the source face box
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
    * Return the box extents for this face centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this face centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /**
    * Convert an AMR abstract box into a face geometry box.  The box indices
    * are cyclically shifted such that the face direction is first.  The face
    * direction runs from the corresponding lower index to the upper index
    * plus one.  All other indices run as in the original box.  The axes
    * are given by X=0, Y=1, and Z=2.
    */
   static hier::Box<DIM> toFaceBox(const hier::Box<DIM>& box, const int axis);

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have face centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const FaceGeometry<DIM>& dst_geometry,
      const FaceGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   FaceGeometry(const FaceGeometry<DIM>&);	// not implemented
   void operator=(const FaceGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "FaceGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "FaceGeometry.C"
#endif
