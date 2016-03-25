//
// File:	OuternodeGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 483 $
// Modified:	$Date: 2005-07-22 16:00:34 -0700 (Fri, 22 Jul 2005) $
// Description:	hier::Box geometry information for outernode centered objects
//

#ifndef included_pdat_OuternodeGeometry
#define included_pdat_OuternodeGeometry

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

template<int DIM> class NodeGeometry;

/*!
  @brief Manages the mapping between the AMR index 
  space and the outernode-centered geometry index space.

  This is a 
  subclass of hier::BoxGeometry<DIM> and it computes intersections between two
  outernode-centered box geometries or between an outernode geometry and 
  a side box geometry.  Outernode data differs from node data in that, given
  a box, an outernode data object represents node-centered data living only 
  on the boundary of the box.  However, for the sides over which outernode
  data is defined, the outernode geometry index space is the same as the side
  geometry index space.  For example, given a three-dimensional box
  [l0:u0,l1:u1,l2:u2], the indices for a three-dimensional outernode data  
  object run as follows:

  - @b X sides (lower and upper) [l1+1:u1-1,l2+1:u2-1]
  - @b Y sides (lower and upper) [l0:u0,l2+1:u2-1]
  - @b Z sides (lower and upper) [l0:u0,l1:u1]

  The reduction in some of the index ranges prevent nodes on edges and
  corners of the box from being repeated by different sides of the box.
  Recall that side data is defined so that the sides associated with a
  given coordinate direction are those whose normal vector lies in that
  direction.  Outernode data matches this conventions.

  @see hier::BoxGeometry
  @see pdat::NodeGeometry
  @see pdat::NodeOverlap
*/

template<int DIM> class OuternodeGeometry : public hier::BoxGeometry<DIM>
{
public:
   /*!
     @brief
     Construct the outernode geometry object given the box and ghost
     cell width.
   */
   OuternodeGeometry(const hier::Box<DIM>& box,
			   const hier::IntVector<DIM>& ghosts);

   /*!
     @brief
     The virtual destructor does nothing interesting.
   */
   virtual ~OuternodeGeometry<DIM>();

   /*!
     @brief
     Return a hier::BoxOverlap<DIM> object to describe the
     overlap in index space between the source outernode
     geometry object (or a node geometry object) and the destination
     object (this).

     Refer to the box geometry class for a detailed
     description of calculateOverlap().
   */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
     @brief
     Return the box extents for this outernode box geometry object.
   */
   const hier::Box<DIM>& getBox() const;

   /*!
     @brief
     Return the ghost cell width for this outernode box geometry object.
   */
   const hier::IntVector<DIM>& getGhosts() const;

private:

   /*!
     @brief
     Compute the overlap
     between the source and destination objects, where the source
     has outernode geometry and the destination node geometry.
   */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const NodeGeometry<DIM>& dst_geometry,
      const OuternodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   /*!
     @brief
     Compute the overlap
     between the source and destination objects, where the source
     has node geometry and the destination outernode geometry.
   */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const OuternodeGeometry<DIM>& dst_geometry,
      const NodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);


   /*!
     @brief
     Compute the overlap
     between the source and destination objects, where the source
     has outernode geometry and the destination outernode geometry.
   */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const OuternodeGeometry<DIM>& dst_geometry,
      const OuternodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   /*! Not implemented */
   OuternodeGeometry(const OuternodeGeometry<DIM>&);
   /*! Not implemented */
   void operator=(const OuternodeGeometry<DIM>&);

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuternodeGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuternodeGeometry.C"
#endif
