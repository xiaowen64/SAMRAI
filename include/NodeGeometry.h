//
// File:	NodeGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Box geometry information for node centered objects
//

#ifndef included_pdat_NodeGeometry
#define included_pdat_NodeGeometry

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
 * Class NodeGeometry<DIM> manages the mapping between the AMR index space
 * and the node-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between node-
 * centered box geometries.  That is, node geometry objects calculate 
 * the node-centered data residing in the intersection of two boxes 
 * defining regions of index space on an AMR patch hierarchy.  For example, 
 * given a three-dimensional box [l0:u0,l1:u1,l2:u2], the indices for a 
 * three-dimensional node data object run as follows:
 * \verbatim
 
      [l0:u0+1,l1:u1+1,l2:u2+1]

 * \endverbatim
 *
 * Note that the intersection between two node-centered boxes can be 
 * complicated since node geometries contain indices on the nodes of 
 * the boxes.  Thus, there may be overlap between two boxes, even though 
 * the boxes do not intersect in the AMR index space.
 *
 * @see hier::BoxGeometry
 * @see pdat::NodeOverlap
 */

template<int DIM> class NodeGeometry : public hier::BoxGeometry<DIM>
{
public:
   /**
    * Construct the node geometry object given the box and ghost cell width.
    */
   NodeGeometry(const hier::Box<DIM>& box, const hier::IntVector<DIM>& ghosts);

   /**
    * The virtual destructor does nothing interesting.
    */
   virtual ~NodeGeometry<DIM>();

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
    * Return the box extents for this node centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /**
    * Return the ghost cell width for this node centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /**
    * Convert an AMR abstract box into a node geometry box.  The lower
    * index is the same, but the upper index is one greater in each
    * dimension.
    */
   static hier::Box<DIM> toNodeBox(const hier::Box<DIM>& box);

private:
   /**
    * Function doOverlap() is the function that computes the overlap
    * between the source and destination objects, where both box geometry
    * objects are guaranteed to have node centered geometry.
    */
   static tbox::Pointer< hier::BoxOverlap<DIM> > doOverlap(
      const NodeGeometry<DIM>& dst_geometry,
      const NodeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   NodeGeometry(const NodeGeometry<DIM>&);	// not implemented
   void operator=(const NodeGeometry<DIM>&);		// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "NodeGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "NodeGeometry.C"
#endif
