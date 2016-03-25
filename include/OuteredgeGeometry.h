//
// File:	OuteredgeGeometry.h
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Release:	$Name:  $
// Revision:	$Revision: 692 $
// Modified:	$Date: 2005-10-28 13:37:46 -0700 (Fri, 28 Oct 2005) $
// Description:	Box geometry information for edge centered objects
//

#ifndef included_pdat_OuteredgeGeometry
#define included_pdat_OuteredgeGeometry

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


template<int DIM> class EdgeGeometry;
   
/*!
 * Class OuteredgeGeometry<DIM> manages the mapping between the AMR index 
 * space and the outeredge-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry and it computes intersections between outeredge-
 * centered box geometries.  Outeredge data resides only on the outer edges
 * of the patch.  
 *
 * Outeredge data is stored in DIM*DIM*2 arrays, containing the data for the 
 * patch boundary sides with each of the possible outward pointing normal 
 * directions. Where an outeredge falls on more than one side (patch edges 
 * and corners), the outeredge belongs to the array associated with the 
 * higher dimensional direction. In each of these arrays, memory allocation 
 * is in column-major ordering (e.g., Fortran style) so that the leftmost 
 * index runs fastest in memory.  
 *
 * The outeredge data is related to the edge data in the following way:
 *
 *    Outeredge box(axis,face_nrml,s) =  EdgeData<DIM>.getBox(axis) 
 *
 * where "axis" corresponds to the box of the standard edge datatype,
 * "face_nrml" is the normal face direction, and "s" indicates the upper
 * or lower face.  Note that when edge_dir = face_dir, there are no outside
 * edges so the data is NULL.
 *
 * A three-dimensional outeredge data object instantiated with 
 * a box [l0:u0,l1:u1,l2:u2] allocates 12 data (i.e., 3x2 pairs) arrays 
 * dimensioned as:
 *
 * \verbatim
 *
 *    a = edge axis
 *    f = face normal dir
 *    s = lower/upper face
 *
 *        (a,f,s) 
 *    0:  (0,0,[0,1])  NULL
 *        (0,1,[0,1])  [l0:u0,l1:l1,l2+1:u2,d],  [l0:u0,u1:u1,l2+1:u2,d]
 *        (0,2,[0,1])  [l0:u0,l1:u1+1,l2:l2,d],  [l0:u0,l1:u1+1,u2:u2,d]
 *        note: trimmed in 2, not trimmed in 1
 *
 *    1:  (1,0,[0,1])  [l0:l0,l1:u1,l2+1:u2,d],  [u0:u0,l1:u1,l2+1:u2,d]
 *        (1,1,[0,1])  NULL
 *        (1,2,[0,1])  [l0:u0+1,l1:u1,l2:l2,d],  [l0:u0+1,l1:u1,u2:u2,d]
 *        note: trimmed in 2, not trimmed in 0
 *
 *    2:  (2,0,[0,1])  [l0:l0,l1+1:u1,l2:u2,d],  [u0:u0,l1+1:u1,l2:u2,d]
 *        (2,1,[0,1])  [l0:u0+1,l1:l1,l2:u2,d],  [l0:u0+1,u1:u1,l2:u2,d]
 *        (2,2,[0,1])  NULL
 *        note: trimmed in 1, not trimmed in 0
 *
 * \endverbatim
 *
 * where 0, 1, and 2 can be thought of as X, Y, Z respectively, and d is the
 * depth of the data.  One- and two-dimensional edge data arrays are managed 
 * similary.  The "a" dimension corresponds with the "axis" of standard 
 * EdgeData<DIM>.
 *
 * Note that the intersection between two outeredge-centered boxes can be 
 * complicated, since edge geometries contain indices on the edges of 
 * boxes.  Thus, there may be overlap between two boxes, even though 
 * the boxes do not intersect in the AMR index space.
 *
 * @see hier::BoxGeometry
 * @see pdat::OuteredgeOverlap
 */

template<int DIM> class OuteredgeGeometry : public hier::BoxGeometry<DIM>
{
public:
   /*!
    * Construct the edge geometry object given the box and ghost cell width.
    */
   OuteredgeGeometry(const hier::Box<DIM>& box, 
                     const hier::IntVector<DIM>& ghosts);

   /*!
    * The virtual destructor does nothing interesting.
    */
   virtual ~OuteredgeGeometry<DIM>();

   /*!
    * Compute the overlap in index space between the source edge box
    * geometry object and the destination box geometry.  Refer to the
    * box geometry class for a detailed description of calculateOverlap().
    */
   virtual tbox::Pointer<hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
    * Return the box extents for this edge centered box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * Return the ghost cell width for this edge centered box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /*!
    * Trim databoxes in particular axis, face_nrml directions to avoid 
    * duplicating data on common outer edges.
    */
   static void trimBoxes(hier::Box<DIM>& boxlo, 
                         hier::Box<DIM>& boxup,
                         const int axis,
                         const int face_nrml);

private:
   /*!
    * Compute overlap between a source outeredge geometry and a destination
    * edge geometry.
    */
   static tbox::Pointer<hier::BoxOverlap<DIM> > doOverlap(
      const pdat::EdgeGeometry<DIM>& dst_geometry,
      const OuteredgeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);


   /*!
    * Compute overlap between a source outeredge geometry and a destination
    * outeredge geometry.
    */
   static tbox::Pointer<hier::BoxOverlap<DIM> > doOverlap(
      const OuteredgeGeometry<DIM>& dst_geometry,
      const OuteredgeGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset);

   OuteredgeGeometry<DIM>(const OuteredgeGeometry<DIM>&);// not implemented
   void operator=(const OuteredgeGeometry<DIM>&);	// not implemented

   hier::Box<DIM> d_box;
   hier::IntVector<DIM> d_ghosts;

};

}
}
#ifndef DEBUG_NO_INLINE
#include "OuteredgeGeometry.I"
#endif
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeGeometry.C"
#endif

