//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/boxgeometry/OuteredgeGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	$Name:  $
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	Box geometry information for edge centered objects
//

#ifndef included_pdat_OuteredgeGeometry
#define included_pdat_OuteredgeGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "EdgeOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

template<int DIM> class EdgeGeometry;
   
/*!
 * Class OuteredgeGeometry<DIM> manages the mapping between the AMR index 
 * and the outeredge geometry index space.  It is a subclass of 
 * hier::BoxGeometry<DIM> and it computes intersections between outeredge
 * box geometries and edge or outeredge box geometries for communication 
 * operations.
 *
 * See header file for OuteredgeData<DIM> class for a more detailed 
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::EdgeGeometry
 * @see pdat::EdgeOverlap
 */

template<int DIM> class OuteredgeGeometry : public hier::BoxGeometry<DIM>
{
public:

   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef EdgeOverlap<DIM> Overlap;


   /*!
    * Convert a given box in the standard cell-centered AMR index space to an
    * outeredge geometry box for the specified axis, face normal, and 
    * lower/upper side.   See OuteredgeData header file for a detailed 
    * description of an outeredge box.
    */
   static hier::Box<DIM> toOuteredgeBox(const hier::Box<DIM>& box,
                                        int axis,
                                        int face_normal,
                                        int side);

   /*!
    * @brief Construct an outeredge geometry object given an AMR index
    * space box and ghost cell width.
    */
   OuteredgeGeometry(const hier::Box<DIM>& box, 
                     const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OuteredgeGeometry<DIM>();

   /*!
    * @brief Compute the overlap in edge-centered index space on the
    * boundaries of the source box geometry and the destination box geometry.
    */
   virtual tbox::Pointer<hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
    * @brief Return the box for this outeredge box geometry object.
    */
   const hier::Box<DIM>& getBox() const;
 
   /*!
    * @brief Return the ghost cell width for this outeredge box geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

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

