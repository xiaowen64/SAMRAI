//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/boxgeometry/SideGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	hier::Box geometry information for side centered objects
//

#ifndef included_pdat_SideGeometry
#define included_pdat_SideGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "SideOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * Class SideGeometry<DIM> manages the mapping between the AMR index space
 * and the side-centered geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between side-
 * centered box geometries for communication operations.
 *
 * See header file for SideData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::SideOverlap
 */

template<int DIM> class SideGeometry : public hier::BoxGeometry<DIM>
{
public:

   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef SideOverlap<DIM> Overlap;


  /*!
    * @brief Convert an AMR index box space box into an side geometry box.
    * An side geometry box extends the given AMR index box space box
    * by one in upper dimension for the side normal coordinate direction.
    */
   static hier::Box<DIM> toSideBox(const hier::Box<DIM>& box,
                                   int side_normal);

   /*!
    * @brief Construct the side geometry object given an AMR index
    * space box, ghost cell width and directions vector indicating
    * which coordinate directions are allocated.
    */
   SideGeometry(const hier::Box<DIM>& box,
                const hier::IntVector<DIM>& ghosts,
                const hier::IntVector<DIM>& directions); 

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~SideGeometry<DIM>();

   /*!
    * @brief Compute the overlap in side-centered index space between
    * the source box geometry and the destination box geometry.
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
    * @brief Return the box for this side centered box geometry
    * object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * @brief Return the ghost cell width for this side centered box
    * geometry object.
    */
   const hier::IntVector<DIM>& getGhosts() const;

   /*!
    * Return constant reference to vector describing which coordinate
    * directions managed by this side geometry object.  
    * 
    * A vector entry of zero indicates that this object will not perform 
    * operations involving the corresponding coordinate direction.  
    * A non-zero value indicates otherwise.
    */
   const hier::IntVector<DIM>& getDirectionVector() const;

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
