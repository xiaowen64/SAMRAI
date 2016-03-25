//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/boxgeometry/OutersideGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	hier::Box geometry information for outerside centered objects
//

#ifndef included_pdat_OutersideGeometry
#define included_pdat_OutersideGeometry

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxGeometry.h"
#include "BoxOverlap.h"
#include "SideOverlap.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI {
    namespace pdat {

template<int DIM> class SideGeometry;

/*!
 * Class OutersideGeometry<DIM> manages the mapping between the AMR index
 * and the outerside geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between outerside
 * box geometries and side or outerside box geometries for communication
 * operations.
 *
 * See header file for OutersideData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::SideGeometry
 * @see pdat::SideOverlap
 */

template<int DIM> class OutersideGeometry : public hier::BoxGeometry<DIM>
{
public:

   /*!
    * The BoxOverlap implemenation for this geometry.
    */
   typedef SideOverlap<DIM> Overlap;


   /*!
    * @brief Construct an outerside geometry object given an AMR index
    * space box and ghost cell width.
    */
   OutersideGeometry(const hier::Box<DIM>& box,
                     const hier::IntVector<DIM>& ghosts);

   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OutersideGeometry<DIM>();

   /*!
    * @brief Compute the overlap in side-centered index space on the
    * boundaries of the source box geometry and the destination box geometry.
    */
   virtual tbox::Pointer< hier::BoxOverlap<DIM> > calculateOverlap(
      const hier::BoxGeometry<DIM>& dst_geometry,
      const hier::BoxGeometry<DIM>& src_geometry,
      const hier::Box<DIM>& src_mask,
      const bool overwrite_interior,
      const hier::IntVector<DIM>& src_offset,
      const bool retry) const;

   /*!
    * @brief Return the box for this outerside box geometry object.
    */
   const hier::Box<DIM>& getBox() const;

   /*!
    * @brief Return the ghost cell width for this outerside box geometry object.
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
