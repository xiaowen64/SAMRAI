//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/OuterfaceGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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

/*!
 * Class OuterfaceGeometry<DIM> manages the mapping between the AMR index
 * and the outerface geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between outerface
 * box geometries and face or outerface box geometries for communication
 * operations.
 *
 * See header file for OuterfaceData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::FaceGeometry
 * @see pdat::FaceOverlap
 */

template<int DIM> class OuterfaceGeometry : public hier::BoxGeometry<DIM>
{
public:
   /*!
    * @brief Construct an outerface geometry object given an AMR index
    * space box and ghost cell width.
    */
   OuterfaceGeometry(const hier::Box<DIM>& box,
                     const hier::IntVector<DIM>& ghosts);
 
   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OuterfaceGeometry<DIM>();
 
   /*!
    * @brief Compute the overlap in face-centered index space on the
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
    * @brief Return the box for this outerface box geometry object.
    */
   const hier::Box<DIM>& getBox() const;
 
   /*!
    * @brief Return the ghost cell width for this outerface box geometry object.
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
