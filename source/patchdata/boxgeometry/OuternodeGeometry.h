//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/patchdata/boxgeometry/OuternodeGeometry.h $
// Package:	SAMRAI patch data geometry
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
 * Class OuternodeGeometry<DIM> manages the mapping between the AMR index
 * and the outernode geometry index space.  It is a subclass of
 * hier::BoxGeometry<DIM> and it computes intersections between outernode
 * box geometries and node or outernode box geometries for communication
 * operations.
 *
 * See header file for OuternodeData<DIM> class for a more detailed
 * description of the data layout.
 *
 * @see hier::BoxGeometry
 * @see pdat::NodeGeometry
 * @see pdat::NodeOverlap
*/

template<int DIM> class OuternodeGeometry : public hier::BoxGeometry<DIM>
{
public:
   /*!
    * @brief Construct an outernode geometry object given an AMR index
    * space box and ghost cell width.
    */
   OuternodeGeometry(const hier::Box<DIM>& box,
                     const hier::IntVector<DIM>& ghosts);
 
   /*!
    * @brief The virtual destructor does nothing interesting.
    */
   virtual ~OuternodeGeometry<DIM>();
 
   /*!
    * @brief Compute the overlap in node-centered index space on the
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
    * @brief Return the box for this outernode box geometry object.
    */
   const hier::Box<DIM>& getBox() const;
 
   /*!
    * @brief Return the ghost cell width for this outernode box geometry object.
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
