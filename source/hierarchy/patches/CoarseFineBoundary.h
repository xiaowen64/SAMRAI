//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/patches/CoarseFineBoundary.h $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1782 $
// Modified:	$LastChangedDate: 2007-12-17 13:04:51 -0800 (Mon, 17 Dec 2007) $
// Description:	For describing coarse-fine boundary interfaces
//

#ifndef included_hier_CoarseFineBoundary
#define included_hier_CoarseFineBoundary

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Array
#include "tbox/Array.h"
#endif
#ifndef included_hier_BoundaryBox
#include "BoundaryBox.h"
#endif
#ifndef included_hier_PatchHierarchy
#include "PatchHierarchy.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif

namespace SAMRAI {
    namespace hier {

/*!
 *  @brief Utility class to construct and maintain a description of the coarse-fine
 *  boundary between a patch level and some coarser level.
 *
 *  A coarse-fine boundary box is a BoundaryBox object, but it is generated
 *  differently than a typical boundary box maintained by a patch geometry object.
 *  The boundary type and location identifiers for regular boundary boxes apply
 *  to coarse-fine boundary boxes.  However, a boundary box serving as a coarse-fine
 *  boundary box describes part of the boundary of a given patch with its next
 *  coarser AMR hierarchy level.  It does not intersect any other patch on the same
 *  level nor does it lie on a physical domain boundary, except where the physical
 *  boundary is periodic and the appropriate continuation of that boundary is part
 *  of a coarser patch level.
 *
 *  The coarse-fine boundary is created from two adjacent hierarchy levels (typically),
 *  but the description lives on (refers to the index space of) the finer level.
 *  Since the coarse-fine boundary describes the boundary to the next coarser level,
 *  the coarsest level (i.e., level zero in an AMR hierarchy) has no coarse-fine
 *  boundary.
 *
 *  Each CoarseFineBoundary object corresponds to one level,
 *  so to represent a hierarchy, you would need an array or list of
 *  such objects.
 */

template<int DIM> class CoarseFineBoundary : public tbox::DescribedClass
{
public:

   /*!
    * @brief Construct a CoarseFineBoundary<DIM> object with no boundary boxes.
    */
   CoarseFineBoundary();

   /*!
    * @brief Construct a CoarseFineBoundary<DIM> object for the specified
    * level in the given patch hierarchy.
    *
    * @param hierarchy       Patch hierarchy in which the patch level resides.
    * @param ln              Level number of level of computed coarse-fine boundary.
    * @param max_ghost_width Max ghost width for which to generate boundary
    *                        boxes.  The ghost width determines the extent
    *                        of the boundary boxes along the level domain boundary,
    *                        similar to regular domain boundary boxes.  Note that
    *                        as in the case of regular boundary boxes, each box
    *                        will always be one cell wide in the direction
    *                        perpendicular to the patch boundary.
    *
    * Note that if level number is zero, the coarse-fine boundary will be empty.
    */
   CoarseFineBoundary<DIM>(
      const PatchHierarchy<DIM>& hierarchy,
      int ln,
      const IntVector<DIM>& max_ghost_width);

   /*!
    * @brief Construct a CoarseFineBoundary<DIM> object for the specified
    * level in the given patch hierarchy.
    *
    * @param hierarchy       Patch hierarchy in which the patch level resides.
    * @param ln              Level number of level of computed coarse-fine boundary.
    * @param max_ghost_width Max ghost width for which to generate boundary
    *                        boxes.  The ghost width determines the extent
    *                        of the boundary boxes along the level domain boundary,
    *                        similar to regular domain boundary boxes.  Note that
    *                        as in the case of regular boundary boxes, each box
    *                        will always be one cell wide in the direction 
    *                        perpendicular to the patch boundary.
    *
    * Note that if level number is zero, the coarse-fine boundary will be empty.
    */
   void computeFromHierarchy(
      const PatchHierarchy<DIM>& hierarchy,
      int ln,
      const IntVector<DIM>& max_ghost_width);

   /*!
    * @brief Construct a CoarseFineBoundary<DIM> object for the specified
    * level based on a given level which is assumed to be the coarsest level
    * (i.e., level zero) in some patch hierarchy.
    *
    * @param level           Patch level of computed coarse-fine boundary.
    * @param level0          Coarsest patch level in hierarchy used to
    *                        compute coarse-fine boundary.
    * @param max_ghost_width Max ghost width for which to generate boundary
    *                        boxes.  The ghost width determines the extent
    *                        of the boundary boxes along the level domain boundary,
    *                        similar to regular domain boundary boxes.  Note that
    *                        as in the case of regular boundary boxes, each box
    *                        will always be one cell wide in the direction
    *                        perpendicular to the patch boundary.
    *
    * Note that if level and level0 are the same, the coarse-fine boundary
    * will be empty.
    */
   void computeFromLevel(
      const PatchLevel<DIM>& level,
      const PatchLevel<DIM>& level0,
      const IntVector<DIM>& max_ghost_width);

   /*!
    * @brief Clear all boundary data.
    */
   void clear();

   //@{
   /*!
    * @name Functions to get the computed coarse-fine boundaries.
    */

   /*!
    * @brief Get an array of boundary boxes of a given type
    * for a specified patch.
    *
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param pn            Patch number
    * @param boundary_type Boundary box type (see BoundaryBox class).    
    */
   const tbox::Array< BoundaryBox<DIM> >& getBoundaries(
      int pn,
      int boundary_type) const;

   /*!
    * @brief Get an array of node boundary boxes for a specified patch
    *        (see BoundaryBox class).
    *
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param pn Patch number
    */
   const tbox::Array< BoundaryBox<DIM> >& getNodeBoundaries(int pn) const;

   /*!
    * @brief Get an array of edge boundary boxes for a specified patch
    *        (see BoundaryBox class).
    *
    * Note that edge boxes are only meaningful if problem dimension is > 1.
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param pn Patch number
    */
   const tbox::Array< BoundaryBox<DIM> >& getEdgeBoundaries(int pn) const;

   /*!
    * @brief Get an array of face boundary boxes for a specified patch
    *        (see BoundaryBox class).
    *
    * Note that face boxes are only meaningful if problem dimension is > 2.
    * The specified patch must exist in the level used to compute
    * the internal state or it is an error.
    *
    * @param pn Patch number
    */
   const tbox::Array< BoundaryBox<DIM> >& getFaceBoundaries(int pn) const;

   //@}

   /*!
    * @brief Print out class data (mostly for debugging).
    */
   virtual void printClassData( std::ostream &os ) const;

private:

   /*!
    * @brief Take a set of boxes representing some domain and
    * append to it the immediate periodic images of the boxes.
    *
    * If there is no periodic directions in the grid,
    * there will be no change.
    *
    * The image boxes help form a virtual domain with which
    * to trick the grid geometry object to compute the coarse-fine
    * boundary instead of the physical boundary.
    *
    * @param boxes Box array to append to.  This function will
    *        append the periodic image boxes to this array.
    * @param shifts Periodic shifts.
    */
   void addPeriodicImageBoxes(
      BoxArray<DIM>& boxes,
      const tbox::Array<tbox::List<IntVector<DIM> > >& shifts);

   /*!
    * @brief Number of patches on the level for which coarse-fine
    * boundary has been computed.
    *
    * This is set to >= 0 when the boundary boxes are generated.
    * Otherwise, it is set to -1.  We do not use the size of
    * d_boundary_boxes to determine if boundary has been generated
    * because it is possible to have no patch on a level.
    */
   int d_npatches;

   /*!
    * @brief Patch boundary boxes describing the coarse-fine boundary.
    *
    * The outer array, sized by DIM times the number of patches on the level,
    * representing for each patch, the DIM types of boundary boxes.
    * The inner array, sized by the number of BoundaryBox<DIM> of
    * a given type, for a given patch.  So, the array of BoundaryBox<DIM>
    * of type i for patch number pn is d_boundary_boxes[pn*DIM+(i-1)].
    * The reason for this is due to the way the boundary boxes are
    * computed in GridGeometry<DIM>::computeBoundaryGeometry.
    */
   tbox::Array< tbox::Array< BoundaryBox<DIM> > > d_boundary_boxes;

};

}
}

#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarseFineBoundary.C"
#endif
