//
// File:	CoarseFineBoundary.C
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	For describing coarse-fine boundary interfaces
//

#ifndef included_hier_CoarseFineBoundary_C
#define included_hier_CoarseFineBoundary_C

#include "CoarseFineBoundary.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

namespace SAMRAI {
    namespace hier {

template<int DIM>  CoarseFineBoundary<DIM>::CoarseFineBoundary()
{
   d_npatches = -1;
}

template<int DIM>  CoarseFineBoundary<DIM>::CoarseFineBoundary(
   const PatchHierarchy<DIM>& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(max_ghost_width > IntVector<DIM>(-1));
#endif
   d_npatches = -1;
   computeFromHierarchy(hierarchy, ln, max_ghost_width);
}

template<int DIM> void CoarseFineBoundary<DIM>::computeFromHierarchy(
   const PatchHierarchy<DIM>& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(max_ghost_width > IntVector<DIM>(-1));
#endif
   const hier::PatchLevel<DIM>& level =
      dynamic_cast<const hier::PatchLevel<DIM>&> (*hierarchy.getPatchLevel(ln));
   const hier::PatchLevel<DIM> &level0 =
      dynamic_cast<const hier::PatchLevel<DIM>&> (*hierarchy.getPatchLevel(0));
   computeFromLevel(level,
                    level0,
                    max_ghost_width);
}

/*
************************************************************************
* Use grid_geometry.computeBoundaryGeometry function,                  *
* setting up the arguments in a way that will generate                 *
* the coarse-fine boundary (instead of the domain boundary).           *
************************************************************************
*/
template<int DIM> void CoarseFineBoundary<DIM>::computeFromLevel(
   const PatchLevel<DIM>& level,
   const PatchLevel<DIM>& level0,
   const IntVector<DIM>& max_ghost_width) 
{
   clear();

   d_npatches = level.getNumberOfPatches();

   const hier::IntVector<DIM>& ratio = level.getRatio();

   tbox::Pointer< hier::GridGeometry<DIM> > grid_geometry = level0.getGridGeometry();

   /*
    * Get the domain's periodic shift.
    */
   const hier::IntVector<DIM> periodic_shift(grid_geometry->getPeriodicShift(ratio));

   bool is_periodic = false;
   for (int i=0; i<DIM; ++i) {
      is_periodic = is_periodic || periodic_shift(i);
   }

   /*
    * The periodic adjusted domain is a (modified) copy of the
    * domain (for the periodic case).  The periodic adjusted
    * level is a (modified) copy of the level boxes (for the
    * periodic case).
    */
   BoxArray<DIM> periodic_adjusted_phys_domain(level0.getBoxes());
   periodic_adjusted_phys_domain.refine(ratio);
   BoxArray<DIM> periodic_adjusted_level_domain(level.getBoxes());

   /*
    * Add to periodic_adjusted_domain the continuation of the
    * level across periodic boundaries.
    */
   if ( is_periodic ) {
      bool is_level_zero = true;
      addPeriodicImageBoxes(periodic_adjusted_phys_domain,
                            ratio,
                            level0,
                            is_level_zero);

      is_level_zero = false;
      addPeriodicImageBoxes(periodic_adjusted_level_domain,
                            ratio,
                            level,
                            is_level_zero);
   }

   /*
    * Here we add some boxes outside of non-periodic boundaries to the
    * adjusted level.  For each patch that touches a regular boundary,
    * grow the patch box and any periodic images of the patch box by
    * the max ghost width.  Remove intersections with the periodic
    * adjusted physical domain.  Add what remains to the adjusted level.
    *
    * This will ensure that ensuing call to create boundary boxes will not
    * create boundary boxes at the locations where the level touches a
    * non-periodic physical boundary, but only where there is a coarse-fine
    * interface in the domain interior (A periodic boundary is considered
    * part of the domain interior for this purpose).
    */

   BoxList<DIM> adjusted_level_domain_list(periodic_adjusted_level_domain);
   const tbox::Array< tbox::List< IntVector<DIM> > >& shifts =
      level.getShiftsForLevel();

   int num_patches = level.getNumberOfPatches();
   for (int p = 0; p < num_patches; p++) {
      if (level.patchTouchesRegularBoundary(p)) {
         const Box<DIM>& patch_box = level.getBoxForPatch(p);

         for (typename tbox::List< IntVector<DIM> >::Iterator sh(shifts[p]);
              sh; sh++) {
            BoxList<DIM> shift_boxes(Box<DIM>::shift(patch_box, sh()));
            shift_boxes.grow(max_ghost_width);
            shift_boxes.removeIntersections(periodic_adjusted_phys_domain);
            adjusted_level_domain_list.unionBoxes(shift_boxes);
         }

         BoxList<DIM> no_shift_boxes(patch_box);
         no_shift_boxes.grow(max_ghost_width);
         no_shift_boxes.removeIntersections(periodic_adjusted_phys_domain);
         adjusted_level_domain_list.unionBoxes(no_shift_boxes);
      }
   }

   BoxArray<DIM> adjusted_level_domain(adjusted_level_domain_list);

   /*
    * Allocate enough room for DIM types of boundary boxes
    * on d_npatches patches.
    */
   d_boundary_boxes.resizeArray(d_npatches*DIM);

   /*
    * Call GridGeometry::computeBoundaryGeometry with arguments contrived
    * such that they give the coarse-fine boundaries instead of the domain
    * boundaries.  The basic algorithm used by
    * GridGeometry::computeBoundaryGeometry is
    * 1. grow boxes by ghost width
    * 2. remove intersection with domain
    * 3. reorganize and classify resulting boxes
    *
    * This is how we get GridGeometry::computeBoundaryGeometry to
    * compute the coarse-fine boundary instead of the physical boundary.
    *
    * Since we handle the periodic boundaries ourselves, do not treat
    * them differently from regular boundaries.  State that all boundaries
    * are non-periodic boundaries.
    *
    * Send the periodic-adjusted level boxes as the domain for the
    * remove-intersection-with-domain operation.  This causes that
    * operation to remove non-coarse-fine (that is, fine-fine) boxes
    * along the periodic boundaries, leaving the coarse-fine boundary
    * boxes.
    *
    * Send the periodic-adjusted domain for the limit-domain intersect
    * operation.  This removes the boundaries that are on the non-periodic
    * boundaries, which is what we want because there is no possibility
    * of a coarse-fine boundary there.
    */
   bool do_all_patches = true;
   IntVector<DIM> use_periodic_shift(0);
   grid_geometry->computeBoundaryBoxesOnLevel(
      d_boundary_boxes.getPointer(),
      level,
      use_periodic_shift,
      max_ghost_width,
      adjusted_level_domain,
      do_all_patches);

}

template<int DIM> void CoarseFineBoundary<DIM>::addPeriodicImageBoxes(
   BoxArray<DIM>& boxes,
   const IntVector<DIM>& shift_refine_ratio,
   const PatchLevel<DIM>& level,
   bool is_level_zero)
{
   int current_size = boxes.getNumberOfBoxes();
   const tbox::Array<tbox::List<IntVector<DIM> > >& shifts =
      level.getShiftsForLevel();

   int ip;
   /*
    * Count number of boxes that must be added to boxes.
    * And resize boxes accordingly before adding the
    * periodic images to it.
    */
   int new_size = current_size;
   for ( ip=0; ip<current_size; ++ip ) {
      new_size += shifts[ip].getNumberItems();
   }
   boxes.resizeBoxArray( new_size );

   /*
    * For all the possible shifts of all patches,
    * compute the shifted box and add it to boxes.
    * This completes the addition of images boxes.
    */
   const int old_size = current_size;

   IntVector<DIM> shift_multiplier(1);
   if (is_level_zero) {
      shift_multiplier = shift_refine_ratio;
   }

   for ( ip=0; ip<old_size; ++ip ) {
      const Box<DIM>& unshifted_box = boxes.getBox(ip);
      const tbox::List< IntVector<DIM> >& shifts_list = shifts[ip];
      if ( ! shifts_list.isEmpty() ) {
         typename tbox::List< IntVector<DIM> >::Iterator sh;
         for ( sh = shifts_list.listStart(); sh; sh++ ) {
            Box<DIM> shifted_box(unshifted_box);
            shifted_box.shift( (*sh) * shift_multiplier );
            boxes(current_size++) = shifted_box;
         }
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( current_size == new_size );
#endif
}

template<int DIM> void CoarseFineBoundary<DIM>::clear() {
   d_npatches = -1;
   d_boundary_boxes.setNull();
}



template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getNodeBoundaries(int pn) const
{
  return getBoundaries( pn, DIM );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getEdgeBoundaries(int pn) const
{
  return getBoundaries( pn, DIM-1 );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getFaceBoundaries(int pn) const
{
  return getBoundaries( pn, DIM-2 );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getBoundaries(int pn,
                                          int boundary_type) const
{
   if ( d_npatches < 0 ) {
      TBOX_ERROR("The boundary boxes have not been computed.");
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( pn >= 0 && pn < d_npatches );
   assert( boundary_type >= 0 && boundary_type <= DIM );
#endif
   return d_boundary_boxes[ pn*DIM + (boundary_type-1) ];
}

template<int DIM> void CoarseFineBoundary<DIM>::printClassData( ostream &os ) const {
   os << "\nCoarseFineBoundary<DIM>::printClassData...";
   os << "\n	number of patches: " << d_npatches;
   int pn, btype;
   for ( pn=0; pn<d_npatches; ++pn ) {
      os << "\n		patch " << pn << '/' << d_npatches;
      for ( btype=0; btype<DIM; ++btype ) {
         os << "\n			type " << btype;
         const tbox::Array< BoundaryBox<DIM> >
            &array_of_boxes = d_boundary_boxes[pn*DIM+btype];
         int num_boxes = array_of_boxes.getSize();
         int bn;
         for ( bn=0; bn<num_boxes; ++bn ) {
            os << "\n				box "
               << bn << "/" << num_boxes << ":";
            os << array_of_boxes[bn].getBox();
         }
      }
   }
   os << "\n";
}


}
}

#endif
