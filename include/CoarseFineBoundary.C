//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/patches/CoarseFineBoundary.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1888 $
// Modified:	$LastChangedDate: 2008-01-22 16:24:44 -0800 (Tue, 22 Jan 2008) $
// Description:	For describing coarse-fine boundary interfaces
//

#ifndef included_hier_CoarseFineBoundary_C
#define included_hier_CoarseFineBoundary_C

#include "CoarseFineBoundary.h"


namespace SAMRAI {
    namespace hier {

template<int DIM>  CoarseFineBoundary<DIM>::CoarseFineBoundary()
{
   d_npatches.resizeArray(1);
   d_npatches[0] = -1;
   d_nblocks = 1;
   d_boundary_boxes.resizeArray(1);
   d_mblk_hierarchy.setNull();
}

template<int DIM>  CoarseFineBoundary<DIM>::CoarseFineBoundary(
   const PatchHierarchy<DIM>& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_ghost_width > IntVector<DIM>(-1));
#endif
   d_mblk_hierarchy.setNull();
   d_npatches.resizeArray(1);
   d_npatches[0] = -1;
   d_nblocks = 1;
   computeFromHierarchy(hierarchy, ln, max_ghost_width);
}

template<int DIM>  CoarseFineBoundary<DIM>::CoarseFineBoundary(
   const tbox::Pointer< MultiblockPatchHierarchy<DIM> >& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
: d_mblk_hierarchy(hierarchy)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_ghost_width > IntVector<DIM>(-1));
#endif

   d_nblocks = hierarchy->getNumberBlocks();
   d_npatches.resizeArray(d_nblocks);

   for (int i = 0; i < d_nblocks; i++) {
      d_npatches[i] = -1;
   }

   computeFromHierarchy(*hierarchy, ln, max_ghost_width);
}

template<int DIM> CoarseFineBoundary<DIM>::~CoarseFineBoundary()
{
}


template<int DIM> void CoarseFineBoundary<DIM>::computeFromHierarchy(
   const PatchHierarchy<DIM>& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_ghost_width > IntVector<DIM>(-1));
#endif
   d_boundary_boxes.resizeArray(1);

   const hier::PatchLevel<DIM>& level =
      dynamic_cast<const hier::PatchLevel<DIM>&> (*hierarchy.getPatchLevel(ln));
   const hier::PatchLevel<DIM> &level0 =
      dynamic_cast<const hier::PatchLevel<DIM>&> (*hierarchy.getPatchLevel(0));
   computeFromLevel(level,
                    level0,
                    max_ghost_width);
}

template<int DIM> void CoarseFineBoundary<DIM>::computeFromHierarchy(
   const MultiblockPatchHierarchy<DIM>& hierarchy,
   int ln,
   const IntVector<DIM>& max_ghost_width)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(max_ghost_width > IntVector<DIM>(-1));
#endif
   d_boundary_boxes.resizeArray(d_nblocks);

   const tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level =
      hierarchy.getPatchLevel(ln);
   const tbox::Pointer<hier::MultiblockPatchLevel<DIM> > level0 =
      hierarchy.getPatchLevel(0);
   computeFromLevel(*level,
                    *level0,
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

   d_npatches[0] = level.getNumberOfPatches();

   const hier::IntVector<DIM>& ratio = level.getRatio();

   tbox::Pointer< hier::GridGeometry<DIM> > grid_geometry = level.getGridGeometry();

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
   BoxArray<DIM> periodic_adjusted_level_domain(level.getBoxes());

   /*
    * Add to periodic_adjusted_domain the continuation of the
    * level across periodic boundaries.
    */
   if ( is_periodic ) {
      addPeriodicImageBoxes(periodic_adjusted_phys_domain,
                            level0.getShiftsForLevel());

      addPeriodicImageBoxes(periodic_adjusted_level_domain,
                            level.getShiftsForLevel());
   }

   periodic_adjusted_phys_domain.refine(ratio);

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
   d_boundary_boxes[0].resizeArray(d_npatches[0]*DIM);

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
      d_boundary_boxes[0].getPointer(),
      level,
      use_periodic_shift,
      max_ghost_width,
      adjusted_level_domain,
      do_all_patches);

}

template<int DIM> void CoarseFineBoundary<DIM>::computeFromLevel(
   const MultiblockPatchLevel<DIM>& level,
   const MultiblockPatchLevel<DIM>& level0,
   const IntVector<DIM>& max_ghost_width)
{
   for (int i = 0; i < d_nblocks; i++) {
      tbox::Pointer< hier::PatchLevel<DIM> > block_level =
         level.getPatchLevelForBlock(i);
      tbox::Pointer< hier::PatchLevel<DIM> > block_level0 =
         level0.getPatchLevelForBlock(i); 

      clear(i);

      d_npatches[i] = block_level->getNumberOfPatches();

      const hier::IntVector<DIM>& ratio = block_level->getRatio();

      tbox::Pointer< hier::GridGeometry<DIM> > grid_geometry =
         block_level->getGridGeometry();

      BoxArray<DIM> phys_domain(block_level0->getBoxes());
      BoxArray<DIM> level_domain(block_level->getBoxes());

      phys_domain.refine(ratio);

      /*
       * Create a pseudo-domain -- the union of the current level's boxes
       * on all blocks in terms of the current block's index space.
       */

      BoxList<DIM> pseudo_domain(phys_domain);

      for (typename tbox::List< typename MultiblockPatchHierarchy<DIM>::Neighbor>::Iterator
           ni(d_mblk_hierarchy->getNeighbors(i)); ni; ni++) {

         BoxList<DIM> neighbor_domain(ni().d_translated_domain);
         neighbor_domain.refine(ratio);

         pseudo_domain.unionBoxes(neighbor_domain);

      } 

      /*
       * In reduced connectivity case, add a box at the singularity point
       * to the pseudo-domain.
       */

      if (d_mblk_hierarchy->reducedConnectivityExists(i)) {
         BoxList<DIM> sing_boxes(d_mblk_hierarchy->getSingularityBoxList(i));
         sing_boxes.refine(ratio);

         sing_boxes.grow(max_ghost_width);

         sing_boxes.removeIntersections(pseudo_domain);
         pseudo_domain.unionBoxes(sing_boxes);
      }

      /*
       * Make a list containing the level boxes for the current block,
       * then add more boxes as a buffer around physical domain boundaries.
       * This prevents physical boundaries from being identified as coarse-fine
       * boundaries.
       */ 

      BoxList<DIM> adjusted_level_domain_list(level_domain);

      int num_patches = block_level->getNumberOfPatches();
      for (int p = 0; p < num_patches; p++) {
         if (block_level->patchTouchesRegularBoundary(p)) {
            const Box<DIM>& patch_box = block_level->getBoxForPatch(p);

            BoxList<DIM> no_shift_boxes(patch_box);
            no_shift_boxes.grow(max_ghost_width);
            no_shift_boxes.removeIntersections(pseudo_domain);
            adjusted_level_domain_list.unionBoxes(no_shift_boxes);

         }
      }

      /*
       * Add buffer of boxes that exist on the current level across
       * block boundaries from the current block.  This prevents block
       * boundaries from being identified as coarse-fine boundaries when they
       * are not.
       */
      for (typename tbox::List<typename MultiblockPatchHierarchy<DIM>::Neighbor>::Iterator
           ni(d_mblk_hierarchy->getNeighbors(i)); ni; ni++) {

         tbox::Pointer< hier::PatchLevel<DIM> > neighbor_level =
            level.getPatchLevelForBlock(ni().d_id);

         hier::BoxArray<DIM> neighbor_boxes(neighbor_level->getBoxes());

         if (neighbor_boxes.size()) {
            d_mblk_hierarchy->translateBoxArray(neighbor_boxes,
                                                ratio,
                                                i,
                                                ni().d_id);

            BoxList<DIM> neighbor_boxes_to_add(phys_domain);
            neighbor_boxes_to_add.grow(max_ghost_width);

            neighbor_boxes_to_add.intersectBoxes(neighbor_boxes);

            adjusted_level_domain_list.unionBoxes(neighbor_boxes_to_add); 
         }
      }

      BoxArray<DIM> adjusted_level_domain(adjusted_level_domain_list);

      d_boundary_boxes[i].resizeArray(d_npatches[i]*DIM);

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
       * Send the adjusted level boxes as the domain for the
       * remove-intersection-with-domain operation.  This causes that
       * operation to remove non-coarse-fine (that is, fine-fine) boxes
       * along the periodic boundaries, leaving the coarse-fine boundary
       * boxes.
       *
       * Send the adjusted domain for the limit-domain intersect
       * operation.  This removes the boundaries that are on the physical
       * boundaries, which is what we want because there is no possibility
       * of a coarse-fine boundary there.
       */
      bool do_all_patches = true;
      IntVector<DIM> use_periodic_shift(0);
      grid_geometry->computeBoundaryBoxesOnLevel(
         d_boundary_boxes[i].getPointer(),
         *block_level,
         use_periodic_shift,
         max_ghost_width,
         adjusted_level_domain,
         do_all_patches);
   }

}

template<int DIM> void CoarseFineBoundary<DIM>::addPeriodicImageBoxes(
   BoxArray<DIM>& boxes,
   const tbox::Array<tbox::List<IntVector<DIM> > >& shifts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( shifts.size() == boxes.getNumberOfBoxes() );
#endif
   int current_size = boxes.getNumberOfBoxes();

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

   for ( ip=0; ip<old_size; ++ip ) {
      const Box<DIM>& unshifted_box = boxes.getBox(ip);
      const tbox::List< IntVector<DIM> >& shifts_list = shifts[ip];
      if ( ! shifts_list.isEmpty() ) {
         typename tbox::List< IntVector<DIM> >::Iterator sh;
         for ( sh = shifts_list.listStart(); sh; sh++ ) {
            Box<DIM> shifted_box(unshifted_box);
            shifted_box.shift( (*sh) );
            boxes(current_size++) = shifted_box;
         }
      }
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( current_size == new_size );
#endif
}

template<int DIM> void CoarseFineBoundary<DIM>::clear(const int block_number) {
   d_npatches[block_number] = -1;
   d_boundary_boxes[block_number].setNull();
}



template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getNodeBoundaries(int patch_num,
                                              int block_num) const
{
  return getBoundaries( patch_num, DIM, block_num );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getEdgeBoundaries(int patch_num,
                                              int block_num) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( DIM < 2 ) {
      TBOX_ERROR("CoarseFineBoundary::getEdgeBoundaries():  There is\n"
                 <<"no edge boundaries in " << DIM << "d.\n");
   }
#endif
  return getBoundaries( patch_num, DIM-1, block_num );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getFaceBoundaries(int patch_num,
                                              int block_num) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( DIM < 3 ) {
      TBOX_ERROR("CoarseFineBoundary::getFaceBoundaries():  There is\n"
                 <<"no face boundaries in " << DIM << "d.\n");
   }
#endif
   return getBoundaries( patch_num, DIM-2, block_num );
}

template<int DIM> const tbox::Array< BoundaryBox<DIM> >&
   CoarseFineBoundary<DIM>::getBoundaries(int patch_num,
                                          int boundary_type,
                                          int block_num) const
{
   if ( d_npatches[block_num] < 0 ) {
      TBOX_ERROR("The boundary boxes have not been computed.");
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( patch_num >= 0 && patch_num < d_npatches[block_num] );
   TBOX_ASSERT( boundary_type >= 0 && boundary_type <= DIM );
#endif
   return d_boundary_boxes[block_num][ patch_num*DIM + (boundary_type-1) ];
}

template<int DIM> void CoarseFineBoundary<DIM>::printClassData( std::ostream &os ) const {
   os << "\nCoarseFineBoundary<DIM>::printClassData...";
   os << "\n	number of patches: " << d_npatches[0];
   int pn, btype;
   for ( pn=0; pn<d_npatches[0]; ++pn ) {
      os << "\n		patch " << pn << '/' << d_npatches[0];
      for ( btype=0; btype<DIM; ++btype ) {
         os << "\n			type " << btype;
         const tbox::Array< BoundaryBox<DIM> >
            &array_of_boxes = d_boundary_boxes[0][pn*DIM+btype];
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
