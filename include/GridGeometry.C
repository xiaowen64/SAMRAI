//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/source/hierarchy/patches/GridGeometry.C $
// Package:	SAMRAI hierarchy package
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1846 $
// Modified:	$LastChangedDate: 2008-01-11 09:51:05 -0800 (Fri, 11 Jan 2008) $
// Description: Base class for geometry management in AMR hierarchy
//

#ifndef included_hier_GridGeometry_C
#define included_hier_GridGeometry_C

#include "SAMRAI_config.h"

#include <stdlib.h>


#include "GridGeometry.h"
#include "BoundaryLookupTable.h"
#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "PatchLevel.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#ifdef DEBUG_NO_INLINE
#include "GridGeometry.I"
#endif
namespace SAMRAI {
    namespace hier {

/*
*************************************************************************
*                                                                       *
* Constructor initializes basic data members to default state.          *
*                                                                       *
*************************************************************************
*/

template<int DIM>  GridGeometry<DIM>::GridGeometry(const std::string &object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
#endif
   d_object_name = object_name;
   d_periodic_shift = IntVector<DIM>(0);
   d_max_data_ghost_width = IntVector<DIM>(-1);

}

/*
*************************************************************************
*                                                                       *
* Empty destructor.                                                     *
*                                                                       *
*************************************************************************
*/

template<int DIM>  GridGeometry<DIM>::~GridGeometry()
{
}

/*
*************************************************************************
*                                                                       *
* Compute boundary boxes for all patches in patch level.  The domain    *
* array describes the interior of the level index space.  Note that     *
* boundaries is assumed to be an array of DIM * #patches Arrays of     *
* BoundaryBoxes.                                                        *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::computeBoundaryBoxesOnLevel(
   tbox::Array<hier::BoundaryBox<DIM> > boundaries[],
   const PatchLevel<DIM>& level,
   const IntVector<DIM>& periodic_shift,
   const IntVector<DIM>& ghost_width,
   const BoxArray<DIM>& domain,
   bool do_all_patches) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(ghost_width >= hier::IntVector<DIM>(0));

   int num_per_dirs = 0;
   for (int i = 0; i < DIM; i++) {
      if (periodic_shift(i)) {
         num_per_dirs++;
      }
   }
#endif

   for (typename PatchLevel<DIM>::Iterator ip(&level); ip; ip++) {
      tbox::Pointer< Patch<DIM> > patch = level.getPatch(ip());
      const int n = DIM * ip();

      if (patch->getPatchGeometry()->getTouchesRegularBoundary() ||
          do_all_patches) {

         Box<DIM> box(patch->getBox());

         /* 
          * patch_boundaries is an array of DIM BoxLists.
          * patch_boundaries[DIM-1] will store boundary boxes of the node type.
          * If DIM > 1, patch_boundaries[DIM-2] will store boundary boxes of
          * the edge type, and if DIM > 2, patch_boundaries[DIM-3] will store
          * boundary boxes of the face type.
          */

         getBoundaryBoxes(&boundaries[n], box, domain, ghost_width,
                          periodic_shift);

         for (int j = 0; j < DIM; j++) {

#ifdef DEBUG_CHECK_ASSERTIONS
            for (int k = 0; k < boundaries[n+j].getSize(); k++) {
               TBOX_ASSERT(checkBoundaryBox(boundaries[n+j][k], *patch,
                                       domain, num_per_dirs, ghost_width));
            }
#endif
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* For each patch in the level, use box intersection operation to        *
* determine what kind of boundaries, if any the patch touches.  Call    *
* Patch<DIM> functions to set flags that store this information once it *
* is found.                                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::findPatchesTouchingBoundaries(
   tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_regular_bdry,
   tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_periodic_bdry,
   const PatchLevel<DIM>& level,
   const IntVector<DIM>& periodic_shift,
   const BoxArray<DIM>& domain) const
{
   touches_regular_bdry.resizeArray(level.getNumberOfPatches());
   touches_periodic_bdry.resizeArray(level.getNumberOfPatches());

   for (int np = 0; np < level.getNumberOfPatches(); np++) {
      touches_regular_bdry[np].resizeArray(DIM);
      touches_periodic_bdry[np].resizeArray(DIM);
      for (int i = 0; i < DIM; i++) {
         touches_regular_bdry[np][i].resizeArray(2);
         touches_periodic_bdry[np][i].resizeArray(2);
         for (int side = 0; side < 2; side++) {
            touches_regular_bdry[np][i][side] = false;
            touches_periodic_bdry[np][i][side] = false;
         }
      }
   }

   BoxList<DIM> domain_interior(domain);

   for (typename PatchLevel<DIM>::Iterator ip(&level); ip; ip++) {
      tbox::Pointer< Patch<DIM> > patch = level.getPatch(ip());
      Box<DIM> box(patch->getBox());

      /*
       * Create a list of boxes inside a layer of one cell outside the patch.
       * Remove the intersections with the domain's interior, so that only
       * boxes outside the physical domain remain in the list.
       */
      BoxList<DIM> bdry_list(Box<DIM>::grow(box, IntVector<DIM>(1)));
      bdry_list.removeIntersections(domain_interior);

      bool touches_any_boundary = false;
      if (bdry_list.getNumberItems() > 0) {
         touches_any_boundary = true;
      }

      if (!touches_any_boundary) {
         for (int nd = 0; nd < DIM; nd++) {
            for (int side = 0; side < 2; side++) {
               touches_regular_bdry[ip()][nd][side] = false;
               touches_periodic_bdry[ip()][nd][side] = false;
            }
         } 
      } else {
         bool bdry_located = false;
         for (int nd = 0; nd < DIM; nd++) {
            BoxList<DIM> lower_list(bdry_list);
            BoxList<DIM> upper_list(bdry_list);

            Box<DIM> test_box(box);

            test_box.growLower(nd, 1);
            lower_list.intersectBoxes(test_box);

            test_box = box; 
            test_box.growUpper(nd, 1);
            upper_list.intersectBoxes(test_box);

            if (lower_list.size()) {
               if (periodic_shift(nd)) {
                  touches_regular_bdry[ip()][nd][0] = false; 
                  touches_periodic_bdry[ip()][nd][0] = true; 
               } else {
                  touches_regular_bdry[ip()][nd][0] = true;
                  touches_periodic_bdry[ip()][nd][0] = false;
               }
               bdry_located = true;
            }

            if (upper_list.size()) {
               if (periodic_shift(nd)) {
                  touches_regular_bdry[ip()][nd][1] = false;
                  touches_periodic_bdry[ip()][nd][1] = true;
               } else {
                  touches_regular_bdry[ip()][nd][1] = true;
                  touches_periodic_bdry[ip()][nd][1] = false;
               }
               bdry_located = true;
            }
         }

         /*
          * By this point, bdry_located will have been set to true almost 
          * every time whenever touches_any_boundary is true.  The only way
          * it will not be true is if the domain is not a parallelpiped, and
          * the patch touches the boundary only at a location such as the
          * concave corner of an L-shaped domain.
          */
         if (!bdry_located) {
            for (int nd = 0; nd < DIM; nd++) {
               touches_periodic_bdry[ip()][nd][0] = false;
               touches_periodic_bdry[ip()][nd][1] = false;

               bool lower_side = false;
               bool upper_side = false;
               for (typename BoxList<DIM>::Iterator bl(bdry_list); bl; bl++) {
                  if (bl().lower()(nd) < box.lower(nd)) {
                     lower_side = true;
                  }
                  if (bl().upper()(nd) > box.upper(nd)) {
                     upper_side = true;
                  }
                  if (lower_side && upper_side) {
                     break;
                  }
               }

               if (lower_side) { 
                  touches_regular_bdry[ip()][nd][0] = true;
               } else {
                  touches_regular_bdry[ip()][nd][0] = false;
               }
               if (upper_side) { 
                  touches_regular_bdry[ip()][nd][1] = true;
               } else {
                  touches_regular_bdry[ip()][nd][1] = false;
               }
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Set geometry data for each patch on level.                            *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void hier::GridGeometry<DIM>::setGeometryOnPatches(
   hier::PatchLevel<DIM>& level,
   const hier::IntVector<DIM>& ratio_to_level_zero,
   tbox::Array< tbox::Array< tbox::Array<bool> > >& touches_regular_bdry,
   tbox::Array<tbox::Array< tbox::Array<bool> > >& touches_periodic_bdry,
   bool defer_boundary_box_creation)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != hier::IntVector<DIM>(0));
   if (DIM > 1) {
      for (int i = 0; i < DIM; i++) {
      TBOX_ASSERT( (ratio_to_level_zero(i)*ratio_to_level_zero((i+1)%DIM) > 0)
              || (ratio_to_level_zero(i) == 1)
              || (ratio_to_level_zero((i+1)%DIM) == 1) );
      }
   }
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(&level); ip; ip++) {
      tbox::Pointer<hier::Patch<DIM> > patch = level.getPatch(ip());
      setGeometryDataOnPatch(*patch, ratio_to_level_zero,
                             touches_regular_bdry[ip()],
                             touches_periodic_bdry[ip()]);
   }

   if (!defer_boundary_box_creation) {
      setBoundaryBoxes(level);
   } 
}

/*
*************************************************************************
*                                                                       *
* Set boundary boxes for each patch on level.                           *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void hier::GridGeometry<DIM>::setBoundaryBoxes(
   hier::PatchLevel<DIM>& level)
{
   tbox::Array<hier::BoundaryBox<DIM> >* boundaries =
      new tbox::Array< hier::BoundaryBox<DIM> >[DIM *
                                                level.getNumberOfPatches()];

   hier::BoxArray<DIM> domain(level.getPhysicalDomain());
   computeBoundaryBoxesOnLevel(
      boundaries,
      level,
      getPeriodicShift(),
      computeMaxGhostWidth(level.getPatchDescriptor()),
      domain);

   for (typename PatchLevel<DIM>::Iterator ip(&level); ip; ip++) {
      tbox::Pointer< Patch<DIM> > patch = level.getPatch(ip());
      const int n = DIM * ip();

      patch->getPatchGeometry()->setBoundaryBoxesOnPatch(&boundaries[n]);
   }

   delete[] boundaries;
}

/*
*************************************************************************
*                                                                       *
* Compute list of  potential periodic shift vectors for each patch      *
* on level.                                                             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::computeShiftsForLevel(
   tbox::Array< tbox::List< IntVector<DIM> > >& shifts,
   const PatchLevel<DIM>& level,
   const BoxArray<DIM>& physical_domain) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(shifts.getSize() == level.getNumberOfPatches());
#endif

   IntVector<DIM> periodic_shift = getPeriodicShift(level.getRatio());

   bool is_periodic = false;
   for (int i = 0; i < DIM; i++) {
      if (periodic_shift(i)) {
         is_periodic = true;
         break;
      }
   }

   int npatches = level.getNumberOfPatches(); 

   if (is_periodic) {
      for (int j = 0; j < npatches; j++) {
         shifts[j].clearItems();
         if (level.patchTouchesPeriodicBoundary(j)) {
            computeShiftsForPatch(shifts[j],
                                  level.getBoxForPatch(j),
                                  physical_domain,
                                  periodic_shift);
         }
      }
   } else {
      for (int k = 0; k < npatches; k++) {
         shifts[k].clearItems();
      }
   }
   
}

/*
*************************************************************************
*                                                                       *
* Compute the valid periodic shifts for the given box.                  *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::computeShiftsForPatch(
   tbox::List< IntVector<DIM> >& shifts, 
   const Box<DIM>& box, 
   const BoxArray<DIM>& domain, 
   const IntVector<DIM>& periodic_shift) const
{

   int num_periodic_dirs = 0;
   int i;

   for (i = 0; i < DIM; i++) {
      if (periodic_shift(i) != 0) {
         num_periodic_dirs++;
      }
   }

   if (num_periodic_dirs > 0) {

      BoundaryLookupTable<DIM>* blut =
         BoundaryLookupTable<DIM>::getLookupTable();

      const tbox::Array<int>& location_index_max =
         blut->getMaxLocationIndices();

      for (int d = 0; d < num_periodic_dirs; d++) {

         const int codim = d+1;

         for (int loc = 0; loc < location_index_max[d]; loc++) {

            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            bool need_to_test = true;
            for (int k = 0; k < dirs.size(); k++) {
               if (periodic_shift(dirs[k]) == 0) {
                  need_to_test = false;
                  break;
               }
            }

            if (need_to_test) {

               Box<DIM> border(box);
               IntVector<DIM> border_shift(0);

               tbox::Array<bool> is_upper(codim);
               for (int j = 0; j < codim; j++) {
                  if (blut->isUpper(loc, codim, j)) {
                     border.lower(dirs[j]) = box.upper(dirs[j]);
                     border.upper(dirs[j]) = box.upper(dirs[j]);
                     border_shift(dirs[j]) = 1;
                     is_upper[j] = true;
                  } else {
                     border.lower(dirs[j]) = box.lower(dirs[j]);
                     border.upper(dirs[j]) = box.lower(dirs[j]);
                     border_shift(dirs[j]) = -1;
                     is_upper[j] = false;
                  }
               }

               border.shift(border_shift);
               BoxList<DIM> border_list(border);

               border_list.removeIntersections(domain);

               if (border_list.size() > 0) {

                  const BoxList<DIM> domain_list(domain);
                  const Box<DIM> domain_bound_box =
                     domain_list.getBoundingBox();

                  if (codim == 1) {

                     IntVector<DIM> new_shift(0);
                     if (is_upper[0]) {
                        new_shift(dirs[0]) =
                           -domain_bound_box.numberCells(dirs[0]);
                     } else {
                        new_shift(dirs[0]) =
                           domain_bound_box.numberCells(dirs[0]);
                     }
                     shifts.addItem(new_shift);

                  } else {

                     bool shift_to_add = true;
                     for (int c = 0; c < codim; c++) {

                        if (is_upper[c]) {
                           if (border.upper(dirs[c]) <=
                               domain_bound_box.upper(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        } else {
                           if (border.lower(dirs[c]) >=
                               domain_bound_box.lower(dirs[c])) {
                              shift_to_add = false;
                              break;
                           }
                        }

                     }

                     if (shift_to_add) {
                        IntVector<DIM> new_shift(0);
                        for (int b = 0; b < codim; b++) {
                           if (is_upper[b]) {
                              new_shift(dirs[b]) =
                                 -domain_bound_box.numberCells(dirs[b]);
                           } else {
                              new_shift(dirs[b]) =
                                 domain_bound_box.numberCells(dirs[b]);
                           }
                        }
                        shifts.addItem(new_shift);
                     }
                  }
               }
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Decompose patch boundary region into pieces depending on spatial dim. *
* Boxes are extended along the boundary to the edge of the ghost layer  *
* if necessary.                                                         *  
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::getBoundaryBoxes(
   tbox::Array< BoundaryBox<DIM> > boundaries[DIM], 
   const Box<DIM>& box, 
   const BoxArray<DIM>& domain_boxes, 
   const IntVector<DIM>& ghosts, 
   const IntVector<DIM> &periodic_shift) const
{

   const Index<DIM> ifirst = box.lower();
   const Index<DIM> ilast  = box.upper();

   int num_per_dirs = 0;
   for (int d = 0; d < DIM; d++) {
      num_per_dirs += (periodic_shift(d) ? 1 : 0);
   }

   if (num_per_dirs == DIM) {

      for (int k = 0; k < DIM; k++) {
         boundaries[k].resizeArray(0);
      }

   } else {
      BoundaryLookupTable<DIM>* blut =
         BoundaryLookupTable<DIM>::getLookupTable();

      const tbox::Array<int>& location_index_max =
         blut->getMaxLocationIndices();
      tbox::Array< BoxList<DIM> > codim_boxlist(DIM);

      for (int d = 0; d < DIM-num_per_dirs; d++) {

         int codim = d+1;
   
         boundaries[d].resizeArray(location_index_max[d]);
         int bdry_array_size = location_index_max[d];
         int num_bboxes = 0; 

         for (int loc = 0; loc < location_index_max[d]; loc++) {
            const tbox::Array<int>& dirs = blut->getDirections(loc, codim);

            tbox::Array<bool> periodic_dir(codim);
            bool all_is_per = true;
            for (int p = 0; p < codim; p++) {
               if (periodic_shift(dirs[p]) == 0) {
                  periodic_dir[p] = false;
                  all_is_per = false;
               } else {
                  periodic_dir[p] = true;
               }
            }

            if (!all_is_per) {
               Box<DIM> border(box);
               IntVector<DIM> border_shift(0);

               for (int i = 0; i < codim; i++) {
                  if (blut->isUpper(loc, codim, i)) {
                     border.lower(dirs[i]) = box.upper(dirs[i]);
                     border.upper(dirs[i]) = box.upper(dirs[i]);
                     border_shift(dirs[i]) = 1;
                  } else {
                     border.lower(dirs[i]) = box.lower(dirs[i]);
                     border.upper(dirs[i]) = box.lower(dirs[i]);
                     border_shift(dirs[i]) = -1;
                  }
               }

               // grow in non-dirs directions
               for (int j = 0; j < DIM; j++) {
                  bool dir_used = false;
                  for (int du = 0; du < codim; du++) {
                     if (dirs[du] == j) {
                        dir_used = true;
                        break;
                     }
                  }
                  if (!dir_used) {
                     border.upper(j) = ilast(j) + ghosts(j);
                     border.lower(j) = ifirst(j) - ghosts(j);
                  }
               }

               BoxList<DIM> domain_list(domain_boxes);
               if (num_per_dirs != 0) {
                  domain_list.grow(periodic_shift);
               }

               /*
                * Intersect border_list with domain, then shift so that
                * true boundary boxes are outside domain.  Then remove
                * intersections with the domain.
                */

               BoxList<DIM> border_list(border);
               border_list.intersectBoxes(domain_list);
 
               border_list.shift(border_shift);

               border_list.removeIntersections(domain_list);

               if (border_list.size() > 0) {
                  for (int bd = 0; bd < d; bd++) {
                     border_list.removeIntersections(codim_boxlist[bd]);

                     if (border_list.size() == 0) {
                        break;
                     }
                  }
               }

               if (border_list.size() > 0) {
                  border_list.coalesceBoxes();
                  for (typename BoxList<DIM>::Iterator bl(border_list);
                       bl; bl++) {
                     if (num_bboxes == bdry_array_size) {
                        boundaries[d].resizeArray(
                           bdry_array_size+location_index_max[d]);
                        bdry_array_size = boundaries[d].size();
                     }

                     BoundaryBox<DIM> boundary_box(bl(), codim, loc);

                     boundaries[d][num_bboxes] = boundary_box;

                     num_bboxes++; 
                  }

                  codim_boxlist[d].unionBoxes(border_list);
               }
            }

            if (loc+1 == location_index_max[d]) {
               boundaries[d].resizeArray(num_bboxes);
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Compute physical domain for index space related to reference domain   *
* by specified ratio.  If any entry of ratio is negative, the reference *
* domain will be coarsened.  Otherwise, it will be refined.             *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::computePhysicalDomain(
   BoxArray<DIM>& domain,
   const IntVector<DIM>& ratio_to_level_zero) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < DIM; i++) {
      TBOX_ASSERT( ratio_to_level_zero(i) != 0 );
   }
   if (DIM > 1) {
      for (i = 0; i < DIM; i++) {
	 TBOX_ASSERT( (ratio_to_level_zero(i)*ratio_to_level_zero((i+1)%DIM) > 0)
		 || (ratio_to_level_zero(i) == 1)
		 || (ratio_to_level_zero((i+1)%DIM) == 1) );
      }
   }
#endif

   domain = d_physical_domain;

   bool coarsen = false;
   IntVector<DIM> tmp_rat = ratio_to_level_zero;
   for (int id = 0; id < DIM; id++) {
      if ( ratio_to_level_zero(id) < 0 ) coarsen = true;
      tmp_rat(id) = abs(ratio_to_level_zero(id)); 
   }

   if ( coarsen ) {
      domain.coarsen(tmp_rat);
   } else {
      domain.refine(tmp_rat);
   }

}

/*
*************************************************************************
*                                                                       *
* Calculates the maximum ghost width for all the variables associated   *
* with the patch descriptor.  This must only be called after all of the *
* variables have been registered with the VariableDatabase.  If a       *
* variable is added that changes the maximum ghost width, then an       *
* assertion failure will result. 
*                                                                       *
*************************************************************************
*/

template<int DIM>
hier::IntVector<DIM> hier::GridGeometry<DIM>::computeMaxGhostWidth(
   const tbox::Pointer<hier::PatchDescriptor<DIM> > descriptor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( !(descriptor.isNull()) );
#endif

   /*
    * If additional variables are added after the Geometry is created and
    * increase the ghost width, an error will result
    */

   hier::IntVector<DIM> ghost_width = descriptor->getMaxGhostWidth();
   
   if ( (hier::IntVector<DIM>::max(d_max_data_ghost_width, ghost_width) 
        != d_max_data_ghost_width) &&
        (d_max_data_ghost_width != hier::IntVector<DIM>(-1)) ) {

      TBOX_ERROR("Error in hier::GridGeometry<DIM> object with name = "
                 << d_object_name << ": in computeMaxGhostWidth():  "
              << "Cannot add variables and increase maximum ghost "
              << "width after creating the GridGeometry!");
   }

   d_max_data_ghost_width = ghost_width;

   return (ghost_width);
}


/*
*************************************************************************
*                                                                       *
* Set physical domain data member from input box array and determine    *
* whether domain is a single box.                                       *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void GridGeometry<DIM>::setPhysicalDomain(const BoxArray<DIM>& domain)
{
   BoxList<DIM> domain_boxes(domain);
   BoxList<DIM> bounding_box(domain_boxes.getBoundingBox());

   bounding_box.removeIntersections(domain_boxes);
   if (bounding_box.getNumberItems() == 0) {
      d_domain_is_single_box = true;
      d_physical_domain = BoxArray<DIM>(domain_boxes.getBoundingBox()); 
   } else {
      d_domain_is_single_box = false;
      d_physical_domain = domain;
   } 
}

/*
*************************************************************************
*                                                                       *
* The argument is an IntVector of length DIM.  It is set to 1          *
* for periodic directions and 0 for all other directions.  In the       *
* periodic directions, the coarse-level shift is calculated and stored  *
* in the IntVector d_periodic_shift. The shift is the number of cells   *
* in each periodic direction and is zero in all other directions.       *
*                                                                       *
*************************************************************************
*/

template<int DIM> void GridGeometry<DIM>::initializePeriodicShift(
   const IntVector<DIM>& directions)
{

   d_periodic_shift = directions;

   int id;
   /*
    * Check incoming array and reset values if necessary.
    */
   for (id = 0; id < DIM; id++) {
      d_periodic_shift(id) = ((d_periodic_shift(id) == 0) ? 0 : 1);
   }
  
   /*
    *  Check if the physical domain is valid for the specified
    *  periodic conditions.  If so, compute the shift in each
    *  dimension based on the the number of cells.
    */
   if ( checkPeriodicValidity( d_physical_domain) ) {

      BoxList<DIM> domain_box_list(d_physical_domain);
      Box<DIM> bounding_box = domain_box_list.getBoundingBox();

      for (id = 0; id < DIM; id++) {
        d_periodic_shift(id) *= bounding_box.numberCells(id);
      }
     
   } else {
      TBOX_ERROR("Error in hier::GridGeometry<DIM> object with name = "
                 << d_object_name << ": in intializePeriodicShift():  "
              << "Domain is not periodic for one (or more) of the dimensions "
              << "specified in the geometry input file!");
   }

}

/*
*************************************************************************
*                                                                       *
* This returns an IntVector of length DIM that is set to the width of  *
* the domain in periodic directions and 0 in all other directions.      *
* the argument contains the refinement ratio relative to the coarsest   *
* level, which is multiplied by d_periodic_shift to get the return      *
* vector.                                                               *
*                                                                       *
*************************************************************************
*/

template<int DIM> IntVector<DIM> GridGeometry<DIM>::getPeriodicShift(
   const IntVector<DIM>& ratio_to_level_zero) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio vector must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   int k;
   for (k = 0; k < DIM; k++) {
      TBOX_ASSERT( ratio_to_level_zero(k) != 0 );
   }
   if (DIM > 1) {
      for (k = 0; k < DIM; k++) {
	 TBOX_ASSERT( (ratio_to_level_zero(k)*ratio_to_level_zero((k+1)%DIM) > 0)
		 || (ratio_to_level_zero(k) == 1)
		 || (ratio_to_level_zero((k+1)%DIM) == 1) );
      }
   }
#endif

   IntVector<DIM> periodic_shift;
   for (int i = 0; i < DIM; i++) {
      if (ratio_to_level_zero(i) > 0) {
         periodic_shift(i) = d_periodic_shift(i)*ratio_to_level_zero(i);
      } else {
        int abs_ratio = abs(ratio_to_level_zero(i));
        periodic_shift(i) = d_periodic_shift(i)/abs_ratio;
      }
   }
   return (periodic_shift);
}


/*
*************************************************************************
*                                                                       *
* This checks if the periodic directions given to the constructor are   *
* valid for the domain.  Periodic directions are valid if the domain    * 
* has exactly two physical boundaries normal to the periodic direction. *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool GridGeometry<DIM>::checkPeriodicValidity(
   const BoxArray<DIM>& domain)
{
   bool is_valid = true;

   IntVector<DIM> valid_direction(1);
   IntVector<DIM> grow_direction(1);

   /*  
    *  Compute the bounding box of a "duplicate" domain + 1 
    *  cell and set the min and max indices of this grown box.
    */
   BoxList<DIM> dup_domain(domain);

   Box<DIM> domain_box = dup_domain.getBoundingBox();
   domain_box.grow(grow_direction);
   int i;
   Index<DIM> min_index(0), max_index(0);
   for (i = 0; i < DIM; i++) {
      //set min/max of the bounding box
      min_index(i) = domain_box.lower(i);
      max_index(i) = domain_box.upper(i);
   }

   /*
    *  Next, for each dimension, grow another "duplicate" domain 
    *  by 1.  Remove the intersections with the original domain,  
    *  and loop through the remaining box list, checking if the 
    *  upper index of the box matches the bounding box max or the 
    *  lower index of the box matches the bounding box min.  If 
    *  not, this dimension is not a valid periodic dimension.
    */
   for (i = 0; i < DIM; i++) {
      BoxList<DIM> dup_domain2(domain);
      IntVector<DIM> grow_one(0);
      grow_one(i) = 1;
      dup_domain2.grow(grow_one);
      dup_domain2.removeIntersections(domain);

      typename BoxList<DIM>::Iterator n;
      for (n = dup_domain2.listStart(); n; n++) {
         Box<DIM> this_box = n();
         Index<DIM> box_lower = this_box.lower();
         Index<DIM> box_upper = this_box.upper();
         if (d_periodic_shift(i) != 0) {
            if ( !((box_lower(i)==min_index(i)) ||
                   (box_upper(i)==max_index(i))) ) { 
               valid_direction(i) = 0;
            }
         }
      }
   }

   for (i = 0; i < DIM; i++) {
      if ((valid_direction(i)==0) &&
          (d_periodic_shift(i) !=0)) {
         is_valid=false;
     }
   }

   return is_valid;
}

/*
*************************************************************************
*                                                                       *
* Perform an error check on a recently-constructed boundary box to      *
* make sure that it is the proper size, is adjacent to a patch, and is  *
* outside the physical domain.                                          *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool GridGeometry<DIM>::checkBoundaryBox(
   const BoundaryBox<DIM>& boundary_box,
   const Patch<DIM>& patch,
   const BoxArray<DIM>& domain,
   const int num_per_dirs,
   const IntVector<DIM> &max_data_ghost_width ) const
{
   bool return_val = true;

   Box<DIM> bbox = boundary_box.getBox();

   /*
    * Test to see that the box is of size 1 in at least 1 direction.
    */
   IntVector<DIM> box_size;

   for (int i = 0; i < DIM; i++) {
      box_size(i) = bbox.numberCells(i);
   }

   if (box_size.min() != 1) {
      return_val = false;
   }

   /*
    * Quick and dirty test to see that boundary box is adjacent to patch
    * boundary, or a patch boundary extended through the ghost region.
    */
   Box<DIM> patch_box = patch.getBox();

   Box<DIM> grow_patch_box(patch_box);

   grow_patch_box.grow(1);

   if (grow_patch_box != (grow_patch_box+bbox)) {
      bool valid_box = false;
      grow_patch_box = patch_box;
      for (int j = 0; j < DIM; j++) {
         if (num_per_dirs == 0) {

            for (int k = 1; k < DIM; k++) { 

               grow_patch_box.grow((j+k)%DIM,
                                   max_data_ghost_width((j+k)%DIM));

            }
 
         } else {

            for (int k = 1; k < DIM; k++) {

               grow_patch_box.grow((j+k)%DIM,
                                   2*max_data_ghost_width((j+k)%DIM));

            }

         } 
         grow_patch_box.grow(j, 1);
         if (grow_patch_box == (grow_patch_box+bbox)) {
            valid_box = true;
         }
         grow_patch_box = patch_box;
      }
      if (!valid_box) {
         return_val = false;
      }
   }

   /*
    * check that the boundary box is outside the physical domain.
    */
   BoxList<DIM> domain_list(domain);
   BoxList<DIM> bbox_list(bbox);

   domain_list.intersectBoxes(bbox_list);
   
   if (domain_list.getNumberItems()) {
      return_val = false;
   }

   return (return_val);
}

 
/*
*************************************************************************
*                                                                       *
* Print object data to the specified output stream.                     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void hier::GridGeometry<DIM>::printClassData(std::ostream& stream) const
{

   stream << "\nhier::GridGeometry<DIM>::printClassData..." << std::endl;
   stream << "hier::GridGeometry<DIM>: this = "
          << (hier::GridGeometry<DIM>*)this << std::endl;
   stream << "d_object_name = " << d_object_name << std::endl;

   const int n = d_physical_domain.getNumberOfBoxes();
   stream << "Number of boxes describing physical domain = " << n << std::endl;
   stream << "Boxes describing physical domain..." << std::endl;
   d_physical_domain.print(stream);


   stream << "\nd_periodic_shift = " << d_periodic_shift << std::endl;

   stream << "d_max_data_ghost_width = " << d_max_data_ghost_width << std::endl;


}

}
}

#endif
