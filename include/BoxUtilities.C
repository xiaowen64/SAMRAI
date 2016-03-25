//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/hierarchy/boxes/BoxUtilities.C $
// Package:	SAMRAI hierarchy
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Routines for processing boxes within a domain of index space.
//

#ifndef included_hier_BoxUtilities_C
#define included_hier_BoxUtilities_C

#include "BoxUtilities.h"


#include <stdlib.h>

#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"

namespace SAMRAI {
   namespace hier {

/*
*************************************************************************
*                                                                       *
* This static private member function is called by findBadCutPoints(),  *
* and the findBadCutPointsForDirection() member functions.  It sets bad *
* cut points near the lower and upper ends of the border box in the     *
* given coordinate direction.                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::findBadCutPointsForBorderAndDirection(
   const int id,
   tbox::Array<bool>& bad_cuts,
   const Box<DIM>& box,
   const Box<DIM>& border,
   const int bad_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT( (0 <= id) && (id < DIM) );
   TBOX_ASSERT(bad_cuts.getSize() == box.numberCells(id));
   TBOX_ASSERT(bad_interval >= 0);
#endif

   if (bad_interval > 0) {

      const int ilo = box.lower(id);
      const int ihi = box.upper(id);

      int iclo, ichi, ic;

      /*
       * Set bad cut points near lower end of border box.
       */
      int mark = border.lower(id);
      if ( mark > (ilo-bad_interval) ) {

         iclo = 
            tbox::MathUtilities<int>::Max(ilo, (mark-bad_interval+1)) - ilo;
         ichi = 
            tbox::MathUtilities<int>::Min(ihi, (mark-1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

         iclo = 
            tbox::MathUtilities<int>::Max(ilo, (mark+1)) - ilo;
         ichi = 
            tbox::MathUtilities<int>::Min(ihi, (mark+bad_interval-1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

      }

      /*
       * Set bad cut points near upper end of border box.
       */
      mark = border.upper(id) + 1;
      if ( mark < (ihi+bad_interval+1) ) {

         iclo = 
            tbox::MathUtilities<int>::Max(ilo, (mark-bad_interval+1)) - ilo;
         ichi = 
            tbox::MathUtilities<int>::Min(ihi, (mark-1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

         iclo = 
            tbox::MathUtilities<int>::Max(ilo, (mark+1)) - ilo;
         ichi = 
            tbox::MathUtilities<int>::Min(ihi, (mark+bad_interval-1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

      }

   }
}

/*
*************************************************************************
*                                                                       *
* Check min size, cut factor, and physical domain constraints for       *
* given box.  If a patch is generated from a box that violates any      *
* of these constraints, then some other routine (e.g., ghost cell       *
* filling, or inter-patch communication) may fail.  Thus, an error      *
* message will be generated describing the violation and the program    *
* will abort.                                                           *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::checkBoxConstraints(
   const Box<DIM>& box,
   const IntVector<DIM>& min_size,
   const IntVector<DIM>& cut_factor,
   const IntVector<DIM>& bad_interval,
   const BoxArray<DIM>& physical_boxes)
{
   int id;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_size > IntVector<DIM>(0));
   TBOX_ASSERT(cut_factor > IntVector<DIM>(0));
   TBOX_ASSERT(bad_interval >= IntVector<DIM>(0)); 
#endif
  
   /*
    * Test box against minimum size constraint.
    */
   tbox::Array<bool> min_is_bad(DIM);
   bool min_violation = false;
   for (id = 0; id < DIM; id++) {
      if ( box.numberCells(id) < min_size(id) ) {
         min_is_bad[id] = true;
         min_violation = true;
      } else {
         min_is_bad[id] = false;
      }
   }

   if (min_violation) {
      tbox::perr << "\nBox = " << box << " -- minimum size = " << min_size << std::endl;
      for (id = 0; id < DIM; id++) { 
         if (min_is_bad[id]) {
            tbox::perr << "min size violated in direction " << id << std::endl;
         }
      }
      TBOX_ERROR("BoxUtilities<DIM>::checkBoxConstraints() error:\n"
                 << "  Box violates minimum size restriction" << std::endl);
   }

   /*
    * Test box against cut factor constraint.
    */
   tbox::Array<bool> factor_is_bad(DIM);
   bool factor_violation = false;
   for (id = 0; id < DIM; id++) {
      if ( (box.numberCells(id) % cut_factor(id)) != 0 ) {
         factor_is_bad[id] = true;
         factor_violation = true;
      } else {
         factor_is_bad[id] = false;
      } 
   }

   if (factor_violation) {
      tbox::perr << "\nBox = " << box << " -- cut factor = " << cut_factor << std::endl;
      for (id = 0; id < DIM; id++) { 
         if (factor_is_bad[id]) {
            tbox::perr << "factor bad in direction " << id << std::endl;
         }
      }
      TBOX_ERROR("BoxUtilities<DIM>::checkBoxConstraints() error:\n"
                 << "  Box violates cut factor restriction" << std::endl);
   }

   if (physical_boxes.getNumberOfBoxes() > 0) {

      tbox::Array<bool> cut_is_bad(DIM);
      for (id = 0; id < DIM; id++) { cut_is_bad[id] = false; }

      bool bad_cut_violation = false;

      /*
       * Test box for bad cut point violation.
       */

      Box<DIM> test_border = box;
      test_border.grow(bad_interval);

      BoxList<DIM> border_boxes(test_border);
      border_boxes.removeIntersections(BoxList<DIM>(physical_boxes));

      if (!border_boxes.isEmpty()) {

         /*
          * Test individual box faces in each direction for bad cuts.
          */

         id = 0;
         while ( (id < DIM) && !bad_cut_violation ) {

            int blo = box.lower(id);
            int bhi = box.upper(id);
            int bad = bad_interval(id);

            /*
             * Test lower box face in single direction.
             */

            Box<DIM> test_box = box;
            test_box.grow(bad_interval);

            test_box.upper(id) = box.lower(id) - 1;

            BoxList<DIM> test_boxes(test_box);
            test_boxes.intersectBoxes(border_boxes);
            test_boxes.simplifyBoxes();

            typename BoxList<DIM>::Iterator tb = test_boxes.listStart();
            while (!bad_cut_violation && tb) {
               if ( (tb().lower(id) > (blo-bad))
                  || (tb().upper(id) < (blo-1)) ) {
                  bad_cut_violation = true;
                  cut_is_bad[id] = true;
               }
               tb++;
            }

            if (!bad_cut_violation) {
   
               /*
                * Test upper box face in single direction.
                */ 
   
               test_box = box;
               test_box.grow(bad_interval);
            
               test_box.lower(id) = box.upper(id) + 1;
      
               test_boxes = BoxList<DIM>(test_box);
               test_boxes.intersectBoxes(border_boxes);
               test_boxes.simplifyBoxes();
   
               tb = test_boxes.listStart();
               while (!bad_cut_violation && tb) {
                  if ( (tb().lower(id) > (bhi+1))
                     || (tb().upper(id) < (bhi+bad)) ) {
                     bad_cut_violation = true;
                     cut_is_bad[id] = true;
                  }
                  tb++;
               }
         
            }
   
            id++;
         }
   
      }

      if (bad_cut_violation) {
         
         tbox::perr << "Box violates bad cut restriction in directions...";
         for (id = 0; id < DIM; id++) {
            if (cut_is_bad[id]) tbox::perr << "\n" << id;
         }
         tbox::perr << "\nBox = " << box << " -- bad cut interval = " 
                                 << bad_interval << std::endl;
         tbox::perr << "Physical domain boxes ... " << std::endl;
         for (int ib = 0; ib < physical_boxes.getNumberOfBoxes(); ib++) {
            tbox::perr << "Box # " << ib << " -- " 
                             << physical_boxes.getBox(ib) << std::endl;
         }
         TBOX_ERROR("BoxUtilities<DIM>::checkBoxConstraints() error:\n"
                    << "  Box violates bad cut restriction" << std::endl);
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Replace each box in the list that is too large with a list of         *
* nonoverlapping smaller boxes whose union covers the same region of    *
* index space as the original box.  The resulting boxes will obey the   *
* minimum size, and cut factor restrictions if the original box does.   *
* However, the maximum size restriction may be sacrified if the box     *
* cannot be chopped at appropriate points.                              *
*                                                                       *
* For each box in the list, we perform the following operations         *
*                                                                       *
*    (1) Determine a set of cut points for each coordinate direction.   *
*        The ideal cuts satisfy all min, max, and factor restrictions   *
*        assuming the box does too.                                     * 
*                                                                       *
*    (2) If step (1) finds that the box may be chopped, we determine    *
*        the bad cut points for the box and adjust the original cut     *
*        points if necessary.  Note that this operation uses the        *
*        physical domain and the bad interval information.              *
*                                                                       *
*    (3) The box is chopped if this is still possible after (1) and (2).*
*                                                                       *
*    (4) If the box is chopped, set the box list to the resulting       *
*        boxes.  Otherwise, put the original box on the list.           * 
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::chopBoxes(
   BoxList<DIM>& boxes,
   const IntVector<DIM>& max_size,
   const IntVector<DIM>& min_size,
   const IntVector<DIM>& cut_factor,
   const IntVector<DIM>& bad_interval,
   const BoxArray<DIM>& physical_boxes)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_size > IntVector<DIM>(0));
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > IntVector<DIM>(0));
   TBOX_ASSERT(bad_interval >= IntVector<DIM>(0));
   TBOX_ASSERT(physical_boxes.getNumberOfBoxes() > 0);
#endif

   BoxList<DIM> in_boxes(boxes);
   boxes.clearItems();

   while ( !in_boxes.isEmpty() ) {

      Box<DIM> box = in_boxes.getFirstItem();
      in_boxes.removeFirstItem();

      BoxList<DIM> tmp_boxes;

      tbox::Array< tbox::List<int> > cut_points(DIM);
      bool chop_box = findBestCutPointsGivenMax(cut_points, 
                                                box, 
                                                max_size, 
                                                min_size, 
                                                cut_factor);

      if (chop_box) {

         for (int id = 0; id < DIM; id++) {
            
            if (cut_points[id].getNumberItems() > 0) {

               tbox::Array<bool> bad_cut_points;

               findBadCutPointsForDirection(id,
                                            bad_cut_points,
                                            box,
                                            physical_boxes,
                                            bad_interval);
               fixBadCutPointsForDirection(id,
                                           cut_points[id],
                                           bad_cut_points,
                                           box,
                                           min_size(id),
                                           cut_factor(id));  

            }

         }

         chopBox(tmp_boxes, 
                 box,
                 cut_points);

         boxes.catenateItems(tmp_boxes);

      } else {

         boxes.appendItem(box);

      }
  
   }

}

/*
*************************************************************************
*                                                                       *
* Chop given box into a collection of boxes according to the collection *
* of cut points specified along each coordinate direction.   This box   *
* list is formed from the resulting boxes.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::chopBox(
   BoxList<DIM>& boxes,
   const Box<DIM>& box,
   const tbox::Array< tbox::List<int> > cut_points)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(cut_points.getSize() == DIM);
#endif

   if (!box.empty()) {

      boxes.clearItems();
      boxes.appendItem(box);

      BoxList<DIM> tmp_boxes;
      for (int id = 0; id < DIM; id++) {

         tmp_boxes.clearItems();

         while ( !boxes.isEmpty() ) {

            Box<DIM> chop_box = boxes.getFirstItem();
            boxes.removeFirstItem();

            if (cut_points[id].getNumberItems() > 0) {
   
               Index<DIM> ilo = chop_box.lower();
               Index<DIM> ihi = chop_box.upper();
               Index<DIM> boxhi = chop_box.upper();

               tbox::List<int>::Iterator cut = cut_points[id].listStart();
#ifdef DEBUG_CHECK_ASSERTIONS
               int last_cut = tbox::MathUtilities<int>::getMin();
#endif
               while (cut) {
                  int cut_val = cut();
#ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(last_cut <= cut_val);
                  last_cut = cut_val;
#endif
                  ihi(id) = cut_val - 1;
                  if ( (ilo(id) < cut_val) && (ihi(id) <= boxhi(id)) ) {
                     Box<DIM> new_box(ilo, ihi);
                     tmp_boxes.appendItem(new_box);
                     ilo(id) = cut_val;
                  }
                  cut++;
               }
   
               ihi(id) = chop_box.upper(id);
               Box<DIM> last_box(ilo, ihi);
               tmp_boxes.appendItem(last_box);
   
            } else {
                tmp_boxes.appendItem(chop_box);
            }

         }

         boxes = tmp_boxes;

      }

   }
}

/*
*************************************************************************
*                                                                       *
* Test each box in this box list for its intersection with the physical *
* domain boundary when it is grown by the given ghost width.  If the    *
* ghost box lies entirely within the domain, or if all of its ghost     *
* cells intersect the domain boundary appropriately, then the box will  *
* not be changed.  Otherwise, the box is removed from the list and is   *
* replaced by a new box formed by growing the original box to boundary. *
* This process eliminates domain boundary intersections which are       *
* deemed unacceptable.  Intersections that are disallowed are those in  *
* which a portion of the domain boundary is parallel to a box face and  *
* lies strictly in the interior of the ghost cell layer adjacent to     *
* that face.  In other words, we eliminate ghost cell regions residing  *
* outside of the domain and which are narrower than the ghost width.    *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::extendBoxesToDomainBoundary(
   BoxList<DIM>& boxes,
   const BoxList<DIM>& domain,
   const IntVector<DIM>& ext_ghosts)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!domain.isEmpty());
   TBOX_ASSERT(ext_ghosts >= IntVector<DIM>(0));
#endif

   bool out_val = false;

   BoxList<DIM> out_boxes;

   while ( !boxes.isEmpty() ) {

      Box<DIM> try_box = boxes.getFirstItem();
      boxes.removeFirstItem();

      out_val = (extendBoxToDomainBoundary(try_box, domain, ext_ghosts) ||
                 out_val);

      out_boxes.appendItem(try_box);

   }

   boxes = out_boxes;

   return(out_val);
}

template<int DIM> bool BoxUtilities<DIM>::extendBoxToDomainBoundary(
   Box<DIM>& box,
   const BoxList<DIM>& domain,
   const IntVector<DIM>& ext_ghosts)
{
   int id;
#ifdef DEBUG_CHECK_ASSERTIONSa
   TBOX_ASSERT(!domain.isEmpty());
   TBOX_ASSERT(ext_ghosts >= IntVector<DIM>(0));
#endif

   bool out_val = false;

   if (!box.empty()) {

      Box<DIM> test_ghost_box = box;
      test_ghost_box.grow(ext_ghosts);

      BoxList<DIM> outside_domain(test_ghost_box);
      outside_domain.removeIntersections(domain);

      if (!outside_domain.isEmpty()) {

         for (id = 0; id < DIM; id++) {
            BoxList<DIM> outside_boxes;
            typename BoxList<DIM>::Iterator lb;

            // Test whether lower end of ghost box extends outside domain
            Box<DIM> test_region = test_ghost_box;
            test_region.upper(id) = box.lower(id)-1;

            outside_boxes = outside_domain;
            outside_boxes.intersectBoxes(test_region);
 
            int box_lo = box.lower(id);
            for (lb=outside_boxes.listStart(); lb; lb++) {
               box_lo = tbox::MathUtilities<int>::Min(box_lo, lb().upper(id)+1);
            }

            // Test whether upper end of ghost box extends outside domain
            test_region = test_ghost_box;
            test_region.lower(id) = box.upper(id)+1;

            outside_boxes = outside_domain;
            outside_boxes.intersectBoxes(test_region);

            int box_hi = box.upper(id);
            for (lb=outside_boxes.listStart(); lb; lb++) {
               box_hi = tbox::MathUtilities<int>::Max(box_hi, lb().lower(id)-1);
            }

            if (!out_val) {
               out_val = ( (box.lower(id) != box_lo) ||
                           (box.upper(id) != box_hi) ); 
            }

            // Adjust box dimensions as necessary
            box.lower(id) = box_lo;
            box.upper(id) = box_hi;

         }

      }

   }

   return(out_val);
}

/*
*************************************************************************
*                                                                       *
* Grow each box in the list that is smaller than the specified minimum  *
* size.  Each box that is grown must remain within the union of the     *
* boxes of the given domain.  If the specified domain is an empty box   *
* list, then each box will be grown to be as large as the minimum size  *
* with no particular restrictions applied.  Note that this operation    *
* may produce overlap regions among boxes on the list in either case.   *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::growBoxesWithinDomain(
   BoxList<DIM>& boxes,
   const BoxList<DIM>& domain,
   const IntVector<DIM>& min_size)
{
   int id;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_size > IntVector<DIM>(0)); 
#endif

   if (!boxes.isEmpty()) {

      BoxList<DIM> out_boxes;

      BoxList<DIM> outside_domain;
      if (domain.isEmpty()) {
         Box<DIM> big_box(boxes.getBoundingBox());
         big_box.grow(min_size);
         outside_domain = BoxList<DIM>(big_box);
         outside_domain.grow(IntVector<DIM>(1));
         outside_domain.removeIntersections(big_box);
      } else {
         outside_domain = domain;
         outside_domain.grow(IntVector<DIM>(1));
         outside_domain.removeIntersections(domain);
      }

      while ( !boxes.isEmpty() ) {

         Box<DIM> try_box = boxes.getFirstItem();
         boxes.removeFirstItem();

         for (id = 0; id < DIM; id++) {

            int grow = min_size(id) - try_box.numberCells(id);
 
            if (grow > 0) {
   
               BoxList<DIM> outside_boxes;
               typename BoxList<DIM>::Iterator lb;
               Box<DIM> test_region;

               // How far may box be grown within domain in lower direction?
               test_region = try_box;
               test_region.lower(id) -= grow;
               test_region.upper(id) = try_box.lower(id) - 1;
 
               outside_boxes = outside_domain;
               outside_boxes.intersectBoxes(test_region);
               
               int grow_lo = try_box.lower(id) - grow;
               for (lb=outside_boxes.listStart(); lb; lb++) {
                  grow_lo = 
                    tbox::MathUtilities<int>::Max(grow_lo, lb().upper(id)+1);
               }
 
               // How far may box be grown within domain in upper direction?
               test_region = try_box;
               test_region.upper(id) += grow;
               test_region.lower(id) = try_box.upper(id) + 1;
 
               outside_boxes = outside_domain;
               outside_boxes.intersectBoxes(test_region);
    
               int grow_up = try_box.upper(id) + grow;
               for (lb=outside_boxes.listStart(); lb; lb++) {
                  grow_up = 
                     tbox::MathUtilities<int>::Min(grow_up, lb().lower(id)-1);
               }
    
               // Adjust box dimensions as necessary
               if ( (grow_up - grow_lo + 1) < min_size(id) ) {
                  try_box.lower(id) = grow_lo;
                  try_box.upper(id) = grow_up;
               } else {
                  int left = try_box.lower(id) - grow_lo;
                  int right = grow_up - try_box.upper(id);
                  int grow_half = grow/2;
    
                  if (left < right) {
                     try_box.lower(id) -= ( (left < grow_half) ? left 
                                                               : grow_half );
                     try_box.upper(id) = try_box.lower(id) + min_size(id) - 1;
                  } else {
                     try_box.upper(id) += ( (right < grow_half) ? right 
                                                                : grow_half );
                     try_box.lower(id) = try_box.upper(id) - min_size(id) + 1;  
                  }
               }

             }

         }

         out_boxes.appendItem(try_box);

      }

      boxes = out_boxes;

   }
}

/*
*************************************************************************
*                                                                       *
* Determine whether this box can be chopped according to specified      *
* max, min, and factor constraints.  If the box may be chopped along    *
* any face, true is returned.  Otherwise, false is returned.  For those *
* directions along which the box may be chopped, the cut points are     *
* computed.  The procedure is as follows:                               *
*                                                                       *
*    (1) Determine which directions chopping is allowed.                *
*    (2) For each direction to chop, determine list of cut points.      *
*                                                                       *
* Important note: By convention, each integer cut point that is         *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::findBestCutPointsGivenMax(
   tbox::Array< tbox::List<int> >& cut_points,
   const Box<DIM>& box,
   const IntVector<DIM>& max_size,
   const IntVector<DIM>& min_size,
   const IntVector<DIM>& cut_factor)
{
   int id;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_size > IntVector<DIM>(0));
   TBOX_ASSERT(min_size <= max_size);
   TBOX_ASSERT(cut_factor > IntVector<DIM>(0));
#endif

   bool chop_ok = false;

   cut_points.resizeArray(DIM);

   for (id = 0; id < DIM; id++) {
      if (findBestCutPointsForDirectionGivenMax(id,
                                                cut_points[id],
                                                box,
                                                max_size(id),
                                                min_size(id),
                                                cut_factor(id))) {
         chop_ok = true;
      }
   }

   return(chop_ok);

}

/*
*************************************************************************
*                                                                       *
* Determine whether this box can be chopped according to specified      *
* max, min, and factor constraints along given coordinate direction.    *
* If the box may be chopped, true is returned; otherwise, false is      *
* returned.  The procedure for determining the cuts is as follows:      *
*                                                                       *
*    (1) Adjust min and max values so that they are integer             *
*        multiples of the cut factor.                                   *
*    (2) Determine number of boxes, min and max box widths.             *
*    (3) Determine list of cut points.                                  *
*                                                                       *
* Important note: By convention, each integer cut point that is         *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::findBestCutPointsForDirectionGivenMax(
   const int idir,
   tbox::List<int>& cut_points,
   const Box<DIM>& box, 
   const int max_size,
   const int min_size,
   const int cut_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > 0);
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > 0);
#endif

   cut_points.clearItems();

   bool chop_ok = ( ((box.numberCells(idir) % cut_factor)
                    || (box.numberCells(idir) <= max_size)
                    || (box.numberCells(idir) < 2*min_size))
                    ? false : true );

   if (chop_ok) {

      int min = min_size; 
      int max = max_size; 

      chop_ok = false;

      int len = box.numberCells(idir);

      if (min % cut_factor) min = (min/cut_factor+1)*cut_factor;
      if (max % cut_factor) max = (max/cut_factor)*cut_factor;

      /* make sure that max >= min.  In the case that
       * max equals min, max is increased only if the len is
       * not divisible by max.  This choice ensures that we
       * choose cut points that satisfy the min constraint
       * but possibly at the expense of breaking the max constraint.
       */
      if ( (max < min) || ((max == min) && ((len % max) != 0)) ) {
         max = tbox::MathUtilities<int>::Min(2*min, len/2);
      }

      int num_boxes = 1;
      int max_width = min;
      int num_wide_boxes = num_boxes;
      int min_width = min;

      num_boxes = (len-1)/max + 1;
      int len_remaining = len - num_boxes * min;

      if (len_remaining > 0) {
         int len_mult = len_remaining/cut_factor;
         num_wide_boxes = len_mult % num_boxes;
         if ( num_wide_boxes != 0 ) {
            max_width += (len_mult/num_boxes + 1) * cut_factor;
            min_width = max_width - cut_factor;
         } else {
            max_width += (len_mult/num_boxes) * cut_factor;
            num_wide_boxes = num_boxes;
            min_width = 0;
         }
      }

      if (num_boxes > 1) {
         int mark = box.lower(idir);
         int wide_count = 0;
         for (int ic = 0; ic < num_boxes-1; ic++) {
            int width = ( (wide_count < num_wide_boxes)
                          ? max_width : min_width );
            mark += width;
            cut_points.appendItem(mark);
            wide_count++;
         }

         chop_ok = true;
      }

   }

   return(chop_ok);

}

/*
*************************************************************************
*                                                                       *
* Determine whether this box may be chopped according to requested      *
* number of cuts along each side.  If the box may be chopped along any  *
* coordinate direction, true is returned.  Otherwise, false is          *
* returned.  For those directions along which the box may be chopped,   *
* the cut points are computed.  The procedure is as follows:            *
*                                                                       *
*    (1) Determine for which directions shopping is allowed.            *
*    (2) For each direction to chop, determine list of cut points.      *
*                                                                       *
* Important note: By convention, each integer cut point that is         *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::findBestCutPointsGivenNumber(
   tbox::Array< tbox::List<int> >& cut_points,
   const Box<DIM>& box,
   const IntVector<DIM>& number_boxes,
   const IntVector<DIM>& min_size,
   const IntVector<DIM>& cut_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > IntVector<DIM>(0));
   TBOX_ASSERT(number_boxes > IntVector<DIM>(0));
   TBOX_ASSERT(cut_factor > IntVector<DIM>(0));
#endif

   int id;

   cut_points.resizeArray(DIM);

   tbox::Array<bool> chop_dir(DIM);
   for (id = 0; id < DIM; id++) {
      cut_points[id].clearItems();
      chop_dir[id] = ( ((number_boxes(id) <= 1)
                       || (box.numberCells(id) % cut_factor(id))
                       || (box.numberCells(id) < 2*min_size(id))
                       || (box.numberCells(id) < 
                           (number_boxes(id)*min_size(id))))
                        ? false : true );
   }

   bool chop_ok = false;

   for (id = 0; id < DIM; id++) {

      if (chop_dir[id]) {

         if (findBestCutPointsForDirectionGivenNumber(id,
                                                      cut_points[id],
                                                      box,
                                                      number_boxes(id),
                                                      min_size(id),
                                                      cut_factor(id))) {
            chop_ok = true;
         }

      }

   }

   return(chop_ok);

}

/*
*************************************************************************
*                                                                       *
* Determine whether this box may be chopped according to requested      *
* number of cuts along given direction.  If the box may be chopped,     *
* true is returned; otherwise, false is returned.  The procedure for    *
* determining the cuts is as follows:                                   *
*                                                                       *
*    (1) Adjust min value so that it is an integer multiple of          *
*        the cut factor.                                                *
*    (2) Determine number of boxes, min and max box widths.             *
*    (3) Determine list of cut points.                                  *
*                                                                       *
* Important note: By convention, each integer cut point that is         *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::findBestCutPointsForDirectionGivenNumber(
   const int idir,
   tbox::List<int>& cut_points,
   const Box<DIM>& box,
   const int num_boxes,
   const int min_size,
   const int cut_factor)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(min_size > 0);
   TBOX_ASSERT(num_boxes > 0);
   TBOX_ASSERT(cut_factor > 0);
#endif

   cut_points.clearItems();

   bool chop_ok = ( ((num_boxes <= 1)
                    || (box.numberCells(idir) % cut_factor)
                    || (box.numberCells(idir) < 2*min_size)
                    || (box.numberCells(idir) < num_boxes*min_size))
                    ? false : true );

   if (chop_ok) {

      chop_ok = false;

      int len = box.numberCells(idir);
      int min = min_size;

      if (min % cut_factor) min = (min/cut_factor+1)*cut_factor;

      int max_width = min;
      int num_wide_boxes = num_boxes;
      int min_width = min;

      int len_remaining = len - num_boxes * min;

      if (len_remaining > 0) {
         int len_mult = len_remaining/cut_factor;
         num_wide_boxes = len_mult % num_boxes;
         if ( num_wide_boxes != 0 ) {
            max_width += (len_mult/num_boxes + 1) * cut_factor;
            min_width = max_width - cut_factor;
         } else {
            max_width += (len_mult/num_boxes) * cut_factor;
            num_wide_boxes = num_boxes;
            min_width = 0;
         }
      }

      if (num_boxes > 1) {
         int mark = box.lower(idir);
         int wide_count = 0;
         for (int ic = 0; ic < num_boxes-1; ic++) {
            int width = ( (wide_count < num_wide_boxes)
                          ? max_width : min_width );
            mark += width;
            cut_points.appendItem(mark);
            wide_count++;
         }

         chop_ok = true;
      }

   }

   return(chop_ok);

}

/*
*************************************************************************
*                                                                       *
* Return true if the box may have bad cut points, potentially.          *
* Otherwise, return false.  Information about which directions may      *
* have bad cut points is returned in the integer vector.  An entry of   *
* zero indicates that there are no bad cut points for the box along     *
* that coordinate direction.  An entry of one indicates that there      *
* may be a bad cut point along that direction. 
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::checkBoxForBadCutPoints(
   IntVector<DIM>& bad_cut_information,
   const Box<DIM>& box,
   const BoxArray<DIM>& physical_boxes,
   const IntVector<DIM>& bad_interval)
{
   bool found_bad = false;

   int id;
   
   bad_cut_information = IntVector<DIM>(0);
   for (id = 0; id < DIM; id++) {
      if ( checkBoxForBadCutPointsInDirection(id,
                                              box,
                                              physical_boxes,
                                              bad_interval) ) {
         bad_cut_information(id) = 1; 
         found_bad = true; 
      }
   }

   return(found_bad);
}

/*
*************************************************************************
*                                                                       *
* Return true if the box may have bad cut points along the given        *
* coordinate direction, potentially.  Otherwise, return false.
*                                                                       *
*************************************************************************
*/

template<int DIM> bool BoxUtilities<DIM>::checkBoxForBadCutPointsInDirection(
   const int id,
   const Box<DIM>& box,
   const BoxArray<DIM>& physical_boxes,
   const IntVector<DIM>& bad_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_interval >= IntVector<DIM>(0));
#endif

   bool found_bad = false;

   if (physical_boxes.getNumberOfBoxes() > 0) {

      BoxList<DIM> level_interior(physical_boxes);

      int bad = bad_interval(id);

      int id2 = 0;
      while ( (id2 < DIM) && !found_bad ) {
         if (id2 != id) {

            int blo = box.lower(id);
            int bhi = box.upper(id);

            /*
             * Test lower box face in direction id2.
             */

            Box<DIM> border = box;
            border.grow(bad_interval);
            border.upper(id2) = box.lower(id2)-1;

            BoxList<DIM> border_boxes(border);
            border_boxes.removeIntersections(level_interior);
            border_boxes.simplifyBoxes();

            typename BoxList<DIM>::Iterator bb = border_boxes.listStart(); 
            while (!found_bad && bb) {
               found_bad = ( (bb().lower(id) > (blo-bad))
                          || (bb().upper(id) < (bhi+bad)) );
               bb++;
            }

            if (!found_bad) {

               /*
                * Test upper box face in direction id2.
                */

               border = box;
               border.grow(bad_interval);
               border.lower(id2) = box.upper(id2)+1;

               border_boxes.clearItems();
               border_boxes.appendItem(border);
               border_boxes.removeIntersections(level_interior);
               border_boxes.simplifyBoxes();

               bb = border_boxes.listStart();
               while (!found_bad && bb) {
                  found_bad = ( (bb().lower(id) > (blo-bad))
                             || (bb().upper(id) < (bhi+bad)) );
                  bb++;
               }
             
            }
     
         }
         id2++;
      }

   }
   return(found_bad);
}
            
/*
*************************************************************************
*                                                                       *
* Determine bad cut points for box based on the specified physical      *
* domain and bad interval.   The cut information is returned as an      *
* array (size = DIM) of arrays (size = number of cells along edge      *
* of the box) of boolean values.  A value of false indicates a          *
* good cut point, a true value indicates that the cut is bad.           *
*                                                                       *
* Important notes: By convention, each integer cut point that is        *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::findBadCutPoints(
   tbox::Array< tbox::Array<bool> >& bad_cuts,
   const Box<DIM>& box, 
   const BoxArray<DIM>& physical_boxes, 
   const IntVector<DIM>& bad_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_cuts.getSize() == DIM);
#endif

   for (int id = 0; id < DIM; id++) {
      findBadCutPointsForDirection(id,
                                   bad_cuts[id],
                                   box,
                                   physical_boxes,
                                   bad_interval);
   }
}

/*
*************************************************************************
*                                                                       *
* Determine bad cut points for box for given coordinate direction       *
* based on the specified physical domain and bad interval.  The cut     *
* information is returned as an array of integer values (size = number  *
* of cells along edge of box.  A value of zero (0) indicates a good     *
* cut point, a non-zero value indicates that the cut is bad.  The       *
* process works as follows:                                             *
*                                                                       *
*    (1) Initialize all cut points to zero (zero = good).               *
*    (2) Determine bad cut points based on domain configuration.        *
*                                                                       *
* Important notes: By convention, each integer cut point that is        *
* computed corresponds to the cell index to the right of cut point.     *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::findBadCutPointsForDirection(
   const int id,
   tbox::Array<bool>& bad_cuts,
   const Box<DIM>& box,
   const BoxArray<DIM>& physical_boxes,
   const IntVector<DIM>& bad_interval)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_interval >= IntVector<DIM>(0));
#endif

   int ic;

   /*
    * Initialize all bad cut points to false; i.e., all are good.
    */
   const int ncells = box.numberCells(id);
   bad_cuts.resizeArray(ncells);
   for (ic = 0; ic < ncells; ic++) {
      bad_cuts[ic] = false;
   }

   /*
    * Determine whether box intersects physical boundary in such a way
    * that a bad cut point may result when the box is grown by the bad
    * interval.  To determine bad cut points for direction i, the box 
    * must be intersected against the domain exterior in each directions 
    * j not equal to i.  First, we check the lower end of box in direction j.  
    * Then, we check the upper end of the box in direction j.  In each case, 
    * the bad cut points are generated by considering each region that 
    * intersects the domain exterior.
    */

   BoxList<DIM> level_interior(physical_boxes);
   Box<DIM> level_bounding_box = level_interior.getBoundingBox();

   for (int id2 = 0; id2 < DIM; id2++) {

      if ( ((DIM == 1 ) && id2 == id) || ((DIM != 1) && (id2 != id)) ) { 
         /*
          * Test lower box face in direction id2.
          */

         Box<DIM> border = box;
         border.grow(bad_interval);
         border.upper(id2) = box.lower(id2)-1; 

         /*
          * limit the width of the border box to the width of the 
          * domain to ensure that bad cut points near the boundary
          * of the box are not missed.
          */
         border.upper(id) = level_bounding_box.upper(id);
         border.lower(id) = level_bounding_box.lower(id);

         BoxList<DIM> border_boxes(border);

	 if ( DIM > 1) {
	    /* 
	     * only remove the level interior if the dimensionality of 
	     * the problem is greater than 1.
	     */
	    border_boxes.removeIntersections(level_interior);
	 }

         if (!border_boxes.isEmpty()) {
            border_boxes.simplifyBoxes();

            for (typename BoxList<DIM>::Iterator bbox(border_boxes); bbox; bbox++) {
                findBadCutPointsForBorderAndDirection(id,
                                                      bad_cuts,
                                                      box,
                                                      bbox(),
                                                      bad_interval(id));
            }
         }

         /*
          * Test upper box face in direction id2.
          */

         border = box;
         border.grow(bad_interval);
         border.lower(id2) = box.upper(id2)+1;

         /*
          * limit the width of the border box to the width of the 
          * domain to ensure that bad cut points near the boundary
          * of the box are not missed.
          */
         border.upper(id) = level_bounding_box.upper(id);
         border.lower(id) = level_bounding_box.lower(id);

         border_boxes.clearItems();
         border_boxes.appendItem(border);

	 if (DIM > 1) {
	    /* 
	     * only remove the level interior if the dimensionality of 
	     * the problem is greater than 1.
	     */
	    border_boxes.removeIntersections(level_interior);
	 }
      
         if (!border_boxes.isEmpty()) {
            border_boxes.simplifyBoxes();
            for (typename BoxList<DIM>::Iterator bbox(border_boxes); bbox; bbox++) {
                findBadCutPointsForBorderAndDirection(id,
                                                      bad_cuts,
                                                      box,
                                                      bbox(),
                                                      bad_interval(id));
            }
         }

      }
   }
}

/*
*************************************************************************
*                                                                       *
* Adjust cut points if they coincide with bad cut points.               *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::fixBadCutPoints(
   tbox::Array< tbox::List<int> >& cuts,
   const tbox::Array< tbox::Array<bool> >& bad_cuts,
   const Box<DIM>& box,
   const IntVector<DIM>& min_size,
   const IntVector<DIM>& cut_factor)
{
   int id;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(cuts.getSize() == DIM);
   TBOX_ASSERT(bad_cuts.getSize() == DIM);
   bool bad_cuts_ok = true;
   for (id = 0; id < DIM; id++) {
      bad_cuts_ok = bad_cuts_ok && 
                    (bad_cuts[id].getSize() == box.numberCells(id)); 
   }
   TBOX_ASSERT(bad_cuts_ok);
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > IntVector<DIM>(0));
   TBOX_ASSERT(cut_factor > IntVector<DIM>(0));
#endif

   for (id = 0; id < DIM; id++) {
      fixBadCutPointsForDirection(id,
                                  cuts[id],
                                  bad_cuts[id],
                                  box,
                                  min_size(id),
                                  cut_factor(id));
   }
}

/*
*************************************************************************
*                                                                       *
* For specified coordinate direction, adjust cut points if they         *
* coincide with bad cut points.  This routine processes cut points      *
* from the beginning of the list and the end of the list simultaneously.*
* When a bad cut is found when processing from the list beginning,      *
* a good cut point is searched for by moving toward the lower end of    *
* the box.  The opposite holds when processing from list end.  The      *
* In either case, a new cut point will be inserted in the list if one   *
* is found.  Otherwise, there will be one less cut point along the box  *
* side.  This routine may be made more robust in the future.            *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::fixBadCutPointsForDirection(
   const int id,
   tbox::List<int>& cuts, 
   const tbox::Array<bool>& bad_cuts, 
   const Box<DIM>& box,
   const int min_in,
   const int fact)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::List<int>::Iterator cut = cuts.listStart();
   TBOX_ASSERT(bad_cuts.getSize() == box.numberCells(id));
   bool cuts_strictly_increase = true;
   if (cut) {
      int prev = cut();
      cut++;
      while (cut && cuts_strictly_increase) {
         if (cut() <= prev) cuts_strictly_increase = false; 
         prev = cut();
         cut++;
      }
   }
   TBOX_ASSERT(cuts_strictly_increase);
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_in > 0);
   TBOX_ASSERT(fact > 0);
   bool cuts_satisfy_factor = true;
   cut = cuts.listStart();
   while (cut && cuts_satisfy_factor) {
      if (((cut() - box.lower(id))%fact) != 0) cuts_satisfy_factor = false; 
      cut++;
   }
   TBOX_ASSERT(cuts_satisfy_factor);
#endif

   /*
    * Do a quick check to see whether there are any bad cut points for
    * the box in the specified coordinate direction.  If not, we are done.
    */
   bool bad_point_exists = false;
   const int ncells = box.numberCells(id);
   for (int ic = 0; ic < ncells; ic++) {
      if (bad_cuts[ic]) bad_point_exists = true;
   }
   
   if (bad_point_exists) { 

      tbox::List<int>::Iterator cutlo = cuts.listStart();

      if (cutlo) {

         int min = min_in;

         if (min % fact) min = (min/fact+1)*fact;

         const int offset = box.lower(id);
         const int ilo = box.lower(id);
         const int ihi = box.upper(id)+1;

         tbox::List<int>::Iterator cuthi = cuts.listEnd();

         while ( cutlo && cuthi && (cutlo() <= cuthi()) ) {
   
            int bad_cut_val, below, above, try_cut;
            tbox::List<int>::Iterator tmplo;
            tbox::List<int>::Iterator tmphi;
    
            if (cutlo == cuthi) {
   
               if (bad_cuts[cutlo()-offset]) {
   
                  bool found_good_cut = false;
   
                  bad_cut_val = cutlo();
                  tmplo = cutlo;
                  tmphi = cutlo;
                  tmplo--;
                  tmphi++;
                  cuts.removeItem(cutlo);
   
                  below = ( tmplo ? tmplo() : ilo );
                  
                  try_cut = bad_cut_val - fact;
                  while( (try_cut >= (below+min))
                         && bad_cuts[try_cut-offset] ) {
                     try_cut -= fact;
                  }

                  if ( try_cut >= (below+min) ) {
                     found_good_cut = true;           
                     if (tmplo) {
                        cuts.addItemAfter(tmplo, try_cut);
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cuts.addItem(try_cut);
                        cutlo = cuts.listStart();
                     }
                     cutlo++;
                  } else {
                     cutlo = tmphi;
                  }
   
                  if (!found_good_cut) {
                    above = ( tmphi ? tmphi() : ihi );
   
                     try_cut = bad_cut_val + fact;
                     while( (try_cut <= (above-min))
                            && bad_cuts[try_cut-offset] ) {
                        try_cut += fact;
                     }
   
                     if ( try_cut <= (above-min) ) {
                        if (tmphi) {
                           cuts.addItemBefore(tmphi, try_cut);
                           cuthi = tmphi;
                           cuthi--;
                        } else {
                           cuts.appendItem(try_cut);
                           cuthi = cuts.listEnd();
                        }
                        cuthi--;
                     } else {
                        cuthi = tmplo;
                     } 
                  }
   
               } else {
                  cutlo++;
                  cuthi--;
               }
   
            } else {
   
               if (bad_cuts[cutlo()-offset]) {
   
                  bad_cut_val = cutlo(); 
                  tmplo = cutlo;
                  tmplo--;
                  cuts.removeItem(cutlo);
   
                  below = ( tmplo ? tmplo() : ilo );
                  
                  try_cut = bad_cut_val - fact;
                  while( (try_cut >= (below+min)) 
                         && bad_cuts[try_cut-offset] ) {
                     try_cut -= fact;
                  }

                  if ( try_cut >= (below+min) ) {
                     if (tmplo) {
                        cuts.addItemAfter(tmplo, try_cut);
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cuts.addItem(try_cut);
                        cutlo = cuts.listStart();
                     }
                     cutlo++;
                  } else {
                     if (tmplo) {
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cutlo = cuts.listStart();
                     }
                  }
   
               } else {
                  cutlo++;
               } 
      
               if (bad_cuts[cuthi()-offset]) {
   
                  bad_cut_val = cuthi(); 
                  tmphi = cuthi;
                  tmphi++; 
                  cuts.removeItem(cuthi);
   
                  above = ( tmphi ? tmphi() : ihi );
   
                  try_cut = bad_cut_val + fact;
                  while( (try_cut <= (above-min)) 
                         && bad_cuts[try_cut-offset] ) {
                     try_cut += fact;
                  }
   
                  if ( try_cut <= (above-min) ) {
                     if (tmphi) {
                        cuts.addItemBefore(tmphi, try_cut);
                        cuthi = tmphi;
                        cuthi--;
                     } else {
                        cuts.appendItem(try_cut);
                        cuthi = cuts.listEnd();
                     }
                     cuthi--;
                  } else {
                     if (tmphi) {
                        cuthi = tmphi;
                        cuthi--;
                     } else {
                        cuthi = cuts.listEnd();
                     }
                  }
   
               } else {
                  cuthi--;
               }
   
            }
   
         }
   
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Decompose each box in this box array into a list of non overlapping   *
* boxes.  Moreover, the regions of index space formed by composing the  *
* union of boxes on each box list are mutually disjoint.                *
*                                                                       *
*************************************************************************
*/

template<int DIM> void BoxUtilities<DIM>::makeNonOverlappingBoxLists(
   tbox::Array< BoxList<DIM> >& box_list_array,
   const BoxArray<DIM>& boxes)
{
   const int nb = boxes.getNumberOfBoxes();

   for (int i = 0; i < box_list_array.getSize(); i++) {
      box_list_array[i].clearItems();
   }

   box_list_array.resizeArray(nb);

   // Copy boxes into a list to preserve the original box array.
   BoxList<DIM> box_list(boxes);

   // Remove portion of index space represented by array box from list.
   // Keep unique pieces on box list.
   for (int ib = 0; ib < nb; ib++) {
      Box<DIM> remove = boxes.getBox(ib);

      for (typename tbox::List<Box<DIM> >::Iterator l(box_list); l; l++) {
         Box<DIM> intersection = remove * l();
         if (intersection == l()) {
            box_list_array[ib].appendItem(l());
         }
      }
      box_list_array[ib].coalesceBoxes();

      box_list.removeIntersections(remove);
   }
}


}
}

#endif
