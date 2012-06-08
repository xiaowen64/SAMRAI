/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for params, tagging, init for gridding.
 *
 ************************************************************************/

#ifndef included_mesh_TagAndInitializeStrategy_C
#define included_mesh_TagAndInitializeStrategy_C

#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>

namespace SAMRAI {
namespace mesh {

TagAndInitializeStrategy::TagAndInitializeStrategy(
   const tbox::Dimension& dim,
   const std::string& object_name):
   d_dim(dim),
   d_object_name(object_name)
{
}

TagAndInitializeStrategy::~TagAndInitializeStrategy()
{
}

/*
 *************************************************************************
 *
 * Sets refine boxes for case where refine region is specified by the
 * user.  The bool return value specifies whether or not the refine
 * boxes have been reset from the last time the method was called
 * (true = they are reset, false = they have NOT changed).
 *
 * Note that if any method which invokes tagging is performed there
 * is always potential that the boxes have changed so this method will
 * always return true in this case.
 *
 *************************************************************************
 */
bool
TagAndInitializeStrategy::getUserSuppliedRefineBoxes(
   hier::BoxContainer& refine_boxes,
   const int level_num,
   const double time)
{
   TBOX_ASSERT(level_num >= 0);
   TBOX_ASSERT(time >= 0.);

   /*
    * The cycle counter and boolean array specifying whether times
    * are used are initialially set based on inputs. There could be
    * circumstances where a set of refine boxes is requested for a
    * level number greater than the number of entries the user
    * supplied in input.  If this occurs, resize the cycle counter
    * and boolean times arrays and set appropriate defaults to
    * avoid logic errors elsewhere.
    */
   if (level_num >= d_refine_boxes_cycle_counter.getSize()) {

      d_refine_boxes_times.resizeArray(level_num + 1);
      d_refine_boxes_cycles.resizeArray(level_num + 1);
      d_refine_boxes_use_times.resizeArray(level_num + 1);
      d_refine_boxes_cycle_counter.resizeArray(level_num + 1);
      d_refine_boxes_old_seq_num.resizeArray(level_num + 1);

      d_refine_boxes_times[level_num].resizeArray(1);
      d_refine_boxes_cycles[level_num].resizeArray(1);
      d_refine_boxes_times[level_num][0] = 0.;
      d_refine_boxes_cycles[level_num][0] = 0;
      d_refine_boxes_use_times[level_num] = false;
      d_refine_boxes_cycle_counter[level_num] = 0;
      d_refine_boxes_old_seq_num[level_num] = -1;

   }

   /*
    * Increment step counter.
    */
   d_refine_boxes_cycle_counter[level_num]++;

   /*
    * Determine which sequence entry in the refine_boxes box array
    * to use.
    */
   int seq_num = 0;
   if (d_refine_boxes_use_times[level_num]) {

      /*
       * If we are using times, the user has supplied an array
       * times = {time1, time2, ...} that corresponds with the
       * refine box array boxes = {boxarr1, boxarr2, ...}. Pick
       * the appropriate sequence number based on the specified
       * time.
       */
      for (int i = 0; i < d_refine_boxes_times[level_num].getSize();
           i++) {
         if (time > d_refine_boxes_times[level_num][i]) seq_num = i;
      }
   } else {

      /*
       * If we are using steps, the user has supplied an array
       * cycles = {cycle1, cycle2, ...} that corresponds with the
       * refine box array boxes = {boxarr1, boxarr2, ...}.  Pick
       * the appropriate seq number based on the counter.
       */
      for (int i = 0; i < d_refine_boxes_cycles[level_num].getSize();
           i++) {
         if (d_refine_boxes_cycle_counter[level_num] >
             d_refine_boxes_cycles[level_num][i])
            seq_num = i;
      }
   }

   /*
    * Print some warnings if the user-supplied entries will not
    * generate any refined boxes for the level.
    */
   hier::BoxContainer empty_boxes;
   if ((d_refine_boxes.getSize() <= level_num) ||
       (d_refine_boxes[level_num][seq_num].size() == 0)) {

      TBOX_WARNING(
         d_object_name << ": getRefineBoxes\n"
                       << "No refine boxes specified for level "
                       << level_num);
      refine_boxes = empty_boxes;

   } else if (d_refine_boxes[level_num].getSize() <= seq_num) {

      if (d_refine_boxes_use_times[level_num]) {
         TBOX_WARNING(
            d_object_name << ": getRefineBoxes\n"
                          << "No refine boxes specified for time sequence "
                          << seq_num << " on level " << level_num
                          << ".\n No refinement will be performed.");
      } else {
         TBOX_WARNING(
            d_object_name << ": getRefineBoxes\n"
                          << "No refine boxes specified for step sequence "
                          << seq_num << " on level " << level_num
                          << ".\n No refinement will be performed.");
      }
      refine_boxes = empty_boxes;

   } else {

      refine_boxes = d_refine_boxes[level_num][seq_num];

   }

   /*
    * If the user has requested their own particular set of refine
    * boxes (i.e. by calling resetRefineBoxes()), overwrite any previously
    * determined refine boxes with their requested set.
    */
   bool use_reset = false;
   if (d_refine_boxes_reset.getSize() > level_num) {
      if (d_refine_boxes_reset[level_num]) {
         use_reset = true;
         refine_boxes = d_reset_refine_boxes[level_num];
      }
   }

   /*
    * If we have not moved to a new sequence number, or otherwise reset
    * boxes from the last time this method was called, then we return
    * "false", indicating boxes have NOT been reset.  If we have reset
    * the boxes, return "true".
    */
   bool modified_refine_boxes;
   if (!use_reset && d_refine_boxes_old_seq_num[level_num] == seq_num) {
      modified_refine_boxes = false;
   } else {
      d_refine_boxes_old_seq_num[level_num] = seq_num;
      modified_refine_boxes = true;
   }

   /*
    * If one of the tagging methods (e.g. gradient detector or
    * Richardson extrapolation) is used, the boxes may be modified
    * even if they have not been modified by user input.
    */
   if (!refineUserBoxInputOnly()) modified_refine_boxes = true;

   return modified_refine_boxes;

}

/*
 *************************************************************************
 *
 * Resets refine boxes for specified level.
 *
 *************************************************************************
 */

void
TagAndInitializeStrategy::resetRefineBoxes(
   const hier::BoxContainer& refine_boxes,
   const int level_num)
{
   TBOX_ASSERT(level_num >= 0);

   int i = d_reset_refine_boxes.getSize();
   if (i <= level_num) {
      d_reset_refine_boxes.resizeArray(level_num + 1);
      d_refine_boxes_reset.resizeArray(level_num + 1);
      for ( ; i < d_reset_refine_boxes.getSize(); ++i) {
         d_refine_boxes_reset[i] = false;
      }
   }

   d_refine_boxes_reset[level_num] = true;
   d_reset_refine_boxes[level_num] = refine_boxes;

}

}
}
#endif
