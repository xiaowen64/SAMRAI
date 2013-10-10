/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Implementation of TreeLoadBalancer::TransitSet.
 *
 ************************************************************************/

#ifndef included_mesh_BoxTransitSet_C
#define included_mesh_BoxTransitSet_C

#include "SAMRAI/mesh/BoxTransitSet.h"
#include "SAMRAI/mesh/BalanceBoxBreaker.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {



/*
*************************************************************************
*************************************************************************
*/
BoxTransitSet::BoxTransitSet() :
   d_set(),
   d_sumload(0),
   d_pparams(0),
   d_bbb(),
   d_tree_degree(2 /* should be replaced by LocalId generator */),
   d_allow_box_breaking(true),
   d_print_steps(false),
   d_print_pop_steps(false),
   d_print_swap_steps(false),
   d_print_break_steps(false)
{
   setTimers();
}

/*
 *************************************************************************
 *
 * This method adjusts the load in this BoxTransitSet by
 * moving work between it (main_bin) and a holding_bin.  It tries to bring
 * main_bin's load to the specified ideal_load.
 *
 * The high_load and low_load define an acceptable range around the
 * ideal_load.  As soon as the main load falls in this range, no
 * further change is tried, even if it may bring the load closer to
 * the ideal.
 *
 * This method makes a best effort and returns the amount of load
 * moved.  It can move BoxInTransits between given sets and, if needed,
 * break some BoxInTransits up to move part of the work.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 * Return amount of load moved from main_bin to hold_bin.  Negative
 * amount means load moved from hold_bin to main_bin.
 *
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoad(
   BoxTransitSet& hold_bin,
   hier::SequentialLocalIdGenerator &id_generator,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load )
{
   BoxTransitSet& main_bin(*this);

   if (d_print_steps) {
      tbox::plog << "  adjustLoad attempting to bring main load from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " or within [" << low_load << ", " << high_load << "]."
                 << std::endl;
   }
   TBOX_ASSERT( low_load <= ideal_load );
   TBOX_ASSERT( high_load >= ideal_load );


   LoadType actual_transfer = 0;

   if ((main_bin.empty() && ideal_load <= 0 ) ||
       (hold_bin.empty() && main_bin.getSumLoad() < ideal_load )) {
      return actual_transfer;
   }

   t_adjust_load->start();

   actual_transfer = adjustLoadByPopping(
      hold_bin,
      ideal_load,
      low_load,
      high_load );

   if (d_print_steps) {
      double balance_penalty = computeBalancePenalty(
         main_bin,
         hold_bin,
         (main_bin.getSumLoad() - ideal_load));
      tbox::plog << "  Balance penalty after adjustLoadByPopping = "
                 << balance_penalty
                 << ", needs " << (ideal_load-main_bin.getSumLoad())
                 << " more with " << main_bin.size() << " main_bin and "
                 << hold_bin.size() << " hold_bin Boxes remaining."
                 << "\n  main_bin now has " << main_bin.getSumLoad()
                 << " in " << main_bin.size() << " boxes."
                 << std::endl;
   }

   /*
    * The algorithm cycles through a do-loop.  Each time around, we
    * try to swap some BoxTransitSet::BoxInTransit between main_bin and hold_bin
    * until we have main_bin's load in [low_load,high_load] or we
    * cannot improve the actual_transfer any further.  Then, we try
    * breaking up a BoxTransitSet::BoxInTransit to improve the results.  If we break
    * some BoxTransitSet::BoxInTransit, we generate some more swapping options that
    * were not there before, so we loop back to try swapping again.
    *
    * If a break phase does not break any Box (and does not generate
    * more swap options), the loop will stop making changes.  We break
    * the loop at that point (and whenever we get main_bin's load in
    * the correct range).  We also break out if there is no improvement,
    * which can happen when the swp and break steps undo each other's
    * work (due to round-off errors).
    */
   do {

      const LoadType old_distance_to_ideal = ideal_load - main_bin.getSumLoad();

      /*
       * Try to balance load through swapping.
       */
      LoadType swap_transfer = adjustLoadBySwapping(
         hold_bin,
         ideal_load,
         low_load,
         high_load);

      actual_transfer += swap_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            main_bin,
            hold_bin,
            (main_bin.getSumLoad() - ideal_load));
         tbox::plog << "  Balance penalty after adjustLoadBySwapping = "
                    << balance_penalty
                    << ", needs " << (ideal_load-main_bin.getSumLoad())
                    << " more with " << main_bin.size() << " main_bin and "
                    << hold_bin.size() << " hold_bin Boxes remaining."
                    << "\n  main_bin now has " << main_bin.getSumLoad()
                    << " in " << main_bin.size() << " boxes."
                    << std::endl;
      }

      // Skip breaking if already in range.
      if (main_bin.getSumLoad() <= high_load && main_bin.getSumLoad() >= low_load ) break;

      /*
       * Skip breaking if adding/subtracting the min load overshoots the range and worsens distance to range.
       */
      if ( tbox::MathUtilities<double>::Abs(main_bin.getSumLoad() - 0.5*(high_load+low_load)) <= 0.5*d_pparams->getMinLoad() ) {
         break;
      }


      if ( d_allow_box_breaking ) {
         /*
          * Assuming that we did the best we could, swapping
          * some BoxTransitSet::BoxInTransit without breaking any, we now break up a Box
          * in the overloaded side for partial transfer to the
          * underloaded side.
          */
         LoadType brk_transfer = adjustLoadByBreaking(
            hold_bin,
            id_generator,
            ideal_load,
            low_load,
            high_load );
         actual_transfer += brk_transfer;

         if (d_print_steps) {
            double balance_penalty = computeBalancePenalty(
               main_bin,
               hold_bin,
               (main_bin.getSumLoad() - ideal_load));
            tbox::plog << "  Balance penalty after adjustLoadByBreaking = "
                       << balance_penalty
                       << ", needs " << (ideal_load-main_bin.getSumLoad())
                       << " more with " << main_bin.size() << " main_bin and "
                       << hold_bin.size() << " hold_bin Boxes remaining."
                       << "\n  main_bin now has " << main_bin.getSumLoad()
                       << " in " << main_bin.size() << " boxes."
                       << std::endl;
         }
         if (brk_transfer == 0) {
            /*
             * If no box can be broken to improve the actual_transfer,
             * there is nothing further we can do.  The swap phase, tried
             * before the break phase, also generated no transfer, so
             * there's no point trying again.  Break out now to save
             * retrying the swap phase.
             */
            if (d_print_steps) {
               tbox::plog << "  adjustLoad stopping due to unsuccessful break."
                          << std::endl;
            }
            break;
         }
      }


      LoadType improvement =
         tbox::MathUtilities<double>::Abs( old_distance_to_ideal
                                           - (ideal_load - main_bin.getSumLoad()) );
      if ( improvement < d_pparams->getLoadComparisonTol() ) {
         break;
      }

      /*
       * Now that we have broken up a Box, redo this loop to
       * see if swapping can produce a better result.
       */
   } while ( ( main_bin.getSumLoad() >= high_load ) ||
             ( main_bin.getSumLoad() <= low_load ) );

   if ( d_print_steps ) {
      const LoadType point_miss = main_bin.getSumLoad() - ideal_load;
      const LoadType range_miss =
         main_bin.getSumLoad() > high_load ? main_bin.getSumLoad() - high_load :
         main_bin.getSumLoad() < low_load ? low_load - main_bin.getSumLoad() : 0;
      tbox::plog << "  adjustLoad point_miss=" << point_miss
                 << "  range_miss="
                 << (range_miss > 0 ? " ":"") // Add space if missed range
                 << (range_miss > 0.5*d_pparams->getMinBoxSize().getProduct() ? " ":"") // Add space if missed range by a lot
                 << range_miss
                 << "  " << main_bin.getSumLoad() << '/'
                 << ideal_load << " [" << low_load << ',' << high_load << ']'
                 << std::endl;
   }

   t_adjust_load->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 * Attempt bring main_bin to within a specific load range by moving
 * one box to/from it from/to hold_bin.  This method is allowed to break
 * the box and move parts of it.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadByBreaking(
   BoxTransitSet& hold_bin,
   hier::SequentialLocalIdGenerator &id_generator,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load )
{
   LoadType actual_transfer = 0;

   if (getSumLoad() > high_load) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      actual_transfer = -hold_bin.adjustLoadByBreaking(
         *this,
         id_generator,
         hold_bin.getSumLoad()-(ideal_load-getSumLoad()),
         hold_bin.getSumLoad()-(high_load-getSumLoad()),
         hold_bin.getSumLoad()-(low_load-getSumLoad()) );
      return actual_transfer;
   }

   BoxTransitSet& main_bin(*this);

   TBOX_ASSERT(low_load <= ideal_load);
   TBOX_ASSERT(ideal_load <= high_load);
   TBOX_ASSERT(main_bin.getSumLoad() <= high_load);

   TBOX_ASSERT(main_bin.size() + hold_bin.size() > 0);

   t_shift_loads_by_breaking->start();

   const LoadType ideal_transfer = ideal_load - main_bin.getSumLoad();
   const LoadType high_transfer = high_load - main_bin.getSumLoad();
   const LoadType low_transfer = low_load - main_bin.getSumLoad();

   if (d_print_steps) {
      tbox::plog << "    adjustLoadByBreaking asked to break off "
                 << ideal_transfer << " [" << low_transfer << ','
                 << high_transfer << "] from one of " << hold_bin.size()
                 << " Boxes to add to set of " << main_bin.size()
                 << " Boxes."
                 << std::endl;
   }


   // Data for the best cutting results so far:
   std::vector<hier::Box> breakoff;
   std::vector<hier::Box> leftover;
   double breakoff_amt = 0.0;
   BoxTransitSet::BoxInTransit breakbox(d_pparams->getMinBoxSize().getDim());

   int break_acceptance_flags[3] = {0,0,0};
   int &found_breakage = break_acceptance_flags[2];

   /*
    * Find best box to break.  Loop in reverse because smaller boxes
    * are cheaper to analyze for bad cuts.
    */
   for (BoxTransitSet::reverse_iterator si = hold_bin.rbegin(); si != hold_bin.rend(); ++si) {

      /*
       * Skip boxes smaller than ideal_transfer.  If we called
       * adjustLoadBySwapping before entering this method, there
       * should not be any such boxes.
       */
      if ( si->d_boxload < ideal_transfer ) {
         continue;
      }

      const BoxTransitSet::BoxInTransit& candidate = *si;

      if (d_print_steps) {
         tbox::plog << "    Considering break candidate " << candidate
                    << std::endl;
      }

      std::vector<hier::Box> trial_breakoff;
      std::vector<hier::Box> trial_leftover;
      double trial_breakoff_amt;

      d_bbb.breakOffLoad(
         trial_breakoff,
         trial_leftover,
         trial_breakoff_amt,
         candidate.d_box,
         ideal_transfer,
         low_transfer,
         high_transfer );

      if (!trial_breakoff.empty()) {

         const bool accept_break = d_bbb.evaluateBreak(
            break_acceptance_flags, breakoff_amt, trial_breakoff_amt,
            ideal_transfer, low_transfer, high_transfer );
         if (d_print_break_steps) {
            tbox::plog << "      Break evaluation:"
                       << "  " << break_acceptance_flags[0]
                       << "  " << break_acceptance_flags[1]
                       << "  " << break_acceptance_flags[2]
                       << std::endl;
         }

         if (d_print_break_steps) {
            tbox::plog << "    Potential to replace " << candidate << " with "
                       << trial_breakoff.size() << " breakoff Boxes and "
                       << trial_leftover.size() << " leftover Boxes."
                       << "  break amount = " << trial_breakoff_amt
                       << "  in-range imp = " << break_acceptance_flags[0]
                       << "  balance imp = " << break_acceptance_flags[1]
                       << "  overal imp = " << break_acceptance_flags[2]
                       << "  accept_break = " << accept_break
                       << std::endl;
         }

         if (accept_break) {
            breakbox = candidate;
            breakoff_amt = trial_breakoff_amt;
            breakoff.swap(trial_breakoff);
            leftover.swap(trial_leftover);
            if ( break_acceptance_flags[0] == 1 ) {
               // We are in the [low,high] range.  That is sufficient.
               break;
            }
         }

      } else {
         if (d_print_break_steps) {
            tbox::plog << "    Break step could not break " << ideal_transfer
                       << " from hold_bin box " << candidate
                       << std::endl;
         }
      }

   }


   if ( found_breakage == 1 ) {
      /*
       * Remove the chosen candidate.  Put its breakoff parts
       * in main_bin and its leftover parts back into hold_bin.
       */
      hold_bin.erase(breakbox);
      for (std::vector<hier::Box>::const_iterator bi = breakoff.begin();
           bi != breakoff.end();
           ++bi) {
         BoxTransitSet::BoxInTransit give_box_in_transit(
            breakbox,
            *bi,
            breakbox.getOwnerRank(),
            id_generator.nextValue());
         give_box_in_transit.d_boxload = static_cast<int>(computeLoad(
                                                             give_box_in_transit.d_orig_box,
                                                             give_box_in_transit.getBox()));
         main_bin.insert(give_box_in_transit);
         actual_transfer += give_box_in_transit.d_boxload;
         if (d_print_break_steps) {
            tbox::plog << "    Breakoff box " << *bi << bi->numberCells()
                       << '|' << bi->size()
                       << " -> " << give_box_in_transit
                       << std::endl;
         }
      }
      for (std::vector<hier::Box>::const_iterator bi = leftover.begin();
           bi != leftover.end();
           ++bi) {
         BoxTransitSet::BoxInTransit keep_box_in_transit(
            breakbox,
            *bi,
            breakbox.getOwnerRank(),
            id_generator.nextValue());
         keep_box_in_transit.d_boxload = static_cast<int>(computeLoad(
                                                             keep_box_in_transit.d_orig_box,
                                                             keep_box_in_transit.getBox()));
         hold_bin.insert(keep_box_in_transit);
         if (d_print_break_steps) {
            tbox::plog << "    Leftover box " << *bi << bi->numberCells()
                       << '|' << bi->size()
                       << " -> " << keep_box_in_transit
                       << std::endl;
         }
      }
   }

   t_shift_loads_by_breaking->stop();
   return actual_transfer;
}



/*
 *************************************************************************
 * Attempt to adjust the load of a main_bin by swapping boxes with
 * a hold_bin.
 *
 * Transfering a BoxTransitSet::BoxInTransit from one BoxTransitSet to another
 * is considered a degenerate "swap" (a BoxTransitSet::BoxInTransit is
 * swapped for nothing) handled by this function.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of main_bin
 * ideal_transfer < 0 means to raise the load of hold_bin
 * The iterative do loop may overshoot the ideal_transfer
 * and may have to swap to shift some of the load
 * back.
 *
 * Return amount of load transfered.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadBySwapping(
   BoxTransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load )
{
   TBOX_ASSERT( high_load >= ideal_load );
   TBOX_ASSERT( low_load <= ideal_load );

   t_adjust_load_by_swapping->start();

   BoxTransitSet& main_bin(*this);

   if (d_print_steps) {
      tbox::plog << "  Attempting to bring main_bin from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by swapping."
                 << std::endl;
   }

   bool found_swap;

   LoadType actual_transfer = 0;

   do {

      /*
       * Ammount we seek to transfer from hi to lo
       * (the "ideal" for this particular iteration).
       * Unlike ideal_transfer and actual_transfer, this quantity is positive.
       */
      LoadType rem_transfer = main_bin.getSumLoad() - ideal_load;
      LoadType low_transfer = main_bin.getSumLoad() - high_load;
      LoadType high_transfer = main_bin.getSumLoad() - low_load;
      if (d_print_swap_steps) {
         tbox::plog << "    Swap progress: " << main_bin.getSumLoad()
                    << " / " << ideal_load << " remaining transfer = "
                    << rem_transfer << " [" << low_transfer << ','
                    << high_transfer << ']' << std::endl;
      }

      LoadType swap_transfer;
      found_swap = swapLoadPair(
         main_bin,
         hold_bin,
         swap_transfer,
         rem_transfer,
         low_transfer,
         high_transfer);
      swap_transfer = -swap_transfer;

      if (found_swap) {
         actual_transfer += swap_transfer;
      }

   } while (found_swap &&
            (main_bin.getSumLoad() < low_load || main_bin.getSumLoad() > high_load ));

   if (d_print_swap_steps) {
      tbox::plog << "  Final balance for adjustLoadBySwapping: "
                 << main_bin.getSumLoad() << " / " << ideal_load
                 << "  Off by " << (main_bin.getSumLoad()-ideal_load)
                 << std::endl;
   }

   t_adjust_load_by_swapping->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 * Attempt to adjust the load of a main_bin by popping the biggest boxes
 * from a source bin and moving them to a destination bin.
 *
 * This method should give results similar to adjustLoadBySwapping,
 * but when the boxes are small in comparison to the load changed, it
 * should be faster.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of main_bin
 * ideal_transfer < 0 means to raise the load of hold_bin
 *
 * Return amount of load transfered.
 *************************************************************************
 */
BoxTransitSet::LoadType
BoxTransitSet::adjustLoadByPopping(
   BoxTransitSet& hold_bin,
   LoadType ideal_load,
   LoadType low_load,
   LoadType high_load )
{
   TBOX_ASSERT( high_load >= ideal_load );
   TBOX_ASSERT( low_load <= ideal_load );

   t_adjust_load_by_popping->start();

   BoxTransitSet& main_bin(*this);

   /*
    * Logic in this method assumes positive transfer from hold_bin
    * (the source) to main_bin (the destination).  When transfering
    * the other way, switch the roles of main_bin and hold_bin.
    */
   BoxTransitSet *src = &hold_bin;
   BoxTransitSet *dst = &main_bin;
   LoadType dst_ideal_load = ideal_load;
   LoadType dst_low_load = low_load;
   LoadType dst_high_load = high_load;

   if ( main_bin.getSumLoad() > ideal_load ) {

      dst_ideal_load = hold_bin.getSumLoad() + ( main_bin.getSumLoad() - ideal_load );
      dst_low_load = hold_bin.getSumLoad() + ( main_bin.getSumLoad() - high_load );
      dst_high_load = hold_bin.getSumLoad() + ( main_bin.getSumLoad() - low_load );

      src = &main_bin;
      dst = &hold_bin;

   }

   if (d_print_pop_steps) {
      tbox::plog << "  Attempting to bring main_bin from "
                 << main_bin.getSumLoad() << " to " << ideal_load
                 << " [" << low_load << ',' << high_load
                 << "] by popping."
                 << std::endl;
   }

   LoadType actual_transfer = 0;
   int acceptance_flags[3] = {0,0,0};

   size_t num_boxes_popped = 0;

   while ( !src->empty() ) {

      const BoxTransitSet::BoxInTransit &candidate_box = *src->begin();

      bool improved = d_bbb.evaluateBreak(
         acceptance_flags,
         dst->getSumLoad(),
         dst->getSumLoad() + candidate_box.d_boxload,
         dst_ideal_load,
         dst_low_load,
         dst_high_load );

      if ( improved ) {

         if ( d_print_pop_steps ) {
            tbox::plog << "    adjustLoadByPopping pop #" << num_boxes_popped
                       << ", " << candidate_box;
         }

         actual_transfer += candidate_box.d_boxload;
         dst->insert(candidate_box);
         src->erase(src->begin());
         ++num_boxes_popped;

         if ( d_print_pop_steps ) {
            tbox::plog << ", main_bin load is " << main_bin.getSumLoad() << '\n';
         }
      }
      if ( ( dst->getSumLoad() >= dst_low_load && dst->getSumLoad() <= high_load ) ||
           !improved ) {
         break;
      }
   }

   if (d_print_pop_steps) {
      tbox::plog << "  Final result in adjustLoadByPopping: "
                 << main_bin.getSumLoad() << " / " << ideal_load
                 << "  Off by " << (main_bin.getSumLoad()-ideal_load)
                 << ".  " << num_boxes_popped << " boxes popped."
                 << std::endl;
   }

   t_adjust_load_by_popping->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 * Find a BoxTransitSet::BoxInTransit in src and a BoxTransitSet::BoxInTransit in dst which when
 * swapped results in shifting close to ideal_shift from src to dst.
 * Make the swap.  Return whether a swap pair was found.
 *************************************************************************
 */
bool
BoxTransitSet::swapLoadPair(
   BoxTransitSet& src,
   BoxTransitSet& dst,
   LoadType& actual_transfer,
   LoadType ideal_transfer,
   LoadType low_transfer,
   LoadType high_transfer ) const
{
   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      bool rval = swapLoadPair(
         dst,
         src,
         actual_transfer,
         -ideal_transfer,
         -high_transfer,
         -low_transfer);
      actual_transfer = -actual_transfer;
      return rval;
   }

   t_find_swap_pair->start();

   if (d_print_swap_steps) {
      tbox::plog << "    swapLoadPair looking for transfer of "
                 << ideal_transfer
                 << " between " << src.size() << "-box src and "
                 << dst.size() << "-box dst." << std::endl;
      tbox::plog << "      src (" << src.size() << "):" << std::endl;
      if ( src.size() < 10 ) {
         for (BoxTransitSet::iterator si = src.begin(); si != src.end(); ++si) {
            tbox::plog << "        " << *si << std::endl;
         }
      }
      tbox::plog << "      dst (" << dst.size() << "):" << std::endl;
      if ( dst.size() < 10 ) {
         for (BoxTransitSet::iterator si = dst.begin(); si != dst.end(); ++si) {
            tbox::plog << "        " << *si << std::endl;
         }
      }
   }

   /*
    * Look for two swap options.  The "high side" option would
    * transfer at least ideal_transfer.  The "low side" option would
    * transfer up to ideal_transfer.
    *
    * Each option is defined by a box from src and a box from dst,
    * designated by the iterators src_hiside, dst_hiside, src_loside
    * and dst_loside.  src_hiside points to the box in the src for the
    * high-side transfer, and similarly for dst_hiside.  src_loside
    * points to the box in the src for the low-side transfer, and
    * similarly for dst_loside.
    *
    * Note that in the degenerate case, the dst box does not exist,
    * and the swap degenerates to moving a box from the src to the
    * dst.
    *
    * Compute the balance_penalty if high and low were swapped.  Keep
    * looking until we find the pair giving the lowest balance_penalty
    * on swapping.
    *
    * isrc and idst point to the current best pair to swap.  new_balance_penalty
    * is the balance_penalty if we swap them.
    *
    * src_test and dst_test are trial pairs to check to see if we can improve on
    * new_balance_penalty.
    *
    * We will look for two "best" pairs:
    *
    * TODO: This method was originally written to compute the best
    * hiside and loside options separately and compare them at the
    * end.  That separation may not be needded anymore.  It may be
    * possible to simplify this method by keeping only the best option
    * at any time.
    */

   // Initialization indicating no swap pair found yet.
   BoxTransitSet::iterator src_hiside = src.end();
   BoxTransitSet::iterator dst_hiside = dst.end();
   BoxTransitSet::iterator src_loside = src.end();
   BoxTransitSet::iterator dst_loside = dst.end();

   // A dummy BoxTransitSet::BoxInTransit for set searches.
   hier::Box dummy_box(d_pparams->getMinBoxSize().getDim());
   BoxTransitSet::BoxInTransit dummy_search_target(d_pparams->getMinBoxSize().getDim());

   // Difference between swap results and ideal, >= 0
   LoadType hiside_transfer = 0.0;
   LoadType loside_transfer = 0.0;


   int loside_acceptance_flags[3] = {0,0,0};
   int hiside_acceptance_flags[3] = {0,0,0};

   if (dst.empty()) {
      /*
       * There is no dst BoxTransitSet::BoxInTransit, so the swap would
       * degnerate to moving a box from src to dst.  Find
       * the best src BoxTransitSet::BoxInTransit to move.
       */
      dummy_search_target = BoxTransitSet::BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
      dummy_search_target.d_boxload = ideal_transfer;
      const BoxTransitSet::iterator src_test = src.lower_bound(dummy_search_target);

      if (d_print_swap_steps) {
         tbox::plog << "  swapLoadPair with empty dst: ";
      }

      if (src_test != src.begin()) {
         BoxTransitSet::iterator src_test1 = src_test;
         --src_test1;
         if ( d_bbb.evaluateBreak( hiside_acceptance_flags, hiside_transfer, src_test1->d_boxload,
                             ideal_transfer, low_transfer, high_transfer ) ) {
            src_hiside = src_test1;
            hiside_transfer = src_hiside->d_boxload;
            if (d_print_swap_steps) {
               tbox::plog << "  hi src: " << (*src_hiside)
                          << " with transfer " << src_hiside->d_boxload
                          << ", off by " << hiside_transfer-ideal_transfer
                          << ", acceptance_flags=" << hiside_acceptance_flags[0]
                          << ',' << hiside_acceptance_flags[1]
                          << ',' << hiside_acceptance_flags[2];
            }
         }
      }
      if (src_test != src.end()) {
         if ( d_bbb.evaluateBreak( loside_acceptance_flags, loside_transfer, src_test->d_boxload,
                             ideal_transfer, low_transfer, high_transfer ) ) {
            src_loside = src_test;
            loside_transfer = src_loside->d_boxload;
            if (d_print_swap_steps) {
               tbox::plog << "  lo src: " << (*src_loside)
                          << " with transfer " << src_loside->d_boxload
                          << ", off by " << loside_transfer-ideal_transfer
                          << ", acceptance_flags=" << loside_acceptance_flags[0]
                          << ',' << loside_acceptance_flags[1]
                          << ',' << loside_acceptance_flags[2];
            }
         }
      }
      if (d_print_swap_steps) {
         tbox::plog << std::endl;
      }

   } else {

      /*
       * Start search through src beginning with the box whose load
       * exceeds the biggest dst box by at least ideal_transfer.
       */
      dummy_search_target = *dst.begin();
      dummy_search_target.d_boxload += ideal_transfer;
      BoxTransitSet::iterator src_beg = src.lower_bound(dummy_search_target);

      for (BoxTransitSet::iterator src_test = src_beg; src_test != src.end(); ++src_test) {

         /*
          * Set dst_test pointing to where we should start looking in dst.
          * Look for a load less than the load of src_test by
          * ideal_transfer.
          */
         dummy_search_target = BoxTransitSet::BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
         dummy_search_target.d_boxload = tbox::MathUtilities<LoadType>::Max(
               src_test->d_boxload - ideal_transfer,
               0);
         BoxTransitSet::iterator dst_test = dst.lower_bound(dummy_search_target);

         if (dst_test != dst.end()) {

            /*
             * lower_bound returned dst_test that would transfer >=
             * ideal_transfer when swapped with src_test.  Check
             * transfererence between src_test and dst_test for the
             * high-side transfer.  Also check the next smaller box in
             * dst for the low-side transfer.
             */

            d_bbb.evaluateBreak( hiside_acceptance_flags, hiside_transfer,
                           src_test->d_boxload - dst_test->d_boxload,
                           ideal_transfer, low_transfer, high_transfer );

            if ( hiside_acceptance_flags[2] == 1 ) {
               src_hiside = src_test;
               dst_hiside = dst_test;
               hiside_transfer = src_hiside->d_boxload - dst_hiside->d_boxload;
               if (d_print_swap_steps) {
                  tbox::plog << "    new hi-swap pair: " << (*src_hiside)
                             << " & " << (*dst_hiside) << " with transfer "
                             << hiside_transfer
                             << " missing by " << hiside_transfer-ideal_transfer
                             << std::endl;
               }
            }

            if (dst_test != dst.begin()) {
               --dst_test; // Now, src_test and dst_test transferer by *less* than ideal_transfer.

               d_bbb.evaluateBreak( loside_acceptance_flags, loside_transfer,
                              src_test->d_boxload - dst_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( loside_acceptance_flags[2] == 1 ) {
                  src_loside = src_test;
                  dst_loside = dst_test;
                  loside_transfer = src_loside->d_boxload - dst_loside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap pair: " << (*src_loside)
                                << " & " << (*dst_loside) << " with transfer "
                                << loside_transfer
                                << " missing by " << loside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
            }

         } else {

            /*
             * The ideal dst to swap is smaller than the smallest dst
             * box.  So the only choice is swapping src_test for nothing.
             * Chech this against the current high- and low-side choices.
             */
            if (src_test->d_boxload > ideal_transfer) {
               // Moving src_test to src is moving too much--hiside.

               d_bbb.evaluateBreak( hiside_acceptance_flags, hiside_transfer,
                              src_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( hiside_acceptance_flags[2] == 1 ) {
                  src_hiside = src_test;
                  dst_hiside = dst.end();
                  hiside_transfer = src_hiside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new hi-swap source: " << (*src_hiside)
                                << " & " << "no dst" << " with transfer "
                                << (src_hiside->d_boxload)
                                << " missing by " << hiside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
            } else {
               // Moving src_test to src is moving (just right or) too little--loside.

               d_bbb.evaluateBreak( loside_acceptance_flags, loside_transfer,
                              src_test->d_boxload,
                              ideal_transfer, low_transfer, high_transfer );

               if ( loside_acceptance_flags[2] == 1 ) {
                  src_loside = src_test;
                  dst_loside = dst.end();
                  loside_transfer = src_loside->d_boxload;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap source: " << (*src_loside)
                                << " & " << "no dst" << " with transfer "
                                << (src_loside->d_boxload)
                                << " missing by " << loside_transfer-ideal_transfer
                                << std::endl;
                  }
               }
               /*
                * Break out of the loop early because there is no
                * point checking smaller src boxes.
                */
               break;
            }
         }

         if ( ( low_transfer <= loside_transfer && loside_transfer <= high_transfer ) ||
              ( low_transfer <= hiside_transfer && hiside_transfer <= high_transfer ) ) {
            // Found a transfer satisfying the range.  Stop searching.
            break;
         }

      }

   }

   /*
    * Swapping does not produce new cuts, so it is ok to omit the penalties
    * arising from cutting.
    */
   double current_balance_penalty = static_cast<double>(ideal_transfer);
   double balance_penalty_loside = static_cast<double>(loside_transfer-ideal_transfer);
   double balance_penalty_hiside = static_cast<double>(hiside_transfer-ideal_transfer);

   if (d_print_swap_steps) {
      tbox::plog.setf(std::ios_base::fmtflags(0),std::ios_base::floatfield);
      tbox::plog.precision(8);
      tbox::plog << "    Swap candidates give penalties (unswap,lo,hi): "
                 << current_balance_penalty << " , " << balance_penalty_loside
                 << " , " << balance_penalty_hiside << std::endl;
   }

   bool found_swap = false;
   BoxTransitSet::iterator isrc = src.end();
   BoxTransitSet::iterator idst = dst.end();
   actual_transfer = 0;

   if ( d_bbb.evaluateBreak( hiside_acceptance_flags, actual_transfer, hiside_transfer,
                       ideal_transfer, low_transfer, high_transfer ) ) {
      isrc = src_hiside;
      idst = dst_hiside;
      actual_transfer = hiside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking hiside." << std::endl;
      }
   }

   if ( d_bbb.evaluateBreak( loside_acceptance_flags, actual_transfer, loside_transfer,
                       ideal_transfer, low_transfer, high_transfer ) ) {
      isrc = src_loside;
      idst = dst_loside;
      actual_transfer = loside_transfer;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking loside." << std::endl;
      }
   }


   if (found_swap) {

      // We can improve balance_penalty by swapping isrc with idst.
      if (d_print_swap_steps) {
         tbox::plog << "    Swapping " << actual_transfer << " units using ";
         if (isrc != src.end()) tbox::plog << *isrc;
         else tbox::plog << "X";
         tbox::plog << " <--> ";
         if (idst != dst.end()) tbox::plog << *idst;
         else tbox::plog << "X";
         tbox::plog << std::endl;
      }

      if (isrc != src.end()) {
         dst.insert(*isrc);
         src.erase(isrc);
      }
      if (idst != dst.end()) {
         src.insert(*idst);
         dst.erase(idst);
      }


   } else {
      if (d_print_swap_steps) {
         if ( isrc == src.end() ) {
            tbox::plog << "    Cannot find swap pair for " << ideal_transfer
                       << " units." << std::endl;
         }
         else {
            tbox::plog << "    Keeping original (no swap)." << std::endl;
         }
      }
   }

   t_find_swap_pair->stop();
   return found_swap;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
BoxTransitSet::setTimers()
{
   t_adjust_load = tbox::TimerManager::getManager()->
      getTimer("mesh::BoxTransitSet::adjustLoad()");
   t_adjust_load_by_popping = tbox::TimerManager::getManager()->
      getTimer("mesh::BoxTransitSet::adjustLoadByPopping()");
   t_adjust_load_by_swapping = tbox::TimerManager::getManager()->
      getTimer("mesh::BoxTransitSet::adjustLoadBySwapping()");
   t_shift_loads_by_breaking = tbox::TimerManager::getManager()->
      getTimer("mesh::BoxTransitSet::adjustLoadByBreaking()");
   t_find_swap_pair = tbox::TimerManager::getManager()->
      getTimer("mesh::BoxTransitSet::swapLoadPair()");
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
