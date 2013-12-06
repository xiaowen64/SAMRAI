/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Utilitiy for breaking boxes during partitioning.
 *
 ************************************************************************/

#ifndef included_mesh_BalanceBoxBreaker
#define included_mesh_BalanceBoxBreaker

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/mesh/PartitioningParams.h"

#include <vector>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Utilities for breaking up boxes during partitioning.
 */

class BalanceBoxBreaker {

public:

   BalanceBoxBreaker(
      const PartitioningParams &pparams,
      bool print_break_steps = false );

   BalanceBoxBreaker( const BalanceBoxBreaker &other );


   /*!
    * @brief Break off a given load size from a given Box.
    *
    * Attempt to break off the ideal_load, or at least a load
    * inside the range [low_load, high_load].
    *
    * @param[out] breakoff Boxes broken off (usually just one).
    *
    * @param[out] leftover Remainder of box after breakoff is gone.
    *
    * @param[out] brk_load The load broken off.
    *
    * @param[in] box Box to break.
    *
    * @param[in] ideal_load Ideal load to break.
    *
    * @param[in] low_load
    *
    * @param[in] high_load
    *
    * @param[in] threshold_width Try to avoid making boxes thinner
    * than this width in any direction.
    *
    * @return whether a successful break was made.
    *
    * @pre ideal_load_to_break > 0
    */
   bool
   breakOffLoad(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& box,
      double ideal_load,
      double low_load,
      double high_load ,
      double threshold_width ) const;


   void setPrintBreakSteps( bool print_break_steps ) {
      d_print_break_steps = print_break_steps;
   }


   /*!
    * @brief Break up box bursty against box solid and adds the pieces
    * in the given container.
    *
    * This version differs from that in BoxContainer in that it tries
    * to minimize slivers.
    */
   static void
   burstBox(
      std::vector<hier::Box>& boxes,
      const hier::Box& bursty,
      const hier::Box& solid );

   static double
   computeBalancePenalty(double imbalance)
   {
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   /*!
    * @brief Compute a score that is low for box widths smaller than
    * some threshold_width.
    */
   static double computeWidthScore(
      const hier::IntVector &box_size,
      double threshold_width );

   /*!
    * @brief Compute a combined width score for multiple boxes.
    */
   static double computeWidthScore(
      const std::vector<hier::Box> &boxes,
      double threshold_width );


private:

   struct TrialBreak {
      TrialBreak( const PartitioningParams &pparams,
                  double threshold_width,
                  double ideal_load, double low_load, double high_load );
      //! @brief Break box from whole and store results.
      void breakBox( const hier::Box &box,
                     const hier::Box &whole );
      void swapWithReversedTrial( TrialBreak &reversed, double whole_box_load );
      //! @brief Compute merits vs doing nothing and return improvement flag.
      bool computeMerits();
      //! @brief Whether this improves over another (or degrades or leaves alone).
      int improvesOver( const TrialBreak &other ) const;

      double d_breakoff_load;
      std::vector<hier::Box> d_breakoff;
      std::vector<hier::Box> d_leftover;
      double d_ideal_load, d_low_load, d_high_load;
      int d_width_score;
      int d_balance_penalty;
      //! @brief Flags from comparing this trial vs doing nothing.
      int d_flags[4];
      const PartitioningParams *d_pparams;
      const double d_threshold_width;
   };

   bool
   breakOffLoad_planar(
      TrialBreak &trial,
      const hier::Box& box,
      const std::vector<std::vector<bool> >& bad_cuts ) const;

   bool
   breakOffLoad_cubic(
      TrialBreak &trial,
      const hier::Box& box,
      const std::vector<std::vector<bool> >& bad_cuts ) const;

   void setTimers();

   const PartitioningParams *d_pparams;

   //@{
   //! @name Debugging and diagnostic data

   bool d_print_break_steps;
   boost::shared_ptr<tbox::Timer> t_break_off_load;
   boost::shared_ptr<tbox::Timer> t_find_bad_cuts;

   //@}

};

}
}
#endif
