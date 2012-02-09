/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a BoxLevel.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevelStatistics_C
#define included_hier_BoxLevelStatistics_C

#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
// #include "SAMRAI/hier/BoxLevelStatistics.I"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {


std::string BoxLevelStatistics::s_quantity_names[NUMBER_OF_QUANTITIES];
int BoxLevelStatistics::s_longest_length;

tbox::StartupShutdownManager::Handler
BoxLevelStatistics::s_initialize_finalize_handler(
   BoxLevelStatistics::initializeCallback,
   0,
   0,
   BoxLevelStatistics::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);


/*
************************************************************************
* Constructor.
************************************************************************
*/
BoxLevelStatistics::BoxLevelStatistics(
   const BoxLevel &box_level )
   : d_box_level(box_level)
{
}


/*
************************************************************************
************************************************************************
*/
BoxLevelStatistics::StatisticalQuantities::StatisticalQuantities()
{
   for ( int i=0; i<NUMBER_OF_QUANTITIES; ++i ) {
      d_values[i] = 0;
   }
}


/*
 ***********************************************************************
 * Compute statistics for local process.
 ***********************************************************************
 */

void BoxLevelStatistics::computeLocalBoxLevelStatistics(
   StatisticalQuantities &sq) const
{
   if (!d_box_level.isInitialized()) {
      TBOX_ERROR("BoxLevelStatistics::computeLocalStatistics cannot compute\n"
                 <<"statistics for uninitialized BoxLevel.");
      return;
   }

   d_box_level.cacheGlobalReducedData();

   const tbox::Dimension& dim(d_box_level.getDim());
   tbox::SAMRAI_MPI mpi(d_box_level.getMPI());

   /*
    * Compute per-processor statistics.  Some quantities are readily
    * available while others are computed in the loop following.
    *
    * Aspect ratio uses a generalized formula that goes to 1 when box
    * has same length on all sides (regardless of dimension),
    * degenerates to the rectangular aspect ratio in 2D, and grows
    * appropriately for dimensions higher than 2.
    */

   sq.d_values[HAS_ANY_BOX] = (d_box_level.getLocalNumberOfBoxes() > 0);
   sq.d_values[NUMBER_OF_CELLS] =
      static_cast<double>(d_box_level.getLocalNumberOfCells());
   sq.d_values[NUMBER_OF_BOXES] =
      static_cast<double>(d_box_level.getLocalNumberOfBoxes());
   sq.d_values[LARGEST_BOX_SIZE] = 0;
   sq.d_values[SMALLEST_BOX_SIZE] = tbox::MathUtilities<double>::getMax();
   sq.d_values[LONGEST_BOX_LEN] = 0;
   sq.d_values[SHORTEST_BOX_LEN] =
      (d_box_level.getLocalNumberOfBoxes() == 0 ?
       0 : tbox::MathUtilities<double>::getMax());
   sq.d_values[LARGEST_ASPECT_RATIO] = 0;
   sq.d_values[SMALLEST_ASPECT_RATIO] = tbox::MathUtilities<double>::getMax();
   sq.d_values[SUM_SURFACE_AREA] = 0.;
   sq.d_values[SUM_NORMALIZED_SURFACE_AREA] = 0.;

   const BoxContainer& boxes = d_box_level.getBoxes();

   for (RealBoxConstIterator ni(boxes); ni.isValid(); ++ni) {

      const Box& mapped_box = *ni;
      const IntVector boxdims = mapped_box.numberCells();
      const int boxvol = boxdims.getProduct();
      const int longdim = boxdims.max();
      const int shortdim = boxdims.min();
      double aspect_ratio = 0.0;
      double surfarea = 0.;
      for (int d = 0; d < dim.getValue(); ++d) {
         surfarea += 2 * double(boxvol) / boxdims(d);
         double tmp = static_cast<double>(boxdims(d))/shortdim - 1.0;
         aspect_ratio += tmp*tmp;
      }
      aspect_ratio = 1.0 + sqrt(aspect_ratio);

      sq.d_values[LARGEST_BOX_SIZE] =
         tbox::MathUtilities<double>::Max(sq.d_values[LARGEST_BOX_SIZE],
                                          boxvol);
      sq.d_values[SMALLEST_BOX_SIZE] =
         tbox::MathUtilities<double>::Min(sq.d_values[SMALLEST_BOX_SIZE],
                                          boxvol);

      sq.d_values[LONGEST_BOX_LEN] =
         tbox::MathUtilities<double>::Max(sq.d_values[LONGEST_BOX_LEN],
                                          longdim);
      sq.d_values[SHORTEST_BOX_LEN] =
         tbox::MathUtilities<double>::Min(sq.d_values[SHORTEST_BOX_LEN],
                                          shortdim);

      sq.d_values[LARGEST_ASPECT_RATIO] =
         tbox::MathUtilities<double>::Max(sq.d_values[LARGEST_ASPECT_RATIO],
                                          aspect_ratio);
      sq.d_values[SMALLEST_ASPECT_RATIO] =
         tbox::MathUtilities<double>::Min(sq.d_values[SMALLEST_ASPECT_RATIO],
                                          aspect_ratio);

      sq.d_values[SUM_SURFACE_AREA] += surfarea;

   }

   /*
    * Smallest surface area possible for the number of cells perfectly
    * distributed in mpi.
    */
   double ideal_surfarea =
      2*dim.getValue() *
      pow(double(d_box_level.getGlobalNumberOfCells()) / mpi.getSize(),
          double(dim.getValue() - 1) / dim.getValue());

   sq.d_values[SUM_NORMALIZED_SURFACE_AREA] =
      sq.d_values[SUM_SURFACE_AREA] / ideal_surfarea;

}


/*
 ***********************************************************************
 * Write out local and globally reduced statistics on the boxes.
 ***********************************************************************
 */

void BoxLevelStatistics::printBoxStats(
   std::ostream& co,
   const std::string& border) const
{
   if (!d_box_level.isInitialized()) {
      co << "BoxLevel is unininitialized.\n";
      return;
   }

   const tbox::Dimension& dim(d_box_level.getDim());
   tbox::SAMRAI_MPI mpi(d_box_level.getMPI());

   StatisticalQuantities sq;
   computeLocalBoxLevelStatistics(sq);

   /*
    * Global reduction for statistics in sq.
    */

   StatisticalQuantities sq_min(sq); // Global min of sq.
   StatisticalQuantities sq_max(sq); // Global max of sq.
   StatisticalQuantities sq_sum(sq); // Global sum of sq.

   int rank_of_min[NUMBER_OF_QUANTITIES];
   int rank_of_max[NUMBER_OF_QUANTITIES];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(
         sq_min.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MINLOC,
         rank_of_min);
      mpi.AllReduce(
         sq_max.d_values,
         NUMBER_OF_QUANTITIES,
         MPI_MAXLOC,
         rank_of_max);
      mpi.AllReduce(sq_sum.d_values, NUMBER_OF_QUANTITIES, MPI_SUM);
   } else {
      for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
         rank_of_min[i] = rank_of_max[i] = 0;
      }
   }

   co.unsetf(std::ios::fixed | std::ios::scientific);
   co.precision(3);

   /*
    * Smallest surface area possible for the number of cells perfectly
    * distributed in mpi.
    */
   double ideal_surfarea =
      2*dim.getValue() *
      pow(double(d_box_level.getGlobalNumberOfCells()) / mpi.getSize(),
          double(dim.getValue() - 1) / dim.getValue());

   co << border << "N = " << d_box_level.getGlobalNumberOfBoxes()
      << " (global number of boxes)\n"
      << border << "P = " << mpi.getSize() << " (number of processes)\n"
      << border << "Ideal surface area is " << ideal_surfarea << " for "
      << (d_box_level.getGlobalNumberOfCells()/mpi.getSize()) << " cells\n"
      << border << std::setw(s_longest_length) << std::string()
      << "    local        min               max             sum    sum/N    sum/P\n";

   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      co << border << std::setw(s_longest_length) << std::left
         << s_quantity_names[i]
         << ' ' << std::setw(8) << std::right << sq.d_values[i]
         << ' ' << std::setw(8) << std::right << sq_min.d_values[i] << " @ "
         << std::setw(6) << std::left << rank_of_min[i]
         << ' ' << std::setw(8) << std::right << sq_max.d_values[i] << " @ "
         << std::setw(6) << std::left << rank_of_max[i]
         << ' ' << std::setw(8) << std::right << sq_sum.d_values[i]
         << ' ' << std::setw(8)
         << std::right << sq_sum.d_values[i] / d_box_level.getGlobalNumberOfBoxes()
         << ' ' << std::setw(8)
         << std::right << sq_sum.d_values[i] / mpi.getSize() << '\n';
   }

   return;
}


/*
 ***********************************************************************
 ***********************************************************************
 */

void BoxLevelStatistics::initializeCallback()
{
   s_quantity_names[HAS_ANY_BOX]                 = "has any box";
   s_quantity_names[NUMBER_OF_CELLS]             = "num cells";
   s_quantity_names[NUMBER_OF_BOXES]             = "num boxes";
   s_quantity_names[LARGEST_BOX_SIZE]            = "largest box size";
   s_quantity_names[SMALLEST_BOX_SIZE]           = "smallest box size";
   s_quantity_names[LONGEST_BOX_LEN]             = "longest box length";
   s_quantity_names[SHORTEST_BOX_LEN]            = "shortest box length";
   s_quantity_names[LARGEST_ASPECT_RATIO]        = "largest aspect ratio";
   s_quantity_names[SMALLEST_ASPECT_RATIO]       = "smallest aspect ratio";
   s_quantity_names[SUM_SURFACE_AREA]            = "sum surf area";
   s_quantity_names[SUM_NORMALIZED_SURFACE_AREA] = "sum norm surf area";
   s_longest_length = 0;
   for ( int i=0; i<NUMBER_OF_QUANTITIES; ++i ) {
      s_longest_length = tbox::MathUtilities<int>::Max(
         s_longest_length, static_cast<int>(s_quantity_names[i].length()));
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */

void BoxLevelStatistics::finalizeCallback()
{
   for ( int i=0; i<NUMBER_OF_QUANTITIES; ++i ) {
      s_quantity_names[i].clear();
   }
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
