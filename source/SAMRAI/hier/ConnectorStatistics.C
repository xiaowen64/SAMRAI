/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Statistical characteristics of a Connector.
 *
 ************************************************************************/
#ifndef included_hier_ConnectorStatistics_C
#define included_hier_ConnectorStatistics_C

#include "SAMRAI/hier/ConnectorStatistics.h"

#include "SAMRAI/tbox/MathUtilities.h"

#ifndef SAMRAI_INLINE
// #include "SAMRAI/hier/ConnectorStatistics.I"
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


std::string ConnectorStatistics::s_quantity_names[NUMBER_OF_QUANTITIES];
int ConnectorStatistics::s_longest_length;

tbox::StartupShutdownManager::Handler
ConnectorStatistics::s_initialize_finalize_handler(
   ConnectorStatistics::initializeCallback,
   0,
   0,
   ConnectorStatistics::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);


/*
************************************************************************
* Constructor.
************************************************************************
*/
ConnectorStatistics::ConnectorStatistics(
   const Connector &connector )
   : d_connector(connector)
{
}


/*
************************************************************************
************************************************************************
*/
ConnectorStatistics::StatisticalQuantities::StatisticalQuantities()
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

void ConnectorStatistics::computeLocalConnectorStatistics(
   StatisticalQuantities &sq) const
{
   if (!d_connector.isInitialized()) {
      TBOX_ERROR("ConnectorStatistics::computeLocalStatistics cannot compute\n"
                 <<"statistics for uninitialized Connector.");
      return;
   }

   d_connector.cacheGlobalReducedData();

   const BoxLevel &base(d_connector.getBase());
   const tbox::SAMRAI_MPI &mpi(base.getMPI());
   const tbox::Dimension& dim(base.getDim());

   /*
    * Whether had boxes should be refined or coarsened before box
    * calculus.
    */
   const bool refine_head = d_connector.getHeadCoarserFlag() &&
      ( d_connector.getRatio() != IntVector::getOne(dim) ||
        !d_connector.ratioIsExact() );
   const bool coarsen_head = !d_connector.getHeadCoarserFlag() &&
      ( d_connector.getRatio() != IntVector::getOne(dim) ||
        !d_connector.ratioIsExact() );


   /*
    * Compute per-processor statistics.  Some quantities are readily
    * available while others are computed in the loops following.
    */

   sq.d_values[NUMBER_OF_BASE_BOXES] = base.getLocalNumberOfBoxes();
   sq.d_values[NUMBER_OF_BASE_CELLS] = base.getLocalNumberOfCells();

   sq.d_values[HAS_ANY_NEIGHBOR_SETS] = d_connector.getLocalNumberOfNeighborSets() > 0;
   sq.d_values[NUMBER_OF_NEIGHBOR_SETS] = d_connector.getLocalNumberOfNeighborSets();

   sq.d_values[HAS_ANY_RELATIONSHIPS] = d_connector.getLocalNumberOfRelationships() > 0;
   sq.d_values[NUMBER_OF_RELATIONSHIPS] =  d_connector.getLocalNumberOfRelationships();
   sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS] = tbox::MathUtilities<double>::getMax();
   sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS] = 0;

   sq.d_values[NUMBER_OF_NEIGHBORS] = 0;
   sq.d_values[NUMBER_OF_LOCAL_NEIGHBORS] = 0.;
   sq.d_values[NUMBER_OF_REMOTE_NEIGHBORS] = 0.;
   sq.d_values[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] = 0.;

   hier::BoxSet visible_neighbors; // All neighbors of local base boxes.

   for ( Connector::ConstNeighborhoodIterator nbi=d_connector.begin();
         nbi!=d_connector.end(); ++nbi ) {

      const BoxSet &neighbors = nbi->second;

      visible_neighbors.insert( neighbors.begin(), neighbors.end() );

      sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS] =
         tbox::MathUtilities<double>::Min(sq.d_values[MIN_NUMBER_OF_RELATIONSHIPS],
                                          neighbors.size());

      sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS] =
         tbox::MathUtilities<double>::Max(sq.d_values[MAX_NUMBER_OF_RELATIONSHIPS],
                                          neighbors.size());

      Box base_box = *base.getBoxStrict(nbi->first);
      base_box.grow(d_connector.getConnectorWidth());

      for ( Connector::ConstNeighborIterator ni=neighbors.begin();
            ni!=neighbors.end(); ++ni ) {

         Box neighbor = *ni;
         if ( refine_head ) {
            neighbor.refine(d_connector.getRatio());
         }
         else if ( coarsen_head ) {
            neighbor.coarsen(d_connector.getRatio());
         }
         neighbor *= base_box;
         const int size = neighbor.size();

         sq.d_values[OVERLAP_SIZE] += size;
         if ( neighbor.getOwnerRank() == mpi.getRank() ) {
            sq.d_values[LOCAL_OVERLAP_SIZE] += size;
         }
         else {
            sq.d_values[REMOTE_OVERLAP_SIZE] += size;
         }

      }

   }


   sq.d_values[NUMBER_OF_NEIGHBORS] = visible_neighbors.size();

   std::set<int> remote_neighbor_owners;
   for ( hier::BoxSet::const_iterator bi=visible_neighbors.begin();
         bi!=visible_neighbors.end(); ++bi ) {
      const Box &neighbor = *bi;

      sq.d_values[NUMBER_OF_LOCAL_NEIGHBORS] +=
         neighbor.getOwnerRank() == mpi.getRank();

      sq.d_values[NUMBER_OF_REMOTE_NEIGHBORS] +=
         neighbor.getOwnerRank() != mpi.getRank();

      remote_neighbor_owners.insert(neighbor.getOwnerRank());

   }
   remote_neighbor_owners.erase(mpi.getRank());
   sq.d_values[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] = remote_neighbor_owners.size();

   return;
}


/*
 ***********************************************************************
 * Write out local and globally reduced statistics on the relationships.
 ***********************************************************************
 */

void ConnectorStatistics::printNeighborStats(
   std::ostream& co,
   const std::string& border) const
{
   if (!d_connector.isInitialized()) {
      co << "Connector is unininitialized.\n";
      return;
   }

   const BoxLevel &base(d_connector.getBase());
   const tbox::SAMRAI_MPI &mpi(base.getMPI());
   const tbox::Dimension& dim(base.getDim());

   StatisticalQuantities sq;
   computeLocalConnectorStatistics(sq);

   /*
    * Global reduction for statistics in sq.
    */

   StatisticalQuantities sq_min(sq); // Global min of sq.
   StatisticalQuantities sq_max(sq); // Global max of sq.
   StatisticalQuantities sq_sum(sq); // Global sum of sq.

   int rank_of_min[NUMBER_OF_QUANTITIES];
   int rank_of_max[NUMBER_OF_QUANTITIES];
   if (mpi.getSize() > 1) {
      mpi.AllReduce(sq_min.d_values, NUMBER_OF_QUANTITIES, MPI_MINLOC, rank_of_min);
      mpi.AllReduce(sq_max.d_values, NUMBER_OF_QUANTITIES, MPI_MAXLOC, rank_of_max);
      mpi.AllReduce(sq_sum.d_values, NUMBER_OF_QUANTITIES, MPI_SUM);
   } else {
      for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
         rank_of_min[i] = rank_of_max[i] = 0;
      }
   }

   co.unsetf(std::ios::fixed | std::ios::scientific);
   co.precision(3);

   co << border << "N = " << base.getGlobalNumberOfBoxes() << " (global number of boxes)\n"
      << border << "P = " << mpi.getSize() << " (number of processes)\n"
      << border << std::setw(s_longest_length) << std::string() << "    local        min               max             sum    sum/N    sum/P\n";

   for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
      co << border << std::setw(s_longest_length) << std::left << s_quantity_names[i]
         << ' ' << std::setw(8) << std::right << sq.d_values[i]
         << ' ' << std::setw(8) << std::right << sq_min.d_values[i] << " @ "
         << std::setw(6) << std::left << rank_of_min[i]
         << ' ' << std::setw(8) << std::right << sq_max.d_values[i] << " @ "
         << std::setw(6) << std::left << rank_of_max[i]
         << ' ' << std::setw(8) << std::right << sq_sum.d_values[i]
         << ' ' << std::setw(8) << std::right << sq_sum.d_values[i] / base.getGlobalNumberOfBoxes()
         << ' ' << std::setw(8) << std::right << sq_sum.d_values[i] / mpi.getSize()
         << '\n';
   }

   return;
}


/*
 ***********************************************************************
 ***********************************************************************
 */

void ConnectorStatistics::initializeCallback()
{

   s_quantity_names[NUMBER_OF_BASE_BOXES] = "num base boxes";
   s_quantity_names[NUMBER_OF_BASE_CELLS] = "num base cells";

   s_quantity_names[HAS_ANY_NEIGHBOR_SETS] = "has any neighbor sets";
   s_quantity_names[NUMBER_OF_NEIGHBOR_SETS] = "num neighbor sets";

   s_quantity_names[HAS_ANY_RELATIONSHIPS] = "has any relationships";
   s_quantity_names[NUMBER_OF_RELATIONSHIPS] = "num relationships";
   s_quantity_names[MIN_NUMBER_OF_RELATIONSHIPS] = "min num relationships";
   s_quantity_names[MAX_NUMBER_OF_RELATIONSHIPS] = "max num relationships";

   s_quantity_names[NUMBER_OF_NEIGHBORS] = "num neighbors";
   s_quantity_names[NUMBER_OF_LOCAL_NEIGHBORS] = "num local neighbors";
   s_quantity_names[NUMBER_OF_REMOTE_NEIGHBORS] = "num remote neighbors";
   s_quantity_names[NUMBER_OF_REMOTE_NEIGHBOR_OWNERS] = "num remote neighbor owners";

   s_quantity_names[OVERLAP_SIZE] = "overlap size";
   s_quantity_names[LOCAL_OVERLAP_SIZE] = "local overlap size";
   s_quantity_names[REMOTE_OVERLAP_SIZE] = "remote overlap size";

   s_longest_length = 0;
   for ( int i=0; i<NUMBER_OF_QUANTITIES; ++i ) {
      s_longest_length = tbox::MathUtilities<int>::Max(
         s_longest_length, s_quantity_names[i].length());
   }
}

/*
 ***************************************************************************
 ***************************************************************************
 */

void ConnectorStatistics::finalizeCallback()
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
