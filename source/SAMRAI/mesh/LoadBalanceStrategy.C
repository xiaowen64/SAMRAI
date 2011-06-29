/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for box load balancing routines. 
 *
 ************************************************************************/

#ifndef included_mesh_LoadBalanceStrategy_C
#define included_mesh_LoadBalanceStrategy_C

#include "SAMRAI/mesh/LoadBalanceStrategy.h"

#include <cstdlib>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

// using namespace std;

int LoadBalanceStrategy::s_sequence_number = 0;

/*
 *************************************************************************
 *									*
 * The constructor and destructor for LoadBalanceStrategy do        *
 * nothing that could be considered even remotely interesting.		*
 *									*
 *************************************************************************
 */

LoadBalanceStrategy::LoadBalanceStrategy()
{
}

LoadBalanceStrategy::~LoadBalanceStrategy()
{
}

void LoadBalanceStrategy::loadBalanceMappedBoxLevel(
   hier::MappedBoxLevel& balance_mapped_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::Connector& unbalanced_to_attractor,
   const hier::Connector& attractor_to_unbalanced,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::MappedBoxLevel& domain_mapped_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   NULL_USE(balance_mapped_box_level);
   NULL_USE(balance_to_anchor);
   NULL_USE(anchor_to_balance);
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(unbalanced_to_attractor);
   NULL_USE(attractor_to_unbalanced);
   NULL_USE(min_size);
   NULL_USE(max_size);
   NULL_USE(domain_mapped_box_level);
   NULL_USE(bad_interval);
   NULL_USE(cut_factor);
   NULL_USE(rank_group);
   TBOX_ERROR(
      "Unusable base method LoadBalanceStrategy::loadBalanceMappedBoxLevel used.");
}

/*
 *************************************************************************
 * Report the load balance on processor, primarily
 * for debugging and checking load balance quality.
 *************************************************************************
 */
void LoadBalanceStrategy::markLoadForPostprocessing(
   int rank,
   double load,
   int nbox)
{
   tbox::plog << "Load mark " << s_sequence_number++
              << " proc " << rank
              << " load " << load
              << " nbox " << nbox
              << "\n";
}

/*
 *************************************************************************
 * Gather and report load balance for a single balancing.
 *************************************************************************
 */
void LoadBalanceStrategy::gatherAndReportLoadBalance(
   double local_load,
   const tbox::SAMRAI_MPI& mpi,
   std::ostream& os) const
{
   int nproc = mpi.getSize();
   std::vector<double> workloads(nproc);
   if (mpi.getSize() > 1) {
      mpi.Allgather(&local_load,
         1,
         MPI_DOUBLE,
         &workloads[0],
         1,
         MPI_DOUBLE);
   } else {
      workloads[0] = local_load;
   }
   this->reportLoadBalance(workloads, os);
}

/*
 *************************************************************************
 * Gather and report load balance for multiple balancings.
 *************************************************************************
 */
void LoadBalanceStrategy::gatherAndReportLoadBalance(
   const std::vector<double>& local_loads,
   const tbox::SAMRAI_MPI& mpi,
   std::ostream& os) const
{
   if (mpi.getSize() > 1) {
      int nproc = mpi.getSize();
      std::vector<double> mutable_local_loads(local_loads);
      std::vector<double> global_workloads(nproc * local_loads.size());
      mpi.Allgather(&mutable_local_loads[0],
         static_cast<int>(local_loads.size()),
         MPI_DOUBLE,
         &global_workloads[0],
         static_cast<int>(local_loads.size()),
         MPI_DOUBLE);
      std::vector<double> workloads_at_seq_i(nproc);
      for (size_t i = 0; i < local_loads.size(); ++i) {
         for (int n = 0; n < nproc; ++n) {
            workloads_at_seq_i[n] = global_workloads[i + n * local_loads.size()];
         }
         os << "================ Sequence " << i << " ===============\n";
         this->reportLoadBalance(workloads_at_seq_i, os);
      }
   } else {
      std::vector<double> workloads(1);
      workloads[0] = local_loads[0];
      this->reportLoadBalance(workloads, os);
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void LoadBalanceStrategy::reportLoadBalance(
   const std::vector<double>& workloads,
   std::ostream& os)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int nproc = mpi.getSize();

   const double demarks[] = { 0.50,
                              0.70,
                              0.85,
                              0.92,
                              0.98,
                              1.02,
                              1.08,
                              1.15,
                              1.30,
                              1.50,
                              2.00 };
   const int ndemarks = 11;

   TBOX_ASSERT((int)workloads.size() == nproc);

   RankAndLoad* rank_and_load = new RankAndLoad[nproc];

   double total_load = 0.0;

   for (int i = 0; i < nproc; ++i) {
      rank_and_load[i].rank = i;
      rank_and_load[i].load = workloads[i];
      total_load += workloads[i];
   }
   qsort((void *)rank_and_load,
      nproc,
      sizeof(RankAndLoad),
      qsortRankAndLoadCompareAscending);

   const double avg_load = total_load / nproc;
   const double min_load = rank_and_load[0].load;
   const int r_min_load = rank_and_load[0].rank;
   const double max_load = rank_and_load[nproc - 1].load;
   const int r_max_load = rank_and_load[nproc - 1].rank;

   os << "total/avg loads: "
      << total_load << " / "
      << avg_load << "\n";
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
   os << std::setprecision(2)
      << "min/max loads: "
      << min_load << " @ P" << r_min_load << " / "
      << max_load << " @ P" << r_max_load << "   "
      << "diffs: "
      << min_load - avg_load << " / "
      << max_load - avg_load << "   "
      << std::setprecision(3)
      << "normalized: "
      << (avg_load != 0 ? min_load / avg_load : 0.0) << " / "
      << (avg_load != 0 ? max_load / avg_load : 0.0) << "\n";

   const char bars[] = "----";
   const char space[] = "   ";
   os.setf(std::ios_base::fixed);
   os << bars;
   for (int n = 0; n < ndemarks; ++n) {
      os << std::setw(4) << std::setprecision(2) << demarks[n] << bars;
   }
   os << '\n';

   double population;
   int irank = 0;
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif
   for (int n = 0; n < ndemarks; ++n) {
      double top = demarks[n];
      int old_irank = irank;
      for ( ; irank < nproc; ++irank)
         if (avg_load == 0 ||
             rank_and_load[irank].load / avg_load > top) break;
      int nrank = irank - old_irank;
      population = 100.0 * nrank / nproc;
      os << std::setw(5) << population << space;
   }
   population = 100.0 * (nproc - irank) / nproc;
   os << population << space;
   os << '\n';

   delete[] rank_and_load;
}

/*
 *************************************************************************
 * for use when sorting loads using the C-library qsort
 *************************************************************************
 */
int LoadBalanceStrategy::qsortRankAndLoadCompareDescending(
   const void* v,
   const void* w)
{
   const RankAndLoad* lv = (const RankAndLoad *)v;
   const RankAndLoad* lw = (const RankAndLoad *)w;
   if (lv->load > lw->load) return -1;

   if (lv->load < lw->load) return 1;

   return 0;
}

/*
 *************************************************************************
 * for use when sorting loads using the C-library qsort
 *************************************************************************
 */
int LoadBalanceStrategy::qsortRankAndLoadCompareAscending(
   const void* v,
   const void* w)
{
   const RankAndLoad* lv = (const RankAndLoad *)v;
   const RankAndLoad* lw = (const RankAndLoad *)w;
   if (lv->load < lw->load) return -1;

   if (lv->load > lw->load) return 1;

   return 0;
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
