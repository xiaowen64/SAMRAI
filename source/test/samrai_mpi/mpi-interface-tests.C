/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   $Description 
 *
 ************************************************************************/

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "mpi-interface-tests.h"

using namespace SAMRAI;
using namespace SAMRAI::tbox;

/*
 * Start up MPI, SAMRAI and run various test of the SAMRAI_MPI interfaces.
 */
void mpiInterfaceTests(
   int& argc,
   char **& argv,
   int& fail_count,
   bool runtime_mpi,
   bool mpi_disabled)
{
   (void)runtime_mpi;

#ifndef HAVE_MPI
   if (runtime_mpi) {
      // Run-time MPI test is only for configuring with MPI.
      return;
   }
#endif
   if (mpi_disabled) {
      tbox::SAMRAI_MPI::initMPIDisabled();
   } else {
      tbox::SAMRAI_MPI::init(&argc, &argv);
   }
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   mpiInterfaceTestBcast(fail_count);

   mpiInterfaceTestAllreduce(fail_count);

   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();
   SAMRAI_MPI::finalize();
}

/*
 * Broadcast test: broadcast the rank from bcast_root.
 */
int mpiInterfaceTestBcast(
   int& fail_count)
{
   SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int rank = mpi.getRank();
   int nproc = mpi.getSize();

   if (nproc != 1) {
      int bcast_root = 0;
      int bcast_data = rank;
      mpi.Bcast(&bcast_data, 1, MPI_INT, bcast_root);
      if (bcast_data != bcast_root) {
         perr << "Broadcast test failed." << std::endl;
         ++fail_count;
         return 1;
      }
   }
   return 0;
}

/*
 * Sum all-reduce test: sum number of processes.
 */
int mpiInterfaceTestAllreduce(
   int& fail_count)
{
   SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int nproc = mpi.getSize();

   if (nproc != 1) {
      int one = 1;
      int sum = -1;
      mpi.Allreduce(&one, &sum, 1, MPI_INT, MPI_SUM);
      if (sum != nproc) {
         perr << "Allreduce test failed." << std::endl;
         ++fail_count;
         return 1;
      }
   }
   return 0;
}
