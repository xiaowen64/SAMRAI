//
// File:        main_stats.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 586 $
// Modified:    $Date: 2005-08-23 10:49:46 -0700 (Tue, 23 Aug 2005) $
// Description: Main program to test statistics operations
//

#include "SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "tbox/SAMRAIManager.h"
#include "tbox/Database.h"
#include "tbox/MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Statistic.h"
#include "tbox/Statistician.h"
#include "tbox/Utilities.h"
#include <string>
using namespace std;

using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("stat.log");
  
   string input_filename;
   string restart_dirname;
   int restore_num = 0;

   bool is_from_restart = false;
 
   if ( (argc != 1) && (argc != 3) ) {
        tbox::pout << "USAGE:  " << argv[0] 
           << "<restart dir> <restore number> [options]\n"
           << "  options:\n"
           << "  none at this time"
           << endl;
        tbox::MPI::abort();
        return (-1);
   } else {
      if (argc == 3) {
         restart_dirname = argv[1];
         restore_num = atoi(argv[2]);
 
         is_from_restart = true;
      }
   }

#ifndef HAVE_HDF5
   is_from_restart = false;
#endif

   /*
    * Open the restart file and read information from file into the
    * restart database.
    */
   if (is_from_restart) {
      tbox::RestartManager::getManager()->
         openRestartFile(restart_dirname, restore_num, tbox::MPI::getNodes());
   }

   tbox::Pointer<tbox::Database> restart_db =
        tbox::RestartManager::getManager()->getRootDatabase();

   /*
    * Create tbox::Statistician.
    */

   tbox::Statistician* statistician = tbox::Statistician::createStatistician();

   /*
    * We no longer need the restart file so lets close it. 
    */
   tbox::RestartManager::getManager()->closeRestartFile();

   /*
    * Create statistics.
    */

   tbox::Pointer<tbox::Statistic> procstat1 = 
      statistician->getStatistic("procstat1", "PROC_STAT"); 
   tbox::Pointer<tbox::Statistic> procstat2 = 
      statistician->getStatistic("procstat2", "PROC_STAT"); 
   tbox::Pointer<tbox::Statistic> procstat3 = 
      statistician->getStatistic("procstat3", "PROC_STAT"); 

   tbox::Pointer<tbox::Statistic> patchstat1 = 
      statistician->getStatistic("patchstat1", "PATCH_STAT"); 
   tbox::Pointer<tbox::Statistic> patchstat2 =
      statistician->getStatistic("patchstat2", "PATCH_STAT");
   tbox::Pointer<tbox::Statistic> patchstat3 =
      statistician->getStatistic("patchstat3", "PATCH_STAT");

   int myrank = tbox::MPI::getRank();
   int i, j, ips;

   /*
    * procstat1 has a sequence length of three on each processor.
    * The value of sequence item i is i*(myrank + 1), where myrank
    * is the tbox::MPI process rank.
    */
   double factor = myrank + 1;
   for (ips = 0; ips < 3; ips++) {
      procstat1->recordProcStat(factor*ips);
   }

   /*
    * procstat2 has a sequence length of one on each processor.
    * The value of sequence item i is (i+1)*(myrank + 2), where myrank
    * is the tbox::MPI process rank.
    */ 
   factor = myrank + 2;
   for (ips = 0; ips < 1; ips++) { 
      procstat2->recordProcStat(factor*(ips+1));
   }

   /*
    * patchstat1 has a sequence length of two on each processor.
    * Each processor has exactly two patches. The value of each sequence 
    * item of the patches will be (patch id*2 + seq #), where the patch 
    * id's on each processor are myrank*2 and myrank*2+1.
    */
   for (ips = 0; ips < 2; ips++) {
      for (i = myrank * 2; i <= myrank * 2 + 1; i++) {
         patchstat1->recordPatchStat(i,                // patch number
                                     i * 2 + ips,      // data value
                                     ips);             // seq number
      }
   }

   /*
    * patchstat2 has a sequence length of two on each processor.
    * Each processor will have a variable number of patches equal to 
    * the tbox::MPI process number + 2 (i.e., processor zero
    * will have two patch values, processor one will have three patch
    * value, etc.).  The value for each patch will be twice the patch
    * number plus the sequence number.
    */
   for (ips = 0; ips < 2; ips++) {
      int patch_start_num = 0;
      for (i = 0; i < myrank; i++) {
         patch_start_num += i+2;
      }
      int num_patches = myrank + 2;
      for (i = patch_start_num; i < patch_start_num + num_patches; i++) {
         patchstat2->recordPatchStat(i,              // patch number
                                     2.0 * i + ips,  // data value
                                     ips);           // sequence number 
      }
   }

   // Test #1:
   int tval = statistician->getNumberProcessorStats();
   if (tval != 3) {
      tbox::perr << "FAILED: - Test #1a: tbox::Statistician::getNumberProcessorStats()\n"
           << "incorrect number of processor statistics found" << endl;
   } else {
      tbox::plog << "Test #1a successful" << endl;
   }
   tval = statistician->getNumberPatchStats();
   if (tval != 3) {
      tbox::perr << "FAILED: - Test #1b: tbox::Statistician::getNumberPatchStats()\n"
           << "incorrect number of patch statistics found" << endl;
   } else {
      tbox::plog << "Test #1b successful" << endl;
   }

   // Test #2:
   tbox::Pointer<tbox::Statistic> tstat; 
   if (!statistician->checkStatisticExists(tstat, "procstat2")) {
      tbox::perr << "FAILED: - Test #2a: tbox::Statistician::checkStatisticExists()\n"
           << "procstat2 added to statistician, but not found" << endl;
   } else {
      tbox::plog << "Test #2a successful" << endl;
      if (tstat->getName() != "procstat2") {
         tbox::perr << "FAILED: - Test #2b: tbox::Statistician::checkStatisticExists()\n"
              << "name of procstat2 does not match statistician entry" << endl; 
      } else {
         tbox::plog << "Test #2b successful" << endl;
      }
   }

   if (!statistician->checkStatisticExists(tstat, "patchstat1")) {
      tbox::perr << "FAILED: - Test #2c: tbox::Statistician::checkStatisticExists()\n"
           << "patchstat1 added to statistician, but not found" << endl;
   } else {
      tbox::plog << "Test #2c successful" << endl;
      if (tstat->getName() != "patchstat1") {
         tbox::perr << "FAILED: - Test #2d: tbox::Statistician::checkStatisticExists()\n" 
              << "name of patchstat1 does not match statistician entry" << endl;
      } else {
         tbox::plog << "Test #2c successful" << endl;
      }
   }

   if (statistician->checkStatisticExists(tstat, "dummy")) {
      tbox::perr << "FAILED: - Test #2e: tbox::Statistician::checkStatisticExists()\n"
           << "wrongly found statistic named dummy in statistician" << endl;
   } else {
      tbox::plog << "Test #2e successful" << endl;
   }
   tstat.setNull();

   // Test #3:
   if (statistician->getProcStatId(procstat1->getName()) 
       != procstat1->getInstanceId()) {
      tbox::perr << "FAILED: - Test #3a: tbox::Statistician::getProcStatId()\n"
           << "procstat1 has wrong instance id in statistician" << endl;
   } else {
      tbox::plog << "Test #3a successful" << endl;
   }
   if (statistician->getPatchStatId(patchstat2->getName()) 
       != patchstat2->getInstanceId()) {
      tbox::perr << "FAILED: - Test #3b: tbox::Statistician::getPatchStatId()\n" 
           << "patchstat2 has wrong instance id in statistician" << endl;
   } else {
      tbox::plog << "Test #3b successful" << endl;
   }
   if (statistician->getProcStatId(patchstat1->getName()) != -1) {
      tbox::perr << "FAILED: - Test #3c: tbox::Statistician::getProcStatId()\n" 
           << "patchstat1 is not a processor statistic" << endl;
   } else {
      tbox::plog << "Test #3c successful" << endl;
   }
   if (statistician->getPatchStatId(procstat2->getName()) != -1) {
      tbox::perr << "FAILED: - Test #3d: tbox::Statistician::getPatchStatId()\n" 
           << "procstat2 is not a patch statistic" << endl;
   } else {
      tbox::plog << "Test #3d successful" << endl;
   }

   // Test #4:
   if (procstat1->getType() != "PROC_STAT") {   
      tbox::perr << "FAILED: - Test #4a: tbox::Statistic::getType()\n" 
           << "procstat1 returns incorrect type" << endl;
   } else {
      tbox::plog << "Test #4a successful" << endl;
   }
   if (patchstat1->getType() != "PATCH_STAT") {
      tbox::perr << "FAILED: - Test #4b: tbox::Statistic::getType()\n"
           << "patchstat1 returns incorrect type" << endl;
   } else {
      tbox::plog << "Test #4b successful" << endl;
   }

   // Test #5:
   int procstat1_seqlen;
   int procstat2_seqlen;
   if (is_from_restart) {
      procstat1_seqlen = 6;
      procstat2_seqlen = 2;
   } else {
      procstat1_seqlen = 3;
      procstat2_seqlen = 1;
   }
   if (procstat1->getStatSequenceLength() != procstat1_seqlen) {
      tbox::perr << "FAILED: - Test #5a: tbox::Statistic::getStatSequenceLength()\n"
           << "procstat1 returns incorrect sequence length" << endl;
   } else {
      tbox::plog << "Test #5a successful" << endl;
   }
   if (procstat2->getStatSequenceLength() != procstat2_seqlen) {
      tbox::perr << "FAILED: - Test #5b: tbox::Statistic::getStatSequenceLength()\n"
           << "procstat2 returns incorrect sequence length" << endl;
   } else {
      tbox::plog << "Test #5b successful" << endl;
   }
   if (patchstat1->getStatSequenceLength() != 2) {
      tbox::perr << "FAILED: - Test #5c: tbox::Statistic::getStatSequenceLength()\n"
           << "patchstat1 returns incorrect sequence length" << endl;
   } else {
      tbox::plog << "Test #5c successful" << endl;
   }
   if (patchstat2->getStatSequenceLength() != 2) {
      tbox::perr << "FAILED: - Test #5d: tbox::Statistic::getStatSequenceLength()\n"
           << "patchstat2 returns incorrect sequence length" << endl;
   } else {
      tbox::plog << "Test #5d successful" << endl;
   }

   statistician->printLocalStatData(tbox::plog);

   statistician->finalize(); 

   statistician->printAllGlobalStatData(tbox::plog);

   statistician->finalize(); 

   statistician->printSpreadSheetOutput();

   statistician->finalize(); 

   if (tbox::MPI::getRank() == 0) {

      // Test #6:
      if (statistician->
         getGlobalProcStatSequenceLength(procstat1->getInstanceId()) != 
                                         procstat1_seqlen) {
         tbox::perr << "FAILED: - Test #6a: "
              << "Statistician::getGlobalProcStatSequenceLength()\n"
              << "incorrect sequence length returned for procstat1" << endl;
      } else {
         tbox::plog << "Test #6a successful" << endl;
      }

      if (statistician->
         getGlobalProcStatSequenceLength(procstat2->getInstanceId()) != 
                                         procstat2_seqlen) {
         tbox::perr << "FAILED: - Test #6b: "
              << "Statistician::getGlobalProcStatSequenceLength()\n"
              << "incorrect sequence length returned for procstat2" << endl;
      } else {
         tbox::plog << "Test #6b successful" << endl;
      }

      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         if (statistician->getGlobalProcStatValue(procstat1->getInstanceId(),
                                                  2, i) != (i+1)*2) {
            tbox::perr << "FAILED: - Test #6c: "
                 << "Statistician::getGlobalProcStatValue()\n"
                 << "incorrect global data returned for procstat1" << endl;
         } else {
            tbox::plog << "Test #6c successful" << endl;
         }
      }
      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         if (statistician->getGlobalProcStatValue(procstat2->getInstanceId(),
                                                  0, i) != i+2) {
            tbox::perr << "FAILED: - Test #6d: "
                 << "Statistician::getGlobalProcStatValue()\n"
                 << "incorrect global data returned for procstat2" << endl;
         } else {
            tbox::plog << "Test #6d successful" << endl;
         }
      }

      // Test #7:
      if (statistician->
         getGlobalPatchStatSequenceLength(patchstat1->getInstanceId()) != 2) {
         tbox::perr << "FAILED: - Test #7a: "
              << "Statistician::getGlobalPatchStatSequenceLength()\n"
              << "incorrect sequence length returned for patchstat1" << endl;
      } else {
         tbox::plog << "Test #7a successful" << endl;
      }

      int num_patches = tbox::MPI::getNodes()*2;
      if (statistician->getGlobalPatchStatNumberPatches(
                        patchstat1->getInstanceId(), 0) != num_patches) {
         tbox::perr << "FAILED: - Test #7b: "
              << "Statistician::getGlobalPatchStatNumberPatches()\n"
              << "incorrect num patches returned for patchstat1" << endl;
      } else {
         tbox::plog << "Test #7b successful" << endl;
      }

      for (i = 0; i < tbox::MPI::getNodes(); i++) {
   
         int pnum = 0;
         for (j = 0; j < 1; j++) {
            pnum = i*2 + j;
         }
         if (statistician->getGlobalPatchStatValue(
                patchstat1->getInstanceId(), 1, pnum) != pnum*2+1) {
            tbox::perr << "FAILED: - Test #7e: "
                 << "Statistician::getGlobalPatchStatValue()\n"
                 << "incorrect global data returned for patch " << pnum 
                 << "in patchstat1." << endl;
         } else {
            tbox::plog << "Test #7e successful" << endl;
         }

      }

      num_patches = 0;
      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         num_patches += i+2;
      }
      if (statistician->getGlobalPatchStatNumberPatches(
                        patchstat2->getInstanceId(), 0) != num_patches) {
         tbox::perr << "FAILED: - Test #7c: "
              << "Statistician::getGlobalPatchStatNumberPatches()\n"
              << "incorrect num patches returned for patchstat2" << endl;
      } else {
         tbox::plog << "Test #7c successful" << endl;
      }

      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         int pnum = 0;
         for (j = 0; j < i; j++) {
            pnum += j+2;
         }
         if (statistician->getGlobalPatchStatPatchMapping(
                           patchstat2->getInstanceId(), 1, pnum) != i) {
            tbox::perr << "FAILED: - Test #7d: "
                 << "Statistician::getGlobalPatchStatPatchMapping()\n"
                 << "incorrect mapping for patch " << pnum << "in " 
                 << "patchstat2." << endl; 
         } else {
            tbox::plog << "Test #7d successful" << endl;
         }

      }

      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         int pnum = 0;
         for (j = 0; j < i; j++) {
            pnum += j+2;
         }
         if (statistician->getGlobalPatchStatValue(
                           patchstat2->getInstanceId(), 1, pnum) != pnum*2+1) {
            tbox::perr << "FAILED: - Test #7f: "
                 << "Statistician::getGlobalPatchStatValue()\n"
                 << "incorrect global data returned for patch " << pnum 
                 << "in patchstat2." << endl;
         } else {
            tbox::plog << "Test #7f successful" << endl;
         }

      }

      // Test #8:
      double sum0 = 0.;
      double sum1 = 0.;
      double sum2 = 0.;

      for (i = 0; i < tbox::MPI::getNodes(); i++) {         
         sum1 += (double) i+1;
         sum2 += (double) 2*(i+1);
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatSum(
                        procstat1->getInstanceId(), 0), sum0)) {
         tbox::perr << "FAILED: - Test #8a: "
              << "Statistician::getGlobalProcStatSum()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8a successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatSum(
                        procstat1->getInstanceId(), 1), sum1)) {
         tbox::perr << "FAILED: - Test #8b: "
              << "Statistician::getGlobalProcStatSum()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8b successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatSum(
                        procstat1->getInstanceId(), 2), sum2)) {
         tbox::perr << "FAILED: - Test #8c: "
              << "Statistician::getGlobalProcStatSum()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8c successful" << endl;
      }

      double max0 = 0.;
      double max1 = tbox::MPI::getNodes();
      double max2 = 2*(tbox::MPI::getNodes());
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMax(
                        procstat1->getInstanceId(), 0), max0)) {
         tbox::perr << "FAILED: - Test #8d: "
              << "Statistician::getGlobalProcStatMax()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8d successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMax(
                        procstat1->getInstanceId(), 1), max1)) {
         tbox::perr << "FAILED: - Test #8e: "
              << "Statistician::getGlobalProcStatMax()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8e successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMax(
                        procstat1->getInstanceId(), 2), max2)) {
         tbox::perr << "FAILED: - Test #8f: "
              << "Statistician::getGlobalProcStatMax()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8f successful" << endl;
      }
      if (statistician->getGlobalProcStatMaxProcessorId(
                  procstat1->getInstanceId(), 1) != tbox::MPI::getNodes()-1) {
         tbox::perr << "FAILED: - Test #8g: "
              << "Statistician::getGlobalProcStatMaxId()\n"
              << "returned incorrect ID." << endl;
      } else {
         tbox::plog << "Test #8g successful" << endl;
      }
      if (statistician->getGlobalProcStatMaxProcessorId(
                  procstat1->getInstanceId(), 2) != tbox::MPI::getNodes()-1) {
         tbox::perr << "FAILED: - Test #8h: "
              << "Statistician::getGlobalProcStatMaxId()\n"
              << "returned incorrect ID." << endl;
      } else {
         tbox::plog << "Test #8h successful" << endl;
      }

      double min0 = 0.;
      double min1 = 1.;
      double min2 = 2.;
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMin(
                        procstat1->getInstanceId(), 0), min0)) {
         tbox::perr << "FAILED: - Test #8i: "
              << "Statistician::getGlobalProcStatMin()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8i successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMin(
                        procstat1->getInstanceId(), 1), min1)) {
         tbox::perr << "FAILED: - Test #8j: "
              << "Statistician::getGlobalProcStatMin()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8j successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalProcStatMin(
                        procstat1->getInstanceId(), 2), min2)) {
         tbox::perr << "FAILED: - Test #8k: "
              << "Statistician::getGlobalProcStatMin()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8k successful" << endl;
      }

      if (statistician->getGlobalProcStatMinProcessorId(
                        procstat1->getInstanceId(), 1) != 0) {
          tbox::pout << statistician->getGlobalProcStatMinProcessorId(
                         procstat1->getInstanceId(), 1) << endl;
         tbox::perr << "FAILED: - Test #8l: "
              << "Statistician::getGlobalProcStatMinProcessorId()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8l successful" << endl;
      }
      if (statistician->getGlobalProcStatMinProcessorId(
                        procstat1->getInstanceId(), 2) != 0) {
          tbox::pout << statistician->getGlobalProcStatMinProcessorId(
                         procstat1->getInstanceId(), 2) << endl;
         tbox::perr << "FAILED: - Test #8m: "
              << "Statistician::getGlobalProcStatMinProcessorId()\n"
              << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #8m successful" << endl;
      }

      // Test #9:
      int pnum = 0;
      sum0 = 0.;
      sum1 = 0.;
      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         for (j = 0; j < i+2; j++) {
            sum0 += (double) pnum*2;
            sum1 += (double) (pnum*2 + 1); 
            pnum += 1;
         }
      }

      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatSum(
          patchstat2->getInstanceId(), 0), sum0)) {
          tbox::perr << "FAILED: - Test #9a: "
               << "Statistician::getGlobalPatchStatSum()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9a successful" << endl;
      }
      
      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatSum(
          patchstat2->getInstanceId(), 1), sum1)) {
          tbox::perr << "FAILED: - Test #9b: "
               << "Statistician::getGlobalPatchStatSum()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9b successful" << endl;
      }
 
   
      max0 = 0.;
      max1 = 0.;   
      pnum = -1;
      for (i = 0; i < tbox::MPI::getNodes(); i++) {
          for (j = 0; j < i+2; j++) {
              pnum += 1;
          }
      }

      max0 = (double) pnum*2;
      max1 = (double) (pnum*2 + 1);
      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatMax(
                  patchstat2->getInstanceId(), 0), max0)) {
          tbox::perr << "FAILED: - Test #9c: "
               << "Statistician::getGlobalPatchStatMax()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9c successful" << endl;
      }

      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatMax(
                  patchstat2->getInstanceId(), 1), max1)) {
          tbox::perr << "FAILED: - Test #9d: "
               << "Statistician::getGlobalPatchStatMax()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9d successful" << endl;
      }

      if (statistician->getGlobalPatchStatMaxPatchId(
                  patchstat2->getInstanceId(), 0) != pnum) {
          tbox::perr << "FAILED: - Test #9e: "
               << "Statistician::getGlobalPatchStatMaxPatchId()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9e successful" << endl;
      }
      
      if (statistician->getGlobalPatchStatMaxPatchId(
                  patchstat2->getInstanceId(), 1) != pnum) {
          tbox::perr << "FAILED: - Test #9f: "
               << "Statistician::getGlobalPatchStatMaxPatchId()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9f successful" << endl;
      }

      min0 = 0.;
      min1 = 1.;
      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatMin(
                  patchstat2->getInstanceId(), 0), min0)) {
          tbox::perr << "FAILED: - Test #9g: "
               << "Statistician::getGlobalPatchStatMin()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9g successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->getGlobalPatchStatMin(
                  patchstat2->getInstanceId(), 1), min1)) {
          tbox::perr << "FAILED: - Test #9h: "
               << "Statistician::getGlobalPatchStatMin()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9h successful" << endl;
      }
      
      if (statistician->getGlobalPatchStatMinPatchId(
                  patchstat2->getInstanceId(), 0) != 0) {
          tbox::perr << "FAILED: - Test #9i: "
               << "Statistician::getGlobalPatchStatMinPatchId()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9i successful" << endl;
      }
      if (statistician->getGlobalPatchStatMinPatchId(
                  patchstat2->getInstanceId(), 1) != 0) {
          tbox::perr << "FAILED: - Test #9j: "
               << "Statistician::getGlobalPatchStatMinPatchId()\n"
               << "incorrect value returned." << endl;
      } else {
          tbox::plog << "Test #9j successful" << endl;
      }

      // Test #10:
      pnum = 0;
      for (i = 0; i < tbox::MPI::getNodes(); i++) {
         max0 = 0.;
         max1 = 0.;
         for (j = 0; j < i+2; j++) {
            max0 += (double) pnum*2;
            max1 += (double) (pnum*2 + 1);
            pnum += 1;
         }
         if (!tbox::Utilities::deq(statistician->
                getGlobalPatchStatProcessorSum(
                   patchstat2->getInstanceId(), i, 0), max0)) {
               tbox::perr << "FAILED: - Test #10a: "
                    << "Statistician::getGlobalPatchStatProcessorSum()\n"
                    << "incorrect value returned." << endl;
         } else {
            tbox::plog << "Test #10a successful" << endl;
         }
         if (!tbox::Utilities::deq(statistician->
                getGlobalPatchStatProcessorSum(
                   patchstat2->getInstanceId(), i, 1), max1)) {
               tbox::perr << "FAILED: - Test #10b: "
                    << "Statistician::getGlobalPatchStatProcessorSum()\n"
                    << "incorrect value returned." << endl;
         } else {
            tbox::plog << "Test #10b successful" << endl;
         }
      }
      
      if (!tbox::Utilities::deq(statistician->
              getGlobalPatchStatProcessorSumMax(
                   patchstat2->getInstanceId(), 0), max0)) {
               tbox::perr << "FAILED: - Test #10c: "
                    << "Statistician::getGlobalPatchStatProcessorSumMax()\n"
                    << "incorrect value returned." << endl;
      } else {
        tbox::plog << "Test #10c successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->
              getGlobalPatchStatProcessorSumMax(
                   patchstat2->getInstanceId(), 1), max1)) {
               tbox::perr << "FAILED: - Test #10d: "
                    << "Statistician::getGlobalPatchStatProcessorSumMax()\n"
                    << "incorrect value returned." << endl;
      } else {
        tbox::plog << "Test #10d successful" << endl;
      }

      int nnodes = tbox::MPI::getNodes();
      if (statistician->getGlobalPatchStatProcessorSumMaxId(
                   patchstat2->getInstanceId(), 0) != nnodes-1) {
               tbox::perr << "FAILED: - Test #10e: "
                    << "Statistician::getGlobalPatchStatProcessorSumMaxId()\n"
                    << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #10e successful" << endl;
      }
      if (statistician->getGlobalPatchStatProcessorSumMaxId(
                   patchstat2->getInstanceId(), 1) != nnodes-1) {
               tbox::perr << "FAILED: - Test #10f: "
                    << "Statistician::getGlobalPatchStatProcessorSumMaxId()\n"
                    << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #10f successful" << endl;
      }

      min0 = 2.;
      min1 = 4.;
      if (!tbox::Utilities::deq(statistician->
              getGlobalPatchStatProcessorSumMin(
                   patchstat2->getInstanceId(), 0), min0)) {
               tbox::perr << "FAILED: - Test #10g: "
                    << "Statistician::getGlobalPatchStatProcessorSumMin()\n"
                    << "incorrect value returned." << endl;
      } else {
        tbox::plog << "Test #10g successful" << endl;
      }
      if (!tbox::Utilities::deq(statistician->
              getGlobalPatchStatProcessorSumMin(
                   patchstat2->getInstanceId(), 1), min1)) {
               tbox::perr << "FAILED: - Test #10h: "
                    << "Statistician::getGlobalPatchStatProcessorSumMin()\n"
                    << "incorrect value returned." << endl;
      } else {
        tbox::plog << "Test #10h successful" << endl;
      }

      if (statistician->getGlobalPatchStatProcessorSumMinId(
                   patchstat2->getInstanceId(), 0) != 0) {
               tbox::perr << "FAILED: - Test #10i: "
                    << "Statistician::getGlobalPatchStatProcessorSumMinId()\n"
                    << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #10i successful" << endl;
      }
      if (statistician->getGlobalPatchStatProcessorSumMinId(
                   patchstat2->getInstanceId(), 1) != 0) {
               tbox::perr << "FAILED: - Test #10j: "
                    << "Statistician::getGlobalPatchStatProcessorSumMinId()\n"
                    << "incorrect value returned." << endl;
      } else {
         tbox::plog << "Test #10j successful" << endl;
      }

   }

   /*
    * We're done.  Write the restart file.  
    */
   string restart_write_dirname = "restart";
#ifdef HAVE_HDF5
   int timestep = 0;
   tbox::RestartManager::getManager()->writeRestartFile(
               restart_write_dirname,
               timestep);
#endif

   tbox::pout << "\nPASSED:  statstest" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();
   return(0);
}

