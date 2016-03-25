//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/test/applications/AutoTester.C $
// Package:     SAMRAI applications
// Copyright:   (c) 1997-2002 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2043 $
// Modified:    $LastChangedDate: 2008-03-12 09:14:32 -0700 (Wed, 12 Mar 2008) $
// Description: Class used for auto testing applications 
//

#include "AutoTester.h"

#include "Box.h"
#include "BoxArray.h"
#include "BoxList.h"
#include "PatchLevel.h"
#include "tbox/List.h"
#include "tbox/MathUtilities.h"


AutoTester::AutoTester(const string& object_name,
                       tbox::Pointer<tbox::Database> input_db)
{
   d_object_name    = object_name;
   d_test_fluxes    = false;
   d_test_iter_num  = 10;
   d_output_correct = false;

   d_write_patch_boxes = false;
   d_read_patch_boxes  = false;
   d_test_patch_boxes_at_steps.resizeArray(0);
   d_simplify_test_boxes = false;

   d_test_patch_boxes_step_count = 0;

   getFromInput(input_db);

   if (d_read_patch_boxes) {
      d_io_box_utility =
         new hier::BoxIOUtility<NDIM>(d_test_patch_boxes_filename,
                          hier::BoxIOUtility<NDIM>::READ);
      if (d_output_correct) {
         d_io_box_utility->printBoxes(tbox::pout);
      }
      
   }

   if (d_write_patch_boxes) {
      d_io_box_utility =
         new hier::BoxIOUtility<NDIM>(d_test_patch_boxes_filename,
                          hier::BoxIOUtility<NDIM>::WRITE);
   }

}

AutoTester::~AutoTester()
{
   if (d_write_patch_boxes) {
      d_io_box_utility->writeLevelBoxesDatabase();
   }

   if (d_io_box_utility )
      delete d_io_box_utility;
}

/*
******************************************************************
*                                                                *
*  Method "evalTestData" compares the result of the run with     *
*  the correct result for runs with the TimeRefinementIntegrator *
*  and HyperbolicLevelIntegrator.                                *
*                                                                *
******************************************************************
*/
int AutoTester::evalTestData( 
   int iter, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   const tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > tri, 
   const tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> > hli, 
   const tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > ga)
{

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified 
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) { 

      /*
       * set precision of output stream.
       */
      tbox::pout.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      tbox::Pointer<hier::PatchLevel<NDIM> > level = 
	 hierarchy->getPatchLevel(nlevels);

      /*
       * Test 0: Time Refinement Integrator
       */
      double time = tri->getIntegratorTime();
      if (d_correct_result.getSize() > 0) {
         if (d_output_correct) {
            tbox::pout << "Test 0: Time Refinement Integrator "
                 << "\n   computed result: " << time;
            
            tbox::pout  << "\n   specified result = " 
                  << d_correct_result[0];
         }
         tbox::pout << endl;
      
         if (tbox::MathUtilities<double>::equalEps(time,d_correct_result[0]) ) {
            tbox::pout << "Test 0: Time Refinement check successful" << endl;
         } else {
            tbox::perr << "Test 0 FAILED: Check Time Refinement Integrator" << endl;
            num_failures++;
         }
      }
      

      /*
       * Test 1: Time Refinement Integrator
       */
      double dt = tri->getLevelDtMax(nlevels);
      if (d_correct_result.getSize() > 1) {
         if (d_output_correct) {
            tbox::pout << "Test 1: Time Refinement Integrator "
                 << "\n   computed result: " << dt;
            tbox::pout  << "\n   specified result = " 
                  << d_correct_result[1];
         }
         tbox::pout << endl;
      
         if (tbox::MathUtilities<double>::equalEps(dt,d_correct_result[1]) ) {
            tbox::pout << "Test 1: Time Refinement check successful" << endl;
         } else {
            tbox::perr << "Test 1 FAILED: Check Time Refinement Integrator" << endl;
            num_failures++;
         }
      }

 
      /*
       * Test 2: Hyperbolic Level Integrator
       */
      dt = hli->getLevelDt(level,time,false);
      if (d_correct_result.getSize() > 2) {
         if (d_output_correct) {
            tbox::pout << "Test 2: Hyperbolic Level Integrator "
                 << "\n   computed result: " << dt;
            
            tbox::pout  << "\n   specified result = " 
                  << d_correct_result[2];
	 }
	 tbox::pout << endl;
      
         if (tbox::MathUtilities<double>::equalEps(dt,d_correct_result[2]) ) {
            tbox::pout << "Test 2: Hyperbolic Level Int check successful" << endl;
         } else {
            tbox::perr << "Test 2 FAILED: Check Hyperbolic Level Integrator" << endl;
            num_failures++;
         }
      }
      

      /*
       * Test 3: Gridding Algorithm
       */
      int n = ga->getMaxLevels();
      if (d_output_correct) {
         tbox::pout << "Test 3: Gridding Algorithm "
              << "\n   computed result: " << n;
         tbox::pout << "\n   correct result = " << nlevels+1;
         tbox::pout << endl; 
      }
      if (n == (nlevels+1) ) {
         tbox::pout << "Test 3: Gridding Algorithm check successful" << endl;
      } else {
         tbox::perr << "Test 3 FAILED: Check Gridding Algorithm" << endl;
         num_failures++;
      }

   }

   if ( (d_test_patch_boxes_at_steps.getSize() > d_test_patch_boxes_step_count)     && (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] == iter) ) 
   {

      int num_levels = hierarchy->getNumberOfLevels();

      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_io_box_utility->printBoxes(tbox::pout);
         }

         for (int ln = 0; ln < num_levels; ln++) {
            hier::BoxArray<NDIM> test_boxes;
            d_io_box_utility->getLevelBoxes(test_boxes,
                                            ln,
                                            d_test_patch_boxes_step_count);

            num_failures +=checkHierarchyBoxes(hierarchy,
                                               ln,
                                               test_boxes,
                                               iter);
         }


      }

      if (d_write_patch_boxes) {

         for (int ln = 0; ln < num_levels; ln++) {
	    tbox::Pointer<hier::PatchLevel<NDIM> > level 
	       = hierarchy->getPatchLevel(ln);
            d_io_box_utility->
               putLevelBoxes(level->getBoxes(),
			     ln,
                             d_test_patch_boxes_step_count);
         }

         if (d_output_correct) {
            d_io_box_utility->printBoxes(tbox::pout);
         }
      }

      d_test_patch_boxes_step_count++;

   }

   return(num_failures);
}

/*
******************************************************************
*                                                                *
*  Method "evalTestData" compares the result of the run with     *
*  the correct result for runs with the MethodOfLinesIntegrator. *
*                                                                *
******************************************************************
*/
int AutoTester::evalTestData(
   int iter,
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   double time,
   const tbox::Pointer<algs::MethodOfLinesIntegrator<NDIM> > mol,
   const tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > ga)
{

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) {

      /*
       * set precision of output stream.
       */
      tbox::pout.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      tbox::Pointer<hier::PatchLevel<NDIM> > level =
         hierarchy->getPatchLevel(nlevels);

      /*
       * Test 0: Time test
       */
      if (d_output_correct) {
         tbox::pout << "Test 0: Simulation Time: "
              << "\n   computed result: " << time;
         if (d_correct_result.getSize() > 0) {
            tbox::pout  << "\n   specified result = "
                  << d_correct_result[0];
         }
         tbox::pout << endl;
      }
      if (tbox::MathUtilities<double>::equalEps(time,d_correct_result[0]) ) {
         tbox::pout << "Test 0: Simulation Time check successful" << endl;
      } else {
         tbox::perr << "Test 0 FAILED: Simulation time incorrect" << endl;
         num_failures++;
      }

      /*
       * Test 1: MethodOfLinesIntegrator
       */
      double dt = mol->getTimestep(hierarchy, time);
      if (d_output_correct) {
         tbox::pout << "Test 1: Method of Lines Integrator "
              << "\n   computed result: " << dt;
         if (d_correct_result.getSize() > 1) {
            tbox::pout  << "\n   specified result = "
                  << d_correct_result[1];
         }
         tbox::pout << endl;
      }
      if (tbox::MathUtilities<double>::equalEps(dt,d_correct_result[1]) ) {
         tbox::pout << "Test 1: MOL Int check successful" << endl;
      } else {
         tbox::perr << "Test 1 FAILED: Check Method of Lines Integrator" << endl;
         num_failures++;
      }

      /*
       * Test 2: Gridding Algorithm
       */
      int n = ga->getMaxLevels();
      if (d_output_correct) {
         tbox::pout << "Test 2: Gridding Algorithm "
              << "\n   computed result: " << n;
         tbox::pout << "\n   correct result = " << nlevels+1;
         tbox::pout << endl;
      }
      if (n == (nlevels+1) ) {
         tbox::pout << "Test 2: Gridding Alg check successful" << endl;
      } else {
         tbox::perr << "Test 2 FAILED: Check Gridding Algorithm" << endl;
         num_failures++;
      }

   }

   if ( (d_test_patch_boxes_at_steps.getSize() > 0) &&
        (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] ==
         iter) ) {

      int num_levels = hierarchy->getNumberOfLevels();

      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_io_box_utility->printBoxes(tbox::pout);
         }

         for (int ln = 0; ln < num_levels; ln++) {
            hier::BoxArray<NDIM> test_boxes;
            d_io_box_utility->getLevelBoxes(test_boxes,
                                            ln,
                                            d_test_patch_boxes_step_count);

            num_failures += checkHierarchyBoxes(hierarchy,
                                                ln,
                                                test_boxes,
                                                iter);
         }


      }

      if (d_write_patch_boxes) {

         if (d_output_correct) {
            d_io_box_utility->printBoxes(tbox::pout);
         }

         for (int ln = 0; ln < num_levels; ln++) {
            tbox::Pointer<hier::PatchLevel<NDIM> > level 
	       = hierarchy->getPatchLevel(ln);
            d_io_box_utility->
               putLevelBoxes(level->getBoxes(),
                             ln,
                             d_test_patch_boxes_step_count);
         }

      }

      d_test_patch_boxes_step_count++;

   }

   return(num_failures);
}

/*
******************************************************************
*                                                                *
*  Get test parameters from input.                               *
*                                                                *
******************************************************************
*/

void AutoTester::getFromInput(tbox::Pointer<tbox::Database> input_db)
{
   tbox::Pointer<tbox::Database> tester_db = 
      input_db->getDatabase(d_object_name);

   /* 
    * Read testing parameters from testing_db 
    */
   if (tester_db->keyExists("test_fluxes")) {
      d_test_fluxes = tester_db->getBool("test_fluxes");
   }

   if (tester_db->keyExists("test_iter_num")) {
      d_test_iter_num = tester_db->getInteger("test_iter_num");
   }

   if (tester_db->keyExists("write_patch_boxes")) {
      d_write_patch_boxes = tester_db->getBool("write_patch_boxes");
   }
   if (tester_db->keyExists("read_patch_boxes")) {
      d_read_patch_boxes = tester_db->getBool("read_patch_boxes");
   }
   if (d_read_patch_boxes && d_write_patch_boxes) {
      tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
           << "Cannot 'read_patch_boxes' and 'write_patch_boxes' \n"
           << "at the same time." << endl;
   }
   if (d_read_patch_boxes || d_write_patch_boxes) {
      if (!tester_db->keyExists("test_patch_boxes_at_steps")) {
         tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
              << "Must provide 'test_patch_boxes_at_steps' data." << endl;
      } else {
         d_test_patch_boxes_at_steps = 
            tester_db->getIntegerArray("test_patch_boxes_at_steps");
      }
      if (!tester_db->keyExists("test_patch_boxes_filename")) {
         tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
              << "Must provide 'test_patch_boxes_filename' data." << endl;
      } else {
         d_test_patch_boxes_filename =
            tester_db->getString("test_patch_boxes_filename");
      }
      if (tester_db->keyExists("simplify_test_boxes")) {
         d_simplify_test_boxes = tester_db->getBool("simplify_test_boxes");
      }
   }

   if (d_test_fluxes) {

      /* 
       * Read expected result for flux test... 
       * Fluxes not verified in this routine.  Rather, we let it
       * write the result and do a "diff" within the script
       */

     tbox::pout << "Do a diff on the resulting *.dat file to verify result."  
          << endl;

   } else {

      /* 
       * Read correct_result array for timestep test... 
       */
      if (tester_db->keyExists("correct_result")) {
         d_correct_result = tester_db->getDoubleArray("correct_result");
      } else {
         TBOX_WARNING("main.C: TESTING is on but no `correct_result' array"
               << "is given in input file." << endl);
      }

      /* Specify whether to output "correct_result" result */
 
      if (tester_db->keyExists("output_correct")) {
         d_output_correct = tester_db->getBool("output_correct");
      }

   }

}

static void getBoxListsSortedBySize(tbox::Array<hier::BoxList<NDIM> >& sorted_boxarrays,
                                    hier::BoxList<NDIM>& in_boxlist) 
{
   in_boxlist.sortDescendingBoxSizes(); 
   
   tbox::List<int> num_box_sizes;
   int check_size = 0;
   if (in_boxlist.getNumberOfItems() > 0) {
      check_size = in_boxlist.getFirstItem().size();
   }
   num_box_sizes.appendItem(0);
   for (hier::BoxList<NDIM>::Iterator lin(in_boxlist); lin; lin++) {
      if (lin().size() == check_size) {
         (num_box_sizes.getLastItem())++;
      } else {
         check_size = lin().size();
         num_box_sizes.appendItem(1);
      }
   }

   sorted_boxarrays.resizeArray(num_box_sizes.getNumberOfItems());

   hier::BoxList<NDIM>::Iterator lin2s(in_boxlist);
   int isort = 0;
   for (tbox::List<int>::Iterator lbs(num_box_sizes); lbs; lbs++) {
      int num_boxes = lbs();
      for (int ib = 0; ib < num_boxes; ib++) {
         hier::Box<NDIM> check_box = lin2s();
         sorted_boxarrays[isort].appendItem(check_box);
         lin2s++;
      }
      isort++;
   } 
}

int AutoTester::checkHierarchyBoxes(
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   const hier::BoxArray<NDIM>& test_boxes,
   int iter)
{
 
   hier::BoxList<NDIM> master_boxes(test_boxes);

   tbox::Pointer<hier::PatchLevel<NDIM> > level 
      = hierarchy->getPatchLevel(level_number);
   hier::BoxList<NDIM> hier_boxes(level->getBoxes());

   if (d_simplify_test_boxes) {
      master_boxes.simplifyBoxes();
      hier_boxes.simplifyBoxes();
   }

   tbox::Array<hier::BoxList<NDIM> > sorted_master_boxes;
   getBoxListsSortedBySize(sorted_master_boxes, master_boxes);

   tbox::Array<hier::BoxList<NDIM> > sorted_boxes;
   getBoxListsSortedBySize(sorted_boxes, hier_boxes);

   int i;

   /*
    * Check to make sure sorted_boxes and sorted_master_boxes are
    * identical.  If not, write an error message.
    */

   bool total_match = 
      (hier_boxes.getNumberOfItems() == master_boxes.getNumberOfItems());

   bool bin_match = 
      (sorted_boxes.getSize() == sorted_master_boxes.getSize());
   
   bool bin_count_match = true;
   if (total_match && bin_match) {
      for (i = 0; i < sorted_boxes.getSize(); i++) {
         bin_count_match = bin_count_match &&
            (sorted_boxes[i].getNumberOfItems() ==
             sorted_master_boxes[i].getNumberOfItems());
      }
   } 

   int num_failures = 0;
   
   if (total_match && bin_match && bin_count_match) {
      tbox::pout << "Test 4: Level " << level_number
           << " Boxes check successful for step " << iter
           << endl << endl;
   } else {
      tbox::perr << "Test 4: FAILED: Level " << level_number
           << " hier::Box configuration doesn't match at step " << iter 
           << endl << endl;
      num_failures++;
   }
  
   if (d_output_correct) {

      tbox::pout << "-------------------------------------------------------" 
           << endl;

      if (!total_match) {
         tbox::pout << "TOTAL NUMBER OF BOXES DOES NOT MATCH " 
              << "ON LEVEL: " << level_number << endl;
      }
      tbox::pout << "Total number of boxes: " << endl;
      tbox::pout << "   hier::PatchLevel boxes -- " << hier_boxes.getNumberOfItems()
           << "   Test boxes -- " << master_boxes.getNumberOfItems() 
           << endl;

      if (!bin_match) {
         tbox::pout << "NUMBER OF BOX SIZES DOES NOT MATCH " 
              << "ON LEVEL: " << level_number << endl;
      }
      tbox::pout << "Number of box sizes: " << endl;
      tbox::pout << "   hier::PatchLevel boxes -- " << sorted_boxes.getSize()
           << "   Test boxes -- " << sorted_master_boxes.getSize()
           << endl;

      if (!bin_count_match) {
         tbox::pout << "BOX SIZE DISTRIBUTION DOES NOT MATCH!" << endl;
      }
      for (i = 0; i < sorted_boxes.getSize(); i++) {
         tbox::pout << "PatchLevel has " 
              << sorted_boxes[i].getNumberOfItems()
              << " boxes";
         if (sorted_boxes[i].getNumberOfItems() > 0) {
            tbox::pout << " of size = " 
                 << sorted_boxes[i].getFirstItem().size();
         }
         tbox::pout << endl;
      }
      for (i = 0; i < sorted_master_boxes.getSize(); i++) {
         tbox::pout << "Test boxes have "
              << sorted_master_boxes[i].getNumberOfItems()
             << " boxes";
         if (sorted_master_boxes[i].getNumberOfItems() > 0) {
            tbox::pout << " of size = " 
                 << sorted_master_boxes[i].getFirstItem().size();
         }
         tbox::pout << endl;
      }

      tbox::pout << "-------------------------------------------------------"
           << endl << endl;

   }

   return (num_failures);
}
