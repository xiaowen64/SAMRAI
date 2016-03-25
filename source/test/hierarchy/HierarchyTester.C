//
// File:        HierarchyTester.C
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 414 $
// Modified:    $Date: 2005-06-01 16:23:06 -0700 (Wed, 01 Jun 2005) $
// Description: Manager class for patch hierarchy refine/coarsen tests.
//

#include "HierarchyTester.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "tbox/Utilities.h"
#include "BergerRigoutsos.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "GriddingAlgorithm.h"
#include "GridGeometry.h"
#include "LoadBalancer.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "StandardTagAndInitialize.h"

using namespace geom;

namespace SAMRAI {

/*
*************************************************************************
*									*
* The constructor initializes object state.                             *
* The destructor deallocates object storage.                            *
*									*
*************************************************************************
*/

HierarchyTester::HierarchyTester(
   const string& object_name,
   Pointer<Database> hier_test_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!hier_test_db.isNull());
#endif

   d_object_name = object_name;

   d_do_refine_test  = false;
   d_do_coarsen_test = false;

   d_ratio = IntVector<NDIM>(0);

   d_initial_patch_hierarchy.setNull();
   d_test_patch_hierarchy.setNull();

   if (hier_test_db->keyExists("do_refine_test")) {
      d_do_refine_test = hier_test_db->getBool("do_refine_test");
      if (d_do_refine_test) {
         tbox::plog << "\nPerforming hierarchy refine test..." << endl;
         if (hier_test_db->keyExists("ratio")) {
            int* tmp_ratio = d_ratio;
            hier_test_db->getIntegerArray("ratio", tmp_ratio, NDIM);
            tbox::plog << "with ratio = " << d_ratio << endl;
         } else {
            TBOX_ERROR(
            "HierarchyTester input error: no 'ratio' found in input" << endl);
         }
      }
   }

   if (!d_do_refine_test) {
      if (hier_test_db->keyExists("do_coarsen_test")) {
         d_do_coarsen_test = hier_test_db->getBool("do_coarsen_test");
      }
      if (d_do_coarsen_test) {
         tbox::plog << "\nPerforming hierarchy coarsen test..." << endl;
         if (hier_test_db->keyExists("ratio")) {
            int* tmp_ratio = d_ratio;
            hier_test_db->getIntegerArray("ratio", tmp_ratio, NDIM);
            tbox::plog << "with ratio = " << d_ratio << endl;
         } else {
            TBOX_ERROR(
            "HierarchyTester input error: no 'ratio' found in input" << endl);
         }
      }
   }

   if (!d_do_refine_test && !d_do_coarsen_test) {
      TBOX_ERROR(
      "HierarchyTester input error: no test specified in input" << endl);
   }

}

HierarchyTester::~HierarchyTester()
{
   d_initial_patch_hierarchy.setNull();
   d_test_patch_hierarchy.setNull();
}

/*
*************************************************************************
*                                                                       *
* Create and configure gridding objects used to build the hierarchy.    *
* Then, create hierarchy and set up dummy data for patch descriptor.    *
*                                                                       *
*************************************************************************
*/

void HierarchyTester::setupInitialHierarchy(
   Pointer<Database> main_input_db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!main_input_db.isNull());
#endif

   Pointer<CartesianGridGeometry<NDIM> > grid_geometry =
      new CartesianGridGeometry<NDIM>(
          "CartesianGridGeometry",
          main_input_db->getDatabase("CartesianGridGeometry"));

   d_initial_patch_hierarchy =
      new PatchHierarchy<NDIM>("InitialPatchHierarchy",
                          grid_geometry);

   Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

   Pointer<LoadBalancer<NDIM> > load_balancer = 
      new LoadBalancer<NDIM>("LoadBalancer", 
                       main_input_db->getDatabase("LoadBalancer"));

   Pointer<StandardTagAndInitialize<NDIM> > dummy_error_detector =
      new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
          this,
          main_input_db->getDatabase("StandardTagAndInitialize"));

   Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
      new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                            main_input_db->getDatabase("GriddingAlgorithm"),
                            dummy_error_detector,
                            box_generator,
                            load_balancer);

   gridding_algorithm->makeCoarsestLevel(d_initial_patch_hierarchy, 
                                         0.0);  // dummy time

   for (int ln = 0; gridding_algorithm->levelCanBeRefined(ln); ln++) {
       gridding_algorithm->makeFinerLevel(d_initial_patch_hierarchy,
                                          0.0,   // dummy time
                                          true,  // indicates initial time 
                                          0);    // dummy tag buffer
   }

}

void HierarchyTester::runHierarchyTestAndVerify()
{
   if (d_do_refine_test) {
      d_test_patch_hierarchy = 
         d_initial_patch_hierarchy->makeRefinedPatchHierarchy(
            "FinePatchHierarchy",
            d_ratio,
            false);
   }

   if (d_do_coarsen_test) {
      d_test_patch_hierarchy = 
         d_initial_patch_hierarchy->makeCoarsenedPatchHierarchy(
            "CoarsePatchHierarchy",
            d_ratio,
            false);
   }

   /*
    **************************************************************
    * Tests 0a-0c check global data in grid geometry.            *
    **************************************************************
    */

   Pointer<GridGeometry<NDIM> > init_geometry = 
      d_initial_patch_hierarchy->getGridGeometry();
   Pointer<GridGeometry<NDIM> > test_geometry = 
      d_test_patch_hierarchy->getGridGeometry();

   // Test #0a:
   if ( init_geometry->getPeriodicShift() != 
        test_geometry->getPeriodicShift() ) {
      tbox::perr << "FAILED: - Test #0a: initial hierarchy has periodic shift " 
           << init_geometry->getPeriodicShift() << " and \n" 
           << "test hierarchy has periodic shift "
           << test_geometry->getPeriodicShift() << endl;
   }

   const BoxArray<NDIM>& init_phys_domain = init_geometry->getPhysicalDomain();
   const BoxArray<NDIM>& test_phys_domain = test_geometry->getPhysicalDomain();

   const int npdboxes = init_phys_domain.getNumberOfBoxes();

   // Test #0b:
   if (d_do_refine_test) {
      for (int ib = 0; ib < npdboxes; ib++) {
         if ( Box<NDIM>::refine(init_phys_domain.getBox(ib), d_ratio) !=
              test_phys_domain.getBox(ib) ) {
            tbox::perr << "FAILED: - Test #0b: test hierarchy physical domain"
                 << " box with array index " << ib
                 << " is not a proper refinement of initial hierarchy"
                 << " physical domain box with same index" << endl;
         }
      }
   }
   if (d_do_coarsen_test) {
      for (int ib = 0; ib < npdboxes; ib++) {
         if ( Box<NDIM>::coarsen(init_phys_domain.getBox(ib), d_ratio) !=
              test_phys_domain.getBox(ib) ) {
            tbox::perr << "FAILED: - Test #0b: test hierarchy physical domain"
                 << " box with array index " << ib
                 << " is not a proper coarsening of initial hierarchy"
                 << " physical domain box with same index" << endl;
         }
      }
   }

   // Test #0c:
   if ( init_geometry->getDomainIsSingleBox() !=
        test_geometry->getDomainIsSingleBox() ) {
      tbox::perr << "FAILED: - Test #0c: initial and test hierarchy do not match"
           << " for geom->getDomainIsSingleBox()" << endl;
   }

   /*
    **************************************************************
    * Tests 1-8 check global data for levels in hierarchy.       *
    **************************************************************
    */

   const int nlevels = d_initial_patch_hierarchy->getNumberLevels();

   // Test #1:
   if (d_test_patch_hierarchy->getNumberLevels() != nlevels) {
      tbox::perr << "FAILED: - Test #1: initial hierarchy has " 
           << nlevels << " levels and \n" 
           << "test hierarchy has " 
           << d_test_patch_hierarchy->getNumberLevels() << "levels" << endl;
   }

   for (int ln = 0; ln < nlevels; ln++) {
      Pointer<PatchLevel<NDIM> > init_level = 
         d_initial_patch_hierarchy->getPatchLevel(ln); 
      Pointer<PatchLevel<NDIM> > test_level = 
         d_test_patch_hierarchy->getPatchLevel(ln); 

      // Test #2:
      if (init_level->getLevelNumber() !=
          test_level->getLevelNumber()) {
         tbox::perr << "FAILED: - Test #2: for level number " << ln
              << " initial hierarchy level number is "
              << init_level->getLevelNumber() 
              << "\n and test hierarchy level number is " 
              << test_level->getLevelNumber() << endl;
      }

      // Test #3:
      if (init_level->getNextCoarserHierarchyLevelNumber() !=
          test_level->getNextCoarserHierarchyLevelNumber()) {
         tbox::perr << "FAILED: - Test #3: for level number " << ln
              << " initial hierarchy next coarser level number is "
              << init_level->getNextCoarserHierarchyLevelNumber() 
              << "\n and test hierarchy next coarser level number is " 
              << test_level->getNextCoarserHierarchyLevelNumber() << endl;
      }

      // Test #4:
      if (init_level->inHierarchy() !=
          test_level->inHierarchy()) {
         tbox::perr << "FAILED: - Test #4: for level number " << ln
              << " initial hierarchy level in hierarchy is "
              << init_level->inHierarchy()
              << "\n and test hierarchy level in hierarchy is "
              << test_level->inHierarchy() << endl;
      }

      // Test #5:
      if (init_level->getNumberOfPatches() !=
          test_level->getNumberOfPatches()) {
         tbox::perr << "FAILED: - Test #5: for level number " << ln
              << " initial hierarchy number of patches is "
              << init_level->getNumberOfPatches()
              << "\n and test hierarchy number of patches is "
              << test_level->getNumberOfPatches() << endl;
      }

      // Test #6:
      if (init_level->getRatio() !=
          test_level->getRatio()) {
         tbox::perr << "FAILED: - Test #6: for level number " << ln
              << " initial hierarchy ratio to level zero is "
              << init_level->getRatio()
              << "\n and test hierarchy ratio to level zero is "
              << test_level->getRatio() << endl;
      }

      // Test #7:
      if (init_level->getRatioToCoarserLevel() !=
          test_level->getRatioToCoarserLevel()) {
         tbox::perr << "FAILED: - Test #7: for level number " << ln
              << " initial hierarchy ratio to coarser level is "
              << init_level->getRatioToCoarserLevel()
              << "\n and test hierarchy ratio to coarser level is "
              << test_level->getRatioToCoarserLevel() << endl;
      }

      const BoxArray<NDIM>& init_domain = init_level->getPhysicalDomain(); 
      const BoxArray<NDIM>& test_domain = test_level->getPhysicalDomain(); 

      const int nboxes = init_domain.getNumberOfBoxes();

      // Test #8:
      if (d_do_refine_test) { 
         for (int ib = 0; ib < nboxes; ib++) {
            if ( Box<NDIM>::refine(init_domain.getBox(ib), d_ratio) !=
                 test_domain.getBox(ib) ) {
               tbox::perr << "FAILED: - Test #8: for level number " << ln
                    << " refined domain box with array index " << ib 
                    << " is not a proper refinement of initial domain "
                    << "box with same index" << endl;
            }
         }
      }
      if (d_do_coarsen_test) {
         for (int ib = 0; ib < nboxes; ib++) {
            if ( Box<NDIM>::coarsen(init_domain.getBox(ib), d_ratio) !=
                 test_domain.getBox(ib) ) {
               tbox::perr << "FAILED: - Test #8: for level number " << ln
                    << " coarsened domain box with array index " << ib 
                    << " is not a proper coarsening of initial domain "
                    << "box with same index" << endl;
            }
         }
      }

     /*
      **************************************************************
      *  Tests 9-13 check global data for patches on each level.   *
      **************************************************************
      */

      const BoxArray<NDIM>& init_boxes = init_level->getBoxes(); 
      const BoxArray<NDIM>& test_boxes = test_level->getBoxes(); 

      const int npatches = init_level->getNumberOfPatches();
      for (int ip = 0; ip < npatches; ip++) {

         // Test #9:
         if (d_do_refine_test) {
            if ( Box<NDIM>::refine(init_boxes.getBox(ip), d_ratio) !=
                 test_boxes.getBox(ip) ) {
               tbox::perr << "FAILED: - Test #9: for level number " << ln
                    << " refined patch box with array index " << ip
                    << " is not a proper refinement of initial domain "
                    << "box with same index" << endl;
            }
         }
         if (d_do_coarsen_test) {
            if ( Box<NDIM>::coarsen(init_boxes.getBox(ip), d_ratio) !=
                 test_boxes.getBox(ip) ) {
               tbox::perr << "FAILED: - Test #9: for level number " << ln
                    << " coarsened patch box with array index " << ip
                    << " is not a proper coarsening of initial domain "
                    << "box with same index" << endl;
            }
         }
   
         // Test #10:
         if ( (init_level->getShiftsForPatch(ip)).getNumberItems() !=
              (test_level->getShiftsForPatch(ip)).getNumberItems() ) {
            tbox::perr << "FAILED: - Test #10: for level number " << ln
                 << " initial and test level have different number of "
                 << "shifts for patch number " << ip << endl;
         }

         // Test #11:
         if ( init_level->getMappingForPatch(ip) !=
              test_level->getMappingForPatch(ip) ) {
            tbox::perr << "FAILED: - Test #11: for level number " << ln
                 << " initial and test level have different processor "
                 << "mapping for patch number " << ip << endl;
         }

         // Test #12:
         if ( init_level->patchTouchesRegularBoundary(ip) !=
              test_level->patchTouchesRegularBoundary(ip) ) {
            tbox::perr << "FAILED: - Test #12: for level number " << ln
                 << " initial and test level do not match for "
                 << "patchTouchesRegularBoundary() "
                 << "for patch number " << ip << endl;
         }

         // Test #13:
         if ( init_level->patchTouchesPeriodicBoundary(ip) !=
              test_level->patchTouchesPeriodicBoundary(ip) ) {
            tbox::perr << "FAILED: - Test #12: for level number " << ln
                 << " initial and test level do not match for "
                 << "patchTouchesPeriodicBoundary() "
                 << "for patch number " << ip << endl;
         }

      }

     /*
      **************************************************************
      *  Tests 14-20 check local data for patches on each level.   *
      **************************************************************
      */
     for (PatchLevel<NDIM>::Iterator tip(init_level); tip; tip++) {
         Pointer<Patch<NDIM> > init_patch = init_level->getPatch(tip()); 
         Pointer<Patch<NDIM> > test_patch = test_level->getPatch(tip()); 

         // Test #14:
         if (d_do_refine_test) {
            if ( Box<NDIM>::refine(init_patch->getBox(), d_ratio) !=
                 test_patch->getBox() ) {
               tbox::perr << "FAILED: - Test #14: for level number " << ln
                    << " box for test level patch " << tip()
                    << " is not a proper refinement of box "
                    << "for initial level patch with same number" << endl;
            }
         }
         if (d_do_coarsen_test) {
            if ( Box<NDIM>::coarsen(init_patch->getBox(), d_ratio) !=
                 test_patch->getBox() ) {
               tbox::perr << "FAILED: - Test #14: for level number " << ln
                    << " box for test level patch " << tip()
                    << " is not a proper coarsening of box "
                    << "for initial level patch with same number" << endl;
            }
         }

         // Test #15:
         if ( init_patch->getPatchNumber() != 
              test_patch->getPatchNumber() ) {
            tbox::perr << "FAILED: - Test #15: for level number " << ln
                 << " initial and test level patches have different patch "
                 << "numbers for patch with index " << tip() << endl;
         }

         // Test #16:
         if ( init_patch->getPatchLevelNumber() != 
              test_patch->getPatchLevelNumber()) {
            tbox::perr << "FAILED: - Test #16: for level number " << ln
                 << " initial and test level patches have different patch "
                 << "level numbers for patch number " << tip() << endl;
         }

         // Test #17:
         if ( init_patch->inHierarchy() != 
              test_patch->inHierarchy()) {
            tbox::perr << "FAILED: - Test #17: for level number " << ln
                 << " initial and test level do not match for "
                 << "inHierarchy() "
                 << "for patch number " << tip() << endl;
         }

         // Test #18:
         if ( init_patch->getPatchGeometry()->getTouchesRegularBoundary() != 
              test_patch->getPatchGeometry()->getTouchesRegularBoundary()) {
            tbox::perr << "FAILED: - Test #18: for level number " << ln
                 << " initial and test level do not match for "
                 << "getTouchesRegularBoundary() "
                 << "for patch number " << tip() << endl;
         }

         // Test #19:
         if ( init_patch->getPatchGeometry()->getTouchesPeriodicBoundary() != 
              test_patch->getPatchGeometry()->getTouchesPeriodicBoundary()) {
            tbox::perr << "FAILED: - Test #19: for level number " << ln
                 << " initial and test level do not match for "
                 << "getTouchesPeriodicBoundary() "
                 << "for patch number " << tip() << endl;
         }

         /*
          **************************************************************
          * Tests 20a-20c check patch geometry data.                   *
          **************************************************************
          */

         Pointer<PatchGeometry<NDIM> > init_patch_geom =
            init_patch->getPatchGeometry();
         Pointer<PatchGeometry<NDIM> > test_patch_geom =
            test_patch->getPatchGeometry();

         // Test #20a:
         if ( init_patch_geom->getRatio() != 
              test_patch_geom->getRatio()) {
            tbox::perr << "FAILED: - Test #20a: for level number " << ln
                 << " patch geometry ratio data does not match "
                 << "for patch number " << tip() << endl;
         }

         // Test #20b:
         if ( init_patch_geom->intersectsPhysicalBoundary() != 
              test_patch_geom->intersectsPhysicalBoundary()) {
            tbox::perr << "FAILED: - Test #20b: for level number " << ln
                 << " intersectsPhysicalBoundary() does not match "
                 << "for patch number " << tip() << endl;
         }

         // Test #20c:
         for (int id = 1; id <= NDIM; id++) {
            if ( (init_patch_geom->getCodimensionBoundary(id)).getSize() !=
                 (test_patch_geom->getCodimensionBoundary(id)).getSize() ) {
               tbox::perr << "FAILED: - Test #20c: for level number " << ln
                    << " number of codimension " << id 
                    << " boundary boxes does not match "
                    << "for patch number " << tip() << endl;
            }
         }

      }
        
   }
}


}
