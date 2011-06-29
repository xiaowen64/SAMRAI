/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Manager class for patch hierarchy refine/coarsen tests. 
 *
 ************************************************************************/

#include "HierarchyTester.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/hier/GridGeometry.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"

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
   const std::string& object_name,
   const tbox::Dimension& dim,
   Pointer<Database> hier_test_db):
   d_dim(dim),
   d_ratio(dim, 0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!hier_test_db.isNull());
#endif

   d_object_name = object_name;

   d_do_refine_test = false;
   d_do_coarsen_test = false;

   d_initial_patch_hierarchy.setNull();
   d_test_patch_hierarchy.setNull();

   if (hier_test_db->keyExists("do_refine_test")) {
      d_do_refine_test = hier_test_db->getBool("do_refine_test");
      if (d_do_refine_test) {
         tbox::plog << "\nPerforming hierarchy refine test..." << std::endl;
         if (hier_test_db->keyExists("ratio")) {
            int* tmp_ratio = &d_ratio[0];
            hier_test_db->getIntegerArray("ratio", tmp_ratio, d_dim.getValue());
            tbox::plog << "with ratio = " << d_ratio << std::endl;
         } else {
            TBOX_ERROR(
               "HierarchyTester input error: no 'ratio' found in input"
               << std::endl);
         }
      }
   }

   if (!d_do_refine_test) {
      if (hier_test_db->keyExists("do_coarsen_test")) {
         d_do_coarsen_test = hier_test_db->getBool("do_coarsen_test");
      }
      if (d_do_coarsen_test) {
         tbox::plog << "\nPerforming hierarchy coarsen test..." << std::endl;
         if (hier_test_db->keyExists("ratio")) {
            int* tmp_ratio = &d_ratio[0];
            hier_test_db->getIntegerArray("ratio", tmp_ratio, d_dim.getValue());
            tbox::plog << "with ratio = " << d_ratio << std::endl;
         } else {
            TBOX_ERROR(
               "HierarchyTester input error: no 'ratio' found in input"
               << std::endl);
         }
      }
   }

   if (!d_do_refine_test && !d_do_coarsen_test) {
      TBOX_ERROR(
         "HierarchyTester input error: no test specified in input" << std::endl);
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
   TBOX_ASSERT(!main_input_db.isNull());
#endif
   Pointer<CartesianGridGeometry> grid_geometry(
      new CartesianGridGeometry(
         d_dim,
         "CartesianGridGeometry",
         main_input_db->getDatabase("CartesianGridGeometry")));

   d_initial_patch_hierarchy =
      new PatchHierarchy("InitialPatchHierarchy",
                         grid_geometry,
                         main_input_db->getDatabase("PatchHierarchy"));

   Pointer<BergerRigoutsos> box_generator(new BergerRigoutsos(d_dim));

   Pointer<TreeLoadBalancer> load_balancer(
      new TreeLoadBalancer(d_dim,
         "TreeLoadBalancer",
         main_input_db->getDatabase("TreeLoadBalancer")));
   load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

   Pointer<StandardTagAndInitialize> dummy_error_detector(
      new StandardTagAndInitialize(d_dim,
         "StandardTagAndInitialize",
         this,
         main_input_db->getDatabase("StandardTagAndInitialize")));

   d_gridding_algorithm =
      new GriddingAlgorithm(d_initial_patch_hierarchy,
         "GriddingAlgorithm",
         main_input_db->getDatabase("GriddingAlgorithm"),
         dummy_error_detector,
         box_generator,
         load_balancer);

   d_gridding_algorithm->makeCoarsestLevel(
      0.0);                                       // dummy time

   for (int ln = 0; d_initial_patch_hierarchy->levelCanBeRefined(ln); ln++) {
      d_gridding_algorithm->makeFinerLevel(
         0.0,                                      // dummy time
         true,                                     // indicates initial time
         0);                                       // dummy tag buffer

   }

   tbox::plog << "Initial hierarchy:\n";
   d_initial_patch_hierarchy->recursivePrint(
      tbox::plog, "   ", 3);

}

int HierarchyTester::runHierarchyTestAndVerify()
{
   int fail_count = 0;

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

   Pointer<GridGeometry> init_geometry =
      d_initial_patch_hierarchy->getGridGeometry();
   Pointer<GridGeometry> test_geometry =
      d_test_patch_hierarchy->getGridGeometry();

   hier::IntVector one_vector(d_dim, 1);

   // Test #0a:
   if (init_geometry->getPeriodicShift(one_vector) !=
       test_geometry->getPeriodicShift(one_vector)) {
      fail_count++;
      tbox::perr << "FAILED: - Test #0a: initial hierarchy has periodic shift "
                 << init_geometry->getPeriodicShift(one_vector) << " and \n"
                 << "test hierarchy has periodic shift "
                 << test_geometry->getPeriodicShift(one_vector) << std::endl;
   }

   const BoxList& init_phys_domain = init_geometry->getPhysicalDomain(0);
   const BoxList& test_phys_domain = test_geometry->getPhysicalDomain(0);

   const int npdboxes = init_phys_domain.getNumberOfBoxes();

   // Test #0b:
   hier::BoxList::Iterator ipditr(init_phys_domain);
   hier::BoxList::Iterator tpditr(test_phys_domain);
   if (d_do_refine_test) {
      for (int ib = 0; ib < npdboxes; ib++, ipditr++, tpditr++) {
         if (Box::refine(*ipditr, d_ratio) != *tpditr) {
            fail_count++;
            tbox::perr << "FAILED: - Test #0b: test hierarchy physical domain"
                       << " box with array index " << ib
                       << " is not a proper refinement of initial hierarchy"
                       << " physical domain box with same index" << std::endl;
         }
      }
   }
   if (d_do_coarsen_test) {
      for (int ib = 0; ib < npdboxes; ib++, ipditr++, tpditr++) {
         if (Box::coarsen(*ipditr, d_ratio) != *tpditr) {
            fail_count++;
            tbox::perr << "FAILED: - Test #0b: test hierarchy physical domain"
                       << " box with array index " << ib
                       << " is not a proper coarsening of initial hierarchy"
                       << " physical domain box with same index" << std::endl;
         }
      }
   }

   // Test #0c:
   if (init_geometry->getDomainIsSingleBox(0) !=
       test_geometry->getDomainIsSingleBox(0)) {
      fail_count++;
      tbox::perr
      << "FAILED: - Test #0c: initial and test hierarchy do not match"
      << " for geom->getDomainIsSingleBox()" << std::endl;
   }

   /*
    **************************************************************
    * Tests 1-8 check global data for levels in hierarchy.       *
    **************************************************************
    */

   const int nlevels = d_initial_patch_hierarchy->getNumberOfLevels();

   // Test #1:
   if (d_test_patch_hierarchy->getNumberOfLevels() != nlevels) {
      fail_count++;
      tbox::perr << "FAILED: - Test #1: initial hierarchy has "
                 << nlevels << " levels and \n"
                 << "test hierarchy has "
                 << d_test_patch_hierarchy->getNumberOfLevels() << "levels"
                 << std::endl;
   }

   for (int ln = 0; ln < nlevels; ln++) {
      Pointer<PatchLevel> init_level =
         d_initial_patch_hierarchy->getPatchLevel(ln);
      Pointer<PatchLevel> test_level =
         d_test_patch_hierarchy->getPatchLevel(ln);

      // Test #2:
      if (init_level->getLevelNumber() !=
          test_level->getLevelNumber()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #2: for level number " << ln
                    << " initial hierarchy level number is "
                    << init_level->getLevelNumber()
                    << "\n and test hierarchy level number is "
                    << test_level->getLevelNumber() << std::endl;
      }

      // Test #4:
      if (init_level->inHierarchy() !=
          test_level->inHierarchy()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #4: for level number " << ln
                    << " initial hierarchy level in hierarchy is "
                    << init_level->inHierarchy()
                    << "\n and test hierarchy level in hierarchy is "
                    << test_level->inHierarchy() << std::endl;
      }

      // Test #5:
      if (init_level->getGlobalNumberOfPatches() !=
          test_level->getGlobalNumberOfPatches()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #5: for level number " << ln
                    << " initial hierarchy number of patches is "
                    << init_level->getGlobalNumberOfPatches()
                    << "\n and test hierarchy number of patches is "
                    << test_level->getGlobalNumberOfPatches() << std::endl;
      }

      // Test #6:
      if (init_level->getRatioToLevelZero() !=
          test_level->getRatioToLevelZero()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #6: for level number " << ln
                    << " initial hierarchy ratio to level zero is "
                    << init_level->getRatioToLevelZero()
                    << "\n and test hierarchy ratio to level zero is "
                    << test_level->getRatioToLevelZero() << std::endl;
      }

      // Test #7:
      if (init_level->getRatioToCoarserLevel() !=
          test_level->getRatioToCoarserLevel()) {
         fail_count++;
         tbox::perr << "FAILED: - Test #7: for level number " << ln
                    << " initial hierarchy ratio to coarser level is "
                    << init_level->getRatioToCoarserLevel()
                    << "\n and test hierarchy ratio to coarser level is "
                    << test_level->getRatioToCoarserLevel() << std::endl;
      }

      const BoxList& init_domain = init_level->getPhysicalDomain(BlockId::zero());
      const BoxList& test_domain = test_level->getPhysicalDomain(BlockId::zero());

      const int nboxes = init_domain.getNumberOfBoxes();

      // Test #8:
      hier::BoxList::Iterator iditr(init_domain);
      hier::BoxList::Iterator tditr(test_domain);
      if (d_do_refine_test) {
         for (int ib = 0; ib < nboxes; ib++, iditr++, tditr++) {
            if (Box::refine(*iditr, d_ratio) != *tditr) {
               fail_count++;
               tbox::perr << "FAILED: - Test #8: for level number " << ln
                          << " refined domain box with array index " << ib
                          << " is not a proper refinement of initial domain "
                          << "box with same index" << std::endl;
            }
         }
      }
      if (d_do_coarsen_test) {
         for (int ib = 0; ib < nboxes; ib++, iditr++, tditr++) {
            if (Box::coarsen(*iditr, d_ratio) != *tditr) {
               fail_count++;
               tbox::perr << "FAILED: - Test #8: for level number " << ln
                          << " coarsened domain box with array index " << ib
                          << " is not a proper coarsening of initial domain "
                          << "box with same index" << std::endl;
            }
         }
      }

      /*
       **************************************************************
       *  Tests 9-13 check global data for patches on each level.   *
       **************************************************************
       */

      IntVector init_connector_width =
         d_initial_patch_hierarchy->getRequiredConnectorWidth(ln,ln);
      IntVector test_connector_width =
         d_test_patch_hierarchy->getRequiredConnectorWidth(ln,ln);
      if (d_do_refine_test) {
         test_connector_width = init_connector_width*d_ratio;
      }
      else {
         init_connector_width = test_connector_width*d_ratio;
      }
      const Connector &init_connector =
         init_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findOrCreateConnector(
            d_initial_patch_hierarchy->getDomainMappedBoxLevel(),
            init_connector_width,
            true /* exact width only */ );
      const Connector &test_connector =
         test_level->getMappedBoxLevel()->getPersistentOverlapConnectors().findConnector(
            d_test_patch_hierarchy->getDomainMappedBoxLevel(),
            test_connector_width,
            true /* exact width only */ );

      for (hier::PatchLevel::Iterator ip(test_level); ip; ip++) {
         const MappedBoxId& mapped_box_id = ip->getMappedBox().getId();
         // Test #9:
         if (d_do_refine_test) {
            if (Box::refine(init_level->getBoxForPatch(mapped_box_id), d_ratio) !=
                test_level->getBoxForPatch(mapped_box_id)) {
               fail_count++;
               tbox::perr << "FAILED: - Test #9: for level number " << ln
                          << " refined patch box with array index " << mapped_box_id 
                          << " is not a proper refinement of initial domain "
                          << "box with same index" << std::endl;
            }
         }

         if (d_do_coarsen_test) {
            if (Box::coarsen(init_level->getBoxForPatch(mapped_box_id), d_ratio) !=
                test_level->getBoxForPatch(mapped_box_id)) {
               fail_count++;
               tbox::perr << "FAILED: - Test #9: for level number " << ln
                          << " coarsened patch box with array index " << mapped_box_id
                          << " is not a proper coarsening of initial domain "
                          << "box with same index" << std::endl;
            }
         }

         // Test #10:
         if (init_connector.getNeighborSet(mapped_box_id) != test_connector.getNeighborSet(mapped_box_id)) {
            fail_count++;
            tbox::perr << "FAILED: - Test #10: for level number " << ln
                       << " initial and test level have different number of "
                       << "domain neighbors for patch number " << mapped_box_id << std::endl;
         }

         // Test #11:
         if (init_level->getMappingForPatch(mapped_box_id) !=
             test_level->getMappingForPatch(mapped_box_id)) {
            fail_count++;
            tbox::perr << "FAILED: - Test #11: for level number " << ln
                       << " initial and test level have different processor "
                       << "mapping for patch number " << mapped_box_id << std::endl;
         }

         // Test #12:
         if (init_level->patchTouchesRegularBoundary(mapped_box_id) !=
             test_level->patchTouchesRegularBoundary(mapped_box_id)) {
            fail_count++;
            tbox::perr << "FAILED: - Test #12: for level number " << ln
                       << " initial and test level do not match for "
                       << "patchTouchesRegularBoundary() "
                       << "for patch number " << mapped_box_id << std::endl;
         }

         // Test #13:
         if (init_level->patchTouchesPeriodicBoundary(mapped_box_id) !=
             test_level->patchTouchesPeriodicBoundary(mapped_box_id)) {
            fail_count++;
            tbox::perr << "FAILED: - Test #12: for level number " << ln
                       << " initial and test level do not match for "
                       << "patchTouchesPeriodicBoundary() "
                       << "for patch number " << mapped_box_id << std::endl;
         }
      }

      /*
       **************************************************************
       *  Tests 14-20 check local data for patches on each level.   *
       **************************************************************
       */
      for (PatchLevel::Iterator tip(test_level); tip; tip++) {
         const MappedBoxId& mapped_box_id = tip->getMappedBox().getId();
         Pointer<Patch> test_patch = test_level->getPatch(mapped_box_id);
         Pointer<Patch> init_patch = init_level->getPatch(mapped_box_id);

         // Test #14:
         if (d_do_refine_test) {
            if (Box::refine(init_patch->getBox(), d_ratio) !=
                test_patch->getBox()) {
               fail_count++;
               tbox::perr << "FAILED: - Test #14: for level number " << ln
                          << " box for test level patch " << mapped_box_id
                          << " is not a proper refinement of box "
                          << "for initial level patch with same number"
                          << std::endl;
            }
         }
         if (d_do_coarsen_test) {
            if (Box::coarsen(init_patch->getBox(), d_ratio) !=
                test_patch->getBox()) {
               fail_count++;
               tbox::perr << "FAILED: - Test #14: for level number " << ln
                          << " box for test level patch " << mapped_box_id
                          << " is not a proper coarsening of box "
                          << "for initial level patch with same number"
                          << std::endl;
            }
         }

         // Test #15:
         if (init_patch->getLocalId() !=
             test_patch->getLocalId()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #15: for level number " << ln
                       <<
            " initial and test level patches have different patch "
                       << "numbers for patch with index " << tip->getLocalId() << std::endl;
         }

         // Test #16:
         if (init_patch->getPatchLevelNumber() !=
             test_patch->getPatchLevelNumber()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #16: for level number " << ln
                       <<
            " initial and test level patches have different patch "
                       << "level numbers for patch number " << tip->getLocalId()
                       << std::endl;
         }

         // Test #17:
         if (init_patch->inHierarchy() !=
             test_patch->inHierarchy()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #17: for level number " << ln
                       << " initial and test level do not match for "
                       << "inHierarchy() "
                       << "for patch number " << tip->getLocalId() << std::endl;
         }

         // Test #18:
         if (init_patch->getPatchGeometry()->getTouchesRegularBoundary() !=
             test_patch->getPatchGeometry()->getTouchesRegularBoundary()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #18: for level number " << ln
                       << " initial and test level do not match for "
                       << "getTouchesRegularBoundary() "
                       << "for patch number " << tip->getLocalId() << std::endl;
         }

         // Test #19:
         if (init_patch->getPatchGeometry()->getTouchesPeriodicBoundary() !=
             test_patch->getPatchGeometry()->getTouchesPeriodicBoundary()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #19: for level number " << ln
                       << " initial and test level do not match for "
                       << "getTouchesPeriodicBoundary() "
                       << "for patch number " << tip->getLocalId() << std::endl;
         }

         /*
          **************************************************************
          * Tests 20a-20c check patch geometry data.                   *
          **************************************************************
          */

         Pointer<PatchGeometry> init_patch_geom =
            init_patch->getPatchGeometry();
         Pointer<PatchGeometry> test_patch_geom =
            test_patch->getPatchGeometry();

         // Test #20a:
         if (init_patch_geom->getRatio() !=
             test_patch_geom->getRatio()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #20a: for level number " << ln
                       << " patch geometry ratio data does not match "
                       << "for patch number " << tip->getLocalId() << std::endl;
         }

         // Test #20b:
         if (init_patch_geom->intersectsPhysicalBoundary() !=
             test_patch_geom->intersectsPhysicalBoundary()) {
            fail_count++;
            tbox::perr << "FAILED: - Test #20b: for level number " << ln
                       << " intersectsPhysicalBoundary() does not match "
                       << "for patch number " << tip->getLocalId() << std::endl;
         }

         // Test #20c:
         for (int id = 1; id <= d_dim.getValue(); id++) {
            if ((init_patch_geom->getCodimensionBoundaries(id)).getSize() !=
                (test_patch_geom->getCodimensionBoundaries(id)).getSize()) {
               fail_count++;
               tbox::perr << "FAILED: - Test #20c: for level number " << ln
                          << " number of codimension " << id
                          << " boundary boxes does not match "
                          << "for patch number " << tip->getLocalId() << std::endl;
            }
         }
      }
   }

   return fail_count;
}

}
