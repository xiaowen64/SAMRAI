//
// File:        NodeMultiblockTest.C
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 1.5 $
// Modified:    $Date: 2004/02/11 23:46:08 $
// Description: AMR communication tests for node-centered patch data
//

#include "NodeMultiblockTest.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#include "BoundaryBox.h"
#include "SkeletonPatchGeometry.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "NodeVariable.h"
#include "MultiblockTester.h"
#include "tbox/Utilities.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"

using namespace SAMRAI;

NodeMultiblockTest::NodeMultiblockTest(
   const string& object_name, 
   tbox::Pointer<tbox::Database> main_input_db, 
   bool do_refine, 
   bool do_coarsen,
   const string& refine_option)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!object_name.empty());
   assert(!main_input_db.isNull());
   assert(!refine_option.empty());
#endif

   d_object_name = object_name;

   d_refine_option = refine_option;

   d_finest_level_number = main_input_db->
                           getDatabase("GriddingAlgorithm")->
                           getInteger("max_levels") - 1;

   int num_blocks = main_input_db->getDatabase("Multiblock")->
                                   getInteger("num_blocks");

   d_skel_grid_geometry =
      new tbox::Pointer<geom::SkeletonGridGeometry<NDIM> >[num_blocks];

   char geom_name[32];

   for (int g = 0; g < num_blocks; g++) {

      sprintf(geom_name, "SkeletonGridGeometry%d", g);

      if (main_input_db->keyExists(geom_name)) {
         d_skel_grid_geometry[g] = new geom::SkeletonGridGeometry<NDIM>(
                                      geom_name,
                                      main_input_db->getDatabase(geom_name));

      } else {
         break;
      }
   }

   setGridGeometrySize(num_blocks);

   for (int b = 0; b < num_blocks; b++) {

      setGridGeometry(d_skel_grid_geometry[b],b);
   }

   readTestInput(main_input_db->getDatabase("NodeMultiblockTest"));
}

NodeMultiblockTest::~NodeMultiblockTest()
{
   delete[] d_skel_grid_geometry;
}

void NodeMultiblockTest::readTestInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void NodeMultiblockTest::registerVariables(MultiblockTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(commtest != (MultiblockTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] =
         new pdat::NodeVariable<NDIM,double>(d_variable_src_name[i],
                                             d_variable_depth[i]);

      commtest->registerVariable(d_variables[i],
                                 d_variables[i],
                                 d_variable_src_ghosts[i], 
                                 d_variable_dst_ghosts[i], 
                                 d_skel_grid_geometry[0], 
                                 d_variable_refine_op[i]); 

   }

}


void NodeMultiblockTest::initializeDataOnPatch(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   int block_number,
   char src_or_dst)
{
   (void) hierarchy;

   if (   (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") 
       || ( (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") 
           && (level_number < d_finest_level_number) ) ) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box<NDIM> dbox = node_data->getGhostBox(); 

         node_data->fillAll((double)block_number);

      }
   }
}

void NodeMultiblockTest::tagCellsToRefine(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number, 
   int tag_index)
{
   (void) hierarchy;

   /*
    * Base class sets tags in box array for each level.
    */
   tagCellsInInputBoxes(patch, level_number, tag_index); 
   
}

void NodeMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw_to_fill) const
{
   (void) time;

   tbox::Pointer< geom::SkeletonPatchGeometry<NDIM> >
      pgeom = patch.getPatchGeometry();

   const tbox::Array<hier::BoundaryBox<NDIM> > node_bdry = 
      pgeom->getCodimensionBoundary(NDIM);
   const int num_node_bdry_boxes = node_bdry.getSize();

#if (NDIM > 1)
   const tbox::Array<hier::BoundaryBox<NDIM> > edge_bdry = 
      pgeom->getCodimensionBoundary(NDIM - 1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();
#endif

#if (NDIM == 3)
   const tbox::Array<hier::BoundaryBox<NDIM> > face_bdry = 
      pgeom->getCodimensionBoundary(NDIM - 2);
   const int num_face_bdry_boxes = face_bdry.getSize();
#endif

   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());

      /*
       * Set node boundary data.
       */
      for (int nb = 0; nb < num_node_bdry_boxes; nb++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[nb],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         hier::Box<NDIM> patch_node_box =
            pdat::NodeGeometry<NDIM>::toNodeBox(patch.getBox()); 
         if (!node_bdry[nb].getIsMultiblockSingularity()) {
            for (pdat::NodeIterator<NDIM> ni(fill_box); ni; ni++) {
               if (!patch_node_box.contains(ni())) {
                  for (int d = 0; d < node_data->getDepth(); d++) { 
                     (*node_data)(ni(),d) =
                        (double)(node_bdry[nb].getLocationIndex()+100);
                  }
               }
            }
         }
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int eb = 0; eb < num_edge_bdry_boxes; eb++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[eb],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         hier::Box<NDIM> patch_node_box =
            pdat::NodeGeometry<NDIM>::toNodeBox(patch.getBox());
         hier::Index<NDIM> plower(patch_node_box.lower());
         hier::Index<NDIM> pupper(patch_node_box.upper());

         if (!edge_bdry[eb].getIsMultiblockSingularity()) {
            for (pdat::NodeIterator<NDIM> ni(fill_box); ni; ni++) {
               if (!patch_node_box.contains(ni())) {
                  bool use_index = true;
                  for (int n = 0; n < NDIM; n++) {
                     if (edge_bdry[eb].getBox().numberCells(n) == 1) {
                        if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }

                  if (use_index) {
                     for (int d = 0; d < node_data->getDepth(); d++) {
                        (*node_data)(ni(),d) =
                           (double)(edge_bdry[eb].getLocationIndex()+100);
                     }
                  }
               }
            }
         }
      }
#endif

#if (NDIM == 3)
      /*
       * Set face boundary data.
       */
      for (int fb = 0; fb < num_face_bdry_boxes; fb++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(face_bdry[fb],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         hier::Box<NDIM> patch_node_box =
            pdat::NodeGeometry<NDIM>::toNodeBox(patch.getBox());
         hier::Index<NDIM> plower(patch_node_box.lower());
         hier::Index<NDIM> pupper(patch_node_box.upper());

         if (!face_bdry[fb].getIsMultiblockSingularity()) {
            for (pdat::NodeIterator<NDIM> ni(fill_box); ni; ni++) {
               if (!patch_node_box.contains(ni())) {
                  bool use_index = true;
                  for (int n = 0; n < NDIM; n++) {
                     if (face_bdry[fb].getBox().numberCells(n) == 1) {
                        if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }
                                                                                
                  if (use_index) {
                     for (int d = 0; d < node_data->getDepth(); d++) {
                        (*node_data)(ni(),d) =
                           (double)(face_bdry[fb].getLocationIndex()+100);
                     }
                  }
               }
            }
         }
      }
#endif 

   }

}

void NodeMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch<NDIM>& patch,
   tbox::List<mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
   const hier::Box<NDIM>& fill_box,
   const hier::BoundaryBox<NDIM>& bbox)
{
   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> pbox(pdat::NodeGeometry<NDIM>::toNodeBox(patch.getBox()));

      hier::Index<NDIM> plower(pbox.lower());
      hier::Index<NDIM> pupper(pbox.upper());

      hier::Box<NDIM> sing_fill_box(node_data->getGhostBox() * fill_box);

      int depth = node_data->getDepth();

      for (pdat::NodeIterator<NDIM> ni(sing_fill_box); ni; ni++) {
         bool use_index = true;
         for (int n = 0; n < NDIM; n++) {
            if (bbox.getBox().numberCells(n) == 1) {
               if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                  use_index = false;
                  break;
               }
            }
         }
         if (use_index) {
            for (int d = 0; d < depth; d++) {
               (*node_data)(ni(),d) = 0.0;
            }
         }
      }

      /*
       * If sing_patches is not empty, that means there is enhanced
       * connectivity, and we get data from other blocks
       */

      if (sing_patches.size()) {

         for (tbox::List
              <mblk::MultiblockRefineSchedule<NDIM>::SingularityPatch>::
              Iterator sp(sing_patches); sp; sp++) {
            tbox::Pointer< pdat::NodeData<NDIM,double> > sing_data =
               sp().d_patch->getPatchData(d_variables[i], getDataContext());
            int sing_neighbor_id = sp().d_id;
            for (pdat::NodeIterator<NDIM> ci(sing_fill_box); ci; ci++) {
               bool use_index = true;
               for (int n = 0; n < NDIM; n++) {
                  if (bbox.getBox().numberCells(n) == 1) {
                     if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; d++) {
                     (*node_data)(ci(),d) += sing_neighbor_id;
                  }
               }
            }
         }

         for (pdat::NodeIterator<NDIM> ci(sing_fill_box); ci; ci++) {
            bool use_index = true;
            for (int n = 0; n < NDIM; n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
            }
               }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(ci(),d) /= sing_patches.size();
               }
            }
         }

      /*
       * In cases of reduced connectivity, there are no other blocks
       * from which to acquire data.
       */

      } else {

         for (pdat::NodeIterator<NDIM> ci(sing_fill_box); ci; ci++) {
            bool use_index = true;
            for (int n = 0; n < NDIM; n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(ci(),d) = (double)bbox.getLocationIndex()+200.0;
               }
            }
         }
      }
   }
}

/*
*************************************************************************
*                                                                       *
* Verify results of communication operations.  This test must be        *
* consistent with data initialization and boundary operations above.    *
*                                                                       *
*************************************************************************
*/
void NodeMultiblockTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<mblk::MultiblockPatchHierarchy<NDIM> > hierarchy, 
   int level_number,
   int block_number)
{

   tbox::plog << "\nEntering NodeMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl; 

   hier::IntVector<NDIM> tgcw(0);
   for (int i = 0; i < d_variables.getSize(); i++) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
                                  getGhostCellWidth());
   }
   hier::Box<NDIM> pbox = patch.getBox();

   tbox::Pointer< pdat::NodeData<NDIM,double> > solution =
      new pdat::NodeData<NDIM,double>(pbox, 1, tgcw);

   hier::Box<NDIM> tbox(pbox);
   tbox.grow(tgcw);

   tbox::List<mblk::MultiblockPatchHierarchy<NDIM>::Neighbor>& neighbors =
      hierarchy->getNeighbors(block_number);
   hier::BoxList<NDIM> singularity(
      hierarchy->getSingularityBoxList(block_number));

   hier::IntVector<NDIM> ratio =
      hierarchy->getPatchLevel(level_number)->getRatio();

   singularity.refine(ratio);

   bool test_failed = false;

   for (int i = 0; i < d_variables.getSize(); i++) {

      double correct = (double)block_number;

      tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());
      int depth = node_data->getDepth();

      hier::Box<NDIM> interior_box(pbox);
      interior_box.grow(-1);

      for (pdat::NodeIterator<NDIM> ci(interior_box); ci; ci++) {
         for (int d = 0; d < depth; d++) {
            double result = (*node_data)(ci(),d);

            if (!tbox::Utilities::deq(correct, result)) {
               tbox::perr << "Test FAILED: ...."
                    << " : node index = " << ci() << endl;
               tbox::perr << "    Variable = " << d_variable_src_name[i]
                    << " : depth index = " << d << endl;
               tbox::perr << "    result = " << result
                    << " : correct = " << correct << endl;
               test_failed = true;
            }
         }
      }

      hier::Box<NDIM> gbox = node_data->getGhostBox();

      hier::Box<NDIM> patch_node_box =
         pdat::NodeGeometry<NDIM>::toNodeBox(pbox);

      for (tbox::List<mblk::MultiblockPatchHierarchy<NDIM>::Neighbor>::
           Iterator ne(neighbors); ne; ne++) {

         correct = ne().d_id;

         hier::BoxList<NDIM> neighbor_ghost(ne().d_translated_domain);
         neighbor_ghost.refine(ratio);
         neighbor_ghost.intersectBoxes(gbox);

         for (hier::BoxList<NDIM>::Iterator ng(neighbor_ghost); ng; ng++) {

            hier::Box<NDIM> patch_node_box =
               pdat::NodeGeometry<NDIM>::toNodeBox(pbox);

            for (pdat::NodeIterator<NDIM> ci(ng()); ci; ci++) {
               if (!patch_node_box.contains(ci())) {
                  for (int d = 0; d < depth; d++) {
                     double result = (*node_data)(ci(),d);

                     if (!tbox::Utilities::deq(correct, result)) {
                        tbox::perr << "Test FAILED: ...."
                             << " : node index = " << ci() << endl;
                        tbox::perr << "  Variable = " << d_variable_src_name[i]
                             << " : depth index = " << d << endl;
                        tbox::perr << "    result = " << result
                             << " : correct = " << correct << endl;
                        test_failed = true;
                     }
                  }
               }
            } 
         } 
      } 

      tbox::Pointer< geom::SkeletonPatchGeometry<NDIM> > pgeom =
         patch.getPatchGeometry();

      for (int b = 0; b < NDIM; b++) {
         tbox::Array<hier::BoundaryBox<NDIM> > bdry =
            pgeom->getCodimensionBoundary(b+1);

         for (int k = 0; k < bdry.size(); k++) {
            hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(bdry[k],
                                                     patch.getBox(),
                                                     tgcw);
            fill_box = fill_box * gbox;

            if (bdry[k].getIsMultiblockSingularity()) {
               correct = 0.0;

               int num_sing_neighbors = 0;
               for (tbox::List
                    <mblk::MultiblockPatchHierarchy<NDIM>::Neighbor>::
                    Iterator ns(neighbors); ns; ns++) {
                  if (ns().d_is_singularity) {
                     hier::BoxList<NDIM> neighbor_ghost(
                                           ns().d_translated_domain);
                     neighbor_ghost.refine(ratio);
                     neighbor_ghost.intersectBoxes(fill_box);
                     if (neighbor_ghost.size()) {
                        num_sing_neighbors++;
                        correct += (double)ns().d_id; 
                     }
                  }
               }

               if (num_sing_neighbors == 0) {

                  correct = (double)bdry[k].getLocationIndex()+200.0;

               } else {

                  correct /= (double) num_sing_neighbors;

               }

            } else {
               correct = (double)(bdry[k].getLocationIndex() + 100);
            }

            for (pdat::NodeIterator<NDIM> ci(fill_box); ci; ci++) {

               if (!patch_node_box.contains(ci())) {

                  bool use_index = true;
                  for (int n = 0; n < NDIM; n++) {
                     if (bdry[k].getBox().numberCells(n) == 1) {
                        if (ci()(n) == patch_node_box.lower()(n) ||
                            ci()(n) == patch_node_box.upper()(n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }

                  if (use_index) {
                     for (int d = 0; d < depth; d++) {
                        double result = (*node_data)(ci(),d);

                        if (!tbox::Utilities::deq(correct, result)) {
                           tbox::perr << "Test FAILED: ...."
                                << " : node index = " << ci() << endl;
                           tbox::perr << "  Variable = " 
                                << d_variable_src_name[i]
                                << " : depth index = " << d << endl;
                           tbox::perr << "    result = " << result
                                << " : correct = " << correct << endl;
                           test_failed = true;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (!test_failed) {   
      tbox::plog << "NodeMultiblockTest Successful!" << endl;
   } else {
      tbox::perr << "Multiblock NodeMultiblockTest FAILED: .\n" << endl;
   }

   solution.setNull();   // just to be anal...

   tbox::plog << "\nExiting NodeMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

}
