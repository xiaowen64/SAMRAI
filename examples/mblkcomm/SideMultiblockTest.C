//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/mblkcomm/SideMultiblockTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: AMR communication tests for side-centered patch data
//

#include "SideMultiblockTest.h"


#include "BoundaryBox.h"
#include "BlockPatchGeometry.h"
#include "SideDoubleConstantRefine.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "SideVariable.h"
#include "MultiblockTester.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"

using namespace SAMRAI;

SideMultiblockTest::SideMultiblockTest(
   const string& object_name, 
   tbox::Pointer<tbox::Database> main_input_db, 
   bool do_refine, 
   bool do_coarsen,
   const string& refine_option)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!main_input_db.isNull());
   TBOX_ASSERT(!refine_option.empty());
#endif

   d_object_name = object_name;

   d_refine_option = refine_option;

   d_finest_level_number = main_input_db->
                           getDatabase("GriddingAlgorithm")->
                           getInteger("max_levels") - 1;

   int num_blocks = main_input_db->getDatabase("Multiblock")->
                                   getInteger("num_blocks");

   d_skel_grid_geometry.resizeArray(num_blocks);

   char geom_name[32];

   for (int g = 0; g < num_blocks; g++) {

      sprintf(geom_name, "BlockGridGeometry%d", g);

      if (main_input_db->keyExists(geom_name)) {
         d_skel_grid_geometry[g] = new geom::BlockGridGeometry<NDIM>(
                                      geom_name,
                                      main_input_db->getDatabase(geom_name),
                                      g);

      } else {
         break;
      }
   }

   tbox::Pointer< hier::MultiblockGridGeometry<NDIM> > mblk_geometry =
      new hier::MultiblockGridGeometry<NDIM>(d_skel_grid_geometry);

   setGridGeometry(mblk_geometry);

   readTestInput(main_input_db->getDatabase("SideMultiblockTest"));
}

SideMultiblockTest::~SideMultiblockTest()
{
}

void SideMultiblockTest::readTestInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void SideMultiblockTest::registerVariables(MultiblockTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (MultiblockTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] =
         new pdat::SideVariable<NDIM,double>(d_variable_src_name[i],
                                             d_variable_depth[i]);

      commtest->registerVariable(d_variables[i],
                                 d_variables[i],
                                 d_variable_src_ghosts[i], 
                                 d_variable_dst_ghosts[i], 
                                 d_skel_grid_geometry[0], 
                                 d_variable_refine_op[i]); 

   }

}


void SideMultiblockTest::initializeDataOnPatch(
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

         tbox::Pointer< pdat::SideData<NDIM,double> > side_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box<NDIM> dbox = side_data->getGhostBox(); 

         side_data->fillAll((double)block_number);

      }
   }
}

void SideMultiblockTest::tagCellsToRefine(
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

void SideMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw_to_fill) const
{
   (void) time;

   tbox::Pointer< geom::BlockPatchGeometry<NDIM> >
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

      tbox::Pointer< pdat::SideData<NDIM,double> > side_data =
         patch.getPatchData(d_variables[i], getDataContext());

      /*
       * Set node boundary data.
       */
      for (int nb = 0; nb < num_node_bdry_boxes; nb++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[nb],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         for (int axis = 0; axis < NDIM; axis++) {
            hier::Box<NDIM> patch_side_box =
               pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(),axis); 
            if (!node_bdry[nb].getIsMultiblockSingularity()) {
               for (pdat::SideIterator<NDIM> ni(fill_box, axis); ni; ni++) {
                  if (!patch_side_box.contains(ni())) {
                     for (int d = 0; d < side_data->getDepth(); d++) { 
                        (*side_data)(ni(),d) =
                           (double)(node_bdry[nb].getLocationIndex()+100);
                     }
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

         for (int axis = 0; axis < NDIM; axis++) {
            hier::Box<NDIM> patch_side_box =
               pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);
            hier::Index<NDIM> plower(patch_side_box.lower());
            hier::Index<NDIM> pupper(patch_side_box.upper());

 
            if (!edge_bdry[eb].getIsMultiblockSingularity()) {
               for (pdat::SideIterator<NDIM> ni(fill_box, axis); ni; ni++) {
                  if (!patch_side_box.contains(ni())) {
                     bool use_index = true;
                     for (int n = 0; n < NDIM; n++) {
                        if (axis == n &&
                            edge_bdry[eb].getBox().numberCells(n) == 1) { 
                           if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < side_data->getDepth(); d++) {
                           (*side_data)(ni(),d) =
                              (double)(edge_bdry[eb].getLocationIndex()+100);
                        }
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

         for (int axis = 0; axis < NDIM; axis++) {
            hier::Box<NDIM> patch_side_box =
               pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);
            hier::Index<NDIM> plower(patch_side_box.lower());
            hier::Index<NDIM> pupper(patch_side_box.upper());
                                                                                
                                                                                
            if (!face_bdry[fb].getIsMultiblockSingularity()) {
               for (pdat::SideIterator<NDIM> ni(fill_box, axis); ni; ni++) {
                  if (!patch_side_box.contains(ni())) {
                     bool use_index = true;
                     for (int n = 0; n < NDIM; n++) {
                        if (axis == n &&
                            face_bdry[fb].getBox().numberCells(n) == 1) {
                           if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }
                                                                                
                     if (use_index) {
                        for (int d = 0; d < side_data->getDepth(); d++) {
                           (*side_data)(ni(),d) =
                              (double)(face_bdry[fb].getLocationIndex()+100);
                        }
                     }
                  }
               }
            }
         }
      }
#endif 

   }

}

void SideMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch<NDIM>& patch,
   tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
   const hier::Box<NDIM>& fill_box,
   const hier::BoundaryBox<NDIM>& bbox)
{
   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::SideData<NDIM,double> > side_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> sing_fill_box(side_data->getGhostBox() * fill_box);
                                                                                
      int depth = side_data->getDepth();

      for (int axis = 0; axis < NDIM; axis++) {
         hier::Box<NDIM> pbox =
            pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);

         hier::Index<NDIM> plower(pbox.lower());
         hier::Index<NDIM> pupper(pbox.upper());

         for (pdat::SideIterator<NDIM> ni(sing_fill_box, axis); ni; ni++) {
            bool use_index = true;
            for (int n = 0; n < NDIM; n++) {
               if (axis == n && bbox.getBox().numberCells(n) == 1) {
                  if (ni()(n) == plower(n) || ni()(n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*side_data)(ni(),d) = 0.0;
               }
            }
         }
      }

      /*
       * If sing_patches is not empty, that means there is enhanced
       * connectivity, and we get data from other blocks
       */

      if (sing_patches.size()) {

         for (tbox::List
              <xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>::
              Iterator sp(sing_patches); sp; sp++) {
            tbox::Pointer< pdat::SideData<NDIM,double> > sing_data =
               sp().d_patch->getPatchData(d_variables[i], getDataContext());
            int sing_neighbor_id = sp().d_id;
            for (int axis = 0; axis < NDIM; axis++) {

               hier::Box<NDIM> pbox =
                  pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);
                                                                                
               hier::Index<NDIM> plower(pbox.lower());
               hier::Index<NDIM> pupper(pbox.upper());

               for (pdat::SideIterator<NDIM> ci(sing_fill_box,axis); ci; ci++) {
                  bool use_index = true;
                  for (int n = 0; n < NDIM; n++) {
                     if (axis == n && bbox.getBox().numberCells(n) == 1) {
                        if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }
                  if (use_index) {
                     for (int d = 0; d < depth; d++) {
                        (*side_data)(ci(),d) += sing_neighbor_id;
                     }
                  }
               }
            }
         }

         for (int axis = 0; axis < NDIM; axis++) {

            hier::Box<NDIM> pbox =
               pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);

            hier::Index<NDIM> plower(pbox.lower());
            hier::Index<NDIM> pupper(pbox.upper());

            for (pdat::SideIterator<NDIM> ci(sing_fill_box, axis); ci; ci++) {
               bool use_index = true;
               for (int n = 0; n < NDIM; n++) {
                  if (axis == n && bbox.getBox().numberCells(n) == 1) {
                     if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; d++) {
                     (*side_data)(ci(),d) /= sing_patches.size();
                  }
               }
            }
         }

      /*
       * In cases of reduced connectivity, there are no other blocks
       * from which to acquire data.
       */

      } else {

         for (int axis = 0; axis < NDIM; axis++) {

            hier::Box<NDIM> pbox =
               pdat::SideGeometry<NDIM>::toSideBox(patch.getBox(), axis);

            hier::Index<NDIM> plower(pbox.lower());
            hier::Index<NDIM> pupper(pbox.upper());

            for (pdat::SideIterator<NDIM> ci(sing_fill_box, axis); ci; ci++) {
               bool use_index = true;
               for (int n = 0; n < NDIM; n++) {
                  if (axis == n && bbox.getBox().numberCells(n) == 1) {
                     if (ci()(n) == plower(n) || ci()(n) == pupper(n)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; d++) {
                     (*side_data)(ci(),d) =
                        (double)bbox.getLocationIndex()+200.0;
                  }
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
bool SideMultiblockTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > hierarchy, 
   int level_number,
   int block_number)
{

   tbox::plog << "\nEntering SideMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl; 

   hier::IntVector<NDIM> tgcw(0);
   for (int i = 0; i < d_variables.getSize(); i++) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
                                  getGhostCellWidth());
   }
   hier::Box<NDIM> pbox = patch.getBox();

   tbox::Pointer< pdat::SideData<NDIM,double> > solution =
      new pdat::SideData<NDIM,double>(pbox, 1, tgcw);

   hier::Box<NDIM> tbox(pbox);
   tbox.grow(tgcw);

   tbox::List<hier::MultiblockPatchHierarchy<NDIM>::Neighbor>& neighbors =
      hierarchy->getNeighbors(block_number);
   hier::BoxList<NDIM> singularity(
      hierarchy->getSingularityBoxList(block_number));

   hier::IntVector<NDIM> ratio =
      hierarchy->getPatchLevel(level_number)->getRatio();

   singularity.refine(ratio);

   bool test_failed = false;

   for (int i = 0; i < d_variables.getSize(); i++) {

      double correct = (double)block_number;

      tbox::Pointer< pdat::SideData<NDIM,double> > side_data =
         patch.getPatchData(d_variables[i], getDataContext());
      int depth = side_data->getDepth();

      hier::Box<NDIM> interior_box(pbox);
      interior_box.grow(-1);

      for (int axis = 0; axis < NDIM; axis++) {
         for (pdat::SideIterator<NDIM> ci(interior_box, axis); ci; ci++) {
            for (int d = 0; d < depth; d++) {
               double result = (*side_data)(ci(),d);

               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                       << " : side index = " << ci() << endl;
                  tbox::perr << "    Variable = " << d_variable_src_name[i]
                       << " : depth index = " << d << endl;
                  tbox::perr << "    result = " << result
                       << " : correct = " << correct << endl;
                  test_failed = true;
               }
            }
         }
      }

      hier::Box<NDIM> gbox = side_data->getGhostBox();

      for (int axis = 0; axis < NDIM; axis++) {
         hier::Box<NDIM> patch_side_box =
            pdat::SideGeometry<NDIM>::toSideBox(pbox, axis);

         for (tbox::List<hier::MultiblockPatchHierarchy<NDIM>::Neighbor>::
              Iterator ne(neighbors); ne; ne++) {

            correct = ne().d_id;

            hier::BoxList<NDIM> neighbor_ghost(ne().d_translated_domain);
            neighbor_ghost.refine(ratio);
            neighbor_ghost.intersectBoxes(gbox);

            for (hier::BoxList<NDIM>::Iterator ng(neighbor_ghost); ng; ng++) {

               for (pdat::SideIterator<NDIM> ci(ng(), axis); ci; ci++) {
                  if (!patch_side_box.contains(ci())) {
                     for (int d = 0; d < depth; d++) {
                        double result = (*side_data)(ci(),d);
   
                        if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                           tbox::perr << "Test FAILED: ...."
                                << " : side index = " << ci() << endl;
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

      tbox::Pointer< geom::BlockPatchGeometry<NDIM> > pgeom =
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
                    <hier::MultiblockPatchHierarchy<NDIM>::Neighbor>::
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

            for (int axis = 0; axis < NDIM; axis++) {
               hier::Box<NDIM> patch_side_box =
                  pdat::SideGeometry<NDIM>::toSideBox(pbox, axis);

               for (pdat::SideIterator<NDIM> ci(fill_box, axis); ci; ci++) {

                  if (!patch_side_box.contains(ci())) {
   
                     bool use_index = true;
                     for (int n = 0; n < NDIM; n++) {
                        if (axis == n && bdry[k].getBox().numberCells(n) == 1) {
                           if (ci()(n) == patch_side_box.lower()(n) ||
                               ci()(n) == patch_side_box.upper()(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < depth; d++) {
                           double result = (*side_data)(ci(),d);

                           if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                              tbox::perr << "Test FAILED: ...."
                                   << " : side index = " << ci() << endl;
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
   }

   if (!test_failed) {   
      tbox::plog << "SideMultiblockTest Successful!" << endl;
   } else {
      tbox::perr << "Multiblock SideMultiblockTest FAILED: .\n" << endl;
   }

   solution.setNull();   // just to be anal...

   tbox::plog << "\nExiting SideMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   return(!test_failed);
}

void SideMultiblockTest::postprocessRefine(
   hier::Patch<NDIM>& fine,
   const hier::Patch<NDIM>& coarse,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& fine_box,
   const hier::IntVector<NDIM>& ratio) const
{
   pdat::SideDoubleConstantRefine<NDIM> ref_op;
                                                                                
   for (int i = 0; i < d_variables.getSize(); i++) {
                                                                                
      int id = hier::VariableDatabase<NDIM>::getDatabase()->
               mapVariableAndContextToIndex(d_variables[i], context);
                                                                                
      ref_op.refine(fine, coarse, id, id, fine_box, ratio);
   }
}

