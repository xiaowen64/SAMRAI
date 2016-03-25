//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/mblkcomm/CellMultiblockTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2147 $
// Modified:    $LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: AMR communication tests for cell-centered patch data
//

#include "CellMultiblockTest.h"


#include "BoundaryBox.h"
#include "BlockPatchGeometry.h"
#include "CellDoubleConstantRefine.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "MultiblockTester.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"

using namespace SAMRAI;

CellMultiblockTest::CellMultiblockTest(
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

   readTestInput(main_input_db->getDatabase("CellMultiblockTest"));
}

CellMultiblockTest::~CellMultiblockTest()
{

}

void CellMultiblockTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

void CellMultiblockTest::registerVariables(MultiblockTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (MultiblockTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] =
         new pdat::CellVariable<NDIM,double>(d_variable_src_name[i],
                                             d_variable_depth[i]);

      commtest->registerVariable(d_variables[i],
                                 d_variables[i],
                                 d_variable_src_ghosts[i], 
                                 d_variable_dst_ghosts[i], 
                                 d_skel_grid_geometry[0], 
                                 d_variable_refine_op[i]); 

   }

}


void CellMultiblockTest::initializeDataOnPatch(
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

         tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box<NDIM> dbox = cell_data->getGhostBox(); 

         cell_data->fillAll((double)block_number);

      }
   }
}

void CellMultiblockTest::tagCellsToRefine(
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

void CellMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw_to_fill) const
{
   (void) time;

   tbox::Pointer< geom::BlockPatchGeometry<NDIM> >
      pgeom = patch.getPatchGeometry();

   const tbox::Array<hier::BoundaryBox<NDIM> > node_bdry = 
      pgeom->getCodimensionBoundaries(NDIM);
   const int num_node_bdry_boxes = node_bdry.getSize();

#if (NDIM > 1)
   const tbox::Array<hier::BoundaryBox<NDIM> > edge_bdry = 
      pgeom->getCodimensionBoundaries(NDIM - 1);
   const int num_edge_bdry_boxes = edge_bdry.getSize();
#endif

#if (NDIM == 3)
   const tbox::Array<hier::BoundaryBox<NDIM> > face_bdry = 
      pgeom->getCodimensionBoundaries(NDIM - 2);
   const int num_face_bdry_boxes = face_bdry.getSize();
#endif

   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
         patch.getPatchData(d_variables[i], getDataContext());

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         if (!node_bdry[ni].getIsMultiblockSingularity()) {
            cell_data->fillAll((double)(node_bdry[ni].getLocationIndex()+100),
                               fill_box);
         }
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         if (!edge_bdry[ei].getIsMultiblockSingularity()) {
            cell_data->fillAll((double)(edge_bdry[ei].getLocationIndex()+100),
                               fill_box);
         }
      }
#endif

#if (NDIM == 3)
      /*
       * Set face boundary data.
       */
      for (int fi = 0; fi < num_face_bdry_boxes; fi++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         if (!face_bdry[fi].getIsMultiblockSingularity()) {
            cell_data->fillAll((double)(face_bdry[fi].getLocationIndex()+100),
                               fill_box);
         }
      }
#endif 

   }

}

void CellMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch<NDIM>& patch,
   tbox::List<xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>&
      sing_patches,
   const hier::Box<NDIM>& fill_box,
   const hier::BoundaryBox<NDIM>& bbox)
{
   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> sing_fill_box(cell_data->getGhostBox() * fill_box);
      cell_data->fillAll(0.0, sing_fill_box);

      int depth = cell_data->getDepth();

      /*
       * If sing_patches is not empty, that means there is enhanced
       * connectivity, and we get data from other blocks
       */

      if (sing_patches.size()) {

         for (tbox::List
              <xfer::MultiblockRefineSchedule<NDIM>::SingularityPatch>::
              Iterator sp(sing_patches); sp; sp++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > sing_data =
               sp().d_patch->getPatchData(d_variables[i], getDataContext());
            int sing_neighbor_id = sp().d_id;
            for (pdat::CellIterator<NDIM> ci(sing_fill_box); ci; ci++) {
               for (int d = 0; d < depth; d++) {
                  (*cell_data)(ci(),d) += sing_neighbor_id;
               }
            }
         }

         for (pdat::CellIterator<NDIM> ci(sing_fill_box); ci; ci++) {
            for (int d = 0; d < depth; d++) {
               (*cell_data)(ci(),d) /= sing_patches.size();
            }
         }

      /*
       * In cases of reduced connectivity, there are no other blocks
       * from which to acquire data.
       */

      } else {

         cell_data->fillAll(
            (double)bbox.getLocationIndex()+200.0, fill_box);  
 
      }
   }
}

void CellMultiblockTest::postprocessRefine(
   hier::Patch<NDIM>& fine,
   const hier::Patch<NDIM>& coarse,
   const tbox::Pointer<hier::VariableContext>& context,
   const hier::Box<NDIM>& fine_box,
   const hier::IntVector<NDIM>& ratio) const
{
   pdat::CellDoubleConstantRefine<NDIM> ref_op;

   for (int i = 0; i < d_variables.getSize(); i++) {

      int id = hier::VariableDatabase<NDIM>::getDatabase()->
               mapVariableAndContextToIndex(d_variables[i], context);

      ref_op.refine(fine, coarse, id, id, fine_box, ratio);
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
bool CellMultiblockTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::MultiblockPatchHierarchy<NDIM> > hierarchy, 
   int level_number,
   int block_number)
{

   tbox::plog << "\nEntering CellMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl; 

   hier::IntVector<NDIM> tgcw(0);
   for (int i = 0; i < d_variables.getSize(); i++) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
                                  getGhostCellWidth());
   }
   hier::Box<NDIM> pbox = patch.getBox();

   tbox::Pointer< pdat::CellData<NDIM,double> > solution =
      new pdat::CellData<NDIM,double>(pbox, 1, tgcw);

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

      tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
         patch.getPatchData(d_variables[i], getDataContext());
      int depth = cell_data->getDepth();

      for (pdat::CellIterator<NDIM> ci(pbox); ci; ci++) {
         for (int d = 0; d < depth; d++) {
            double result = (*cell_data)(ci(),d);

            if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
               tbox::perr << "Test FAILED: ...."
                    << " : cell index = " << ci() << endl;
               tbox::perr << "    Variable = " << d_variable_src_name[i]
                    << " : depth index = " << d << endl;
               tbox::perr << "    result = " << result
                    << " : correct = " << correct << endl;
               test_failed = true;
            }
         }
      }

      hier::Box<NDIM> gbox = cell_data->getGhostBox();

      for (tbox::List<hier::MultiblockPatchHierarchy<NDIM>::Neighbor>::
           Iterator ne(neighbors); ne; ne++) {

         correct = ne().d_id;

         hier::BoxList<NDIM> neighbor_ghost(ne().d_translated_domain);
         neighbor_ghost.refine(ratio);
         neighbor_ghost.intersectBoxes(gbox);

         for (hier::BoxList<NDIM>::Iterator ng(neighbor_ghost); ng; ng++) {

            for (pdat::CellIterator<NDIM> ci(ng()); ci; ci++) {
               for (int d = 0; d < depth; d++) {
                  double result = (*cell_data)(ci(),d);

                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     tbox::perr << "Test FAILED: ...."
                          << " : cell index = " << ci() << endl;
                     tbox::perr << "    Variable = " << d_variable_src_name[i]
                          << " : depth index = " << d << endl;
                     tbox::perr << "    result = " << result
                          << " : correct = " << correct << endl;
                     test_failed = true;
                  }
               }
            } 
         } 
      } 

      tbox::Pointer< geom::BlockPatchGeometry<NDIM> > pgeom =
         patch.getPatchGeometry();

      for (int b = 0; b < NDIM; b++) {
         tbox::Array<hier::BoundaryBox<NDIM> > bdry =
            pgeom->getCodimensionBoundaries(b+1);

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

            for (pdat::CellIterator<NDIM> ci(fill_box); ci; ci++) {
               for (int d = 0; d < depth; d++) {
                  double result = (*cell_data)(ci(),d);

                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     tbox::perr << "Test FAILED: ...."
                          << " : cell index = " << ci() << endl;
                     tbox::perr << "    Variable = " << d_variable_src_name[i]
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

   if (!test_failed) {   
      tbox::plog << "CellMultiblockTest Successful!" << endl;
   } else {
      tbox::perr << "Multiblock CellMultiblockTest FAILED: \n" << endl;
   }

   solution.setNull();   // just to be anal...

   tbox::plog << "\nExiting CellMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   return(!test_failed);
}
