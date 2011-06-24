/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data 
 *
 ************************************************************************/

#include "NodeMultiblockTest.h"

#include "SAMRAI/geom/SAMRAITransferOperatorRegistry.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/NodeVariable.h"

#include "MultiblockTester.h"

using namespace SAMRAI;

NodeMultiblockTest::NodeMultiblockTest(
   const string& object_name,
   const tbox::Dimension& dim,
   tbox::Pointer<tbox::Database> main_input_db,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option):
   PatchMultiblockTestStrategy(dim),
   d_dim(dim)
{
   NULL_USE(do_refine);
   NULL_USE(do_coarsen);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!main_input_db.isNull());
   TBOX_ASSERT(!refine_option.empty());
#endif

   d_object_name = object_name;

   d_refine_option = refine_option;

   d_finest_level_number = main_input_db->
      getDatabase("PatchHierarchy")->
      getInteger("max_levels") - 1;

   int num_blocks = main_input_db->getDatabase("BlockGridGeometry")->
      getInteger("num_blocks");

   char geom_name[32];

   sprintf(geom_name, "BlockGridGeometry");

   if (main_input_db->keyExists(geom_name)) {
      getGridGeometry() = new hier::GridGeometry(
               dim,
               geom_name,
               tbox::Pointer<hier::TransferOperatorRegistry>(
                  new geom::SAMRAITransferOperatorRegistry(dim)),
               main_input_db->getDatabase(geom_name));

   } else {
      TBOX_ERROR("NodeMultiblockTest: could not find entry `"
         << geom_name << "' in input.");
   }

   readTestInput(main_input_db->getDatabase("NodeMultiblockTest"));
}

NodeMultiblockTest::~NodeMultiblockTest()
{
}

void NodeMultiblockTest::readTestInput(
   tbox::Pointer<tbox::Database> db)
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

void NodeMultiblockTest::registerVariables(
   MultiblockTester* commtest)
{
   TBOX_ASSERT(commtest != (MultiblockTester *)NULL);

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] =
         new pdat::NodeVariable<double>(d_dim,
                                        d_variable_src_name[i],
                                        d_variable_depth[i]);

      commtest->registerVariable(d_variables[i],
         d_variables[i],
         d_variable_src_ghosts[i],
         d_variable_dst_ghosts[i],
         getGridGeometry(),
         d_variable_refine_op[i]);

   }

}

void NodeMultiblockTest::initializeDataOnPatch(
   hier::Patch& patch,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   int level_number,
   const hier::BlockId& block_id,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(src_or_dst);

   if ((d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
       || ((d_refine_option == "INTERIOR_FROM_COARSER_LEVEL")
           && (level_number < d_finest_level_number))) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer<pdat::NodeData<double> > node_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box dbox = node_data->getGhostBox();

         node_data->fillAll((double)block_id.getBlockValue());

      }
   }
}

void NodeMultiblockTest::tagCellsToRefine(
   hier::Patch& patch,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   int level_number,
   int tag_index)
{
   (void)hierarchy;

   /*
    * Base class sets tags in box array for each level.
    */
   tagCellsInInputBoxes(patch, level_number, tag_index);

}

void NodeMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw_to_fill) const
{
   (void)time;

   tbox::Pointer<hier::PatchGeometry>
   pgeom = patch.getPatchGeometry();

   const tbox::Array<hier::BoundaryBox> node_bdry =
      pgeom->getCodimensionBoundaries(d_dim.getValue());
   const int num_node_bdry_boxes = node_bdry.getSize();

   tbox::Array<hier::BoundaryBox> edge_bdry;
   int num_edge_bdry_boxes = 0;
   if (d_dim > tbox::Dimension(1)) {
      edge_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 1);
      num_edge_bdry_boxes = edge_bdry.getSize();
   }

   tbox::Array<hier::BoundaryBox> face_bdry;
   int num_face_bdry_boxes = 0;
   if (d_dim == tbox::Dimension(3)) {
      face_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 2);
      num_face_bdry_boxes = face_bdry.getSize();
   }

   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer<pdat::NodeData<double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());

      /*
       * Set node boundary data.
       */
      for (int nb = 0; nb < num_node_bdry_boxes; nb++) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[nb],
               patch.getBox(),
               gcw_to_fill);

         hier::Box patch_node_box =
            pdat::NodeGeometry::toNodeBox(patch.getBox());
         if (!node_bdry[nb].getIsMultiblockSingularity()) {
            for (pdat::NodeIterator ni(fill_box); ni; ni++) {
               if (!patch_node_box.contains(ni())) {
                  for (int d = 0; d < node_data->getDepth(); d++) {
                     (*node_data)(ni(), d) =
                        (double)(node_bdry[nb].getLocationIndex() + 100);
                  }
               }
            }
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int eb = 0; eb < num_edge_bdry_boxes; eb++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[eb],
                  patch.getBox(),
                  gcw_to_fill);

            hier::Box patch_node_box =
               pdat::NodeGeometry::toNodeBox(patch.getBox());
            hier::Index plower(patch_node_box.lower());
            hier::Index pupper(patch_node_box.upper());

            if (!edge_bdry[eb].getIsMultiblockSingularity()) {
               for (pdat::NodeIterator ni(fill_box); ni; ni++) {
                  if (!patch_node_box.contains(ni())) {
                     bool use_index = true;
                     for (int n = 0; n < d_dim.getValue(); n++) {
                        if (edge_bdry[eb].getBox().numberCells(n) == 1) {
                           if (ni() (n) == plower(n) || ni() (n) ==
                               pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < node_data->getDepth(); d++) {
                           (*node_data)(ni(), d) =
                              (double)(edge_bdry[eb].getLocationIndex() + 100);
                        }
                     }
                  }
               }
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fb = 0; fb < num_face_bdry_boxes; fb++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fb],
                  patch.getBox(),
                  gcw_to_fill);

            hier::Box patch_node_box =
               pdat::NodeGeometry::toNodeBox(patch.getBox());
            hier::Index plower(patch_node_box.lower());
            hier::Index pupper(patch_node_box.upper());

            if (!face_bdry[fb].getIsMultiblockSingularity()) {
               for (pdat::NodeIterator ni(fill_box); ni; ni++) {
                  if (!patch_node_box.contains(ni())) {
                     bool use_index = true;
                     for (int n = 0; n < d_dim.getValue(); n++) {
                        if (face_bdry[fb].getBox().numberCells(n) == 1) {
                           if (ni() (n) == plower(n) || ni() (n) ==
                               pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < node_data->getDepth(); d++) {
                           (*node_data)(ni(), d) =
                              (double)(face_bdry[fb].getLocationIndex() + 100);
                        }
                     }
                  }
               }
            }
         }
      }

   }

}

void NodeMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   tbox::List<tbox::Pointer<hier::Patch> >& sing_patches,
   const hier::Box& fill_box,
   const hier::BoundaryBox& bbox)
{
   for (int i = 0; i < d_variables.getSize(); i++) {

      tbox::Pointer<pdat::NodeData<double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box pbox(pdat::NodeGeometry::toNodeBox(patch.getBox()));

      hier::Index plower(pbox.lower());
      hier::Index pupper(pbox.upper());

      hier::Box sing_fill_box(node_data->getGhostBox() * fill_box);

      int depth = node_data->getDepth();

      for (pdat::NodeIterator ni(sing_fill_box); ni; ni++) {
         bool use_index = true;
         for (int n = 0; n < d_dim.getValue(); n++) {
            if (bbox.getBox().numberCells(n) == 1) {
               if (ni() (n) == plower(n) || ni() (n) == pupper(n)) {
                  use_index = false;
                  break;
               }
            }
         }
         if (use_index) {
            for (int d = 0; d < depth; d++) {
               (*node_data)(ni(), d) = 0.0;
            }
         }
      }

      /*
       * If sing_patches is not empty, that means there is enhanced
       * connectivity, and we get data from other blocks
       */

      if (sing_patches.size()) {

         for (tbox::List<tbox::Pointer<hier::Patch> >::Iterator
              sp(sing_patches); sp; sp++) {
            tbox::Pointer<pdat::NodeData<double> > sing_data =
               sp()->getPatchData(d_variables[i], getDataContext());
            tbox::Pointer<hier::PatchGeometry> patch_geom =
               sp()->getPatchGeometry();
            int sing_neighbor_id = 
               sp()->getMappedBox().getBlockId().getBlockValue();
            for (pdat::NodeIterator ci(sing_fill_box); ci; ci++) {
               bool use_index = true;
               for (int n = 0; n < d_dim.getValue(); n++) {
                  if (bbox.getBox().numberCells(n) == 1) {
                     if (ci() (n) == plower(n) || ci() (n) == pupper(n)) {
                        use_index = false;
                        break;
                     }
                  }
               }
               if (use_index) {
                  for (int d = 0; d < depth; d++) {
                     (*node_data)(ci(), d) += sing_neighbor_id;
                  }
               }
            }
         }

         for (pdat::NodeIterator ci(sing_fill_box); ci; ci++) {
            bool use_index = true;
            for (int n = 0; n < d_dim.getValue(); n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if (ci() (n) == plower(n) || ci() (n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(ci(), d) /= sing_patches.size();
               }
            }
         }

         /*
          * In cases of reduced connectivity, there are no other blocks
          * from which to acquire data.
          */

      } else {

         for (pdat::NodeIterator ci(sing_fill_box); ci; ci++) {
            bool use_index = true;
            for (int n = 0; n < d_dim.getValue(); n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if (ci() (n) == plower(n) || ci() (n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(ci(),
                               d) = (double)bbox.getLocationIndex() + 200.0;
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
bool NodeMultiblockTest::verifyResults(
   const hier::Patch& patch,
   const tbox::Pointer<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::BlockId& block_id)
{

   tbox::plog << "\nEntering NodeMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl;

   hier::IntVector tgcw(d_dim, 0);
   for (int i = 0; i < d_variables.getSize(); i++) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
         getGhostCellWidth());
   }
   hier::Box pbox = patch.getBox();

   tbox::Pointer<pdat::NodeData<double> > solution(
      new pdat::NodeData<double>(pbox, 1, tgcw));

   hier::Box tbox(pbox);
   tbox.grow(tgcw);

   const tbox::List<hier::GridGeometry::Neighbor>& neighbors =
      hierarchy->getGridGeometry()->getNeighbors(block_id.getBlockValue());
   hier::BoxList singularity(
      hierarchy->getGridGeometry()->getSingularityBoxList(block_id.getBlockValue()));

   hier::IntVector ratio =
      hierarchy->getPatchLevel(level_number)->getRatioToLevelZero();

   singularity.refine(ratio);

   bool test_failed = false;

   for (int i = 0; i < d_variables.getSize(); i++) {

      double correct = (double)block_id.getBlockValue();

      tbox::Pointer<pdat::NodeData<double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());
      int depth = node_data->getDepth();

      hier::Box interior_box(pbox);
      interior_box.grow(hier::IntVector(d_dim, -1));

      for (pdat::NodeIterator ci(interior_box); ci; ci++) {
         for (int d = 0; d < depth; d++) {
            double result = (*node_data)(ci(), d);

            if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
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

      tbox::Pointer<hier::PatchGeometry> pgeom =
         patch.getPatchGeometry();

      hier::Box gbox = node_data->getGhostBox();

      hier::Box patch_node_box =
         pdat::NodeGeometry::toNodeBox(pbox);

      hier::BoxList sing_node_boxlist;
      for (hier::BoxList::Iterator si(singularity); si; si++) {
         sing_node_boxlist.addItem(pdat::NodeGeometry::toNodeBox(si()));
      }

      for (tbox::List<hier::GridGeometry::Neighbor>::
           Iterator ne(neighbors); ne; ne++) {

         correct = ne().getBlockNumber();

         hier::BoxList neighbor_ghost(ne().getTranslatedDomain());

         hier::BoxList neighbor_node_ghost;
         for (hier::BoxList::Iterator nn(neighbor_ghost); nn; nn++) {
            neighbor_node_ghost.addItem(
               pdat::NodeGeometry::toNodeBox(nn()));
         }

         neighbor_node_ghost.refine(ratio);

         neighbor_node_ghost.intersectBoxes(
            pdat::NodeGeometry::toNodeBox(gbox));

         neighbor_node_ghost.removeIntersections(sing_node_boxlist);

         for (hier::BoxList::Iterator ng(neighbor_node_ghost); ng; ng++) {

            for (hier::BoxIterator ci(ng()); ci; ci++) {
               if (!patch_node_box.contains(ci())) {
                  pdat::NodeIndex ni(ci(), hier::IntVector(d_dim, 0));
                  for (int d = 0; d < depth; d++) {
                     double result = (*node_data)(ni, d);

                     if (!tbox::MathUtilities<double>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : node index = " << ni << endl;
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

      for (int b = 0; b < d_dim.getValue(); b++) {
         tbox::Array<hier::BoundaryBox> bdry =
            pgeom->getCodimensionBoundaries(b + 1);

         for (int k = 0; k < bdry.size(); k++) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bdry[k],
                  patch.getBox(),
                  tgcw);
            fill_box = fill_box * gbox;

            if (bdry[k].getIsMultiblockSingularity()) {
               correct = 0.0;

               int num_sing_neighbors = 0;
               for (tbox::List
                    <hier::GridGeometry::Neighbor>::
                    Iterator ns(neighbors); ns; ns++) {
                  if (ns().isSingularity()) {
                     hier::BoxList neighbor_ghost(
                        ns().getTranslatedDomain());
                     neighbor_ghost.refine(ratio);
                     neighbor_ghost.intersectBoxes(fill_box);
                     if (neighbor_ghost.size()) {
                        num_sing_neighbors++;
                        correct += block_id.getBlockValue();
                     }
                  }
               }

               if (num_sing_neighbors == 0) {

                  correct = (double)bdry[k].getLocationIndex() + 200.0;

               } else {

                  correct /= (double)num_sing_neighbors;

               }

            } else {
               correct = (double)(bdry[k].getLocationIndex() + 100);
            }

            for (pdat::NodeIterator ci(fill_box); ci; ci++) {

               if (!patch_node_box.contains(ci())) {

                  bool use_index = true;
                  for (int n = 0; n < d_dim.getValue(); n++) {
                     if (bdry[k].getBox().numberCells(n) == 1) {
                        if (ci() (n) == patch_node_box.lower() (n) ||
                            ci() (n) == patch_node_box.upper() (n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }

                  if (use_index) {
                     for (int d = 0; d < depth; d++) {
                        double result = (*node_data)(ci(), d);

                        if (!tbox::MathUtilities<double>::equalEps(correct,
                               result)) {
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

   return !test_failed;
}
