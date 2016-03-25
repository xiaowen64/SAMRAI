//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/NodeDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2147 $
// Modified:    $LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: AMR communication tests for node-centered patch data
//

#include "NodeDataTest.h"


#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "NodeGeometry.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "NodeVariable.h"
#include "CommTester.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "Variable.h"
#include "VariableDatabase.h"

namespace SAMRAI {

NodeDataTest::NodeDataTest(
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

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }

   d_refine_option = refine_option;

   d_Acoef = 0.0;
   d_Bcoef = 0.0;
   d_Ccoef = 0.0;
   d_Dcoef = 0.0;

   d_finest_level_number = main_input_db->
                           getDatabase("GriddingAlgorithm")->
                           getInteger("max_levels") - 1;

   d_cart_grid_geometry = new geom::CartesianGridGeometry<NDIM>(
                          "CartesianGridGeometry",
                           main_input_db->getDatabase("CartesianGridGeometry"));

   setGridGeometry(d_cart_grid_geometry);

   readTestInput(main_input_db->getDatabase("NodePatchDataTest"));

}

NodeDataTest::~NodeDataTest()
{
}

void NodeDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   /*
    * Read coeeficients of linear profile to test interpolation.
    */
   if (db->keyExists("Acoef")) {
      d_Acoef = db->getDouble("Acoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Acoeff' found." << endl);
   }
   if (db->keyExists("Dcoef")) {
      d_Dcoef = db->getDouble("Dcoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Dcoef' found." << endl); 
   }
#if (NDIM > 1)
   if (db->keyExists("Bcoef")) {
      d_Bcoef = db->getDouble("Bcoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Bcoef' found." << endl); 
   }
#endif
#if (NDIM > 2)
   if (db->keyExists("Ccoef")) {
      d_Ccoef = db->getDouble("Ccoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Ccoef' found." << endl); 
   }
#endif
   
   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void NodeDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] = new pdat::NodeVariable<NDIM,double>(d_variable_src_name[i],
                                                d_variable_depth[i]);

      if (d_do_refine) {
         commtest->registerVariable(d_variables[i],
                                    d_variables[i],
                                    d_variable_src_ghosts[i], 
                                    d_variable_dst_ghosts[i], 
                                    d_cart_grid_geometry, 
                                    d_variable_refine_op[i]); 
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables[i], 
                                    d_variables[i],
                                    d_variable_src_ghosts[i],
                                    d_variable_dst_ghosts[i],
                                    d_cart_grid_geometry, 
                                    d_variable_coarsen_op[i]);   
      }

   }

}

void NodeDataTest::setLinearData(
   tbox::Pointer< pdat::NodeData<NDIM,double> > data, 
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const pdat::NodeIndex<NDIM> loweri(patch.getBox().lower(), (pdat::NodeIndex<NDIM>::Corner) 0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   for (pdat::NodeIterator<NDIM> ci(sbox); ci; ci++) {

      /*
       * Compute spatial location of node center and
       * set data to linear profile.
       */

      x = lowerx[0] + dx[0]*(ci()(0) - loweri(0));
      y = z = 0.;
#if (NDIM > 1)
      y = lowerx[1] + dx[1]*(ci()(1) - loweri(1));
#endif
#if (NDIM > 2)
      z = lowerx[2] + dx[2]*(ci()(2) - loweri(2));
#endif

      for (int d = 0; d < depth; d++) {
         (*data)(ci(),d) = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
      }

   }

}


void NodeDataTest::initializeDataOnPatch(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   char src_or_dst)
{
   (void) hierarchy;

   if (d_do_refine) {

      if (   (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") 
          || ( (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") 
              && (level_number == 0) ) ) {

         for (int i = 0; i < d_variables.getSize(); i++) {

            tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
               patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = node_data->getBox(); 

            setLinearData(node_data, dbox, patch);

         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box<NDIM> dbox = node_data->getGhostBox();

         setLinearData(node_data, dbox,
                             patch);//, hierarchy, level_number);

      } 
        
   }

}

void NodeDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::NodeData<NDIM,double> >& data,
   const hier::Box<NDIM>& interior,
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   const pdat::NodeIndex<NDIM> loweri(interior.lower(), (pdat::NodeIndex<NDIM>::Corner) 0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   for (pdat::NodeIterator<NDIM> ci(interior); ci; ci++) {

      /*
       * Compute spatial location of edge and
       * compare data to linear profile.
       */

      x = lowerx[0] + dx[0]*(ci()(0) - loweri(0));
      y = z = 0.;
#if (NDIM > 1)
      y = lowerx[1] + dx[1]*(ci()(1) - loweri(1));
#endif
#if (NDIM > 2)
      z = lowerx[2] + dx[2]*(ci()(2) - loweri(2));
#endif

      double value;
      for (int d = 0; d < depth; d++) {
         value = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
         if (!(tbox::MathUtilities<double>::equalEps((*data)(ci(),d), value))) {
            tbox::perr << "FAILED: -- patch interior not properly filled" << endl;
         }
      }

   }

}

void NodeDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) time;

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

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

      tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> patch_interior = node_data->getBox();
      checkPatchInteriorData(node_data, patch_interior, pgeom);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(node_data, fill_box, patch);
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(node_data, fill_box, patch);
      }
#endif

#if (NDIM == 3)
      /*
       * Set face boundary data.
       */
      for (int fi = 0; fi < num_face_bdry_boxes; fi++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(node_data, fill_box, patch);
      }
#endif

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
bool NodeDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering NodeDataTest::verifyResults..." << endl;
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

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox,
                             patch);//, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
            patch.getPatchData(d_variables[i], getDataContext());
         int depth = node_data->getDepth();
         hier::Box<NDIM> dbox = node_data->getGhostBox();

         for (pdat::NodeIterator<NDIM> ci(dbox); ci; ci++) {
            double correct = (*solution)(ci());
            for (int d = 0; d < depth; d++) {
               double result = (*node_data)(ci(),d);
               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...." 
                       << " : node index = " << ci() << endl;
                  tbox::perr << "    hier::Variable<NDIM> = " << d_variable_src_name[i]
                       << " : depth index = " << d << endl;
                  tbox::perr << "    result = " << result
                       << " : correct = " << correct << endl;
                  test_failed = true;
               } 
            }
         }

      }
      if (!test_failed) {   
         tbox::plog << "Node test Successful!" << endl;
      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting NodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);

}

}
