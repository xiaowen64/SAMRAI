//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/OuternodeDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: AMR communication tests for node-centered patch data
//

#include "OuternodeDataTest.h"


#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "NodeData.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "NodeVariable.h"
#include "OuternodeVariable.h"
#include "PatchData.h"
#include "CommTester.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "Variable.h"
#include "VariableDatabase.h"

namespace SAMRAI {

OuternodeDataTest::OuternodeDataTest(
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
   if ( d_do_refine ) {
      TBOX_ERROR("There is no refine test for Outernode data type, because\n"
		 <<"Outernode refinement does not exist at this time.");
      /*
	The refine codes are still kept in this class in case we
	somehow define Outernode refinement in the future.
      */
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

   readTestInput(main_input_db->getDatabase("OuternodePatchDataTest"));

}

OuternodeDataTest::~OuternodeDataTest()
{
}

void OuternodeDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

void OuternodeDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables_src.resizeArray(nvars);
   d_variables_dst.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables_src[i] =
	new pdat::OuternodeVariable<NDIM,double>(d_variable_src_name[i],
				      d_variable_depth[i]);
      d_variables_dst[i] =
	new pdat::NodeVariable<NDIM,double>(d_variable_dst_name[i],
				 d_variable_depth[i]);

      if (d_do_refine) {
         commtest->registerVariable(d_variables_src[i],
                                    d_variables_dst[i],
                                    d_variable_src_ghosts[i], 
                                    d_variable_dst_ghosts[i], 
                                    d_cart_grid_geometry, 
                                    d_variable_refine_op[i]); 
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables_src[i], 
                                    d_variables_dst[i],
                                    d_variable_src_ghosts[i],
                                    d_variable_dst_ghosts[i],
                                    d_cart_grid_geometry, 
                                    d_variable_coarsen_op[i]);   
      }

   }

}


void OuternodeDataTest::setLinearData(
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



void OuternodeDataTest::setLinearData(
   tbox::Pointer< pdat::OuternodeData<NDIM,double> > data, 
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
   TBOX_ASSERT(box == patch.getBox() );
   if ( box != data->getBox() ) {
      TBOX_ERROR("Box is not identical to data box, which is\n"
		 <<"required for testing Outernode communication.");
   }
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const pdat::NodeIndex<NDIM> loweri(patch.getBox().lower(), (pdat::NodeIndex<NDIM>::Corner) 0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   int n, s;
   for ( n=0; n<NDIM; ++n ) {
      for ( s=0; s<2; ++s ) {
	 const hier::Box<NDIM> databox = data->getDataBox(n,s);
	 for (hier::Box<NDIM>::Iterator bi(databox); bi; bi++) {

	    /*
	     * Compute spatial location of node center and
	     * set data to linear profile.
	     */

	    x = lowerx[0] + dx[0]*(bi()(0) - loweri(0));
	    y = z = 0.;
#if (NDIM > 1)
	    y = lowerx[1] + dx[1]*(bi()(1) - loweri(1));
#endif
#if (NDIM > 2)
	    z = lowerx[2] + dx[2]*(bi()(2) - loweri(2));
#endif

	    pdat::NodeIndex<NDIM> ni(bi(), (pdat::NodeIndex<NDIM>::Corner)0);
	    for (int d = 0; d < depth; d++) {
	       (*data)(ni,d) = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
	    }
	 }
      }

   }

}


void OuternodeDataTest::initializeDataOnPatch(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   char src_or_dst)
{
   (void) hierarchy;
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
   variable_db->printClassData();
   tbox::Array<tbox::Pointer<hier::Variable<NDIM> > > &variables =
      src_or_dst == 's' ? d_variables_src : d_variables_dst;

   if (d_do_refine) {

      if (   (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") 
          || ( (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") 
              && (level_number == 0) ) ) {

         for (int i = 0; i < variables.getSize(); i++) {

	    tbox::Pointer<hier::PatchData<NDIM> > data = patch.getPatchData(variables[i],
							 getDataContext());
	    TBOX_ASSERT( !data.isNull() );

	    tbox::Pointer<pdat::OuternodeData<NDIM,double> > onode_data = data;
	    tbox::Pointer<pdat::NodeData<NDIM,double> > node_data = data;

            hier::Box<NDIM> dbox = data->getBox(); 

	    if ( !node_data.isNull() ) {
	       setLinearData(node_data, dbox, patch);
	    }
	    if ( !onode_data.isNull() ) {
	       setLinearData(onode_data, dbox, patch);
	    }

         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < variables.getSize(); i++) {

         tbox::Pointer<hier::PatchData<NDIM> > data = patch.getPatchData(variables[i],
						      getDataContext());
	 TBOX_ASSERT( !data.isNull() );
         tbox::Pointer<pdat::OuternodeData<NDIM,double> > onode_data = data;
         tbox::Pointer<pdat::NodeData<NDIM,double> > node_data = data;

         hier::Box<NDIM> dbox = data->getGhostBox();

	 if ( !node_data.isNull() ) {
	    setLinearData(node_data, dbox, patch);
	 }
	 if ( !onode_data.isNull() ) {
	    setLinearData(onode_data, dbox, patch);
	 }

      } 
        
   }

}

void OuternodeDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::OuternodeData<NDIM,double> >& data,
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

void OuternodeDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   TBOX_ERROR("Only coarsen operations can be done with this test.\n"
	      << "Coarsen operations should not need physical bc.\n");
   (void) time;

   return;
}

/*
*************************************************************************
*                                                                       *
* Verify results of communication operations.  This test must be        *
* consistent with data initialization and boundary operations above.    *
*                                                                       *
*************************************************************************
*/
bool OuternodeDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OuternodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl; 

      hier::IntVector<NDIM> tgcw(0);
      for (int i = 0; i < d_variables_dst.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
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

      for (int i = 0; i < d_variables_dst.getSize(); i++) {

         tbox::Pointer< pdat::NodeData<NDIM,double> > node_data =
            patch.getPatchData(d_variables_dst[i], getDataContext());
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
         tbox::plog << "Outernode test Successful!" << endl;
      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting OuternodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);
}

}
