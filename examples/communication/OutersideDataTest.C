//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/OutersideDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: AMR communication tests for outerside-centered patch data
//

#include "OutersideDataTest.h"


#include "ArrayData.h"
#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CommTester.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "SideVariable.h"
#include "OutersideGeometry.h"
#include "OutersideVariable.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"

namespace SAMRAI {

OutersideDataTest::OutersideDataTest(
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
      TBOX_ERROR("There is no refine test for Outerside data type, because\n"
		 <<"Outerside refinement does not exist at this time.");
      /*
	The refine codes are still kept in this class in case we
	implement Outerside refinement in the future.
      */
   }
   
   d_refine_option = refine_option;

   d_test_direction.resizeArray(0);
   d_use_fine_value_at_interface.resizeArray(0);

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
   readTestInput(main_input_db->getDatabase("OutersidePatchDataTest"));

}

OutersideDataTest::~OutersideDataTest()
{
}

void OutersideDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!db.isNull());
#endif

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));

   tbox::Pointer<tbox::Database> var_data = db->getDatabase("VariableData");
   tbox::Array<string> var_keys = var_data->getAllKeys();
   int nkeys = var_keys.getSize();

   d_test_direction.resizeArray(nkeys);
   d_use_fine_value_at_interface.resizeArray(nkeys);

   for (int i = 0; i < nkeys; i++) {
      tbox::Pointer<tbox::Database> var_db = var_data->getDatabase(var_keys[i]);
     
      if (var_db->keyExists("test_direction")) {
         d_test_direction[i] = var_db->getInteger("test_direction");
      } else {
         d_test_direction[i] = -1;
      } 

      if (var_db->keyExists("use_fine_value_at_interface")) {
         d_use_fine_value_at_interface[i] = 
            var_db->getBool("use_fine_value_at_interface");
      } else {
         d_use_fine_value_at_interface[i] = true;
      }
      
   }

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

}

void OutersideDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables_src.resizeArray(nvars);
   d_variables_dst.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables_src[i] = 
         new pdat::OutersideVariable<NDIM,double>(d_variable_src_name[i],
                                  d_variable_depth[i]);
      d_variables_dst[i] = 
         new pdat::SideVariable<NDIM,double>(d_variable_dst_name[i],
                                  d_variable_depth[i],
                                  d_use_fine_value_at_interface[i],
                                  d_test_direction[i]);

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

void OutersideDataTest::initializeDataOnPatch(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   char src_or_dst)
{
   (void) hierarchy;
   hier::VariableDatabase<NDIM>* variable_db =
      hier::VariableDatabase<NDIM>::getDatabase();
   variable_db->printClassData();
   tbox::Array<tbox::Pointer<hier::Variable<NDIM> > > &variables =
      src_or_dst == 's' ? d_variables_src : d_variables_dst;
   
   if (d_do_refine) {

      if (   (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") 
          || ( (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") 
              && (level_number == 0) ) ) {

         for (int i = 0; i < variables.getSize(); i++) {

	    tbox::Pointer<hier::PatchData<NDIM> > data =
               patch.getPatchData(variables[i], getDataContext());

	    TBOX_ASSERT( !data.isNull() );

            tbox::Pointer<pdat::OutersideData<NDIM,double> > oside_data = data;
	    tbox::Pointer<pdat::SideData<NDIM,double> > side_data = data;

            hier::Box<NDIM> dbox = data->getBox(); 

	    if ( !side_data.isNull() ) {
	       setLinearData(side_data, dbox, patch);
	    }
	    if ( !oside_data.isNull() ) {
	       setLinearData(oside_data, dbox, patch);
	    }
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < variables.getSize(); i++) {

         tbox::Pointer<hier::PatchData<NDIM> > data =
            patch.getPatchData(variables[i], getDataContext());

	 TBOX_ASSERT( !data.isNull() );

         tbox::Pointer<pdat::OutersideData<NDIM,double> > oside_data = data;
         tbox::Pointer<pdat::SideData<NDIM,double> > side_data = data;

         hier::Box<NDIM> dbox = data->getGhostBox();

         if ( !side_data.isNull() ) {
            setLinearData(side_data, dbox, patch);
         }
         if ( !oside_data.isNull() ) {
            setLinearData(oside_data, dbox, patch);
         }
      } 
        
   }

}

void OutersideDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::OutersideData<NDIM,double> >& data,
   const hier::Box<NDIM>& interior,
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x=0., y=0., z=0.;

   const int depth = data->getDepth();

   for (int axis = 0; axis < NDIM; axis++) {
         const pdat::SideIndex<NDIM> loweri(interior.lower(), axis, 0);
         for (pdat::SideIterator<NDIM> si(interior, axis); si; si++) {

            /*
             * Compute spatial location of face and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0]*(si()(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0]*(si()(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
#if (NDIM > 1)
            if (axis == 1) {
               y = lowerx[1] + dx[1]*(si()(1) - loweri(1));
            } else {
               y = lowerx[1] + dx[1]*(si()(1) - loweri(1) + 0.5);
            }
#endif
#if (NDIM > 2)
            if (axis == 2) {
               z = lowerx[2] + dx[2]*(si()(2) - loweri(2));
            } else {
               z = lowerx[2] + dx[2]*(si()(2) - loweri(2) + 0.5);
            }
#endif

            double value;
            for (int d = 0; d < depth; d++) {
               value = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;

               if (!(tbox::MathUtilities<double>::equalEps((*data)(si(),d), value))) {
                  tbox::perr << "FAILED: -- patch interior not properly filled" << endl;
               }
            }
         }
   }
}

void OutersideDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) time;
   TBOX_ERROR("Only coarsen operations can be done with this test.\n"
	      << "Coarsen operations should not need physical bc.\n");
}

void OutersideDataTest::setLinearData(
   tbox::Pointer< pdat::SideData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >
      pgeom = patch.getPatchGeometry();
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   hier::IntVector<NDIM> directions = data->getDirectionVector();

   for (int axis = 0; axis < NDIM; axis++) {
      if (directions(axis)) {
         const pdat::SideIndex<NDIM> loweri(patch.getBox().lower(), axis, 0);
         for (pdat::SideIterator<NDIM> ei(sbox, axis); ei; ei++) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0]*(ei()(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0]*(ei()(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
#if (NDIM > 1)
            if (axis == 1) {
               y = lowerx[1] + dx[1]*(ei()(1) - loweri(1));
            } else {
               y = lowerx[1] + dx[1]*(ei()(1) - loweri(1) + 0.5);
            }
#endif
#if (NDIM > 2)
            if (axis == 2) {
               z = lowerx[2] + dx[2]*(ei()(2) - loweri(2));
            } else {
               z = lowerx[2] + dx[2]*(ei()(2) - loweri(2) + 0.5);
            }
#endif

            for (int d = 0; d < depth; d++) {
               (*data)(ei(),d) = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
            }

         }
      }
   }
}

void OutersideDataTest::setLinearData(
   tbox::Pointer< pdat::OutersideData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >
      pgeom = patch.getPatchGeometry();
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x=0., y=0., z=0.;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   for (int axis = 0; axis < NDIM; axis++) {
      for (int s = 0; s < 2; s++) {
         const hier::Box<NDIM> databox = data->getArrayData(axis, s).getBox();
         const pdat::SideIndex<NDIM> loweri(patch.getBox().lower(), axis, 0);
	 for (hier::Box<NDIM>::Iterator bi(databox); bi; bi++) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0]*(bi()(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0]*(bi()(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
#if (NDIM > 1)
            if (axis == 1) {
               y = lowerx[1] + dx[1]*(bi()(1) - loweri(1));
            } else {
               y = lowerx[1] + dx[1]*(bi()(1) - loweri(1) + 0.5);
            }
#endif
#if (NDIM > 2)
            if (axis == 2) {
               z = lowerx[2] + dx[2]*(bi()(2) - loweri(2));
            } else {
               z = lowerx[2] + dx[2]*(bi()(2) - loweri(2) + 0.5);
            }
#endif

            double value =  d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
            pdat::SideIndex<NDIM> si(bi(), axis, 0);
            for (int d = 0 ; d < depth; d++) {
               (*data)(si,s,d) = value;
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

bool OutersideDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OutersideDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl; 

      hier::IntVector<NDIM> tgcw(0);
      for (int i = 0; i < d_variables_dst.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
                                     getGhostCellWidth());
      }
      hier::Box<NDIM> pbox = patch.getBox();

      tbox::Pointer< pdat::SideData<NDIM,double> > solution =
         new pdat::SideData<NDIM,double>(pbox, 1, tgcw);

      hier::Box<NDIM> tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox, patch);//, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables_dst.getSize(); i++) {
         
         tbox::Pointer< pdat::SideData<NDIM,double> > side_data =
            patch.getPatchData(d_variables_dst[i], getDataContext());
         int depth = side_data->getDepth();
         hier::Box<NDIM> dbox = side_data->getGhostBox();

         hier::IntVector<NDIM> directions = side_data->getDirectionVector();

         for (int id = 0; id < NDIM; id++) {
           if (directions(id)) {
               for (pdat::SideIterator<NDIM> si(dbox, id); si; si++) {
                  double correct = (*solution)(si());
                  for (int d = 0; d < depth; d++) {
                     double result = (*side_data)(si(),d);
                     if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                        tbox::perr << "Test FAILED: ...." 
                             << " : side_data index = " << si() << endl;
                        tbox::perr << "    hier::Variable<NDIM> = "
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
      if (!test_failed) {   
         tbox::plog << "Outerside test Successful!" << endl;
      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting OutersidedataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);

}

}
