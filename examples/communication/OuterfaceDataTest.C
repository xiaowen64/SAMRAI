//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/OuterfaceDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:     $Name:  $
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: AMR communication tests for outerface-centered patch data
//

#include "OuterfaceDataTest.h"


#include "ArrayData.h"
#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CommTester.h"
#include "FaceGeometry.h"
#include "FaceIndex.h"
#include "FaceIterator.h"
#include "FaceVariable.h"
#include "OuterfaceGeometry.h"
#include "OuterfaceVariable.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"

namespace SAMRAI {

OuterfaceDataTest::OuterfaceDataTest(
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

   readTestInput(main_input_db->getDatabase("OuterfacePatchDataTest"));

}

OuterfaceDataTest::~OuterfaceDataTest()
{
}

void OuterfaceDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

void OuterfaceDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables_src.resizeArray(nvars);
   d_variables_dst.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables_src[i] = 
         new pdat::OuterfaceVariable<NDIM,double>(d_variable_src_name[i],
                                                  d_variable_depth[i]);

      d_variables_dst[i] = 
         new pdat::FaceVariable<NDIM,double>(d_variable_dst_name[i],
                                             d_variable_depth[i],
                                             d_use_fine_value_at_interface[i]);

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

void OuterfaceDataTest::initializeDataOnPatch(
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

            tbox::Pointer< hier::PatchData<NDIM> > data =
               patch.getPatchData(variables[i], getDataContext());

	    TBOX_ASSERT( !data.isNull() );
            
	    tbox::Pointer<pdat::OuterfaceData<NDIM,double> > oface_data = data;
	    tbox::Pointer<pdat::FaceData<NDIM,double> > face_data = data;

            hier::Box<NDIM> dbox = data->getBox();

	    if ( !face_data.isNull() ) {
	       setLinearData(face_data, dbox, patch);
	    }
	    if ( !oface_data.isNull() ) {
	       setLinearData(oface_data, dbox, patch);
	    }
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < variables.getSize(); i++) {

         tbox::Pointer<hier::PatchData<NDIM> > data =
            patch.getPatchData(variables[i], getDataContext());
         
	 TBOX_ASSERT( !data.isNull() );
         
         tbox::Pointer<pdat::OuterfaceData<NDIM,double> > oface_data = data;
         tbox::Pointer<pdat::FaceData<NDIM,double> > face_data = data;
         
         
         hier::Box<NDIM> dbox = data->getGhostBox();
         
         if ( !face_data.isNull() ) {
            setLinearData(face_data, dbox, patch);
         }
         if ( !oface_data.isNull() ) {
            setLinearData(oface_data, dbox, patch);
         }
         
      } 
        
   }

}

void OuterfaceDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::OuterfaceData<NDIM,double> >& data,
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
      const pdat::FaceIndex<NDIM> loweri(interior.lower(), axis, 0);
      for (pdat::FaceIterator<NDIM> fi(interior, axis); fi; fi++) {

         /*
          * Compute spatial location of face and
          * set data to linear profile.
          */
         
         if (axis == 0) {
            x = lowerx[0] + dx[0]*(fi()(0) - loweri(0));
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(1) - loweri(1) + 0.5);
#endif
#if (NDIM > 2)            
            z = lowerx[2] + dx[2]*(fi()(2) - loweri(2) + 0.5);
#endif
         }
         else if (axis == 1) {
            x = lowerx[0] + dx[0]*(fi()(NDIM-1) - loweri(NDIM-1) + 0.5);
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(0) - loweri(0));
#endif
#if (NDIM > 2)
            z = lowerx[2] + dx[2]*(fi()(1) - loweri(1) + 0.5);
#endif
         }
         else if (axis == 2) {
            x = lowerx[0] + dx[0]*(fi()(1) - loweri(1) + 0.5);
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(2) - loweri(2) + 0.5);
#endif
#if (NDIM > 2)
            z = lowerx[2] + dx[2]*(fi()(0) - loweri(0));
#endif
         }

         double value;
         for (int d = 0; d < depth; d++) {
            value = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
            if (!(tbox::MathUtilities<double>::equalEps((*data)(fi(),d), value))) {
               tbox::perr << "FAILED: -- patch interior not properly filled" << endl;
            }
         }
      }
   }
}

void OuterfaceDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) time;
}

void OuterfaceDataTest::setLinearData(
   tbox::Pointer< pdat::FaceData<NDIM,double> > data,
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
      const pdat::FaceIndex<NDIM> loweri(patch.getBox().lower(), axis, 0);
      for (pdat::FaceIterator<NDIM> fi(sbox, axis); fi; fi++) {

         /*
          * Compute spatial location of cell center and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0]*(fi()(0) - loweri(0));
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(1) - loweri(1) + 0.5);
#endif
#if (NDIM > 2)            
            z = lowerx[2] + dx[2]*(fi()(2) - loweri(2) + 0.5);
#endif
         }
         else if (axis == 1) {
            x = lowerx[0] + dx[0]*(fi()(NDIM-1) - loweri(NDIM-1) + 0.5);
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(0) - loweri(0));
#endif
#if (NDIM > 2)
            z = lowerx[2] + dx[2]*(fi()(1) - loweri(1) + 0.5);
#endif
         }
         else if (axis == 2) {
            x = lowerx[0] + dx[0]*(fi()(1) - loweri(1) + 0.5);
#if (NDIM > 1)
            y = lowerx[1] + dx[1]*(fi()(2) - loweri(2) + 0.5);
#endif
#if (NDIM > 2)
            z = lowerx[2] + dx[2]*(fi()(0) - loweri(0));
#endif
         }


         for (int d = 0; d < depth; d++) {
            (*data)(fi(),d) = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
         }

      }
   }

}

void OuterfaceDataTest::setLinearData(
   tbox::Pointer< pdat::OuterfaceData<NDIM,double> > data,
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
      const hier::Box<NDIM> facebox =
         pdat::FaceGeometry<NDIM>::toFaceBox(box, axis);
      for (int f = 0; f < 2; f++) {
         const hier::Box<NDIM> databox = data->getArrayData(axis, f).getBox();

         const pdat::FaceIndex<NDIM> loweri(patch.getBox().lower(), axis, 0);
	 for (hier::Box<NDIM>::Iterator bi(databox); bi; bi++) {
            
            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */
            int axis2 =0;
            
            if (axis == 0) {
               x = lowerx[0] + dx[0]*(bi()(0) - loweri(0));
#if (NDIM > 1)
               y = lowerx[1] + dx[1]*(bi()(1) - loweri(1) + 0.5);
#endif
#if (NDIM > 2)            
               z = lowerx[2] + dx[2]*(bi()(2) - loweri(2) + 0.5);
#endif
            }
            else if (axis == 1) {
               axis2 = NDIM - 1;
               x = lowerx[0] + dx[0]*(bi()(NDIM-1) - loweri(NDIM-1) + 0.5);
#if (NDIM > 1)
               y = lowerx[1] + dx[1]*(bi()(0) - loweri(0));
#endif
#if (NDIM > 2)
               z = lowerx[2] + dx[2]*(bi()(1) - loweri(1) + 0.5);
#endif
            }
            else if (axis == 2) {
               axis2 = 1;
               x = lowerx[0] + dx[0]*(bi()(1) - loweri(1) + 0.5);
#if (NDIM > 1)
               y = lowerx[1] + dx[1]*(bi()(2) - loweri(2) + 0.5);
#endif
#if (NDIM > 2)
               z = lowerx[2] + dx[2]*(bi()(0) - loweri(0));
#endif
            }
            double value =  d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;

            for (int d = 0 ; d < depth; d++) {
                data->getArrayData(axis,f)(bi(),d) = value;
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

bool OuterfaceDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OuterfaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl; 

      hier::IntVector<NDIM> tgcw(0);
      for (int i = 0; i < d_variables_dst.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
                                     getGhostCellWidth());
      }
      hier::Box<NDIM> pbox = patch.getBox();

      tbox::Pointer< pdat::FaceData<NDIM,double> > solution =
         new pdat::FaceData<NDIM,double>(pbox, 1, tgcw);

      hier::Box<NDIM> tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox, patch);//, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables_dst.getSize(); i++) {

         tbox::Pointer< pdat::FaceData<NDIM,double> > face_data =
            patch.getPatchData(d_variables_dst[i], getDataContext());
         int depth = face_data->getDepth();
         hier::Box<NDIM> dbox = face_data->getGhostBox();

         for (int id = 0; id < NDIM; id++) {
            for (pdat::FaceIterator<NDIM> fi(dbox, id); fi; fi++) {
               double correct = (*solution)(fi());
               for (int d = 0; d < depth; d++) {
                  double result = (*face_data)(fi(),d);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     tbox::perr << "Test FAILED: ...." 
                          << " : face_data index = " << fi() << endl;
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
      if (!test_failed) {   
         tbox::plog << "Outerface test Successful!" << endl;
      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting OuterfaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);
}

}
