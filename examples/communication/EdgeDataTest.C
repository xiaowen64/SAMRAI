//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/examples/communication/EdgeDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: AMR communication tests for edge-centered patch data
//

#include "EdgeDataTest.h"


#include "ArrayData.h"
#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CommTester.h"
#include "EdgeGeometry.h"
#include "EdgeIndex.h"
#include "EdgeIterator.h"
#include "EdgeVariable.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"

namespace SAMRAI {

static double GETVALUE(const double a, const double b) 
{
   return(a*((double)NDIM) + b);
}

EdgeDataTest::EdgeDataTest(
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

   readTestInput(main_input_db->getDatabase("EdgePatchDataTest"));

}

EdgeDataTest::~EdgeDataTest()
{
}

void EdgeDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

}

void EdgeDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] = 
         new pdat::EdgeVariable<NDIM,double>(d_variable_src_name[i],
                                  d_variable_depth[i],
                                  d_use_fine_value_at_interface[i]);
        //                          d_test_direction[i]);

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

void EdgeDataTest::setConstantData(
   tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   double ndimfact,
   double axfact) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   if (!box.empty()) {

      for (int axis = 0; axis < NDIM; axis++) {
         pdat::ArrayData<NDIM,double>& array = data->getArrayData(axis);
 
         double value = GETVALUE(ndimfact, axis);
         hier::Box<NDIM> edge_box = pdat::EdgeGeometry<NDIM>::toEdgeBox(box, axis);
         array.fillAll(value, edge_box);
      }
   }

}


void EdgeDataTest::setConservativeData(
   tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch,
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
   TBOX_ASSERT(!hierarchy.isNull());
   TBOX_ASSERT( (level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()) );
#endif

   int i,j;
   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

   hier::BoxArray<NDIM> domain = level->getPhysicalDomain();
   int ncells = 0;
   for (i = 0; i < domain.getNumberOfBoxes(); i++) {
      ncells += domain(i).size();
   }

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;
  
   if (level_number == 0) {

      /*
       * Set edge value on level zero as follows:
       *
       *    u0(i,j,k) = (j + k)/ncells
       *    u1(i,j,k) = (i + k)/ncells
       *    u2(i,j,k) = (i + j)/ncells
       */

      for (int axis = 0; axis < NDIM; axis++) { 
         for (pdat::CellIterator<NDIM> ci(sbox); ci; ci++) {
            double value = 0.0;
            for (i = 0; i < NDIM; i++) {
               if (i == axis) {
                  value += (double)(ci()(i));
               }
            }
            value /= ncells;
#if (NDIM == 1)
            for (int edge = 0; edge < 1; edge++) {
#endif
#if (NDIM == 2)
            for (int edge = pdat::EdgeIndex<NDIM>::Lower;
                 edge <= pdat::EdgeIndex<NDIM>::Upper; edge++) {
#endif
#if (NDIM == 3)
            for (int edge = pdat::EdgeIndex<NDIM>::LowerLeft;
                 edge <= pdat::EdgeIndex<NDIM>::UpperRight; edge++) {
#endif
               pdat::EdgeIndex<NDIM> si(ci(), axis, edge);

               for (int d = 0; d < depth; d++) {
                  (*data)(si,d) = value;
               }

            }
         }
      }
      
   } else {
      
      /*
       * Set edge value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is value on the underlying coarse edge, (ci,cj,ck) is
       * the underlying coarse edge index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell edge centers.
       */
      
      hier::IntVector<NDIM> ratio = level->getRatio();
      const int max_ratio = ratio.max();

      tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
      const double* dx = pgeom->getDx();

      int coarse_ncells = ncells;
      double* delta = new double[max_ratio*NDIM];
      for (j = 0; j < NDIM; j++) {
          coarse_ncells /= ratio(j);
          double coarse_dx = dx[j] * ratio(j);
          for (i = 0; i < ratio(j); i++) {
             delta[j*max_ratio+i] = (i + 0.5)*dx[j] - coarse_dx*0.5;
          }
       }

      for (int axis = 0; axis < NDIM; axis++) {
         hier::IntVector<NDIM> ci;
         hier::IntVector<NDIM> del;
         for (pdat::CellIterator<NDIM> fi(sbox); fi; fi++) {
            double value = 0.0;
            for (i = 0; i < NDIM; i++) {
               if (i == axis) {
                  int findx = fi()(i);
                  ci(i) = ( (findx < 0) ? (findx+1)/ratio(i)-1
                            : findx/ratio(i) );
                  del(i) = (int)delta[i*max_ratio + findx-ci(i)*ratio(i)];
                  value += (double)(ci(i));
               }
            }
            value /= coarse_ncells; 
            
            for (j = 0; j < NDIM; j++) {
               if (j == axis) {
                  value += ci(j) * del(j);
               }
            }            

#if (NDIM == 1)
            for (int edge = 0; edge < 1; edge++) {
#endif
#if (NDIM == 2)
            for (int edge = pdat::EdgeIndex<NDIM>::Lower;
                 edge <= pdat::EdgeIndex<NDIM>::Upper; edge++) {
#endif
#if (NDIM == 3)
            for (int edge = pdat::EdgeIndex<NDIM>::LowerLeft;
                 edge <= pdat::EdgeIndex<NDIM>::UpperRight; edge++) {
#endif
               pdat::EdgeIndex<NDIM> si(fi(), axis, edge);

               for (int d = 0; d < depth; d++) {
                  (*data)(si,d) = value;
               }

            }

         }
      }
      delete[] delta;
      
   }
   
}
   

                                 

void EdgeDataTest::initializeDataOnPatch(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
   int level_number,
   char src_or_dst)
{
   d_hierarchy = hierarchy;

   if (d_do_refine) {

      if (   (d_refine_option == "INTERIOR_FROM_SAME_LEVEL") 
          || ( (d_refine_option == "INTERIOR_FROM_COARSER_LEVEL") 
              && (level_number == 0) ) ) {

         for (int i = 0; i < d_variables.getSize(); i++) {

            tbox::Pointer< pdat::EdgeData<NDIM,double> > edge_data =
               patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = edge_data->getBox(); 

            setLinearData(edge_data, dbox,
                          patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::EdgeData<NDIM,double> > edge_data =
            patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = edge_data->getGhostBox();

            setConservativeData(edge_data, dbox,
                                patch, hierarchy, level_number);

      } 
        
   }

}

void EdgeDataTest::setConstantBoundaryData(
   tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
   const hier::BoundaryBox<NDIM>& bbox,
   double ndimfact,
   double axfact) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif
   int lid = bbox.getLocationIndex();
   const hier::IntVector<NDIM>& gcw = data->getGhostCellWidth();
   hier::Box<NDIM> fillbox = bbox.getBox();

   if ( bbox.getBoundaryType() == NDIM ) {

      for (int axis = 0; axis < NDIM; axis++) {
         if (axis == 0) {
            if (lid % 2) {
               fillbox.growLower(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            } else {
               fillbox.growUpper(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            } 
         }
         if (axis == 1) {
            if ( (lid == 0) || (lid == 1) || (lid == 4) || (lid == 5) ) {
               fillbox.growLower(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            } else {
               fillbox.growUpper(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            }
         }
         if (axis == 2) {
            if (lid < 4) {
               fillbox.growLower(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            } else {
               fillbox.growUpper(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
            }
         }
      }

   } else if ( bbox.getBoundaryType() == 1 ) {
      for (int axis = 0; axis < NDIM; axis++) {
         if (lid == 2*axis) {
            fillbox.growLower(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
         } else if (lid == 2*axis+1) {
            fillbox.growUpper(axis, tbox::MathUtilities<int>::Max(gcw(axis)-1,0));
         }
      }
   } else if ( bbox.getBoundaryType() == 2 ) {
      switch(lid) {
         case 0: {
            fillbox.growLower(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            fillbox.growLower(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 1: {
            fillbox.growUpper(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            fillbox.growLower(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 2: {
            fillbox.growLower(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            fillbox.growUpper(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 3: {
            fillbox.growUpper(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            fillbox.growUpper(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 4: {
            fillbox.growLower(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growLower(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 5: {
            fillbox.growLower(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growUpper(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 6: {
            fillbox.growUpper(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growLower(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 7: {
            fillbox.growUpper(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growUpper(2, tbox::MathUtilities<int>::Max(gcw(2)-1,0));
            break;
         }
         case 8: {
            fillbox.growLower(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growLower(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            break;
         }
         case 9: {
            fillbox.growUpper(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growLower(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            break;
         }
         case 10: {
            fillbox.growLower(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growUpper(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            break;
         }
         case 11: {
            fillbox.growUpper(0, tbox::MathUtilities<int>::Max(gcw(0)-1,0));
            fillbox.growUpper(1, tbox::MathUtilities<int>::Max(gcw(1)-1,0));
            break;
         }
      }
   }
   
   for (int id = 0; id < NDIM; id++) {
      double value = GETVALUE(ndimfact, axfact); 
      data->getArrayData(id).fillAll(value,
                                     pdat::EdgeGeometry<NDIM>::toEdgeBox(fillbox,id));
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

bool EdgeDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
#if (NDIM > 1)
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering EdgeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl; 

      hier::IntVector<NDIM> tgcw(0);
      for (int i = 0; i < d_variables.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
                                     getGhostCellWidth());
      }
      hier::Box<NDIM> pbox = patch.getBox();

      tbox::Pointer< pdat::EdgeData<NDIM,double> > solution =
         new pdat::EdgeData<NDIM,double>(pbox, 1, tgcw);

      hier::Box<NDIM> tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_coarsen) {
         setConservativeData(solution, tbox,
                             patch, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::EdgeData<NDIM,double> > edge_data =
            patch.getPatchData(d_variables[i], getDataContext());
         int depth = edge_data->getDepth();
         hier::Box<NDIM> dbox = edge_data->getGhostBox();

         if (d_do_refine) {

             setLinearData(solution, tbox, patch);
         }
   
         for (int id = 0; id < NDIM; id++) {
            for (pdat::EdgeIterator<NDIM> si(dbox, id); si; si++) {
               double correct = (*solution)(si());
               for (int d = 0; d < depth; d++) {
                  double result = (*edge_data)(si(),d);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     test_failed = true;
                     tbox::perr << "Test FAILED: ...." 
                          << " : edge_data index = " << si() << endl;
                     tbox::perr << "    hier::Variable<NDIM> = " << d_variable_src_name[i]
                          << " : depth index = " << d << endl;
                     tbox::perr << "    result = " << result
                          << " : correct = " << correct << endl;
                  }
               }
            }
         }

      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting EdgeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

#endif

   return (!test_failed);
}

void EdgeDataTest::setLinearData(
   tbox::Pointer< pdat::EdgeData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   for (int axis = 0; axis < NDIM; axis++) {
      const pdat::EdgeIndex<NDIM> loweri(patch.getBox().lower(), axis, 0);
      for (pdat::EdgeIterator<NDIM> ei(sbox, axis); ei; ei++) {

         /*
          * Compute spatial location of cell center and
          * set data to linear profile.
          */

         if (axis != 0) {
            x = lowerx[0] + dx[0]*(ei()(0) - loweri(0));
         } else {
            x = lowerx[0] + dx[0]*(ei()(0) - loweri(0) + 0.5);
         }
         y = z = 0.;
#if (NDIM > 1)
         if (axis != 1) {
            y = lowerx[1] + dx[1]*(ei()(1) - loweri(1));
         } else {
            y = lowerx[1] + dx[1]*(ei()(1) - loweri(1) + 0.5);
         }
#endif
#if (NDIM > 2)
         if (axis != 2) {
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

void EdgeDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::EdgeData<NDIM,double> >& data,
   const hier::Box<NDIM>& interior,
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   for (int axis = 0; axis < NDIM; axis++) {
      const pdat::EdgeIndex<NDIM> loweri(interior.lower(), axis, 0);
      for (pdat::EdgeIterator<NDIM> ei(interior, axis); ei; ei++) {

         /*
          * Compute spatial location of cell center and
          * set data to linear profile.
          */

         if (axis != 0) {
            x = lowerx[0] + dx[0]*(ei()(0) - loweri(0));
         } else {
            x = lowerx[0] + dx[0]*(ei()(0) - loweri(0) + 0.5);
         }
         y = z = 0.;
#if (NDIM > 1)
         if (axis != 1) {
            y = lowerx[1] + dx[1]*(ei()(1) - loweri(1));
         } else {
            y = lowerx[1] + dx[1]*(ei()(1) - loweri(1) + 0.5);
         }
#endif
#if (NDIM > 2)
         if (axis != 2) {
            z = lowerx[2] + dx[2]*(ei()(2) - loweri(2));
         } else {
            z = lowerx[2] + dx[2]*(ei()(2) - loweri(2) + 0.5);
         }
#endif

         double value;
         for (int d = 0; d < depth; d++) {
            value = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
            if (!(tbox::MathUtilities<double>::equalEps((*data)(ei(),d), value))) {
               tbox::perr << "FAILED:  -- patch interior not properly filled" << endl;
            }
         }

      }
   }

}


void EdgeDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw) const
{
   (void) time;

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();

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

      tbox::Pointer< pdat::EdgeData<NDIM,double> > edge_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> patch_interior = edge_data->getBox();
      checkPatchInteriorData(edge_data, patch_interior, pgeom);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {
         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(edge_data, fill_box, patch);
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {
         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(edge_data, fill_box, patch);
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

         setLinearData(edge_data, fill_box, patch);
      }
#endif

   }

}

}
