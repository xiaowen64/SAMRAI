//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/CellDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2147 $
// Modified:    $LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: AMR communication tests for cell-centered patch data
//

#include "CellDataTest.h"


#include "BoundaryBox.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "CommTester.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "Variable.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"

namespace SAMRAI {

CellDataTest::CellDataTest(
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

   readTestInput(main_input_db->getDatabase("CellPatchDataTest"));

}

CellDataTest::~CellDataTest()
{
}

void CellDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

void CellDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] = new pdat::CellVariable<NDIM,double>(d_variable_src_name[i],
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

void CellDataTest::setLinearData(
   tbox::Pointer< pdat::CellData<NDIM,double> > data, 
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
   const pdat::CellIndex<NDIM> loweri(patch.getBox().lower());
   const pdat::CellIndex<NDIM> upperi(patch.getBox().upper());
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;

   for (pdat::CellIterator<NDIM> ci(sbox); ci; ci++) {

      /*
       * Compute spatial location of cell center and
       * set data to linear profile.
       */

      x = lowerx[0] + dx[0]*(ci()(0) - loweri(0) + 0.5);
      y = z = 0.;
#if (NDIM > 1)
      y = lowerx[1] + dx[1]*(ci()(1) - loweri(1) + 0.5);
#endif
#if (NDIM > 2)
      z = lowerx[2] + dx[2]*(ci()(2) - loweri(2) + 0.5);
#endif

      for (int d = 0; d < depth; d++) {
         (*data)(ci(),d) = d_Dcoef + d_Acoef*x + d_Bcoef*y + d_Ccoef*z;
      }

   }

}

void CellDataTest::setConservativeData(
   tbox::Pointer< pdat::CellData<NDIM,double> > data,
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

   int i,j,d;
   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

   hier::BoxArray<NDIM> domain = level->getPhysicalDomain();
   int ncells = 0;
   for (i = 0; i < domain.getNumberOfBoxes(); i++) {
      ncells += domain[i].size();
   }

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box; 

   if (level_number == 0) {

      /*
       * Set cell value on level zero to u(i,j,k) = (i + j + k)/ncells.
       */

      for (pdat::CellIterator<NDIM> fi(sbox); fi; fi++) {
         double value = 0.0;
         for (i = 0; i < NDIM; i++) {
            value += (double)(fi()(i));
         }
         value /= ncells;
         for (d = 0; d < depth; d++) {
            (*data)(fi(),d) = value;
         } 
      }

   } else {

       /*
        * Set cell value on level > 0 to 
        *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
        * where u_c is underlying coarse value, (ci,cj,ck) is
        * the underlying coarse cell index, and (del_i,del_j,del_k)
        * is the vector between the coarse and fine cell centers.
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

       hier::IntVector<NDIM> ci;
       hier::IntVector<NDIM> del;
       for (pdat::CellIterator<NDIM> fi(sbox); fi; fi++) {

          double value = 0.0;
          for (i = 0; i < NDIM; i++) {
             int findx = fi()(i);
             ci(i) = ( (findx < 0) ? (findx+1)/ratio(i)-1
                                   : findx/ratio(i) );
             del(i) = (int)delta[i*max_ratio + findx-ci(i)*ratio(i)];
             value += (double)(ci(i));
          }
          value /= coarse_ncells; 
 
          for (j = 0; j < NDIM; j++) {
             value += ci(j) * del(j);   
          }

          for (d = 0; d < depth; d++) {
             (*data)(fi(),d) = value;
          }

       }
       delete[] delta;

   }

}
                                 

void CellDataTest::initializeDataOnPatch(
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

            tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
               patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = cell_data->getBox(); 

            setLinearData(cell_data, dbox, patch);

         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
            patch.getPatchData(d_variables[i], getDataContext());

         hier::Box<NDIM> dbox = cell_data->getGhostBox();

         setConservativeData(cell_data, dbox,
                             patch, hierarchy, level_number);

      } 
        
   }

}

void CellDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::CellData<NDIM,double> >& data,
   const hier::Box<NDIM>& interior,
   const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> >& pgeom) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   const pdat::CellIndex<NDIM> loweri(interior.lower());
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   for (pdat::CellIterator<NDIM> ci(interior); ci; ci++) {

      /*
       * Compute spatial location of cell center and
       * compare data to linear profile.
       */

      x = lowerx[0] + dx[0]*(ci()(0) - loweri(0) + 0.5);
      y = z = 0.;
#if (NDIM > 1)
      y = lowerx[1] + dx[1]*(ci()(1) - loweri(1) + 0.5);
#endif
#if (NDIM > 2)
      z = lowerx[2] + dx[2]*(ci()(2) - loweri(2) + 0.5);
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
   

void CellDataTest::setPhysicalBoundaryConditions(
   hier::Patch<NDIM>& patch,
   const double time,
   const hier::IntVector<NDIM>& gcw_to_fill) const
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

      tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> patch_interior = cell_data->getBox();
      checkPatchInteriorData(cell_data, patch_interior, pgeom);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
                                                  patch.getBox(),
                                                  gcw_to_fill);
                                                      
         setLinearData(cell_data, fill_box, patch);
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                                                  patch.getBox(),
                                                  gcw_to_fill);

         setLinearData(cell_data, fill_box, patch);
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

         setLinearData(cell_data, fill_box, patch);
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
bool CellDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;

   bool test_failed = false;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering CellDataTest::verifyResults..." << endl;
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

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setConservativeData(solution, tbox,
                             patch, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::CellData<NDIM,double> > cell_data =
            patch.getPatchData(d_variables[i], getDataContext());
         int depth = cell_data->getDepth();
         hier::Box<NDIM> dbox = cell_data->getGhostBox();

         for (pdat::CellIterator<NDIM> ci(dbox); ci; ci++) {
            double correct = (*solution)(ci());
            for (int d = 0; d < depth; d++) {
               double result = (*cell_data)(ci(),d);
               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...." 
                       << " : cell index = " << ci() << endl;
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
         tbox::plog << "CellDataTest Successful!" << endl;
      }

      solution.setNull();   // just to be anal...

      tbox::plog << "\nExiting CellDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);

}

}
