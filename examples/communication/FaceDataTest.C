//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/communication/FaceDataTest.C $
// Package:     SAMRAI tests
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2147 $
// Modified:    $LastChangedDate: 2008-04-23 16:48:12 -0700 (Wed, 23 Apr 2008) $
// Description: AMR communication tests for face-centered patch data
//

#include "FaceDataTest.h"


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
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"

namespace SAMRAI {

FaceDataTest::FaceDataTest(
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

   readTestInput(main_input_db->getDatabase("FacePatchDataTest"));

}

FaceDataTest::~FaceDataTest()
{
}

void FaceDataTest::readTestInput(tbox::Pointer<tbox::Database> db)
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

void FaceDataTest::registerVariables(CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester*)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i] = 
         new pdat::FaceVariable<NDIM,double>(d_variable_src_name[i],
                                  d_variable_depth[i],
                                  d_use_fine_value_at_interface[i]);

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


void FaceDataTest::setConservativeData(
   tbox::Pointer< pdat::FaceData<NDIM,double> > data,
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
      ncells += domain[i].size();
   }

   const int depth = data->getDepth();

   const hier::Box<NDIM> sbox = data->getGhostBox() * box;
  
   if (level_number == 0) {

      /*
       * Set face value on level zero as follows:
       *
       *    u0(i,j,k) = (j + k)/ncells
       *    u1(i,j,k) = (i + k)/ncells
       *    u2(i,j,k) = (i + j)/ncells
       */

      for (int axis = 0; axis < NDIM; axis++) { 
         for (pdat::CellIterator<NDIM> ci(sbox); ci; ci++) {
            double value = 0.0;
            for (i = 0; i < NDIM; i++) {
               if (i != axis) {
                  value += (double)(ci()(i));
               }
            }
            value /= ncells;
            for (int face = pdat::FaceIndex<NDIM>::Lower;
                 face <= pdat::FaceIndex<NDIM>::Upper; face++) {
               pdat::FaceIndex<NDIM> si(ci(), axis, face);
               for (int d = 0; d < depth; d++) {
                  (*data)(si,d) = value;
               }
            }
         }
      }

   } else {

      /*
       * Set face value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is value on the underlying coarse face, (ci,cj,ck) is
       * the underlying coarse face index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell face centers.
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
                if (i != axis) {
                   int findx = fi()(i);
                   ci(i) = ( (findx < 0) ? (findx+1)/ratio(i)-1
                             : findx/ratio(i) );
                   del(i) = (int)delta[i*max_ratio + findx-ci(i)*ratio(i)];
                   value += (double)(ci(i));
                }
             }
             value /= coarse_ncells; 
            
             for (j = 0; j < NDIM; j++) {
                if (j != axis) {
                   value += ci(j) * del(j);
                }
             } 

             for (int face = pdat::FaceIndex<NDIM>::Lower;
                  face <= pdat::FaceIndex<NDIM>::Upper; face++) {
                pdat::FaceIndex<NDIM> si(fi(), axis, face);  
                for (int d = 0; d < depth; d++) {
                   (*data)(si,d) = value;
                }
             }
          }
       }
       delete[] delta;

   }

}
                                 

void FaceDataTest::initializeDataOnPatch(
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

            tbox::Pointer< pdat::FaceData<NDIM,double> > face_data =
               patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = face_data->getBox();

            setLinearData(face_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::FaceData<NDIM,double> > face_data =
            patch.getPatchData(d_variables[i], getDataContext());

            hier::Box<NDIM> dbox = face_data->getGhostBox();

            setConservativeData(face_data, dbox,
                                patch, hierarchy, level_number);

      } 
        
   }

}

void FaceDataTest::checkPatchInteriorData(
   const tbox::Pointer< pdat::FaceData<NDIM,double> >& data,
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

void FaceDataTest::setPhysicalBoundaryConditions(
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

      tbox::Pointer< pdat::FaceData<NDIM,double> > face_data =
         patch.getPatchData(d_variables[i], getDataContext());

      hier::Box<NDIM> patch_interior = face_data->getBox();
      checkPatchInteriorData(face_data, patch_interior, pgeom);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
                                                  patch.getBox(),
                                                  gcw);


         setLinearData(face_data, fill_box, patch);
      }

#if (NDIM > 1)
      /*
       * Set edge boundary data.
       */
      for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(face_data, fill_box, patch);
      }
#endif

#if (NDIM > 2)
      /*
       * Set face boundary data.
       */
      for (int fi = 0; fi < num_face_bdry_boxes; fi++) {

         hier::Box<NDIM> fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                                                  patch.getBox(),
                                                  gcw);

         setLinearData(face_data, fill_box, patch);
      }
#endif

   }

}

void FaceDataTest::setLinearData(
   tbox::Pointer< pdat::FaceData<NDIM,double> > data,
   const hier::Box<NDIM>& box,
   hier::Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!data.isNull());
#endif

   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
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

/*
*************************************************************************
*                                                                       *
* Verify results of communication operations.  This test must be        *
* consistent with data initialization and boundary operations above.    *
*                                                                       *
*************************************************************************
*/

bool FaceDataTest::verifyResults(
   hier::Patch<NDIM>& patch, 
   const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy, 
   int level_number)
{
   (void) hierarchy;
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering FaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl; 

      hier::IntVector<NDIM> tgcw(0);
      for (int i = 0; i < d_variables.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
                                     getGhostCellWidth());
      }
      hier::Box<NDIM> pbox = patch.getBox();

      tbox::Pointer< pdat::FaceData<NDIM,double> > solution =
         new pdat::FaceData<NDIM,double>(pbox, 1, tgcw);

      hier::Box<NDIM> tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_coarsen) {
         setConservativeData(solution, tbox,
                             patch, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         tbox::Pointer< pdat::FaceData<NDIM,double> > face_data =
            patch.getPatchData(d_variables[i], getDataContext());
         int depth = face_data->getDepth();
         hier::Box<NDIM> dbox = face_data->getGhostBox();

         if (d_do_refine) {
            setLinearData(solution, tbox, patch);
         }
   
         for (int id = 0; id < NDIM; id++) {
            for (pdat::FaceIterator<NDIM> si(dbox, id); si; si++) {
               double correct = (*solution)(si());
               for (int d = 0; d < depth; d++) {
                  double result = (*face_data)(si(),d);
                  if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                     test_failed = true;
                     tbox::perr << "Test FAILED: ...." 
                          << " : face_data index = " << si() << endl;
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

      tbox::plog << "\nExiting FaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl; 

   }

   return (!test_failed);
}

}
