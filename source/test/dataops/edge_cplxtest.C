//
// File:        edge_cplxtest.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 590 $
// Modified:    $Date: 2005-08-29 10:37:19 -0700 (Mon, 29 Aug 2005) $
// Description: Main program to test edge-centered complex patch data ops
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iomanip>
using namespace std;

#include "tbox/MPI.h"
#include "tbox/PIO.h"
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif

#include "tbox/SAMRAIManager.h"
#include "tbox/Pointer.h"

#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "EdgeData.h"
#include "HierarchyDataOpsComplex.h"
#include "HierarchyEdgeDataOpsComplex.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "EdgeIndex.h"
#include "EdgeIterator.h"
#include "EdgeVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDescriptor.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "tbox/Complex.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"
#include "VariableContext.h"

using namespace SAMRAI;

/* Helper function declarations */
bool complexDataSameAsValue(int desc_id, dcomplex value,
   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy);
bool doubleDataSameAsValue(int desc_id, double value,
   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy);

#define NVARS 4

int main( int argc, char *argv[] ) {
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
   tbox::PIO::logAllNodes("edge_cplxtest.log");

   int ln, iv;

   // Make a dummy hierarchy domain
   double lo[2] = {0.0, 0.0};
   double hi[2] = {1.0, 0.5};

   hier::Box<2> coarse0(hier::Index<2>(0,0), hier::Index<2>(9,2));
   hier::Box<2> coarse1(hier::Index<2>(0,3), hier::Index<2>(9,4));
   hier::Box<2> fine0(hier::Index<2>(4,4), hier::Index<2>(7,7));
   hier::Box<2> fine1(hier::Index<2>(8,4), hier::Index<2>(13,7));
   hier::IntVector<2> ratio(2);

   hier::BoxArray<2> coarse_domain(2);
   hier::BoxArray<2> fine_domain(2);
   coarse_domain(0) = coarse0;
   coarse_domain(1) = coarse1;
   fine_domain(0) = fine0;
   fine_domain(1) = fine1;

   tbox::Pointer<geom::CartesianGridGeometry<2> > geometry =
   new geom::CartesianGridGeometry<2>("CartesianGeometry", lo, hi, coarse_domain);

   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy = 
      new hier::PatchHierarchy<2>("PatchHierarchy", geometry);

   // Note: For these simple tests we allow at most 2 processors.
   const int nproc = tbox::MPI::getNodes();
   assert(nproc < 3);

   const int n_coarse_boxes = coarse_domain.getNumberOfBoxes();
   const int n_fine_boxes = fine_domain.getNumberOfBoxes();
   hier::ProcessorMapping mapping0(n_coarse_boxes);
   hier::ProcessorMapping mapping1(n_fine_boxes);

   int ib;
   for (ib = 0; ib < n_coarse_boxes; ib++) {
      if (nproc > 1) {
         mapping0.setProcessorAssignment(ib, ib);
      } else {
         mapping0.setProcessorAssignment(ib, 0);
      }
   }

   for (ib = 0; ib < n_fine_boxes; ib++) {
      if (nproc > 1) {
         mapping1.setProcessorAssignment(ib, ib);
      } else {
         mapping1.setProcessorAssignment(ib, 0);
      }
   }

   hierarchy->makeNewPatchLevel(0, hier::IntVector<2>(1), coarse_domain, mapping0);
   hierarchy->makeNewPatchLevel(1, ratio, fine_domain, mapping1);

   // Create instance of hier::Variable<NDIM> database
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer< hier::VariableContext > dummy = variable_db->getContext("dummy");
   const hier::IntVector<2> no_ghosts(0);

   // Make some dummy variables and data on the hierarchy
   tbox::Pointer< pdat::EdgeVariable<2,dcomplex> > fvar[NVARS];
   int svindx[NVARS];
   fvar[0] = new pdat::EdgeVariable<2,dcomplex>("fvar0", 1);
   svindx[0] = variable_db->registerVariableAndContext(
               fvar[0], dummy, no_ghosts);
   fvar[1] = new pdat::EdgeVariable<2,dcomplex>("fvar1", 1);
   svindx[1] = variable_db->registerVariableAndContext(
               fvar[1], dummy, no_ghosts);
   fvar[2] = new pdat::EdgeVariable<2,dcomplex>("fvar2", 1);
   svindx[2] = variable_db->registerVariableAndContext(
               fvar[2], dummy, no_ghosts);
   fvar[3] = new pdat::EdgeVariable<2,dcomplex>("fvar3", 1);
   svindx[3] = variable_db->registerVariableAndContext(
               fvar[3], dummy, no_ghosts);

   tbox::Pointer< pdat::EdgeVariable<2,double> >
   swgt = new pdat::EdgeVariable<2,double>("swgt", 1);
   int swgt_id = variable_db->registerVariableAndContext(
                 swgt, dummy, no_ghosts);

   // allocate data on hierarchy
   for (ln = 0; ln < 2; ln++) {
      hierarchy->getPatchLevel(ln)->allocatePatchData(swgt_id);
      for (iv = 0; iv < NVARS; iv++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(svindx[iv]);
      }
   }

   tbox::Pointer<math::HierarchyDataOpsComplex<2> > edge_ops = 
   new math::HierarchyEdgeDataOpsComplex<2>(hierarchy, 0, 1);
   assert(!edge_ops.isNull());

   tbox::Pointer<math::HierarchyDataOpsReal<2,double> > swgt_ops = 
   new math::HierarchyEdgeDataOpsReal<2,double>(hierarchy, 0 , 1);

   tbox::Pointer<hier::Patch<2> > patch;
   tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom;

   // Initialize control volume data for edge-centered components
   hier::Box<2> coarse_fine = fine0+fine1; 
   coarse_fine.coarsen(ratio);
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
         tbox::Pointer< pdat::EdgeData<2,double> > data;
         patch = level->getPatch(ip());
         pgeom = patch->getPatchGeometry();
         const double* dx = pgeom->getDx();
         const double edge_vol = dx[0]*dx[1];
         data = patch->getPatchData(swgt_id);
         data->fillAll(edge_vol);
         pdat::EdgeIndex<2> fi;
         int plo0 = patch->getBox().lower(0);
         int phi0 = patch->getBox().upper(0);
         int plo1 = patch->getBox().lower(1);
         int phi1 = patch->getBox().upper(1);
         int ic;

         if (ln == 0) {
            data->fillAll(0.0, (coarse_fine * patch->getBox()) );

            if (patch->getPatchNumber() == 0) {
               //bottom edge boundaries
               for (ic = plo0; ic <= phi0; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,plo1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 0.5;
               }
               //left and right edge boundaries
               for (ic = plo1; ic <= phi1; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(plo0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 0.5;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(phi0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 0.5;
               }
            } else {
               //top and bottom edge boundaries
               for (ic = plo0; ic <= phi0; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,plo1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) = 0.0;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,phi1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 0.5;
               }
               //left and right edge boundaries
               for (ic = plo1; ic <= phi1; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(plo0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 0.5;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(phi0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 0.5;
               }
            }
         } else {
            if (patch->getPatchNumber() == 0) {
               // top and bottom coarse-fine edge boundaries
               for (ic = plo0; ic <= phi0; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,plo1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 1.5;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,phi1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 1.5;
               }
               //left coarse-fine edge boundaries
               for (ic = plo1; ic <= phi1; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(plo0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 1.5;
               }
            } else {
               // top and bottom coarse-fine edge boundaries
               for (ic = plo0; ic <= phi0; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,plo1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) *= 1.5;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(ic,phi1), pdat::EdgeIndex<2>::X, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 1.5;
               }
               //left and right coarse-fine edge boundaries
               for (ic = plo1; ic <= phi1; ic++) {
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(plo0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Lower);
                  (*data)(fi) = 0.0;
                  fi = pdat::EdgeIndex<2>(hier::Index<2>(phi0,ic), pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Upper);
                  (*data)(fi) *= 1.5;
               }
            }
         }
      }
   }

   // Test #1: Print out control volume data and compute its integral

   // Test #1a: Check control volume data set properly
   // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
   // 0.0025 on fine level
/*   bool vol_test_passed = true;
   for (ln = 0; ln < 2; ln++) {
      for (hier::PatchLevel<2>::Iterator ip(hierarchy->getPatchLevel(ln)); ip; ip++) {
         patch = hierarchy->getPatchLevel(ln)->getPatch(ip());
         tbox::Pointer< pdat::CellData<2,double> > svdata = patch->getPatchData(cwgt_id);

         for (pdat::CellIterator<2> c(svdata->getBox());c && vol_test_passed;c++) {
            pdat::CellIndex<2> cell_index = c();

            if (ln == 0) {
               if ((coarse_fine * patch->getBox()).contains(cell_index)) {
                  if ( !tbox::Utilities::deq((*svdata)(cell_index),0.0) ) {
                     vol_test_passed = false;
                  }
               } else {
                  if ( !tbox::Utilities::deq((*svdata)(cell_index),0.01) ) {
                     vol_test_passed = false;
                  }
               }
            }

            if (ln == 1) {
               if ( !tbox::Utilities::deq((*svdata)(cell_index),0.0025) ) {
                  vol_test_passed = false;
               }
            }
         }
      }
   }
   if (!vol_test_passed) {
      tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << endl;
      cwgt_ops->printData(cwgt_id, tbox::plog);
   }
*/
/*   tbox::plog << "edge control volume data" << endl;
   swgt_ops->printData(swgt_id, tbox::plog);
*/

   // Test #1b: math::HierarchyEdgeDataOpsComplex<2>::sumControlVolumes()
   // Expected: norm = 1.0
   double norm = 
          edge_ops->sumControlVolumes(svindx[0], swgt_id);
   if ( !tbox::Utilities::deq(norm,1.0)) {
      tbox::perr << "FAILED: - Test #1b: math::HierarchyEdgeDataOpsComplex<2>::sumControlVolumes()\n"
           << "Expected value = 1.0 , Computed value = "
           << norm << endl;
   }

   // Test #2: math::HierarchyEdgeDataOpsComplex<2>::numberOfEntries()
   // Expected: num_data_points = 209
   int num_data_points = edge_ops->numberOfEntries(svindx[0]);
   if ( num_data_points != 209 ) {
      tbox::perr << "FAILED: - Test #2: math::HierarchyEdgeDataOpsComplex<2>::numberOfEntries()\n"
           << "Expected value = 209 , Computed value = "
           << num_data_points << endl;
   }

   // Test #3a: math::HierarchyEdgeDataOpsComplex<2>::setToScalar()
   // Expected: v0 = (2.0,1.5)
   dcomplex val0 = dcomplex(2.0,1.5);
   edge_ops->setToScalar(svindx[0], val0);
   if (!complexDataSameAsValue(svindx[0],val0,hierarchy)) {
      tbox::perr << "FAILED: - Test #3a: math::HierarchyEdgeDataOpsComplex<2>::setToScalar()\n"
           << "Expected: v0 = " << val0 << endl;
      edge_ops->printData(svindx[0], tbox::plog);
   }

   // Test #3b: math::HierarchyEdgeDataOpsComplex<2>::setToScalar()
   // Expected: v1 = (4.0, 3.0)
   dcomplex val1(4.0,3.0);
   edge_ops->setToScalar(svindx[1], dcomplex(4.0, 3.0));
   if (!complexDataSameAsValue(svindx[1],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #3b: math::HierarchyEdgeDataOpsComplex<2>::setToScalar()\n"
           << "Expected: v1 = " << val1 << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #4: math::HierarchyEdgeDataOpsComplex<2>::copyData()
   // Expected: v2 = v1 = (4.0, 3.0)
   edge_ops->copyData(svindx[2], svindx[1]);
   if (!complexDataSameAsValue(svindx[2],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #4: math::HierarchyEdgeDataOpsComplex<2>::copyData()\n"
           << "Expected: v2 = " << val1 << endl;
      edge_ops->printData(svindx[2], tbox::plog);
   }

   // Test #5: math::HierarchyEdgeDataOpsComplex<2>::swapData()
   // Expected:  v0 = (4.0, 3.0), v1 = (2.0,1.5)
   edge_ops->swapData(svindx[0], svindx[1]);
   if (!complexDataSameAsValue(svindx[0],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #5a: math::HierarchyEdgeDataOpsComplex<2>::swapData()\n"
           << "Expected: v0 = " << val1 << endl;
      edge_ops->printData(svindx[0], tbox::plog);
   }
   if (!complexDataSameAsValue(svindx[1],val0,hierarchy)) {
      tbox::perr << "FAILED: - Test #5b: math::HierarchyEdgeDataOpsComplex<2>::swapData()\n"
           << "Expected: v1 = " << val0 << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #6: math::HierarchyEdgeDataOpsComplex<2>::scale()
   // Expected:  v2 = 0.25 * v2 = (1.0,0.75)
   edge_ops->scale(svindx[2], 0.25, svindx[2]);
   dcomplex val_scale(1.0,0.75);
   if (!complexDataSameAsValue(svindx[2],val_scale,hierarchy)) {
      tbox::perr << "FAILED: - Test #6: math::HierarchyEdgeDataOpsComplex<2>::scale()\n"
           << "Expected: v2 = " << val_scale << endl;
      edge_ops->printData(svindx[2], tbox::plog);
   }

   // Test #7: math::HierarchyEdgeDataOpsComplex<2>::add()
   // Expected:  v3 = v0 + v1 = (6.0, 4.5) 
   edge_ops->add(svindx[3], svindx[0], svindx[1]);
   dcomplex val_add(6.0,4.5);
   if (!complexDataSameAsValue(svindx[3],val_add,hierarchy)) {
      tbox::perr << "FAILED: - Test #7: math::HierarchyEdgeDataOpsComplex<2>::add()\n"
           << "Expected: v3 = " << val_add << endl;
      edge_ops->printData(svindx[3], tbox::plog);
   }

   // Reset v0: v0 = (0.0,4.5)
   edge_ops->setToScalar(svindx[0], dcomplex(0.0,4.5));

   // Test #8: math::HierarchyEdgeDataOpsComplex<2>::subtract()
   // Expected:  v1 = v3 - v0 = (6.0,0.0)
   edge_ops->subtract(svindx[1], svindx[3], svindx[0]);  
   dcomplex val_sub(6.0,0.0);
   if (!complexDataSameAsValue(svindx[1],val_sub,hierarchy)) {
      tbox::perr << "FAILED: - Test #8: math::HierarchyEdgeDataOpsComplex<2>::subtract()\n"
           << "Expected: v1 = " << val_sub << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #9a: math::HierarchyEdgeDataOpsComplex<2>::addScalar()
   // Expected:  v1 = v1 + (0.0,-4.0) = (6.0,-4.0) 
   edge_ops->addScalar(svindx[1], svindx[1], dcomplex(0.0,-4.0));
   dcomplex val_addScalar(6.0,-4.0);
   if (!complexDataSameAsValue(svindx[1],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9a: math::HierarchyEdgeDataOpsComplex<2>::addScalar()\n"
           << "Expected: v1 = " << val_addScalar << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #9b: math::HierarchyEdgeDataOpsComplex<2>::addScalar()
   // Expected:  v2 = v2 + (0.0,0.25) = (1.0,1.0) 
   edge_ops->addScalar(svindx[2], svindx[2], dcomplex(0.0,0.25)); 
   val_addScalar = dcomplex(1.0,1.0);
   if (!complexDataSameAsValue(svindx[2],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9b: math::HierarchyEdgeDataOpsComplex<2>::addScalar()\n"
           << "Expected: v2 = " << val_addScalar << endl;
      edge_ops->printData(svindx[2], tbox::plog);
   }

   // Test #9c: math::HierarchyEdgeDataOpsComplex<2>::addScalar()
   // Expected:  v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
   edge_ops->addScalar(svindx[2], svindx[2], dcomplex(3.0,-4.0));
   val_addScalar = dcomplex(4.0,-3.0);
   if (!complexDataSameAsValue(svindx[2],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9c: math::HierarchyEdgeDataOpsComplex<2>::addScalar()\n"
           << "Expected: v2 = " << val_addScalar << endl;
      edge_ops->printData(svindx[2], tbox::plog);
   }

   // Rest v3:  v3 = (0.5, 0.0)
   edge_ops->setToScalar(svindx[3], dcomplex(0.5, 0.0));

   // Test #10: math::HierarchyEdgeDataOpsComplex<2>::multiply()
   // Expected:  v1 = v3 * v1 = (3.0,-2.0) 
   edge_ops->multiply(svindx[1], svindx[3], svindx[1]);
   dcomplex val_mult(3.0,-2.0); 
   if (!complexDataSameAsValue(svindx[1],val_mult,hierarchy)) {
      tbox::perr << "FAILED: - Test #10: math::HierarchyEdgeDataOpsComplex<2>::multiply()\n"
           << "Expected: v1 = " << val_mult << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #11: math::HierarchyEdgeDataOpsComplex<2>::divide()
   // Expected:  v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
   edge_ops->divide(svindx[0], svindx[2], svindx[1]);
   dcomplex val_div(1.3846153846154,-0.076923076923077);
   if (!complexDataSameAsValue(svindx[0],val_div,hierarchy)) {
      tbox::perr << "FAILED: - Test #11: math::HierarchyEdgeDataOpsComplex<2>::divide()\n"
           << "Expected: v0 = " << val_div << endl;
      edge_ops->printData(svindx[0], tbox::plog);
   }

   // Test #12: math::HierarchyEdgeDataOpsComplex<2>::reciprocal()
   // Expected:  v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
   edge_ops->reciprocal(svindx[1], svindx[1]); 
   dcomplex val_rec(0.23076923076923, 0.15384615384615);
   if (!complexDataSameAsValue(svindx[1],val_rec,hierarchy)) {
      tbox::perr << "FAILED: - Test #12: math::HierarchyEdgeDataOpsComplex<2>::reciprocal()\n"
           << "Expected: v1 = " << val_rec << endl;
      edge_ops->printData(svindx[1], tbox::plog);
   }

   // Test #13:  Place some bogus values on coarse level
   tbox::Pointer< pdat::EdgeData<2,dcomplex> > sdata;

   // set values

   tbox::Pointer<hier::PatchLevel<2> > level_zero = hierarchy->getPatchLevel(0);
   for (hier::PatchLevel<2>::Iterator ip(level_zero); ip; ip++) {
      patch = level_zero->getPatch(ip());
      sdata = patch->getPatchData(svindx[2]);
      hier::Index<2> index0(2,2);
      hier::Index<2> index1(5,3);
      if (patch->getBox().contains(index0)) {
         (*sdata)(pdat::EdgeIndex<2>(index0, pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Lower), 0) = dcomplex(100.0,-50.0);
      }
      if (patch->getBox().contains(index1)) {
         (*sdata)(pdat::EdgeIndex<2>(index1, pdat::EdgeIndex<2>::Y, pdat::EdgeIndex<2>::Upper), 0) = dcomplex(-1000.0,20.0);
      }
   }

   // check values
   bool bogus_value_test_passed = true;
   for (hier::PatchLevel<2>::Iterator ipp(level_zero); ipp; ipp++) {
      patch = level_zero->getPatch(ipp());
      sdata = patch->getPatchData(svindx[2]);
      pdat::EdgeIndex<2> index0(hier::Index<2>(2,2),pdat::EdgeIndex<2>::Y,pdat::EdgeIndex<2>::Lower);
      pdat::EdgeIndex<2> index1(hier::Index<2>(5,3),pdat::EdgeIndex<2>::Y,pdat::EdgeIndex<2>::Upper);

      // check X axis data
      for (pdat::EdgeIterator<2> c(sdata->getBox(),pdat::EdgeIndex<2>::X);c && bogus_value_test_passed;c++) {
         pdat::EdgeIndex<2> edge_index = c();

         if (!tbox::Utilities::ceq((*sdata)(edge_index),
            dcomplex(4.0,-3.0)) ) {
            bogus_value_test_passed = false;
         }
      }

      // check Y axis data
      for (pdat::EdgeIterator<2> cc(sdata->getBox(),pdat::EdgeIndex<2>::Y);cc && bogus_value_test_passed;cc++) {
         pdat::EdgeIndex<2> edge_index = cc();

         if (edge_index == index0) {
            if (!tbox::Utilities::ceq((*sdata)(edge_index),
                dcomplex(100.0,-50.0)) ) {
               bogus_value_test_passed = false;
            }
         } else {
            if (edge_index == index1) {
               if (!tbox::Utilities::ceq((*sdata)(edge_index),
                  dcomplex(-1000.0,20.0)) ) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (!tbox::Utilities::ceq((*sdata)(edge_index),
                  dcomplex(4.0,-3.0)) ) {
                  bogus_value_test_passed = false;
               }
            }
         }
      }
   }
   if (!bogus_value_test_passed) {
      tbox::perr << "FAILED: - Test #13:  Place some bogus values on coarse level" << endl;
      edge_ops->printData(svindx[2], tbox::plog);
   }

   // Test norms on patch data with svindx[2] on hierarchy with bogus values

   // Test #14: math::HierarchyEdgeDataOpsComplex<2>::L1Norm() - w/o control weight
   // Expected:  bogus_l1_norm = 2217.003379
   double bogus_l1_norm = edge_ops->L1Norm(svindx[2]);
   if ( !tbox::Utilities::deq(bogus_l1_norm,2217.003379) ) {
      tbox::perr << "FAILED: - Test #14: math::HierarchyEdgeDataOpsComplex<2>::L1Norm()"
           << " - w/o control weight\n"
           << "Expected value = 2217.003379, Computed value = "
           << setprecision(12) << bogus_l1_norm << endl;
   }

   // Test #15: math::HierarchyEdgeDataOpsComplex<2>::L1Norm() - w/control weight
   // Expected:  correct_l1_norm = 5.0
   double correct_l1_norm = edge_ops->L1Norm(svindx[2],swgt_id);
   if ( !tbox::Utilities::deq(correct_l1_norm,5.0) ) {
      tbox::perr << "FAILED: - Test #15: math::HierarchyEdgeDataOpsComplex<2>::L1Norm()"
           << " - w/control weight\n"
           << "Expected value = 5.0, Computed value = "
           << correct_l1_norm << endl;
   }

   // Test #16: math::HierarchyEdgeDataOpsComplex<2>::L2Norm()
   // Expected:  l2_norm = 5.0
   double l2_norm = edge_ops->L2Norm(svindx[2],swgt_id);
   if ( !tbox::Utilities::deq(l2_norm,5.0) ) {
      tbox::perr << "FAILED: - Test #16: math::HierarchyEdgeDataOpsComplex<2>::L2Norm()\n"
           << "Expected value = 5.0, Computed value = "
           << l2_norm << endl;
   }

   // Test #17: math::HierarchyEdgeDataOpsComplex<2>::maxNorm() - w/o control weight
   // Expected:  bogus_max_norm = 1000.19998
   double bogus_max_norm = edge_ops->maxNorm(svindx[2]);
   if ( !tbox::Utilities::deq(bogus_max_norm,1000.19998) ) {
      tbox::perr << "FAILED: - Test #17: math::HierarchyEdgeDataOpsComplex<2>::maxNorm()"
           << " - w/o control weight\n"
           << "Expected value = 1000.19998, Computed value = "
           << bogus_max_norm << endl;
   }

   // Test #18: math::HierarchyEdgeDataOpsComplex<2>::maxNorm() - w/control weight
   // Expected:  max_norm = 5.0
   double max_norm = edge_ops->maxNorm(svindx[2],swgt_id);
   if ( !tbox::Utilities::deq(max_norm,5.0) ) {
      tbox::perr << "FAILED: - Test #18: math::HierarchyEdgeDataOpsComplex<2>::maxNorm()"
           << " - w/control weight\n"
           << "Expected value = 5.0, Computed value = "
           << max_norm << endl;
   }

   // Reset data and test sums, axpy's
   edge_ops->setToScalar(svindx[0], dcomplex(1.0, -3.0)); 
   edge_ops->setToScalar(svindx[1], dcomplex(2.5, 3.0)); 
   edge_ops->setToScalar(svindx[2], dcomplex(7.0, 0.0)); 

   // Test #19: math::HierarchyEdgeDataOpsComplex<2>::linearSum()
   // Expected:  v3 = (2.0,5.0)
   edge_ops->linearSum(svindx[3], dcomplex(2.0,0.0), svindx[1], dcomplex(0.0, -1.0), svindx[0]);
   dcomplex val_linearSum(2.0,5.0);
   if (!complexDataSameAsValue(svindx[3],val_linearSum,hierarchy)) {
      tbox::perr << "FAILED: - Test #19: math::HierarchyEdgeDataOpsComplex<2>::linearSum()\n"
           << "Expected: v3 = " << val_linearSum << endl;
      edge_ops->printData(svindx[3], tbox::plog);
   }

   // Test #20: math::HierarchyEdgeDataOpsComplex<2>::axmy()
   // Expected:  v3 = (6.5,12.0)
   edge_ops->axmy(svindx[3], 3.0, svindx[1], svindx[0]);
   dcomplex val_axmy(6.5,12.0);
   if (!complexDataSameAsValue(svindx[3],val_axmy,hierarchy)) {
      tbox::perr << "FAILED: - Test #20: math::HierarchyEdgeDataOpsComplex<2>::axmy()\n"
           << "Expected: v3 = " << val_axmy << endl;
      edge_ops->printData(svindx[3], tbox::plog);
   }

   // Test #21a: math::HierarchyEdgeDataOpsComplex<2>::dot() - (ind2) * (ind1)
   // Expected:  cdot = (17.5,-21.0)
   dcomplex cdot = edge_ops->dot(svindx[2], svindx[1], swgt_id);
   dcomplex ans_2_dot_1(17.5,-21.0);
   if ( !tbox::Utilities::ceq(cdot,ans_2_dot_1) ) {
      tbox::perr << "FAILED: - Test #21a: math::HierarchyEdgeDataOpsComplex<2>::dot() - (ind2) * (ind1)\n"
           << "Expected value = (17.5,-21.0), Computed value = "
           << cdot << endl;
   }

   // Test #21b: math::HierarchyEdgeDataOpsComplex<2>::dot() - (ind1) * (ind2)
   // Expected:  cdot = (17.5,-1.0)
   dcomplex cdot2 = edge_ops->dot(svindx[1], svindx[2], swgt_id); 
   dcomplex ans_1_dot_2(17.5,21.0);
   if ( !tbox::Utilities::ceq(cdot2,ans_1_dot_2) ) {
      tbox::perr << "FAILED: - Test #21b: math::HierarchyEdgeDataOpsComplex<2>::dot() - (ind1) * (ind2)\n"
           << "Expected value = (17.5,21.0), Computed value = "
           << cdot2 << endl;
   }

   // Test #22: math::HierarchyEdgeDataOpsComplex<2>::abs()
   // Expected:  abs(v0) = 5.0
   edge_ops->setToScalar(svindx[0], dcomplex(4.0, -3.0));
   edge_ops->abs(swgt_id, svindx[0]);
   if (!doubleDataSameAsValue(swgt_id,5.0,hierarchy)) {
      tbox::perr << "FAILED: - Test #22: math::HierarchyEdgeDataOpsComplex<2>::abs()\n"
           << "Expected: abs(v0) = 5.0" << endl;
      swgt_ops->printData(swgt_id, tbox::plog);
   }

   // deallocate data on hierarchy
   for (ln = 0; ln < 2; ln++) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(swgt_id);
      for (iv = 0; iv < NVARS; iv++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(svindx[iv]);
      }
   }

   for (iv = 0; iv < NVARS; iv++) {
      fvar[iv].setNull();
   }
   swgt.setNull();

   geometry.setNull();
   hierarchy.setNull();
   edge_ops.setNull(); 
   swgt_ops.setNull(); 

   tbox::pout << "\nPASSED:  edge cplxtest" << endl;

   tbox::SAMRAIManager::shutdown();
   tbox::MPI::finalize();
 
   return(0); 
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
bool
complexDataSameAsValue(int desc_id, dcomplex value,
   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy)
{
   bool test_passed = true;

   int ln;
   tbox::Pointer<hier::Patch<2> > patch;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::EdgeData<2,dcomplex> > svdata = patch->getPatchData(desc_id);

         for (pdat::EdgeIterator<2> c(svdata->getBox(),1); c && test_passed ; c++) {
            pdat::EdgeIndex<2> edge_index = c();
            if ( !tbox::Utilities::ceq((*svdata) (edge_index),value) ) {
               test_passed = false;
            }
         }
      }
   }

   return (test_passed);
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
bool
doubleDataSameAsValue(int desc_id, double value,
   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy)
{
   bool test_passed = true;

   int ln;
   tbox::Pointer<hier::Patch<2> > patch;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::EdgeData<2,double> > svdata = patch->getPatchData(desc_id);

         for (pdat::EdgeIterator<2> c(svdata->getBox(),1); c && test_passed ; c++) {
            pdat::EdgeIndex<2> edge_index = c();
            if ( !tbox::Utilities::deq((*svdata) (edge_index),value) ) {
               test_passed = false;
            }
         }
      }
   }

   return (test_passed);
}

