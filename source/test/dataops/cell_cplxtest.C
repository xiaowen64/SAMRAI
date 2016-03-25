//
// File:        cell_cplxtest.C
// Package:     SAMRAI test
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 676 $
// Modified:    $Date: 2005-10-14 15:23:35 -0700 (Fri, 14 Oct 2005) $
// Description: Main program to test cell-centered complex patch data ops
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
#include "CellData.h"
#include "HierarchyDataOpsComplex.h"
#include "HierarchyCellDataOpsComplex.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyCellDataOpsReal.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
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
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy);
bool doubleDataSameAsValue(int desc_id, double value,
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy);

#define NVARS 4

int main( int argc, char *argv[] ) {
   tbox::MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();
//   tbox::PIO::logOnlyNodeZero("cell_cplxtest.log");
   tbox::PIO::logAllNodes("cell_cplxtest.log");

   int ln, iv;

   // Make a dummy hierarchy domain
   double lo[NDIM];
   double hi[NDIM];

   hier::Index<NDIM> clo0;
   hier::Index<NDIM> chi0;
   hier::Index<NDIM> clo1;
   hier::Index<NDIM> chi1;
   hier::Index<NDIM> flo0;
   hier::Index<NDIM> fhi0;
   hier::Index<NDIM> flo1;
   hier::Index<NDIM> fhi1;

   for (int i = 0; i < NDIM; i++) {
      lo[i] = 0.0;
      clo0(i) = 0;
      flo0(i) = 4;
      fhi0(i) = 7;
      if (i == 1) {
         hi[i] = 0.5;
         chi0(i) = 2;
         clo1(i) = 3;
         chi1(i) = 4;
      } else {
         hi[i] = 1.0;
         chi0(i) = 9;
         clo1(i) = 0;
         chi1(i) = 9;
      }
      if (i == 0) {
         flo1(i) = 8;
         fhi1(i) = 13;
      } else {
         flo1(i) = flo0(i);
         fhi1(i) = fhi0(i);
      }
   }

   hier::Box<NDIM> coarse0(clo0, chi0);
   hier::Box<NDIM> coarse1(clo1, chi1);
   hier::Box<NDIM> fine0(flo0, fhi0);
   hier::Box<NDIM> fine1(flo1, fhi1);
   hier::IntVector<NDIM> ratio(2);

   hier::BoxArray<NDIM> coarse_domain(2);
   hier::BoxArray<NDIM> fine_domain(2);
   coarse_domain(0) = coarse0;
   coarse_domain(1) = coarse1;
   fine_domain(0) = fine0;
   fine_domain(1) = fine1;

   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > geometry =
   new geom::CartesianGridGeometry<NDIM>("CartesianGeometry", lo, hi, coarse_domain);

   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = 
      new hier::PatchHierarchy<NDIM>("PatchHierarchy", geometry);

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

   hierarchy->makeNewPatchLevel(0, hier::IntVector<NDIM>(1), coarse_domain, mapping0);
   hierarchy->makeNewPatchLevel(1, ratio, fine_domain, mapping1);

   // Create instance of hier::Variable<NDIM> database
   hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer< hier::VariableContext > dummy = variable_db->getContext("dummy");
   const hier::IntVector<NDIM> no_ghosts(0);

   // Make some dummy variables and register them with hier::VariableDatabase<NDIM>
   tbox::Pointer< pdat::CellVariable<NDIM,dcomplex> > cvar[NVARS];
   int cvindx[NVARS];
   cvar[0] = new pdat::CellVariable<NDIM,dcomplex>("cvar0", 1);
   cvindx[0] = variable_db->registerVariableAndContext(
      cvar[0], dummy, no_ghosts);
   cvar[1] = new pdat::CellVariable<NDIM,dcomplex>("cvar1", 1);
   cvindx[1] = variable_db->registerVariableAndContext(
      cvar[1], dummy, no_ghosts);
   cvar[2] = new pdat::CellVariable<NDIM,dcomplex>("cvar2", 1);
   cvindx[2] = variable_db->registerVariableAndContext(
      cvar[2], dummy, no_ghosts);
   cvar[3] = new pdat::CellVariable<NDIM,dcomplex>("cvar3", 1);
   cvindx[3] = variable_db->registerVariableAndContext(
      cvar[3], dummy, no_ghosts);

   tbox::Pointer< pdat::CellVariable<NDIM,double> > cwgt;
   cwgt = new pdat::CellVariable<NDIM,double>("cwgt", 1);
   int cwgt_id = variable_db->registerVariableAndContext(
      cwgt, dummy, no_ghosts);

   // allocate data on hierarchy
   for (ln = 0; ln < 2; ln++) {
      hierarchy->getPatchLevel(ln)->allocatePatchData(cwgt_id);
      for (iv = 0; iv < NVARS; iv++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(cvindx[iv]);
      }
   }

   tbox::Pointer<math::HierarchyDataOpsComplex<NDIM> > cell_ops = 
   new math::HierarchyCellDataOpsComplex<NDIM>(hierarchy, 0, 1);
   assert(!cell_ops.isNull());

   tbox::Pointer<math::HierarchyDataOpsReal<NDIM,double> > cwgt_ops = 
   new math::HierarchyCellDataOpsReal<NDIM,double>(hierarchy, 0 , 1);

   tbox::Pointer<hier::Patch<NDIM> > patch;
   tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom;

   // Initialize control volume data for cell-centered components
   hier::Box<NDIM> coarse_fine = fine0+fine1; 
   coarse_fine.coarsen(ratio);
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         pgeom = patch->getPatchGeometry();
         const double* dx = pgeom->getDx();
         double cell_vol = dx[0];
         for (int i = 1; i < NDIM; i++) {
            cell_vol *= dx[i];
         }
         tbox::Pointer< pdat::CellData<NDIM,double> > cvdata =
            patch->getPatchData(cwgt_id);
         cvdata->fillAll(cell_vol);
         if (ln == 0) cvdata->fillAll(0.0, (coarse_fine * patch->getBox()) );
      }
   }

   // Test #1: Print out control volume data and compute its integral

   // Test #1a: Check control volume data set properly
   // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
   // 0.0025 on fine level for 2d.  0.001 and 0.000125 for 3d.
   bool vol_test_passed = true;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = patch->getPatchData(cwgt_id);

         for (pdat::CellIterator<NDIM> c(cvdata->getBox());c && vol_test_passed;c++) {
            pdat::CellIndex<NDIM> cell_index = c();

            if (ln == 0) {
               if ((coarse_fine * patch->getBox()).contains(cell_index)) {
                  if ( !tbox::Utilities::deq((*cvdata)(cell_index),0.0) ) {
                     vol_test_passed = false;
                  } 
               } else {
#if (NDIM == 2)
                  if ( !tbox::Utilities::deq((*cvdata)(cell_index),0.01) ) {
#endif
#if (NDIM == 3)
                  if ( !tbox::Utilities::deq((*cvdata)(cell_index),0.001) ) {
#endif
                     vol_test_passed = false;
                  }
               }
            }

            if (ln == 1) {
#if (NDIM == 2)
               if ( !tbox::Utilities::deq((*cvdata)(cell_index),0.0025) ) {
#endif
#if (NDIM == 3)
               if ( !tbox::Utilities::deq((*cvdata)(cell_index),0.000125) ) {
#endif
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

   // Test #1b: math::HierarchyCellDataOpsComplex<NDIM>::sumControlVolumes()
   // Expected: norm = 0.5
   double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
   if ( !tbox::Utilities::deq(norm,0.5)) {
      tbox::perr << "FAILED: - Test #1b: math::HierarchyCellDataOpsComplex<NDIM>::sumControlVolumes()\n"
           << "Expected value = 0.5 , Computed value = " 
           << norm << endl;
   }

   // Test #2: math::HierarchyCellDataOpsComplex<NDIM>::numberOfEntries()
   // Expected: num_data_points = 90 in 2d, 660 in 3d
   int num_data_points = cell_ops->numberOfEntries(cvindx[0]);
#if (NDIM == 2)
   if ( num_data_points != 90 ) {
      tbox::perr << "FAILED: - Test #2: math::HierarchyCellDataOpsReal<NDIM>::numberOfEntries()\n"
           << "Expected value = 90 , Computed value = "
#endif
#if (NDIM == 3)
   if ( num_data_points != 660 ) {
      tbox::perr << "FAILED: - Test #2: math::HierarchyCellDataOpsReal<NDIM>::numberOfEntries()\n"
           << "Expected value = 660 , Computed value = "
#endif
           << num_data_points << endl;
   }

   // Test #3a: math::HierarchyCellDataOpsComplex<NDIM>::setToScalar()
   // Expected: v0 = (2.0,1.5)
   dcomplex val0 = dcomplex(2.0,1.5);
   cell_ops->setToScalar(cvindx[0], val0);
   if (!complexDataSameAsValue(cvindx[0],val0,hierarchy)) {
      tbox::perr << "FAILED: - Test #3a: math::HierarchyCellDataOpsComplex<NDIM>::setToScalar()\n"
           << "Expected: v0 = " << val0 << endl;
      cell_ops->printData(cvindx[0], tbox::plog);
   }

   // Test #3b: math::HierarchyCellDataOpsComplex<NDIM>::setToScalar()
   // Expected: v1 = (4.0,3.0)
   dcomplex val1(4.0, 3.0);
   cell_ops->setToScalar(cvindx[1], val1);
   if (!complexDataSameAsValue(cvindx[1],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #3b: math::HierarchyCellDataOpsComplex<NDIM>::setToScalar()\n"
           << "Expected: v1 = " << val1 << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #4: math::HierarchyCellDataOpsComplex<NDIM>::copyData()
   // Expected: v2 = v1 = (4.0, 3.0)
   cell_ops->copyData(cvindx[2], cvindx[1]);
   if (!complexDataSameAsValue(cvindx[2],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #4: math::HierarchyCellDataOpsComplex<NDIM>::copyData()\n"
           << "Expected: v2 = v1 = " << val1 << endl;
      cell_ops->printData(cvindx[2], tbox::plog);
   }

   // Test #5: math::HierarchyCellDataOpsComplex<NDIM>::swapData()
   // Expected: v0 = (4.0, 3.0), v1 = (2.0,1.5)
   cell_ops->swapData(cvindx[0], cvindx[1]);
   if (!complexDataSameAsValue(cvindx[0],val1,hierarchy)) {
      tbox::perr << "FAILED: - Test #5a: math::HierarchyCellDataOpsComplex<NDIM>::swapData()\n"
           << "Expected: v0 = " << val1 << endl;
      cell_ops->printData(cvindx[0], tbox::plog);
   }
   if (!complexDataSameAsValue(cvindx[1],val0,hierarchy)) {
      tbox::perr << "FAILED: - Test #5b: math::HierarchyCellDataOpsComplex<NDIM>::swapData()\n"
           << "Expected: v1 = " << val0 << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #6: math::HierarchyCellDataOpsComplex<NDIM>::scale()
   // Expected: v2 = 0.25 * v2 = (1.0,0.75)
   cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
   dcomplex val_scale(1.0,0.75);
   if (!complexDataSameAsValue(cvindx[2],val_scale,hierarchy)) {
      tbox::perr << "FAILED: - Test #6: math::HierarchyCellDataOpsComplex<NDIM>::scale()\n"
           << "Expected: v2 = " << val_scale << endl;
      cell_ops->printData(cvindx[2], tbox::plog);
   }

   // Test #7: math::HierarchyCellDataOpsComplex<NDIM>::add()
   // Expected: v3 = v0 + v1 = (6.0, 4.5) 
   cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
   dcomplex val_add(6.0,4.5);
   if (!complexDataSameAsValue(cvindx[3],val_add,hierarchy)) {
      tbox::perr << "FAILED: - Test #7: math::HierarchyCellDataOpsComplex<NDIM>::add()\n"
           << "Expected: v3 = " << val_add << endl;
      cell_ops->printData(cvindx[3], tbox::plog);
   }

   // Reset v0:  v0 = (0.0,4.5)
   cell_ops->setToScalar(cvindx[0], dcomplex(0.0,4.5));

   // Test #8: math::HierarchyCellDataOpsComplex<NDIM>::subtract()
   // Expected: v1 = v3 - v0 = (6.0,0.0)
   cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);  
   dcomplex val_sub(6.0,0.0);
   if (!complexDataSameAsValue(cvindx[1],val_sub,hierarchy)) {
      tbox::perr << "FAILED: - Test #8: math::HierarchyCellDataOpsComplex<NDIM>::subtract()\n"
           << "Expected: v1 = " << val_sub << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #9a: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()
   // Expected: v1 = v1 + (0.0,-4.0) = (6.0,-4.0) 
   cell_ops->addScalar(cvindx[1], cvindx[1], dcomplex(0.0,-4.0));
   dcomplex val_addScalar(6.0,-4.0);
   if (!complexDataSameAsValue(cvindx[1],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9a: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()\n"
           << "Expected: v1 = " << val_addScalar << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #9b: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()
   // Expected:   v2 = v2 + (0.0,0.25) = (1.0,1.0) 
   cell_ops->addScalar(cvindx[2], cvindx[2], dcomplex(0.0,0.25)); 
   val_addScalar = dcomplex(1.0,1.0);
   if (!complexDataSameAsValue(cvindx[2],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9c: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()\n"
           << "Expected: v2 = " << val_addScalar << endl;
      cell_ops->printData(cvindx[2], tbox::plog);
   }

   // Test #9c: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()
   // Expected:  v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
   cell_ops->addScalar(cvindx[2], cvindx[2], dcomplex(3.0,-4.0));
   val_addScalar = dcomplex(4.0,-3.0);
   if (!complexDataSameAsValue(cvindx[2],val_addScalar,hierarchy)) {
      tbox::perr << "FAILED: - Test #9d: math::HierarchyCellDataOpsComplex<NDIM>::addScalar()\n"
           << "Expected: v2 = " << val_addScalar << endl;
      cell_ops->printData(cvindx[2], tbox::plog);
   }

   // Reset v3: v3 = (0.5, 0.0)
   cell_ops->setToScalar(cvindx[3], dcomplex(0.5, 0.0));

   // Test #10: math::HierarchyCellDataOpsComplex<NDIM>::multiply()
   // Expected:  v1 = v3 * v1 = (3.0,-2.0) 
   cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
   dcomplex val_mult(3.0,-2.0);
   if (!complexDataSameAsValue(cvindx[1],val_mult,hierarchy)) {
      tbox::perr << "FAILED: - Test #10: math::HierarchyCellDataOpsComplex<NDIM>::multiply()\n"
           << "Expected: v1 = " << val_mult << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #11: math::HierarchyCellDataOpsComplex<NDIM>::divide()
   // Expected:  v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
   cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
   dcomplex val_div(1.3846153846154,-0.076923076923077);
   if (!complexDataSameAsValue(cvindx[0],val_div,hierarchy)) {
      tbox::perr << "FAILED: - Test #11: math::HierarchyCellDataOpsComplex<NDIM>::divide()\n"
           << "Expected: v0 = " << val_div << endl;
      cell_ops->printData(cvindx[0], tbox::plog);
   }

   // Test #12: math::HierarchyCellDataOpsComplex<NDIM>::reciprocal()
   // Expected:  v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
   cell_ops->reciprocal(cvindx[1], cvindx[1]); 
   dcomplex val_rec(0.23076923076923, 0.15384615384615);
   if (!complexDataSameAsValue(cvindx[1],val_rec,hierarchy)) {
      tbox::perr << "FAILED: - Test #12: math::HierarchyCellDataOpsComplex<NDIM>::reciprocal()\n"
           << "Expected: v1 = " << val_rec << endl;
      cell_ops->printData(cvindx[1], tbox::plog);
   }

   // Test #13:  Place some bogus values on coarse level
   tbox::Pointer< pdat::CellData<NDIM,dcomplex> > cdata;

   // set values
   tbox::Pointer<hier::PatchLevel<NDIM> > level_zero 
      = hierarchy->getPatchLevel(0);
   for (hier::PatchLevel<NDIM>::Iterator ip(level_zero); ip; ip++) {
      patch = level_zero->getPatch(ip());
      cdata = patch->getPatchData(cvindx[2]);
      hier::Index<NDIM> index0(2);
      hier::Index<NDIM> index1(3);
      index1(0) = 5;
      if (patch->getBox().contains(index0)) {
         (*cdata)(pdat::CellIndex<NDIM>(index0), 0) = dcomplex(100.0,-50.0);
      }
      if (patch->getBox().contains(index1)) {
         (*cdata)(pdat::CellIndex<NDIM>(index1), 0) = dcomplex(-1000.0,20.0);
      }
   }

   // check values
   bool bogus_value_test_passed = true;
   for (hier::PatchLevel<NDIM>::Iterator ipp(level_zero); ipp; ipp++) {
      patch = level_zero->getPatch(ipp());
      cdata = patch->getPatchData(cvindx[2]);
      hier::Index<NDIM> index0(2);
      hier::Index<NDIM> index1(3);
      index1(0) = 5;

      for (pdat::CellIterator<NDIM> c(cdata->getBox());c && bogus_value_test_passed;c++) {
            pdat::CellIndex<NDIM> cell_index = c();

         if (cell_index == index0) {
            if (!tbox::Utilities::ceq((*cdata)(cell_index),
                dcomplex(100.0,-50.0)) ) {
               bogus_value_test_passed = false;
            }
         } else {
            if (cell_index == index1) {
               if (!tbox::Utilities::ceq((*cdata)(cell_index),
                  dcomplex(-1000.0,20.0)) ) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (!tbox::Utilities::ceq((*cdata)(cell_index),
                  dcomplex(4.0,-3.0)) ) {
                  bogus_value_test_passed = false;
               }
            }
         }
      }
   }
   if (!bogus_value_test_passed) {
      tbox::perr << "FAILED: - Test #13:  Place some bogus values on coarse level" << endl;
      cell_ops->printData(cvindx[2], tbox::plog);
   }

   // Test norms on patch data with cvindx[2] on hierarchy with bogus values


   // Test #14: math::HierarchyCellDataOpsComplex<NDIM>::L1Norm() - w/o control weight
   // Expected:  bogus_l1_norm = 1552.00337888 in 2d, 4402.00337888 in 3d.
   double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
#if (NDIM == 2)
   if ( !tbox::Utilities::deq(bogus_l1_norm,1552.00337888) ) {
#endif
#if (NDIM == 3)
   if ( !tbox::Utilities::deq(bogus_l1_norm,4402.00337888) ) {
#endif
      tbox::perr << "FAILED: - Test #14: math::HierarchyCellDataOpsComplex<NDIM>::L1Norm()"
           << " - w/o control weight\n"
#if (NDIM == 2)
           << "Expected value = 1552.00337888 , Computed value = "
#endif
#if (NDIM == 3)
           << "Expected value = 4402.00337888 , Computed value = "
#endif
           << setprecision(12) << bogus_l1_norm << endl;
   }

   // Test #15: math::HierarchyCellDataOpsComplex<NDIM>::L1Norm() - w/control weight
   // Expected:  correct_l1_norm = 2.5
   double correct_l1_norm = cell_ops->L1Norm(cvindx[2],cwgt_id);
   if ( !tbox::Utilities::deq(correct_l1_norm,2.5) ) {
      tbox::perr << "FAILED: - Test #15: math::HierarchyCellDataOpsComplex<NDIM>::L1Norm()"
           << " - w/control weight\n" 
           << "Expected value = 2.5, Computed value = "
           << correct_l1_norm << endl;
   }

   // Test #16: math::HierarchyCellDataOpsComplex<NDIM>::L2Norm() 
   // Expected:  l2_norm = 3.53553390593
   double l2_norm = cell_ops->L2Norm(cvindx[2],cwgt_id);
   if ( !tbox::Utilities::deq(l2_norm,3.53553390593) ) {
      tbox::perr << "FAILED: - Test #16: math::HierarchyCellDataOpsComplex<NDIM>::L2Norm()\n"  
           << "Expected value = 3.53553390593, Computed value = "
           << l2_norm << endl;
   }

   // Test #17: math::HierarchyCellDataOpsComplex<NDIM>::maxNorm() - w/o control weight
   // Expected:  bogus_max_norm = 1000.19998
   double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
   if ( !tbox::Utilities::deq(bogus_max_norm,1000.19998) ) {
      tbox::perr << "FAILED: - Test #17: math::HierarchyCellDataOpsComplex<NDIM>::maxNorm() "
           << "- w/o control weight\n"
           << "Expected value = 1000.19998, Computed value = "
           << bogus_max_norm << endl;
   }

   // Test #18: math::HierarchyCellDataOpsComplex<NDIM>::maxNorm() - w/control weight
   // Expected:  max_norm = 5.0
   double max_norm = cell_ops->maxNorm(cvindx[2],cwgt_id);
   if ( !tbox::Utilities::deq(max_norm,5.0) ) {
      tbox::perr << "FAILED: - Test #18: math::HierarchyCellDataOpsComplex<NDIM>::maxNorm() "
           << "- w/control weight\n" 
           << "Expected value = 5.0, Computed value = "
           << max_norm << endl;
   }

   // Reset data and test sums, axpy's
   cell_ops->setToScalar(cvindx[0], dcomplex(1.0, -3.0)); 
   cell_ops->setToScalar(cvindx[1], dcomplex(2.5, 3.0)); 
   cell_ops->setToScalar(cvindx[2], dcomplex(7.0, 0.0)); 

   // Test #19: math::HierarchyCellDataOpsComplex<NDIM>::linearSum()
   // Expected:  v3 = (2.0,5.0)
   cell_ops->linearSum(cvindx[3], dcomplex(2.0,0.0), cvindx[1], dcomplex(0.0, -1.0), cvindx[0]);
   dcomplex val_linearSum(2.0,5.0);
   if (!complexDataSameAsValue(cvindx[3],val_linearSum,hierarchy)) {
      tbox::perr << "FAILED: - Test #19: math::HierarchyCellDataOpsComplex<NDIM>::linearSum()\n"
           << "Expected: v3 = " << val_linearSum << endl;
      cell_ops->printData(cvindx[3], tbox::plog);
   }
  
   // Test #20: math::HierarchyCellDataOpsComplex<NDIM>::axmy()
   // Expected:  v3 = (6.5,12.0)
   cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
   dcomplex val_axmy(6.5,12.0);
   if (!complexDataSameAsValue(cvindx[3],val_axmy,hierarchy)) {
      tbox::perr << "FAILED: - Test #20: math::HierarchyCellDataOpsComplex<NDIM>::axmy()\n"
           << "Expected: v3 = " << val_axmy << endl;
      cell_ops->printData(cvindx[3], tbox::plog);
   }

   // Test #21a: math::HierarchyCellDataOpsComplex<NDIM>::dot()
   // Expected:  cdot = (8.75,-10.5)
   dcomplex cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
   dcomplex ans_2_dot_1(8.75,-10.5);
   if ( !tbox::Utilities::ceq(cdot,ans_2_dot_1) ) {
      tbox::perr << "FAILED: - Test #21a: math::HierarchyCellDataOpsComplex<NDIM>::dot()\n"
           << "Expected value = (8.75,-10.5), Computed value = "
           << cdot << endl;
   }

   // Test #21b: math::HierarchyCellDataOpsComplex<NDIM>::dot()
   // Expected:  cdot = (8.75,10.5)
   dcomplex cdot2 = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id); 
   dcomplex ans_1_dot_2(8.75,10.5);
   if ( !tbox::Utilities::ceq(cdot2,ans_1_dot_2) ) {
      tbox::perr << "FAILED: - Test #21b: math::HierarchyCellDataOpsComplex<NDIM>::dot()\n"
           << "Expected value = (8.75,10.5), Computed value = "
           << cdot2 << endl;
   }

   // Test #22: math::HierarchyCellDataOpsComplex<NDIM>::abs()
   // Expected:  abs(v0) = 5.0
   cell_ops->setToScalar(cvindx[0], dcomplex(4.0, -3.0));
   cell_ops->abs(cwgt_id, cvindx[0]);
   if (!doubleDataSameAsValue(cwgt_id,5.0,hierarchy)) {
      tbox::perr << "FAILED: - Test #22: math::HierarchyCellDataOpsComplex<NDIM>::abs()\n"
           << "Expected: abs(v0) = 5.0" << endl;
      cwgt_ops->printData(cwgt_id, tbox::plog);
   }

   // deallocate data on hierarchy
   for (ln = 0; ln < 2; ln++) {
      hierarchy->getPatchLevel(ln)->deallocatePatchData(cwgt_id);
      for (iv = 0; iv < NVARS; iv++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(cvindx[iv]);
      }
   }

   for (iv = 0; iv < NVARS; iv++) {
      cvar[iv].setNull();
   }
   cwgt.setNull();

   geometry.setNull();
   hierarchy.setNull();
   cell_ops.setNull(); 
   cwgt_ops.setNull(); 

   tbox::pout << "\nPASSED:  cell cplxtest" << endl;

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
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy)
{
   bool test_passed = true;

   int ln;
   tbox::Pointer<hier::Patch<NDIM> > patch;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::CellData<NDIM,dcomplex> > cvdata = patch->getPatchData(desc_id);

         for (pdat::CellIterator<NDIM> c(cvdata->getBox()); c && test_passed ; c++) {
            pdat::CellIndex<NDIM> cell_index = c();
            if ( !tbox::Utilities::ceq((*cvdata) (cell_index),value) ) {
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
   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy)
{
   bool test_passed = true;

   int ln;
   tbox::Pointer<hier::Patch<NDIM> > patch;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = patch->getPatchData(desc_id);

         for (pdat::CellIterator<NDIM> c(cvdata->getBox()); c && test_passed ; c++) {
            pdat::CellIndex<NDIM> cell_index = c();
            if ( !tbox::Utilities::deq((*cvdata) (cell_index),value) ) {
               test_passed = false;
            }
         }
      }
   }

   return (test_passed);
}

