//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dataops/cell_hiertest.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program to test cell-centered patch data ops
//

#include "SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
using namespace std;

#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"

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
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"
#include "VariableContext.h"

using namespace SAMRAI;

/* Helper function declarations */
bool doubleDataSameAsValue(int desc_id, double value,
			   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy);

#define NVARS 4

int main( int argc, char *argv[] ) {
   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("cell_hiertest.log");

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
      const int nproc = tbox::SAMRAI_MPI::getNodes();
      TBOX_ASSERT(nproc < 3);

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

      // Make some dummy variables and data on the hierarchy
      tbox::Pointer< pdat::CellVariable<NDIM,double> > cvar[NVARS];
      int cvindx[NVARS];
      cvar[0] = new pdat::CellVariable<NDIM,double>("cvar0", 1);
      cvindx[0] = variable_db->registerVariableAndContext(
	 cvar[0], dummy, no_ghosts);
      cvar[1] = new pdat::CellVariable<NDIM,double>("cvar1", 1);
      cvindx[1] = variable_db->registerVariableAndContext(
	 cvar[1], dummy, no_ghosts);
      cvar[2] = new pdat::CellVariable<NDIM,double>("cvar2", 1);
      cvindx[2] = variable_db->registerVariableAndContext(
	 cvar[2], dummy, no_ghosts);
      cvar[3] = new pdat::CellVariable<NDIM,double>("cvar3", 1);
      cvindx[3] = variable_db->registerVariableAndContext(
	 cvar[3], dummy, no_ghosts);

      tbox::Pointer< pdat::CellVariable<NDIM,double> >
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

      tbox::Pointer<math::HierarchyDataOpsReal<NDIM,double> > cell_ops = 
	 new math::HierarchyCellDataOpsReal<NDIM,double> (hierarchy, 0, 1);
      TBOX_ASSERT(!cell_ops.isNull());

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
      // 0.0025 on fine level
      bool vol_test_passed = true;
      for (ln = 0; ln < 2; ln++) {

	 tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
	    patch = level->getPatch(ip());
	    tbox::Pointer< pdat::CellData<NDIM,double> > cvdata =
	       patch->getPatchData(cwgt_id);

	    for (pdat::CellIterator<NDIM> c(cvdata->getBox());
		 c && vol_test_passed;c++) {
	       pdat::CellIndex<NDIM> cell_index = c();

	       if (ln == 0) {
		  if ((coarse_fine * patch->getBox()).contains(cell_index)) {
		     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0) ) {
			vol_test_passed = false;
		     }
		  } else {
#if (NDIM == 2)
		     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.01) ) {
#endif
#if (NDIM == 3)
			if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.001) ) {
#endif
			   vol_test_passed = false;
			}
		     }
		  }

		  if (ln == 1) {
#if (NDIM == 2)
		     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0025) ) {
#endif
#if (NDIM == 3)
			if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.000125) ) {
#endif
			   vol_test_passed = false;
			}
		     }
		  }
	       }
	    }
	    if (!vol_test_passed) {
	       num_failures++;
	       tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << endl;
	       cwgt_ops->printData(cwgt_id, tbox::pout);
	    }

	    // Test #1b: math::HierarchyCellDataOpsReal<NDIM>::sumControlVolumes()
	    // Expected: norm = 0.5
	    double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
	    if ( !tbox::MathUtilities<double>::equalEps(norm,0.5)) {
	       num_failures++;
	       tbox::perr << "FAILED: - Test #1b: math::HierarchyCellDataOpsReal<NDIM>::sumControlVolumes()\n"
			  << "Expected value = 0.5 , Computed value = "
			  << norm << endl;
	    }

	    // Test #2: math::HierarchyCellDataOpsReal<NDIM>::numberOfEntries()
	    // Expected: num_data_points = 90 for 2D, 660 for 3D
	    int num_data_points = cell_ops->numberOfEntries(cvindx[0]);
#if (NDIM == 2)
	    if ( num_data_points != 90 ) {
	       num_failures++;
	       tbox::perr << "FAILED: - Test #2: math::HierarchyCellDataOpsReal<NDIM>::numberOfEntries()\n"
			  << "Expected value = 90 , Computed value = "
#endif
#if (NDIM == 3)
		  if ( num_data_points != 660 ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #2: math::HierarchyCellDataOpsReal<NDIM>::numberOfEntries()\n"
				<< "Expected value = 660 , Computed value = "
#endif
				<< num_data_points << endl;
		  }

	       // Test #3a: math::HierarchyCellDataOpsReal<NDIM>::setToScalar()
	       // Expected: v0 = 2.0
	       double val0 = 2.0;
	       cell_ops->setToScalar(cvindx[0], val0);
	       if (!doubleDataSameAsValue(cvindx[0],val0,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #3a: math::HierarchyCellDataOpsReal<NDIM>::setToScalar()\n"
			     << "Expected: v0 = " << val0 << endl;
		  cell_ops->printData(cvindx[0], tbox::pout);
	       }

	       // Test #3b: math::HierarchyCellDataOpsReal<NDIM>::setToScalar()
	       // Expected: v1 = (4.0)
	       cell_ops->setToScalar(cvindx[1], 4.0);
	       double val1 = 4.0;
	       if (!doubleDataSameAsValue(cvindx[1],val1,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #3b: math::HierarchyCellDataOpsReal<NDIM>::setToScalar()\n"
			     << "Expected: v1 = " << val1 << endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }

	       // Test #4: math::HierarchyCellDataOpsReal<NDIM>::copyData()
	       // Expected: v2 = v1 = (4.0)
	       cell_ops->copyData(cvindx[2], cvindx[1]);
	       if (!doubleDataSameAsValue(cvindx[2],val1,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #4: math::HierarchyCellDataOpsReal<NDIM>::copyData()\n"
			     << "Expected: v2 = " << val1 << endl;
		  cell_ops->printData(cvindx[2], tbox::pout);
	       }

	       // Test #5: math::HierarchyCellDataOpsReal<NDIM>::swapData()
	       // Expected: v0 = (4.0), v1 = (2.0)
	       cell_ops->swapData(cvindx[0], cvindx[1]);
	       if (!doubleDataSameAsValue(cvindx[0],val1,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #5a: math::HierarchyCellDataOpsReal<NDIM>::swapData()\n"
			     << "Expected: v0 = " << val1 << endl;
		  cell_ops->printData(cvindx[0], tbox::pout);
	       }
	       if (!doubleDataSameAsValue(cvindx[1],val0,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #5b: math::HierarchyCellDataOpsReal<NDIM>::swapData()\n"
			     << "Expected: v1 = " << val0 << endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }

	       // Test #6: math::HierarchyCellDataOpsReal<NDIM>::scale()
	       // Expected: v2 = 0.25 * v2 = (1.0)
	       cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
	       double val_scale = 1.0;
	       if (!doubleDataSameAsValue(cvindx[2],val_scale,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #6: math::HierarchyCellDataOpsReal<NDIM>::scale()\n"
			     << "Expected: v2 = " << val_scale << endl;
		  cell_ops->printData(cvindx[2], tbox::pout);
	       }

	       // Test #7: math::HierarchyCellDataOpsReal<NDIM>::add()
	       // Expected: v3 = v0 + v1 = (6.0) 
	       cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
	       double val_add = 6.0;
	       if (!doubleDataSameAsValue(cvindx[3],val_add,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #7: math::HierarchyCellDataOpsReal<NDIM>::add()\n"
			     << "Expected: v3 = " << val_add << endl;
		  cell_ops->printData(cvindx[3], tbox::pout);
	       }

	       // Reset v0: v0 = (0.0)
	       cell_ops->setToScalar(cvindx[0], 0.0);

	       // Test #8: math::HierarchyCellDataOpsReal<NDIM>::subtract()
	       // Expected: v1 = v3 - v0 = (6.0)
	       cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);  
	       double val_sub = 6.0;
	       if (!doubleDataSameAsValue(cvindx[1],val_sub,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #8: math::HierarchyCellDataOpsReal<NDIM>::subtract()\n"
			     << "Expected: v1 = " << val_sub << endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }
  
	       // Test #9a: math::HierarchyCellDataOpsReal<NDIM>::addScalar()
	       // Expected: v1 = v1 + (0.0) = (6.0) 
	       cell_ops->addScalar(cvindx[1], cvindx[1], 0.0);
	       double val_addScalar = 6.0;
	       if (!doubleDataSameAsValue(cvindx[1],val_addScalar,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #9a: math::HierarchyCellDataOpsReal<NDIM>::addScalar()\n"
			     << "Expected: v1 = " << val_addScalar << endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }
  
	       // Test #9b: math::HierarchyCellDataOpsReal<NDIM>::addScalar()
	       // Expected: v2 = v2 + (0.0) = (1.0) 
	       cell_ops->addScalar(cvindx[2], cvindx[2], 0.0); 
	       val_addScalar = 1.0;
	       if (!doubleDataSameAsValue(cvindx[2],val_addScalar,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #9b: math::HierarchyCellDataOpsReal<NDIM>::addScalar()\n"
			     << "Expected: v2 = " << val_addScalar << endl;
		  cell_ops->printData(cvindx[2], tbox::pout);
	       }

	       // Test #9c: math::HierarchyCellDataOpsReal<NDIM>::addScalar()
	       // Expected: v2 = v2 + (3.0) = (4.0)
	       cell_ops->addScalar(cvindx[2], cvindx[2], 3.0);
	       val_addScalar = 4.0;
	       if (!doubleDataSameAsValue(cvindx[2],val_addScalar,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #9c: math::HierarchyCellDataOpsReal<NDIM>::addScalar()\n"
			     << "Expected: v2 = " << val_addScalar << endl;
		  cell_ops->printData(cvindx[2], tbox::pout);
	       }

	       // Reset v3:  v3 = (0.5)
	       cell_ops->setToScalar(cvindx[3], 0.5);

	       // Test #10: math::HierarchyCellDataOpsReal<NDIM>::multiply()
	       // Expected: v1 = v3 * v1 = (3.0) 
	       cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
	       double val_mult = 3.0;
	       if (!doubleDataSameAsValue(cvindx[1],val_mult,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #10: math::HierarchyCellDataOpsReal<NDIM>::multiply()\n"
			     << "Expected: v1 = " << val_mult<< endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }

	       // Test #11: math::HierarchyCellDataOpsReal<NDIM>::divide()
	       // Expected: v0 = v2 / v1 = 1.3333333333
	       cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
	       double val_div = 1.33333333333;
	       if (!doubleDataSameAsValue(cvindx[0],val_div,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #11: math::HierarchyCellDataOpsReal<NDIM>::divide()\n"
			     << "Expected: v0 = " << val_div << endl;
		  cell_ops->printData(cvindx[0], tbox::pout);
	       }

	       // Test #12: math::HierarchyCellDataOpsReal<NDIM>::reciprocal()
	       // Expected:  v1 = 1 / v1 = (0.333333333)
	       cell_ops->reciprocal(cvindx[1], cvindx[1]); 
	       double val_rec = 0.33333333333;
	       if (!doubleDataSameAsValue(cvindx[1],val_rec,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #12: math::HierarchyCellDataOpsReal<NDIM>::reciprocal()\n"
			     << "Expected: v1 = " << val_rec << endl;
		  cell_ops->printData(cvindx[1], tbox::pout);
	       }

	       // Test #13: math::HierarchyCellDataOpsReal<NDIM>::abs()
	       // Expected:  v3 = abs(v2) = 4.0
	       cell_ops->abs(cvindx[3], cvindx[2]);
	       double val_abs = 4.0;
	       if (!doubleDataSameAsValue(cvindx[3],val_abs,hierarchy)) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #13: math::HierarchyCellDataOpsReal<NDIM>::abs()\n"
			     << "Expected: v3 = " << val_abs << endl;
		  cell_ops->printData(cvindx[3], tbox::pout);
	       }

	       // Test #14: Place some bogus values on coarse level 
	       tbox::Pointer< pdat::CellData<NDIM,double> > cdata;

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
		     (*cdata)(pdat::CellIndex<NDIM>(index0), 0) = 100.0;
		  }
		  if (patch->getBox().contains(index1)) {
		     (*cdata)(pdat::CellIndex<NDIM>(index1), 0) = -1000.0;
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

		  for (pdat::CellIterator<NDIM> c(cdata->getBox());
		       c && bogus_value_test_passed;c++) {
		     pdat::CellIndex<NDIM> cell_index = c();

		     if (cell_index == index0) {
			if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index),100.0)) {
			   bogus_value_test_passed = false;
			}
		     } else {
			if (cell_index == index1) {
			   if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index),-1000.0)) {
			      bogus_value_test_passed = false;
			   }
			} else {
			   if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index), 4.0)) {
			      bogus_value_test_passed = false;
			   }
			}
		     }
		  }
	       }
	       if (!bogus_value_test_passed) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #14:  Place some bogus values on coarse level" << endl;
		  cell_ops->printData(cvindx[2], tbox::pout);
	       }


	       // Test #15: math::HierarchyCellDataOpsReal<NDIM>::L1Norm() - w/o control weight
	       // Expected:  bogus_l1_norm = 1452 in 2d, 3732 in 3d
	       double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
#if (NDIM == 2)
	       if ( !tbox::MathUtilities<double>::equalEps(bogus_l1_norm,1452) ) {
		  num_failures++;
		  tbox::perr << "FAILED: - Test #15: math::HierarchyCellDataOpsReal<NDIM>::L1Norm()"
			     << " - w/o control weight\n"
			     << "Expected value = 1452, Computed value = "
#endif
#if (NDIM == 3)
		     if ( !tbox::MathUtilities<double>::equalEps(bogus_l1_norm,3732) ) {
			num_failures++;
			tbox::perr << "FAILED: - Test #15: math::HierarchyCellDataOpsReal<NDIM>::L1Norm()"
				   << " - w/o control weight\n"
				   << "Expected value = 3732, Computed value = "
#endif
				   << setprecision(12) << bogus_l1_norm << endl;
		     }

		  // Test #16: math::HierarchyCellDataOpsReal<NDIM>::L1Norm() - w/control weight
		  // Expected:  correct_l1_norm = 2.0
		  double correct_l1_norm = cell_ops->L1Norm(cvindx[2],cwgt_id);
		  if ( !tbox::MathUtilities<double>::equalEps(correct_l1_norm,2.0) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #16: math::HierarchyCellDataOpsReal<NDIM>::L1Norm()"
				<< " - w/control weight\n"
				<< "Expected value = 2.0, Computed value = "
				<< correct_l1_norm << endl;
		  }

		  // Test #17: math::HierarchyCellDataOpsReal<NDIM>::L2Norm()
		  // Expected:  l2_norm = 2.82842712475
		  double l2_norm = cell_ops->L2Norm(cvindx[2],cwgt_id);
		  if ( !tbox::MathUtilities<double>::equalEps(l2_norm,2.82842712475) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #17: math::HierarchyCellDataOpsReal<NDIM>::L2Norm()\n"
				<< "Expected value = 2.82842712475, Computed value = "
				<< l2_norm << endl;
		  }

		  // Test #18: math::HierarchyCellDataOpsReal<NDIM>::L2Norm() - w/o control weight
		  // Expected:  bogus_max_norm = 1000.0
		  double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
		  if ( !tbox::MathUtilities<double>::equalEps(bogus_max_norm,1000.0) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #18: math::HierarchyCellDataOpsReal<NDIM>::L2Norm()"
				<< " - w/o control weight\n"
				<< "Expected value = 1000.0, Computed value = "
				<< bogus_max_norm << endl;
		  }

		  // Test #19: math::HierarchyCellDataOpsReal<NDIM>::L2Norm() - w/control weight
		  // Expected:  max_norm = 4.0
		  double max_norm = cell_ops->maxNorm(cvindx[2],cwgt_id);
		  if ( !tbox::MathUtilities<double>::equalEps(max_norm,4.0) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #19: math::HierarchyCellDataOpsReal<NDIM>::L2Norm()"
				<< " - w/control weight\n"
				<< "Expected value = 4.0, Computed value = "
				<< max_norm << endl;
		  }

		  // Reset data and test sums, axpy's
		  cell_ops->setToScalar(cvindx[0], 1.00); 
		  cell_ops->setToScalar(cvindx[1], 2.5); 
		  cell_ops->setToScalar(cvindx[2], 7.0); 


		  // Test #20: math::HierarchyCellDataOpsReal<NDIM>::linearSum()
		  // Expected:  v3 = 5.0
		  cell_ops->linearSum(cvindx[3], 2.0, cvindx[1], 0.00, cvindx[0]);
		  double val_linearSum = 5.0;
		  if (!doubleDataSameAsValue(cvindx[3],val_linearSum,hierarchy)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #20: math::HierarchyCellDataOpsReal<NDIM>::linearSum()\n"
				<< "Expected: v3 = " << val_linearSum << endl;
		     cell_ops->printData(cvindx[3], tbox::pout);
		  }

		  // Test #21: math::HierarchyCellDataOpsReal<NDIM>::axmy()
		  // Expected:  v3 = 6.5
		  cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
		  double val_axmy = 6.5;
		  if (!doubleDataSameAsValue(cvindx[3],val_axmy,hierarchy)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #21: math::HierarchyCellDataOpsReal<NDIM>::axmy()\n"
				<< "Expected: v3 = " << val_axmy << endl;
		     cell_ops->printData(cvindx[3], tbox::pout);
		  }

		  // Test #22a: math::HierarchyCellDataOpsReal<NDIM>::dot() - (ind2) * (ind1)
		  // Expected:  cdot = 8.75
		  double cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
		  if ( !tbox::MathUtilities<double>::equalEps(cdot,8.75) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #22a: math::HierarchyCellDataOpsReal<NDIM>::dot() - (ind2) * (ind1)\n"
				<< "Expected Value = 8.75, Computed Value = "
				<< cdot << endl;
		  }

		  // Test #22b: math::HierarchyCellDataOpsReal<NDIM>::dot() - (ind1) * (ind2)
		  // Expected:  cdot = 8.75
		  cdot = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id); 
		  if ( !tbox::MathUtilities<double>::equalEps(cdot,8.75) ) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #22b: math::HierarchyCellDataOpsReal<NDIM>::dot() - (ind1) * (ind2)\n"
				<< "Expected Value = 8.75, Computed Value = "
				<< cdot << endl;
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

		  if (num_failures == 0) {
		     tbox::pout << "\nPASSED:  cell hiertest" << endl;
		  }
	       }

	       tbox::SAMRAIManager::shutdown();
	       tbox::SAMRAI_MPI::finalize();
 
	       return(num_failures); 
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
			   if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
			      test_passed = false;
			   }
			}
		     }
		  }

		  return (test_passed);
	       }
