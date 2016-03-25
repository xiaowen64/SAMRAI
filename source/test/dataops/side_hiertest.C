//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dataops/side_hiertest.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program to test side-centered patch data ops
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
#include "SideData.h"
#include "HierarchyDataOpsComplex.h"
#include "HierarchySideDataOpsComplex.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "SideIndex.h"
#include "SideIterator.h"
#include "SideVariable.h"
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
			   tbox::Pointer<hier::PatchHierarchy<2> > hierarchy);

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

      tbox::PIO::logAllNodes("side_hiertest.log");
   
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
   
      hierarchy->makeNewPatchLevel(0, hier::IntVector<2>(1), coarse_domain, mapping0);
      hierarchy->makeNewPatchLevel(1, ratio, fine_domain, mapping1);
   
      // Create instance of hier::Variable<NDIM> database
      hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
      tbox::Pointer< hier::VariableContext > dummy = variable_db->getContext("dummy");
      const hier::IntVector<2> no_ghosts(0);
   
      // Make some dummy variables and data on the hierarchy
      tbox::Pointer< pdat::SideVariable<2,double> > fvar[NVARS];
      int svindx[NVARS];
      fvar[0] = new pdat::SideVariable<2,double>("fvar0", 1);
      svindx[0] = variable_db->registerVariableAndContext(
	 fvar[0], dummy, no_ghosts);
      fvar[1] = new pdat::SideVariable<2,double>("fvar1", 1);
      svindx[1] = variable_db->registerVariableAndContext(
	 fvar[1], dummy, no_ghosts);
      fvar[2] = new pdat::SideVariable<2,double>("fvar2", 1);
      svindx[2] = variable_db->registerVariableAndContext(
	 fvar[2], dummy, no_ghosts);
      fvar[3] = new pdat::SideVariable<2,double>("fvar3", 1);
      svindx[3] = variable_db->registerVariableAndContext(
	 fvar[3], dummy, no_ghosts);
   
      tbox::Pointer< pdat::SideVariable<2,double> >
	 swgt = new pdat::SideVariable<2,double>("swgt", 1);
      int swgt_id = variable_db->registerVariableAndContext(
	 swgt, dummy, no_ghosts);
   
      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
	 tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
	 level->allocatePatchData(swgt_id);
	 for (iv = 0; iv < NVARS; iv++) {
	    level->allocatePatchData(svindx[iv]);
	 }
      }
      
      tbox::Pointer< math::HierarchyDataOpsReal<2,double> > side_ops = 
	 new math::HierarchySideDataOpsReal<2,double>(hierarchy, 0, 1);
      TBOX_ASSERT(!side_ops.isNull());

      tbox::Pointer<math::HierarchyDataOpsReal<2,double> > swgt_ops = 
	 new math::HierarchySideDataOpsReal<2,double>(hierarchy, 0 , 1);

      tbox::Pointer<hier::Patch<2> > patch;
      tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom;

      // Initialize control volume data for side-centered components
      hier::Box<2> coarse_fine = fine0+fine1; 
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
	 tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
	    tbox::Pointer< pdat::SideData<2,double> > data;
	    patch = level->getPatch(ip());
	    pgeom = patch->getPatchGeometry();
	    const double* dx = pgeom->getDx();
	    const double side_vol = dx[0]*dx[1];
	    data = patch->getPatchData(swgt_id);
	    data->fillAll(side_vol);
	    pdat::SideIndex<2> fi;
	    int plo0 = patch->getBox().lower(0);
	    int phi0 = patch->getBox().upper(0);
	    int plo1 = patch->getBox().lower(1);
	    int phi1 = patch->getBox().upper(1);
	    int ic;

	    if (ln == 0) {
	       data->fillAll(0.0, (coarse_fine * patch->getBox()) );

	       if (patch->getPatchNumber() == 0) {
		  //bottom side boundaries
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,plo1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 0.5;
		  }
		  //left and right side boundaries
		  for (ic = plo1; ic <= phi1; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(plo0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 0.5;
		     fi = pdat::SideIndex<2>(hier::Index<2>(phi0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Upper);
		     (*data)(fi) *= 0.5;
		  }
	       } else {
		  //top and bottom side boundaries
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,plo1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Lower);
		     (*data)(fi) = 0.0;
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,phi1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Upper);
		     (*data)(fi) *= 0.5;
		  }
		  //left and right side boundaries
		  for (ic = plo1; ic <= phi1; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(plo0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 0.5;
		     fi = pdat::SideIndex<2>(hier::Index<2>(phi0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Upper);
		     (*data)(fi) *= 0.5;
		  }
	       }
	    } else {
	       if (patch->getPatchNumber() == 0) {
		  // top and bottom coarse-fine side boundaries
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,plo1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 1.5;
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,phi1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Upper);
		     (*data)(fi) *= 1.5;
		  }
		  //left coarse-fine side boundaries
		  for (ic = plo1; ic <= phi1; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(plo0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 1.5;
		  }
	       } else {
		  // top and bottom coarse-fine side boundaries
		  for (ic = plo0; ic <= phi0; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,plo1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Lower);
		     (*data)(fi) *= 1.5;
		     fi = pdat::SideIndex<2>(hier::Index<2>(ic,phi1), pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Upper);
		     (*data)(fi) *= 1.5;
		  }
		  //left and right coarse-fine side boundaries
		  for (ic = plo1; ic <= phi1; ic++) {
		     fi = pdat::SideIndex<2>(hier::Index<2>(plo0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Lower);
		     (*data)(fi) = 0.0;
		     fi = pdat::SideIndex<2>(hier::Index<2>(phi0,ic), pdat::SideIndex<2>::X, pdat::SideIndex<2>::Upper);
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
     tbox::Pointer< pdat::SideData<2,double> > cvdata = patch->getPatchData(cwgt_id);

     for (pdat::SideIterator<2> c(cvdata->getBox(),1);c && vol_test_passed;c++) {
     pdat::SideIndex<2> side_index = c();

     if (ln == 0) {
     if ((coarse_fine * patch->getBox()).contains(side_index)) {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(side_index),0.0) ) {
     vol_test_passed = false;
     }
     } else {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(side_index),0.01) ) {
     vol_test_passed = false;
     }
     }
     }

     if (ln == 1) {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(side_index),0.0025) ) {
     vol_test_passed = false;
     }
     }
     }
     }
     }
     if (!vol_test_passed) {
     num_failures++;
     tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << endl;
     cwgt_ops->printData(cwgt_id, tbox::plog);
     }
*/
      // Print out control volume data and compute its integral
/*   tbox::plog << "side control volume data" << endl;
     swgt_ops->printData(swgt_id, tbox::plog);
*/

      // Test #1b: math::HierarchySideDataOpsReal<2>::sumControlVolumes()
      // Expected: norm = 1.0
      double norm = 
	 side_ops->sumControlVolumes(svindx[0], swgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(norm,1.0)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #1b: math::HierarchySideDataOpsReal<2>::sumControlVolumes()\n"
		    << "Expected value = 1.0 , Computed value = "
		    << norm << endl;
      }

      // Test #2: math::HierarchySideDataOpsReal<2>::numberOfEntries()
      // Expected: num_data_points = 209
      int num_data_points = side_ops->numberOfEntries(svindx[0]);
      if ( num_data_points != 209 ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #2: math::HierarchySideDataOpsReal<2>::numberOfEntries()\n"
		    << "Expected value = 209 , Computed value = "
		    << num_data_points << endl;
      }

      // Test #3a: math::HierarchySideDataOpsReal<2>::setToScalar()
      // Expected: v0 = 2.0
      double val0 = double(2.0);
      side_ops->setToScalar(svindx[0], val0);
      if (!doubleDataSameAsValue(svindx[0],val0,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3a: math::HierarchySideDataOpsReal<2>::setToScalar()\n"
		    << "Expected: v0 = " << val0 << endl;
	 side_ops->printData(svindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchySideDataOpsReal<2>::setToScalar()
      // Expected: v1 = (4.0)
      side_ops->setToScalar(svindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(svindx[1],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3b: math::HierarchySideDataOpsReal<2>::setToScalar()\n"
		    << "Expected: v1 = " << val1 << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #4: math::HierarchySideDataOpsReal<2>::copyData()
      // Expected:  v2 = v1 = (4.0)
      side_ops->copyData(svindx[2], svindx[1]);
      if (!doubleDataSameAsValue(svindx[2],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #4: math::HierarchySideDataOpsReal<2>::copyData()\n"
		    << "Expected: v2 = " << val1 << endl;
	 side_ops->printData(svindx[2], tbox::plog);
      }

      // Test #5: math::HierarchySideDataOpsReal<2>::swapData()
      // Expected:  v0 = (4.0), v1 = (2.0)
      side_ops->swapData(svindx[0], svindx[1]);
      if (!doubleDataSameAsValue(svindx[0],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #5a: math::HierarchySideDataOpsReal<2>::swapData()\n"
		    << "Expected: v0 = " << val1 << endl;
	 side_ops->printData(svindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(svindx[1],val0,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #5b: math::HierarchySideDataOpsReal<2>::swapData()\n"
		    << "Expected: v1 = " << val0 << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #6: math::HierarchySideDataOpsReal<2>::scale()
      // Expected:  v2 = 0.25 * v2 = (1.0)
      side_ops->scale(svindx[2], 0.25, svindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(svindx[2],val_scale,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #6: math::HierarchySideDataOpsReal<2>::scale()\n"
		    << "Expected: v2 = " << val_scale << endl;
	 side_ops->printData(svindx[2], tbox::plog);
      }

      // Test #7: math::HierarchySideDataOpsReal<2>::add()
      // Expected: v3 = v0 + v1 = (6.0) 
      side_ops->add(svindx[3], svindx[0], svindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(svindx[3],val_add,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #7: math::HierarchySideDataOpsReal<2>::add()\n"
		    << "Expected: v3 = " << val_add << endl;
	 side_ops->printData(svindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      side_ops->setToScalar(svindx[0], 0.0);

      // Test #8: math::HierarchySideDataOpsReal<2>::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      side_ops->subtract(svindx[1], svindx[3], svindx[0]);  
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(svindx[1],val_sub,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #8: math::HierarchySideDataOpsReal<2>::subtract()\n"
		    << "Expected: v1 = " << val_sub << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchySideDataOpsReal<2>::addScalar()
      // Expected:  v1 = v1 + (0.0) = (6.0) 
      side_ops->addScalar(svindx[1], svindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(svindx[1],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9a: math::HierarchySideDataOpsReal<2>::addScalar()\n"
		    << "Expected: v1 = " << val_addScalar << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchySideDataOpsReal<2>::addScalar()
      // Expected:  v2 = v2 + (0.0) = (1.0) 
      side_ops->addScalar(svindx[2], svindx[2], 0.0); 
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(svindx[2],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9b: math::HierarchySideDataOpsReal<2>::addScalar()\n"
		    << "Expected: v2 = " << val_addScalar << endl;
	 side_ops->printData(svindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchySideDataOpsReal<2>::addScalar()
      // Expected:  v2 = v2 + (3.0) = (4.0)
      side_ops->addScalar(svindx[2], svindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(svindx[2],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9c: math::HierarchySideDataOpsReal<2>::addScalar()\n"
		    << "Expected: v2 = " << val_addScalar << endl;
	 side_ops->printData(svindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5)
      side_ops->setToScalar(svindx[3], 0.5);

      // Test #10: math::HierarchySideDataOpsReal<2>::multiply()
      // Expected: v1 = v3 * v1 = (3.0) 
      side_ops->multiply(svindx[1], svindx[3], svindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(svindx[1],val_mult,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #10 math::HierarchySideDataOpsReal<2>::multiply()\n"
		    << "Expected: v1 = " << val_mult << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #11: math::HierarchySideDataOpsReal<2>::divide()
      // Expected: v0 = v2 / v1 =  1.33333333333333
      side_ops->divide(svindx[0], svindx[2], svindx[1]);
      double val_div = 1.333333333333;
      if (!doubleDataSameAsValue(svindx[0],val_div,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #11 math::HierarchySideDataOpsReal<2>::divide()\n"
		    << "Expected: v0 = " << val_div << endl;
	 side_ops->printData(svindx[0], tbox::plog);
      }

      // Test #12: math::HierarchySideDataOpsReal<2>::reciprocal()
      // Expected: v1 = 1 / v1 = (0.333333333)
      side_ops->reciprocal(svindx[1], svindx[1]); 
      double val_rec = 0.33333333333333;
      if (!doubleDataSameAsValue(svindx[1],val_rec,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #12 math::HierarchySideDataOpsReal<2>::reciprocal()\n"
		    << "Expected: v1 = " << val_rec << endl;
	 side_ops->printData(svindx[1], tbox::plog);
      }

      // Test #13: math::HierarchySideDataOpsReal<2>::abs()
      // Expected: v3 = abs(v2) = 4.0
      side_ops->abs(svindx[3], svindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(svindx[3],val_abs,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #13 math::HierarchySideDataOpsReal<2>::abs()\n"
		    << "Expected: v3 = " << val_abs << endl;
	 side_ops->printData(svindx[3], tbox::plog);
      }

   
      // Test #14: Place some bogus values on coarse level 
      tbox::Pointer< pdat::SideData<2,double> > cdata;

      // set values
      tbox::Pointer<hier::PatchLevel<2> > level_zero = hierarchy->getPatchLevel(0);
      for (hier::PatchLevel<2>::Iterator ip(level_zero); ip; ip++) {
	 patch = level_zero->getPatch(ip());
	 cdata = patch->getPatchData(svindx[2]);
	 hier::Index<2> index0(2,2);
	 hier::Index<2> index1(5,3);
	 if (patch->getBox().contains(index0)) {
	    (*cdata)(pdat::SideIndex<2>(index0, pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Lower), 0) = 100.0;
	 }
	 if (patch->getBox().contains(index1)) {
	    (*cdata)(pdat::SideIndex<2>(index1, pdat::SideIndex<2>::Y, pdat::SideIndex<2>::Upper), 0) = -1000.0;
	 }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel<2>::Iterator ipp(level_zero); ipp; ipp++) {
	 patch = level_zero->getPatch(ipp());
	 cdata = patch->getPatchData(svindx[2]);
	 pdat::SideIndex<2> index0(hier::Index<2>(2,2),pdat::SideIndex<2>::Y,pdat::SideIndex<2>::Lower);
	 pdat::SideIndex<2> index1(hier::Index<2>(5,3),pdat::SideIndex<2>::Y,pdat::SideIndex<2>::Upper);

	 // check X axis data
	 for (pdat::SideIterator<2> c(cdata->getBox(),pdat::SideIndex<2>::X);c && bogus_value_test_passed;c++) {
	    pdat::SideIndex<2> side_index = c();

	    if (!tbox::MathUtilities<double>::equalEps((*cdata)(side_index), 4.0)) {
	       bogus_value_test_passed = false;
	    }
	 }

	 // check Y axis data
	 for (pdat::SideIterator<2> cc(cdata->getBox(),pdat::SideIndex<2>::Y);cc && bogus_value_test_passed;cc++) {
            pdat::SideIndex<2> side_index = cc();

	    if (side_index == index0) {
	       if (!tbox::MathUtilities<double>::equalEps((*cdata)(side_index),100.0)) {
		  bogus_value_test_passed = false;
	       }
	    } else {
	       if (side_index == index1) {
		  if (!tbox::MathUtilities<double>::equalEps((*cdata)(side_index),-1000.0)) {
		     bogus_value_test_passed = false;
		  }
	       } else {
		  if (!tbox::MathUtilities<double>::equalEps((*cdata)(side_index), 4.0)) {
		     bogus_value_test_passed = false;
		  }
	       }
	    }
	 }
      }
      if (!bogus_value_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #14:  Place some bogus values on coarse level" << endl;
	 side_ops->printData(svindx[2], tbox::plog);
      }


      // Test #15: math::HierarchySideDataOpsReal<2>::L1Norm() - w/o control weights
      // Expected: bogus_l1_norm = 1984.00
      double bogus_l1_norm = side_ops->L1Norm(svindx[2]);
      if ( !tbox::MathUtilities<double>::equalEps(bogus_l1_norm,1984.00) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #15: math::HierarchySideDataOpsReal<2>::L1Norm()"
		    << " - w/o control weights\n"
		    << "Expected value = 1984.00, Computed value = "
		    << setprecision(12) << bogus_l1_norm << endl;
      }

      // Test #16: math::HierarchySideDataOpsReal<2>::L1Norm() - w/control weights
      // Expected: correct_l1_norm = 4.0
      double correct_l1_norm = side_ops->L1Norm(svindx[2],swgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(correct_l1_norm,4.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #16: math::HierarchySideDataOpsReal<2>::L1Norm()"
		    << " - w/control weights\n"
		    << "Expected value = 4.0, Computed value = "
		    << correct_l1_norm << endl;
      }

      // Test #17: math::HierarchySideDataOpsReal<2>::L2Norm()
      // Expected: l2_norm = 4.0
      double l2_norm = side_ops->L2Norm(svindx[2],swgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(l2_norm,4.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #17: math::HierarchySideDataOpsReal<2>::L2Norm()\n"
		    << "Expected value = 4.0, Computed value = "
		    << l2_norm << endl;
      }

      // Test #18: math::HierarchySideDataOpsReal<2>::maxNorm() - w/o control weights
      // Expected: bogus_max_norm = 1000.0
      double bogus_max_norm = side_ops->maxNorm(svindx[2]);
      if ( !tbox::MathUtilities<double>::equalEps(bogus_max_norm,1000.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #18: math::HierarchySideDataOpsReal<2>::maxNorm()"
		    << " - w/o control weights\n"
		    << "Expected value = 1000.0, Computed value = "
		    << bogus_max_norm << endl;
      }

      // Test #19: math::HierarchySideDataOpsReal<2>::maxNorm() - w/control weights
      // Expected: max_norm = 4.0
      double max_norm = side_ops->maxNorm(svindx[2],swgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(max_norm,4.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #19: math::HierarchySideDataOpsReal<2>::maxNorm()"
		    << " - w/control weights\n"
		    << "Expected value = 4.0, Computed value = "
		    << max_norm << endl;
      }

      // Reset data and test sums, axpy's
      side_ops->setToScalar(svindx[0], 1.0); 
      side_ops->setToScalar(svindx[1], 2.5); 
      side_ops->setToScalar(svindx[2], 7.0); 

      // Test #20: math::HierarchySideDataOpsReal<2>::linearSum()
      // Expected: v3 = 5.0
      side_ops->linearSum(svindx[3], 2.0, svindx[1], 0.0, svindx[0]);
      double val_linearSum = 5.0;
      if (!doubleDataSameAsValue(svindx[3],val_linearSum,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #20: math::HierarchySideDataOpsReal<2>::linearSum()\n"
		    << "Expected: v3 = " << val_linearSum << endl;
	 side_ops->printData(svindx[3], tbox::perr);
      }

      // Test #21: math::HierarchySideDataOpsReal<2>::axmy()
      // Expected: v3 = 6.5
      side_ops->axmy(svindx[3], 3.0, svindx[1], svindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(svindx[3],val_axmy,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #21: math::HierarchySideDataOpsReal<2>::axmy()\n"
		    << "Expected: v3 = " << val_axmy << endl;
	 side_ops->printData(svindx[3], tbox::plog);
      }

      // Test #22a: math::HierarchySideDataOpsReal<2>::dot() - (ind2) * (ind1)
      // Expected: cdot = 17.5
      double cdot = side_ops->dot(svindx[2], svindx[1], swgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(cdot,17.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #22a: math::HierarchySideDataOpsReal<2>::dot() - (ind2) * (ind1)\n"
		    << "Expected Value = 17.5, Computed Value = "
		    << cdot << endl;
      }

      // Test #22b: math::HierarchySideDataOpsReal<2>::dot() - (ind2) * (ind1)
      // Expected: cdot = 17.5
      cdot = side_ops->dot(svindx[1], svindx[2], swgt_id); 
      if ( !tbox::MathUtilities<double>::equalEps(cdot,17.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #22b: math::HierarchySideDataOpsReal<2>::dot() - (ind2) * (ind1)\n"
		    << "Expected Value = 17.5, Computed Value = "
		    << cdot << endl;
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
      side_ops.setNull(); 
      swgt_ops.setNull(); 

      if (num_failures == 0) {
	 tbox::pout << "\nPASSED:  side hiertest" << endl;
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
		      tbox::Pointer<hier::PatchHierarchy<2> > hierarchy)
{
   bool test_passed = true;

   int ln;
   tbox::Pointer<hier::Patch<2> > patch;
   for (ln = 0; ln < 2; ln++) {
      tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
         patch = level->getPatch(ip());
         tbox::Pointer< pdat::SideData<2,double> > cvdata = patch->getPatchData(desc_id);

         for (pdat::SideIterator<2> c(cvdata->getBox(),1); c && test_passed ; c++) {
            pdat::SideIndex<2> side_index = c();
            if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (side_index),value) ) {
               test_passed = false;
            }
         }
      }
   }

   return (test_passed);
}

