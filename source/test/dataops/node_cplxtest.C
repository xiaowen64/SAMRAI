//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dataops/node_cplxtest.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program to test node-centered complex patch data ops
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
#include "NodeData.h"
#include "HierarchyDataOpsComplex.h"
#include "HierarchyNodeDataOpsComplex.h"
#include "HierarchyDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "NodeIndex.h"
#include "NodeIterator.h"
#include "NodeVariable.h"
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
bool complexDataSameAsValue(int desc_id, dcomplex value,
			    tbox::Pointer<hier::PatchHierarchy<2> > hierarchy);
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

// tbox::PIO::logOnlyNodeZero("node_cplxtest.log");
      tbox::PIO::logAllNodes("node_cplxtest.log");

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
      tbox::Pointer< pdat::NodeVariable<2,dcomplex> > nvar[NVARS];
      int nvindx[NVARS];
      nvar[0] = new pdat::NodeVariable<2,dcomplex>("nvar0", 1);
      nvindx[0] = variable_db->registerVariableAndContext(
	 nvar[0], dummy, no_ghosts);
      nvar[1] = new pdat::NodeVariable<2,dcomplex>("nvar1", 1);
      nvindx[1] = variable_db->registerVariableAndContext(
	 nvar[1], dummy, no_ghosts);
      nvar[2] = new pdat::NodeVariable<2,dcomplex>("nvar2", 1);
      nvindx[2] = variable_db->registerVariableAndContext(
	 nvar[2], dummy, no_ghosts);
      nvar[3] = new pdat::NodeVariable<2,dcomplex>("nvar3", 1);
      nvindx[3] = variable_db->registerVariableAndContext(
	 nvar[3], dummy, no_ghosts);

      tbox::Pointer< pdat::NodeVariable<2,double> >
	 nwgt = new pdat::NodeVariable<2,double>("nwgt", 1);
      int nwgt_id = variable_db->registerVariableAndContext(
	 nwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
	 hierarchy->getPatchLevel(ln)->allocatePatchData(nwgt_id);
	 for (iv = 0; iv < NVARS; iv++) {
	    hierarchy->getPatchLevel(ln)->allocatePatchData(nvindx[iv]);
	 }
      }

      tbox::Pointer<math::HierarchyDataOpsComplex<2> > node_ops = 
	 new math::HierarchyNodeDataOpsComplex<2>(hierarchy, 0, 1);
      TBOX_ASSERT(!node_ops.isNull());

      tbox::Pointer<math::HierarchyDataOpsReal<2,double> > nwgt_ops = 
	 new math::HierarchyNodeDataOpsReal<2,double>(hierarchy, 0 , 1);

      tbox::Pointer<hier::Patch<2> > patch;
      tbox::Pointer<geom::CartesianPatchGeometry<2> > pgeom;

      // Initialize control volume data for node-centered components
      hier::Box<2> coarse_fine = fine0+fine1; 
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
	 tbox::Pointer<hier::PatchLevel<2> > level = hierarchy->getPatchLevel(ln);
	 for (hier::PatchLevel<2>::Iterator ip(level); ip; ip++) {
	    tbox::Pointer< pdat::NodeData<2,double> > data;
	    patch = level->getPatch(ip());
	    pgeom = patch->getPatchGeometry();
	    const double* dx = pgeom->getDx();
	    const double node_vol = dx[0]*dx[1];
	    data = patch->getPatchData(nwgt_id);
	    data->fillAll(node_vol);
	    pdat::NodeIndex<2> ni;
	    int plo0 = patch->getBox().lower(0);
	    int phi0 = patch->getBox().upper(0);
	    int plo1 = patch->getBox().lower(1);
	    int phi1 = patch->getBox().upper(1);
	    int ic;

	    if (ln == 0) {
	       data->fillAll(0.0, (coarse_fine * patch->getBox()) );

	       if (patch->getPatchNumber() == 0) {
		  //bottom face boundaries
		  for (ic = plo0; ic < phi0; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,plo1), pdat::NodeIndex<2>::LowerRight);
		     (*data)(ni) *= 0.5;
		  }
		  //left and right face boundaries
		  for (ic = plo1; ic <= phi1; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(plo0,ic), pdat::NodeIndex<2>::UpperLeft);
		     (*data)(ni) *= 0.5;
		     ni = pdat::NodeIndex<2>(hier::Index<2>(phi0,ic), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 0.5;
		  }
		  // corner boundaries
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(plo0,plo1), pdat::NodeIndex<2>::LowerLeft)) *= 0.25;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(phi0,plo1), pdat::NodeIndex<2>::LowerRight)) *= 0.25;
	       } else {
		  //top and bottom face boundaries
		  for (ic = plo0; ic < phi0; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,phi1), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 0.5;
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,plo1), pdat::NodeIndex<2>::LowerRight);
		     (*data)(ni) = 0.0;
		  }
		  //left and right face boundaries
		  for (ic = plo1; ic < phi1; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(plo0,ic), pdat::NodeIndex<2>::UpperLeft);
		     (*data)(ni) *= 0.5;
		     ni = pdat::NodeIndex<2>(hier::Index<2>(phi0,ic), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 0.5;
		  }
		  // corner boundaries
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(plo0,plo1), pdat::NodeIndex<2>::LowerLeft)) = 0.0;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(plo0,phi1), pdat::NodeIndex<2>::UpperLeft)) *= 0.25;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(phi0,plo1), pdat::NodeIndex<2>::LowerRight)) = 0.0;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(phi0,phi1), pdat::NodeIndex<2>::UpperRight)) *= 0.25;
	       }
	    } else {
	       if (patch->getPatchNumber() == 0) {
		  // top and bottom coarse-fine face boundaries
		  for (ic = plo0; ic <= phi0; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,plo1), pdat::NodeIndex<2>::LowerRight);
		     (*data)(ni) *= 1.5;
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,phi1), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 1.5;
		  }
		  //left coarse-fine face boundaries
		  for (ic = plo1; ic < phi1; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(plo0,ic), pdat::NodeIndex<2>::UpperLeft);
		     (*data)(ni) *= 1.5;
		  }
		  // coarse-fine corner boundaries
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(plo0,plo1), pdat::NodeIndex<2>::LowerLeft)) *= 2.25;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(plo0,phi1), pdat::NodeIndex<2>::UpperLeft)) *= 2.25;
	       } else {
		  // top and bottom coarse-fine face boundaries
		  for (ic = plo0; ic < phi0; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,plo1), pdat::NodeIndex<2>::LowerRight);
		     (*data)(ni) *= 1.5;
		     ni = pdat::NodeIndex<2>(hier::Index<2>(ic,phi1), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 1.5;
		  }
		  //right coarse-fine face boundaries
		  for (ic = plo1; ic < phi1; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(phi0,ic), pdat::NodeIndex<2>::UpperRight);
		     (*data)(ni) *= 1.5;
		  }
		  // coarse-fine corner boundaries
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(phi0,plo1), pdat::NodeIndex<2>::LowerRight)) *= 2.25;
		  (*data)(pdat::NodeIndex<2>(hier::Index<2>(phi0,phi1), pdat::NodeIndex<2>::UpperRight)) *= 2.25;
		  //shared left boundaries
		  for (ic = plo1; ic <= phi1+1; ic++) {
		     ni = pdat::NodeIndex<2>(hier::Index<2>(plo0,ic), pdat::NodeIndex<2>::LowerLeft);
		     (*data)(ni) = 0;
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
     tbox::Pointer< pdat::NodeData<2,double> > cvdata = patch->getPatchData(cwgt_id);

     for (pdat::NodeIterator<2> c(cvdata->getBox());c && vol_test_passed;c++) {
     pdat::NodeIndex<2> cell_index = c();

     if (ln == 0) {
     if ((coarse_fine * patch->getBox()).contains(cell_index)) {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0) ) {
     vol_test_passed = false;
     }
     } else {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.01) ) {
     vol_test_passed = false;
     }
     }
     }

     if (ln == 1) {
     if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),0.0025) ) {
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
/*   tbox::plog << "node control volume data" << endl;
     nwgt_ops->printData(nwgt_id, tbox::plog);
*/

      // Test #1b: math::HierarchyNodeDataOpsComplex<2>::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = 
	 node_ops->sumControlVolumes(nvindx[0], nwgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(norm,0.5)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #1: math::HierarchyNodeDataOpsComplex<2>::sumControlVolumes()\n"
		    << "Expected value = 0.5 , Computed value = "
		    << norm << endl;
      }

      // Test #2: math::HierarchyNodeDataOpsComplex<2>::numberOfEntries()
      // Expected: num_data_points = 121
      int num_data_points = node_ops->numberOfEntries(nvindx[0]);
      if ( num_data_points != 121 ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #2: math::HierarchyNodeDataOpsComplex<2>::numberOfEntries()\n"
		    << "Expected value = 121 , Computed value = "
		    << num_data_points << endl;
      }

      // Test #3a: math::HierarchyNodeDataOpsComplex<2>::setToScalar()
      // Expected: v0 = (2.0,1.5)
      dcomplex val0 = dcomplex(2.0,1.5);
      node_ops->setToScalar(nvindx[0], val0);
      if (!complexDataSameAsValue(nvindx[0],val0,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3a: math::HierarchyNodeDataOpsComplex<2>::setToScalar()\n"
		    << "Expected: v0 = " << val0 << endl;
	 node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyNodeDataOpsComplex<2>::setToScalar()
      // Expected:  v1 = (4.0, 3.0)
      dcomplex val1(4.0,3.0);
      node_ops->setToScalar(nvindx[1], dcomplex(4.0, 3.0));
      if (!complexDataSameAsValue(nvindx[1],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3b: math::HierarchyNodeDataOpsComplex<2>::setToScalar()\n"
		    << "Expected: v1 = " << val1 << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyNodeDataOpsComplex<2>::copyData()
      // Expected:   v2 = v1 = (4.0, 3.0)
      node_ops->copyData(nvindx[2], nvindx[1]);
      if (!complexDataSameAsValue(nvindx[2],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #4: math::HierarchyNodeDataOpsComplex<2>::copyData()\n"
		    << "Expected: v2 = " << val1 << endl;
	 node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyNodeDataOpsComplex<2>::swapData()
      // Expected:  v0 = (4.0, 3.0), v1 = (2.0,1.5)
      node_ops->swapData(nvindx[0], nvindx[1]);
      if (!complexDataSameAsValue(nvindx[0],val1,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #5a: math::HierarchyNodeDataOpsComplex<2>::swapData()\n"
		    << "Expected: v0 = " << val1 << endl;
	 node_ops->printData(nvindx[0], tbox::plog);
      }
      if (!complexDataSameAsValue(nvindx[1],val0,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #5b: math::HierarchyNodeDataOpsComplex<2>::swapData()\n"
		    << "Expected: v1 = " << val0 << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyNodeDataOpsComplex<2>::scale()
      // Expected:  v2 = 0.25 * v2 = (1.0,0.75)
      node_ops->scale(nvindx[2], 0.25, nvindx[2]);
      dcomplex val_scale(1.0,0.75);
      if (!complexDataSameAsValue(nvindx[2],val_scale,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #6: math::HierarchyNodeDataOpsComplex<2>::scale()\n"
		    << "Expected: v2 = " << val_scale << endl;
	 node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyNodeDataOpsComplex<2>::add()
      // Expected: v3 = v0 + v1 = (6.0, 4.5) 
      node_ops->add(nvindx[3], nvindx[0], nvindx[1]);
      dcomplex val_add(6.0,4.5);
      if (!complexDataSameAsValue(nvindx[3],val_add,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #7: math::HierarchyNodeDataOpsComplex<2>::add()\n"
		    << "Expected: v3 = " << val_add << endl;
	 node_ops->printData(nvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0,4.5)
      node_ops->setToScalar(nvindx[0], dcomplex(0.0,4.5));

      // Test #8: math::HierarchyNodeDataOpsComplex<2>::subtract()
      // Expected: v1 = v3 - v0 = (6.0,0.0)
      node_ops->subtract(nvindx[1], nvindx[3], nvindx[0]);  
      dcomplex val_sub(6.0,0.0);
      if (!complexDataSameAsValue(nvindx[1],val_sub,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #8: math::HierarchyNodeDataOpsComplex<2>::subtract()\n"
		    << "Expected: v1 = " << val_sub << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyNodeDataOpsComplex<2>::addScalar()
      // Expected: v1 = v1 + (0.0,-4.0) = (6.0,-4.0) 
      node_ops->addScalar(nvindx[1], nvindx[1], dcomplex(0.0,-4.0));
      dcomplex val_addScalar(6.0,-4.0);
      if (!complexDataSameAsValue(nvindx[1],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9a: math::HierarchyNodeDataOpsComplex<2>::addScalar()\n"
		    << "Expected: v1 = " << val_addScalar << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyNodeDataOpsComplex<2>::addScalar()
      // Expected: v2 = v2 + (0.0,0.25) = (1.0,1.0) 
      node_ops->addScalar(nvindx[2], nvindx[2], dcomplex(0.0,0.25)); 
      val_addScalar = dcomplex(1.0,1.0);
      if (!complexDataSameAsValue(nvindx[2],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9b: math::HierarchyNodeDataOpsComplex<2>::addScalar()\n"
		    << "Expected: v2 = " << val_addScalar << endl;
	 node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyNodeDataOpsComplex<2>::addScalar()
      // Expected: v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
      node_ops->addScalar(nvindx[2], nvindx[2], dcomplex(3.0,-4.0));
      val_addScalar = dcomplex(4.0,-3.0);
      if (!complexDataSameAsValue(nvindx[2],val_addScalar,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9c: math::HierarchyNodeDataOpsComplex<2>::addScalar()\n"
		    << "Expected: v2 = " << val_addScalar << endl;
	 node_ops->printData(nvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5, 0.0)
      node_ops->setToScalar(nvindx[3], dcomplex(0.5, 0.0));

      // Test #10: math::HierarchyNodeDataOpsComplex<2>::multiply()
      // Expected: v1 = v3 * v1 = (3.0,-2.0) 
      node_ops->multiply(nvindx[1], nvindx[3], nvindx[1]);
      dcomplex val_mult(3.0,-2.0);
      if (!complexDataSameAsValue(nvindx[1],val_mult,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #10: math::HierarchyNodeDataOpsComplex<2>::multiply()\n"
		    << "Expected: v1 = " << val_mult << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyNodeDataOpsComplex<2>::divide()
      // Expected: v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
      node_ops->divide(nvindx[0], nvindx[2], nvindx[1]);
      dcomplex val_div(1.3846153846154,-0.076923076923077);
      if (!complexDataSameAsValue(nvindx[0],val_div,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #11: math::HierarchyNodeDataOpsComplex<2>::divide()\n"
		    << "Expected: v0 = " << val_div << endl;
	 node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyNodeDataOpsComplex<2>::reciprocal()
      // Expected: v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
      node_ops->reciprocal(nvindx[1], nvindx[1]); 
      dcomplex val_rec(0.23076923076923, 0.15384615384615);
      if (!complexDataSameAsValue(nvindx[1],val_rec,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #12: math::HierarchyNodeDataOpsComplex<2>::reciprocal()\n"
		    << "Expected: v1 = " << val_rec << endl;
	 node_ops->printData(nvindx[1], tbox::plog);
      }
      // Test #13: Place some bogus values on coarse level 
      tbox::Pointer< pdat::NodeData<2,dcomplex> > ndata;

      // set values
      tbox::Pointer<hier::PatchLevel<2> > level_zero = hierarchy->getPatchLevel(0);
      for (hier::PatchLevel<2>::Iterator ip(level_zero); ip; ip++) {
	 patch = level_zero->getPatch(ip());
	 ndata = patch->getPatchData(nvindx[2]);
	 hier::Index<2> index0(2,2);
	 hier::Index<2> index1(5,3);
	 if (patch->getBox().contains(index0)) {
	    (*ndata)(pdat::NodeIndex<2>(index0, pdat::NodeIndex<2>::LowerLeft), 0) = dcomplex(100.0,-50.0);
	 }
	 if (patch->getBox().contains(index1)) {
	    (*ndata)(pdat::NodeIndex<2>(index1, pdat::NodeIndex<2>::UpperRight), 0) = dcomplex(-1000.0,20.0);
	 }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel<2>::Iterator ipp(level_zero); ipp; ipp++) {
	 patch = level_zero->getPatch(ipp());
	 ndata = patch->getPatchData(nvindx[2]);
	 pdat::NodeIndex<2> index0(hier::Index<2>(2,2),pdat::NodeIndex<2>::LowerLeft);
	 pdat::NodeIndex<2> index1(hier::Index<2>(5,3),pdat::NodeIndex<2>::UpperRight);

	 for (pdat::NodeIterator<2> c(ndata->getBox());c && bogus_value_test_passed;c++) {
	    pdat::NodeIndex<2> node_index = c();

	    if (node_index == index0) {
	       if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
							    dcomplex(100.0,-50.0)) ) {
		  bogus_value_test_passed = false;
	       }
	    } else {
	       if (node_index == index1) {
		  if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
							       dcomplex(-1000.0,20.0)) ) {
		     bogus_value_test_passed = false;
		  }
	       } else {
		  if (!tbox::MathUtilities<dcomplex>::equalEps((*ndata)(node_index),
							       dcomplex(4.0,-3.0)) ) {
		     bogus_value_test_passed = false;
		  }
	       }
	    }
	 }
      }
      if (!bogus_value_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #13:  Place some bogus values on coarse level" << endl;
	 node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test norms on patch data with nvindx[2] on hierarchy with bogus values

      // Test #14: math::HierarchyNodeDataOpsComplex<2>::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1787.0034
      double bogus_l1_norm = node_ops->L1Norm(nvindx[2]);
      if ( !tbox::MathUtilities<double>::equalEps(bogus_l1_norm,1787.0034) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #14: math::HierarchyCellDataOpsComplex<2>::L1Norm()"
		    << "Expected value = 1787.0034, Computed value = "
		    << setprecision(12) << bogus_l1_norm << endl;
      }

      // Test #15: math::HierarchyNodeDataOpsComplex<2>::L1Norm() - w/ control weight
      // Expected: l1_norm = 2.5
      double correct_l1_norm = node_ops->L1Norm(nvindx[2],nwgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(correct_l1_norm,2.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #15: math::HierarchyNodeDataOpsComplex<2>::L1Norm()"
		    << " - w/ control weight\n"
		    << "Expected value = 2.5, Computed value = "
		    << correct_l1_norm << endl;
      }

      // Test #16: math::HierarchyNodeDataOpsComplex<2>::L2Norm()
      // Expected: l2_norm = 3.53553390593
      double l2_norm = node_ops->L2Norm(nvindx[2],nwgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(l2_norm,3.53553390593) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #16: math::HierarchyNodeDataOpsComplex<2>::L2Norm()\n"
		    << "Expected value = 3.53553390593, Computed value = "
		    << l2_norm << endl;
      }

      // Test #17: math::HierarchyNodeDataOpsComplex<2>::maxNorm() -w/o control weight
      // Expected: bogus_max_norm = 1000.19998
      double bogus_max_norm = node_ops->maxNorm(nvindx[2]);
      if ( !tbox::MathUtilities<double>::equalEps(bogus_max_norm,1000.19998) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #17: Bogus maxNorm of v2\n"
		    << "Expected value = 1000.19998, Computed value = "
		    << bogus_max_norm << endl;
      }

      // Test #18: math::HierarchyNodeDataOpsComplex<2>::maxNorm() -w/control weight
      // Expected: max_norm = 5.0
      double max_norm = node_ops->maxNorm(nvindx[2],nwgt_id);
      if ( !tbox::MathUtilities<double>::equalEps(max_norm,5.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #18: maxNorm of v2\n"
		    << "Expected value = 5.0, Computed value = "
		    << max_norm << endl;
      }

      // Reset data and test sums, axpy's
      node_ops->setToScalar(nvindx[0], dcomplex(1.0, -3.0)); 
      node_ops->setToScalar(nvindx[1], dcomplex(2.5, 3.0)); 
      node_ops->setToScalar(nvindx[2], dcomplex(7.0, 0.0)); 

      // Test #19: math::HierarchyNodeDataOpsComplex<2>::linearSum()
      // Expected:  v3 = (2.0,5.0)
      node_ops->linearSum(nvindx[3], dcomplex(2.0,0.0), nvindx[1], dcomplex(0.0, -1.0), nvindx[0]);
      dcomplex val_linearSum(2.0,5.0);
      if (!complexDataSameAsValue(nvindx[3],val_linearSum,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #19: math::HierarchyNodeDataOpsComplex<2>::linearSum()\n"
		    << "Expected: v3 = " << val_linearSum << endl;
	 node_ops->printData(nvindx[3], tbox::plog);
      }
  
      // Test #20: math::HierarchyNodeDataOpsComplex<2>::axmy()
      // Expected:  v3 = (6.5,12.0)
      node_ops->axmy(nvindx[3], 3.0, nvindx[1], nvindx[0]);
      dcomplex val_axmy(6.5,12.0);
      if (!complexDataSameAsValue(nvindx[3],val_axmy,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #20: math::HierarchyNodeDataOpsComplex<2>::axmy()\n"
		    << "Expected: v3 = " << val_axmy << endl;
	 node_ops->printData(nvindx[3], tbox::plog);
      }
  
      // Test #21a: math::HierarchyNodeDataOpsComplex<2>::dot()
      // Expected:  cdot = (8.75,-10.5)
      dcomplex cdot = node_ops->dot(nvindx[2], nvindx[1], nwgt_id);
      dcomplex ans_2_dot_1(8.75,-10.5);
      if ( !tbox::MathUtilities<dcomplex>::equalEps(cdot,ans_2_dot_1) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #21a: math::HierarchyNodeDataOpsComplex<2>::dot()\n"
		    << "Expected value = (8.75,-10.5), Computed value = "
		    << cdot << endl;
      }


      // Test #21b: math::HierarchyNodeDataOpsComplex<2>::dot()
      // Expected:  cdot = (8.75,10.5)
      dcomplex cdot2 = node_ops->dot(nvindx[1], nvindx[2], nwgt_id); 
      dcomplex ans_1_dot_2(8.75,10.5);
      if ( !tbox::MathUtilities<dcomplex>::equalEps(cdot2,ans_1_dot_2) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #21b: math::HierarchyNodeDataOpsComplex<2>::dot()\n"
		    << "Expected value = (8.75,10.5), Computed value = "
		    << cdot2 << endl;
      }

      // Test #22: math::HierarchyNodeDataOpsComplex<NDIM>::abs()
      // Expected: abs(v0) = 5.0
      node_ops->setToScalar(nvindx[0], dcomplex(4.0, -3.0));
      node_ops->abs(nwgt_id, nvindx[0]);
      if (!doubleDataSameAsValue(nwgt_id,5.0,hierarchy)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #22: math::HierarchyNodeDataOpsComplex<2>::abs()\n"
		    << "Expected: abs(v0) = 5.0" << endl;
	 nwgt_ops->printData(nwgt_id, tbox::plog);
      }


      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
	 hierarchy->getPatchLevel(ln)->deallocatePatchData(nwgt_id);
	 for (iv = 0; iv < NVARS; iv++) {
	    hierarchy->getPatchLevel(ln)->deallocatePatchData(nvindx[iv]);
	 }
      }

      for (iv = 0; iv < NVARS; iv++) {
	 nvar[iv].setNull();
      }
      nwgt.setNull();

      geometry.setNull();
      hierarchy.setNull();
      node_ops.setNull(); 
      nwgt_ops.setNull(); 

      if (num_failures == 0) {
	 tbox::pout << "\nPASSED:  node cplxtest" << endl;
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
         tbox::Pointer< pdat::NodeData<2,dcomplex> > nvdata = patch->getPatchData(desc_id);

         for (pdat::NodeIterator<2> c(nvdata->getBox()); c && test_passed ; c++) {
            pdat::NodeIndex<2> node_index = c();
            if ( !tbox::MathUtilities<dcomplex>::equalEps((*nvdata) (node_index),value) ) {
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
         tbox::Pointer< pdat::NodeData<2,double> > nvdata = patch->getPatchData(desc_id);

         for (pdat::NodeIterator<2> c(nvdata->getBox()); c && test_passed ; c++) {
            pdat::NodeIndex<2> node_index = c();
            if ( !tbox::MathUtilities<double>::equalEps((*nvdata) (node_index),value) ) {
               test_passed = false;
            }
         }
      }
   }

   return (test_passed);
}
