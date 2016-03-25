//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/dataops/cell_patchtest.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Main program to test cell patch data operations.
//

#include "SAMRAI_config.h"

#include <typeinfo>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
using namespace std;

#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"

#include "tbox/SAMRAIManager.h"

#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "CellVariable.h"
#include "ComponentSelector.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchDataFactory.h"
#include "tbox/Pointer.h"
#include "ProcessorMapping.h"
#include <string>

#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"
#include "VariableDatabase.h"
#include "VariableContext.h"

#include "PatchCellDataOpsReal.h"

using namespace std;
using namespace SAMRAI;

/* Helper function declarations */
bool doubleDataSameAsValue(int desc_id, double value,
			   tbox::Pointer<hier::Patch<NDIM> > patch);

int main( int argc, char *argv[] ) {

   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("cell_patchtest.log");

      /* Make a dummy mesh domain with one patch */
      double lo[NDIM];
      double hi[NDIM];
      for (int i = 0; i < NDIM; i++) {
	 lo[i] = 0.0;
	 if (i == 1) {
	    hi[i] = 0.5;
	 } else {
	    hi[i] = 1.0;
	 }
      }
      hier::Index<NDIM> indxlo(0);
      hier::Index<NDIM> indxhi(9);
      indxhi(1) = 4;
      hier::Box<NDIM> patch_box(indxlo, indxhi); 
      hier::BoxArray<NDIM> grid_domain(1);
      grid_domain(0) = patch_box;
      hier::IntVector<NDIM> ratio(1);

      geom::CartesianGridGeometry<NDIM> geometry("CartesianGeometry",
						 lo, hi, grid_domain);
      hier::ComponentSelector patch_components;

      tbox::Pointer<hier::Patch<NDIM> > tpatch =
	 new hier::Patch<NDIM>(
	    patch_box, 
	    hier::VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor());

      /* Make a variety of data on the patch. */

      /* Make three contexts for patch */
      tbox::Pointer<hier::VariableContext> ghost_width_1_context = 
	 hier::VariableDatabase<NDIM>::getDatabase()->getContext("ghost_width_1");
      tbox::Pointer<hier::VariableContext> ghost_width_2_context = 
	 hier::VariableDatabase<NDIM>::getDatabase()->getContext("ghost_width_2");
      tbox::Pointer<hier::VariableContext> ghost_width_3_context = 
	 hier::VariableDatabase<NDIM>::getDatabase()->getContext("ghost_width_3");

      /* Make ghost cell IntVectors which are used when variables
       * and contexts are registered 
       */
      hier::IntVector<NDIM> nghosts_1(1);
      hier::IntVector<NDIM> nghosts_2(2);
      hier::IntVector<NDIM> nghosts_3(3);

      /* Make cell-centered double variable for patch */
      tbox::Pointer< pdat::CellVariable<NDIM,double> > cell_double_variable = 
	 new pdat::CellVariable<NDIM,double>("cell_double_variable", 1);

      int cdvindx[3];

      /* 
       *Register cell-centered double variable and 3 contexts with
       * hier::VariableDatabase.
       */
      cdvindx[0] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_double_variable, ghost_width_1_context,nghosts_1);
      patch_components.setFlag(cdvindx[0]);

      cdvindx[1] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_double_variable, ghost_width_2_context,nghosts_2);
      patch_components.setFlag(cdvindx[1]);
  
      cdvindx[2] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_double_variable, ghost_width_3_context,nghosts_3);
      patch_components.setFlag(cdvindx[2]);


      /* Make control volume for cell-centered patch variables */
      tbox::Pointer<hier::VariableContext> ghost_width_0_context = 
	 hier::VariableDatabase<NDIM>::getDatabase()->getContext("ghost_width_0");
      hier::IntVector<NDIM> nghosts_0(0);
      tbox::Pointer< pdat::CellVariable<NDIM,double> > cwgt =
	 new pdat::CellVariable<NDIM,double>("cwgt", 1);
      int cwgt_id =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cwgt, ghost_width_0_context,nghosts_0);
      patch_components.setFlag(cwgt_id);

      int ccvindx[3];

      ccvindx[0] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_double_variable, ghost_width_1_context,nghosts_1);
      patch_components.setFlag(ccvindx[0]);

      ccvindx[1] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_double_variable, ghost_width_2_context,nghosts_2);
      patch_components.setFlag(ccvindx[1]);

      // Make two cell-centered int variables for the patch 
      tbox::Pointer< pdat::CellVariable<NDIM,int> > cell_int_variable =
	 new pdat::CellVariable<NDIM,int>("cell_int_variable", 1);

      int civindx[3];

      civindx[0] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_int_variable, ghost_width_1_context,nghosts_1);
      patch_components.setFlag(civindx[0]);

      civindx[1] =
	 hier::VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
	    cell_int_variable, ghost_width_2_context,nghosts_2);
      patch_components.setFlag(civindx[1]);

      // Test #1: Check the state of hier::PatchDescriptor<NDIM>
      int desc_id;
      string var_ctxt_name;
      bool descriptor_test_passed = true;
      bool name_error_indx[6];
      bool factory_error_indx[6];
      for (desc_id = 0; desc_id < 6; desc_id++) {
	 name_error_indx[desc_id] = false;
	 factory_error_indx[desc_id] = false;
      }
      //make strings to be used for comparison during tests

      string cell_double_variable1("cell_double_variable##ghost_width_1");
      string cell_double_variable2("cell_double_variable##ghost_width_2");
      string cell_double_variable3("cell_double_variable##ghost_width_3");
      string cwgt_variable0("cwgt##ghost_width_0");
      string cell_int_variable1("cell_int_variable##ghost_width_1");
      string cell_int_variable2("cell_int_variable##ghost_width_2");

      for (desc_id = 0; desc_id < 6; desc_id++) {
	 var_ctxt_name =
	    hier::VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->
            mapIndexToName(desc_id);

	 // test flag is used to overcome a compiler bug in GCC which 
	 // complained if the the comparison was done inside the if;
	 // must be too complicated for it.
	 bool test;

	 switch (desc_id) {
	    case 0:
	       if (var_ctxt_name != cell_double_variable1) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }

	       test = typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			     ->getPatchDescriptor()
			     ->getPatchDataFactory(desc_id)) == 
		  typeid(pdat::CellDataFactory<NDIM, double >);

	       if( !test ) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;

	    case 1:
	       if (var_ctxt_name != cell_double_variable2) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }

	       test = typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			     ->getPatchDescriptor()
			     ->getPatchDataFactory(desc_id)) == 
		  typeid(pdat::CellDataFactory<NDIM, double >);
	       if (!test) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;
	    
	    case 2:
	       if (var_ctxt_name != cell_double_variable3) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }
	       test = typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			     ->getPatchDescriptor()
			     ->getPatchDataFactory(desc_id)) == 
		  typeid(pdat::CellDataFactory<NDIM, double >);
      
	       if (!test) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;
	    
	    case 3: 
	       if (var_ctxt_name != cwgt_variable0) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }
	       test =  typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			      ->getPatchDescriptor()
			      ->getPatchDataFactory(desc_id)) == 
		  typeid(pdat::CellDataFactory<NDIM, double >);
	       if(!test) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;

	    case 4:
	       if (var_ctxt_name != cell_int_variable1) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }
	       test = typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			     ->getPatchDescriptor()
			     ->getPatchDataFactory(desc_id)) ==  
		  typeid(pdat::CellDataFactory<NDIM, int >);
	       if(!test) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;

	    case 5:
	       if (var_ctxt_name != cell_int_variable2) {
		  descriptor_test_passed = false;
		  name_error_indx[desc_id] = true;
	       }
	       test = typeid(*hier::VariableDatabase<NDIM>::getDatabase()
			     ->getPatchDescriptor()
			     ->getPatchDataFactory(desc_id)) == 
		  typeid(pdat::CellDataFactory<NDIM, int >);
	       if(!test) {
		  descriptor_test_passed = false;
		  factory_error_indx[desc_id] = true;
	       }
	       break;
	 };
      }

      if (!descriptor_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #1: State of PatchDescriptor" << endl;

	 for (desc_id = 0; desc_id < 6; desc_id++) {
	    if (name_error_indx[desc_id] == true) {
	       tbox::plog << "Name for index = " << desc_id << " incorrect" << endl;
	    }
	    if (factory_error_indx[desc_id] == true) {
	       tbox::plog << "Factory for index = " << desc_id << " incorrect" << endl;
	    }
	 }
      }

      // Test #2: Check state of hier::Patch<NDIM> before allocating storage
      if (tpatch->getBox() != patch_box) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #2a: hier::Patch<NDIM> box incorrectly set\n"
		    << "Expected: d_box = " << patch_box  << "\n"
		    << "Set to: d_box = " << tpatch->getBox() << endl;
      }
      if (tpatch->getPatchNumber() != -1) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #2b: hier::Patch<NDIM> number incorrectly set\n"
		    << "Expected: d_patch_number = -1\n"
		    << "Set to: d_patch_number = " << tpatch->getPatchNumber() << endl;
      }
/*   HOW TO GET NUMBER OF COMPONENTS ON PATCH 
     FOR NOW JUST CHECK THAT getPatchData(0) returns NULL

     num_failures++;
     tbox::perr << "FAILED: - Test #2c: Number of components allocated incorrect\n"
     << "Expected: number of components = 0\n"
     << "Got: number of components = 0\n"
*/
      for (desc_id = 0; desc_id < 6; desc_id++) {
	 if (tpatch->checkAllocated(desc_id)) {
	    num_failures++;
	    tbox::perr << "FAILED: - Test #2c." << desc_id 
		       << ": Descriptor slot " << desc_id 
		       << " should not be allocated but is!" << endl;
	 }
      }

      // Allocate all data on patch 
      tpatch->allocatePatchData(patch_components);


      // Test #3: Check state of hier::Patch<NDIM> after allocating storage
      if (tpatch->getBox() != patch_box) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3a: hier::Patch<NDIM> box incorrectly set\n"
		    << "Expected: d_box = " << patch_box  << "\n"
		    << "Set to: d_box = " << tpatch->getBox() << endl;
      }
      if (tpatch->getPatchNumber() != -1) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #3b: hier::Patch<NDIM> number incorrectly set\n"
		    << "Expected: d_patch_number = -1\n"
		    << "Set to: d_patch_number = " << tpatch->getPatchNumber() << endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; desc_id++) {

	 if (!tpatch->checkAllocated(desc_id)) {
	    num_failures++;
	    tbox::perr << "FAILED: - Test #3c.0: Descriptor index " << desc_id 
		       << " should be allocated but isn't!" << endl;
	 } else {

	    string patch_data_name = 
	       typeid(*tpatch->getPatchData(desc_id)).name();

	    hier::IntVector<NDIM> ghost_width = 
	       tpatch->getPatchData(desc_id)->getGhostCellWidth();
   
	    switch (desc_id) {
	       case 0:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++; 
		     tbox::perr << "FAILED: - Test #3c.0.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(1)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.0.b: Ghost width incorrect\n"
				<< "Expected: (1,1)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
      
	       case 1:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.1.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(2)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.1.b: Ghost width incorrect\n"
				<< "Expected: (2,2)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 2:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.2.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(3)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.2.b: Ghost width incorrect\n"
				<< "Expected: (3,3)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 3:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.3.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(0)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.3.b: Ghost width incorrect\n"
				<< "Expected: (0,0)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 4:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, int >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.4.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, int >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(1)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.4.b: Ghost width incorrect\n"
				<< "Expected: (1,1)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 5:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, int >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.5.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, int >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(2)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #3c.5.b: Ghost width incorrect\n"
				<< "Expected: (2,2)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	    };
	 }
      }

      // Initialize control volume data for cell-centered data
      const double* dx = geometry.getDx();
      double cell_vol = dx[0];
      for (int i = 1; i < NDIM; i++) {
	 cell_vol *= dx[i];
      }

      tbox::Pointer< pdat::CellData<NDIM,double> > weight = tpatch->getPatchData(cwgt_id);
      weight->fillAll(cell_vol);

      // Simple tests of cell data operations 

      math::PatchCellDataOpsReal<NDIM,double> cdops_double;

      // Get pointers to patch data objects
      tbox::Pointer< pdat::CellData<NDIM,double> > cddata0 = tpatch->getPatchData(cdvindx[0]);
      tbox::Pointer< pdat::CellData<NDIM,double> > cddata1 = tpatch->getPatchData(cdvindx[1]);
      tbox::Pointer< pdat::CellData<NDIM,double> > cddata2 = tpatch->getPatchData(cdvindx[2]);

      tbox::Pointer< pdat::CellData<NDIM,int> > cidata0 = tpatch->getPatchData(civindx[0]);
      tbox::Pointer< pdat::CellData<NDIM,int> > cidata1 = tpatch->getPatchData(civindx[1]);

      // Test #4a: math::PatchCellDataOpsReal<NDIM>::setToScalar()
      // Expected: cddata0 = 0.0
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());
      double val0 = 0.0;
      if (!doubleDataSameAsValue(cdvindx[0],val0,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #4a: math::PatchCellDataOpsReal<NDIM>::setToScalar()\n"
		    << "Expected: cddata0 = " << val0 << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #4b: math::PatchCellDataOpsReal<NDIM>::setToScalar()
      // Expected: cddata1 = 1.0
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
      double val1 = 1.0;
      if (!doubleDataSameAsValue(cdvindx[1],val1,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #4b: math::PatchCellDataOpsReal<NDIM>::setToScalar()\n"
		    << "Expected: cddata1 = " << val1 << endl;
	 cdops_double.printData(cddata1, cddata1->getGhostBox(), tbox::plog);
      }

      // Test #4c: math::PatchCellDataOpsReal<NDIM>::setToScalar()
      // Expected: cddata2 = 2.0
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());
      double val2 = 2.0;
      if (!doubleDataSameAsValue(cdvindx[2],val2,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #4b: math::PatchCellDataOpsReal<NDIM>::setToScalar()\n"
		    << "Expected: cddata2 = " << val2 << endl;
	 cdops_double.printData(cddata2, cddata2->getGhostBox(), tbox::plog);
      }

      // Test #5: math::PatchCellDataOpsReal<NDIM>::add()
      // Expected: cddata0 =  cddata1 + cddata2
      cdops_double.add(cddata0, cddata1, cddata2, tpatch->getBox());
      double val_add = 3.0;
      if (!doubleDataSameAsValue(cdvindx[0],val_add,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #5: math::PatchCellDataOpsReal<NDIM>::add()\n"
		    << "Expected: cddata0 = " << val_add << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Test #6: math::PatchCellDataOpsReal<NDIM>::subtract() on [(3,1),(5,2)]
      // Expected: cddata0 = cddata0 - cddata2
      hier::Index<NDIM> indx0(1);
      hier::Index<NDIM> indx1(2);
      indx0(0) = 3;
      indx1(0) = 5;
      cdops_double.subtract(cddata0, cddata0, cddata2, 
			    hier::Box<NDIM>(indx0, indx1));
      bool subtract_inbox_test_passed = true;
      hier::Box<NDIM> inbox(indx0, indx1);
      double val_inbox = 1.0;
      double val_not_inbox = 3.0;
      tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = tpatch->getPatchData(cdvindx[0]);

      for (pdat::CellIterator<NDIM> c(cvdata->getBox()); c && 
	      subtract_inbox_test_passed ; c++) {
	 pdat::CellIndex<NDIM> cell_index = c();
    
	 double value; 
	 if ( inbox.contains(cell_index) ) { 
	    value = val_inbox;
	 } else {
	    value = val_not_inbox;
	 }

	 if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
	    subtract_inbox_test_passed = false;
	 }
      }

      if (!subtract_inbox_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #6: math::PatchCellDataOpsReal<NDIM>::subtract() on [(3,1),(5,2)]\n"
		    << "Expected: cddata0 = 1.0 in [(3,1),(5,2)]\n"
		    << "          cddata0 = 3.0 outside box\n"  << endl;
	 cdops_double.printData(cddata0,  tpatch->getBox(), tbox::plog);
      }
 
      // Test #7: math::PatchCellDataOpsReal<NDIM>::scale()
      // Expected: cddata0 = 0.4 * cddata2
      cdops_double.scale(cddata0, 0.4, cddata2, tpatch->getBox());
      double val_scale = 0.8;
      if (!doubleDataSameAsValue(cdvindx[0],val_scale,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #7: math::PatchCellDataOpsReal<NDIM>::scale()\n"
		    << "Expected: cddata0 = " << val_scale << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }
 
      // Test #8: math::PatchCellDataOpsReal<NDIM>::multiply()
      // Expected: cddata0 = cddata0 * cddata2
      cdops_double.multiply(cddata0, cddata0, cddata2, tpatch->getBox());
      double val_mult = 1.6;
      if (!doubleDataSameAsValue(cdvindx[0],val_mult,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #8: math::PatchCellDataOpsReal<NDIM>::multiply()\n"
		    << "Expected: cddata0 = " << val_mult << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }
 
      // Test #9: math::PatchCellDataOpsReal<NDIM>::divide() in box [(3,1),(5,2)]
      // Expected: cddata0 = cddata0/cddata2
      cdops_double.divide(cddata0, cddata0, cddata2,
			  hier::Box<NDIM>(indx0, indx1));
      bool divide_inbox_test_passed = true;
      val_inbox = 0.8;
      val_not_inbox = 1.6;
      cvdata = tpatch->getPatchData(cdvindx[0]);

      for (pdat::CellIterator<NDIM> cc(cvdata->getBox()); cc && 
	      divide_inbox_test_passed ; cc++) {
	 pdat::CellIndex<NDIM> cell_index = cc();
    
	 double value; 
	 if ( inbox.contains(cell_index) ) { 
	    value = val_inbox;
	 } else {
	    value = val_not_inbox;
	 }

	 if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
	    divide_inbox_test_passed = false;
	 }
      }

      if (!divide_inbox_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #9: math::PatchCellDataOpsReal<NDIM>::divide() on [(3,1),(5,2)]\n"
		    << "Expected: cddata0 = 1.0 in [(3,1),(5,2)]\n"
		    << "          cddata0 = 3.0 outside box\n" << endl;
	 cdops_double.printData(cddata0,  tpatch->getBox(), tbox::plog);
      }
 
      // Test #10: math::PatchCellDataOpsReal<NDIM>::reciprocal()
      // Expected: cddata0 = 1/cddata2
      cdops_double.reciprocal(cddata0, cddata2, tpatch->getBox());
      double val_rec = 0.5;
      if (!doubleDataSameAsValue(cdvindx[0],val_rec,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #10: math::PatchCellDataOpsReal<NDIM>::reciprocal()\n"
		    << "Expected: cddata0 = " << val_rec << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      // Reset cddata1 and cddata2 
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());
 
      // Test #11: math::PatchCellDataOpsReal<NDIM>::linearSum()
      // Expected: cddata0 = 10*cddata1 + 20*cddata2
      cdops_double.linearSum(cddata0,10.0,cddata1,20.0,cddata2, tpatch->getBox());
      double val_linearSum = 50.0;
      if (!doubleDataSameAsValue(cdvindx[0],val_linearSum,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #11: math::PatchCellDataOpsReal<NDIM>::linearSum()\n"
		    << "Expected: cddata0 = " << val_linearSum << endl;
	 cdops_double.printData(cddata0, cddata0->getGhostBox(), tbox::plog);
      }

      cdops_double.setRandomValues(cddata0,1.0, 0.001, hier::Box<NDIM>(indx0, indx1));
      tbox::plog << "\ncddata0 = random " << endl;
      cdops_double.printData(cddata0, hier::Box<NDIM>(indx0, indx1) , tbox::plog );
 
      cdops_double.setRandomValues(cddata0,1.0, 0.001, hier::Box<NDIM>(indx0, indx1));
      tbox::plog << "\ncddata0 = random " << endl;
      cdops_double.printData(cddata0, hier::Box<NDIM>(indx0, indx1) , tbox::plog);

      // Reset cddata0 to 0.0 
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());

      // Test #12: math::PatchCellDataOpsReal<NDIM>::linearSum() on box = [(3,1),(5,2)]
      // Expected: cddata0 = 10*cddata1 + 20*cddata2 on [(3,1),(5,2)]
      cdops_double.linearSum(cddata0,10.0,cddata1,20.0,cddata2,
			     hier::Box<NDIM>(indx0,indx1));
      bool restricted_linSum_test_passed = true;
      val_inbox = 50.0;
      val_not_inbox = 0.0;
      cvdata = tpatch->getPatchData(cdvindx[0]);

      for (pdat::CellIterator<NDIM> cci(cvdata->getBox()); cci && 
	      restricted_linSum_test_passed ; cci++) {
	 pdat::CellIndex<NDIM> cell_index = cci();
    
	 double value; 
	 if ( inbox.contains(cell_index) ) { 
	    value = val_inbox;
	 } else {
	    value = val_not_inbox;
	 }

	 if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
	    restricted_linSum_test_passed = false;
	 }
      }

      if (!restricted_linSum_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #12: math::PatchCellDataOpsReal<NDIM>::linearSum()\n"
		    << "Expected: cddata0 = " << val_linearSum 
		    << " on box = [(3,1),(5,2)]" << endl;
	 cdops_double.printData(cddata0,  tpatch->getBox(), tbox::plog);
      }
 
// set individual data points and check min/max routines
      hier::Index<NDIM> newindx0(indx0);
      hier::Index<NDIM> newindx1(indx0);
      newindx0(1) = 2;
      cdops_double.setToScalar(cddata1, 0.0003, hier::Box<NDIM>(indx0, newindx0));
      newindx0(0) = 1;
      cdops_double.setToScalar(cddata1, 12345.0, hier::Box<NDIM>(newindx0, newindx0));
      newindx0(0) = 5;
      newindx0(1) = 3;
      newindx1(0) = 5;
      newindx1(1) = 4; 
      cdops_double.setToScalar(cddata1, 21.0, hier::Box<NDIM>(newindx0, newindx1));

      // Test #13: math::PatchCellDataOpsReal<NDIM>::setToScalar() on box
      // Expected: cddata1 = 0.0003 in [(3,1),(3,2)]
      //           cddata1 = 12345.0 in [(1,2),(1,2)]
      //           cddata1 = 21.0 in [(5,3),(5,4)]
      //           cddata1 = 1.0 everywhere else
      bool setToScalar_onBox_test_passed = true;
      newindx0 = indx0;
      newindx0(1) = 2;
      hier::Box<NDIM> box1(indx0,newindx0);
      newindx0(0) = 1;
      hier::Box<NDIM> box2(newindx0,newindx0);
      newindx0(0) = 5;
      newindx0(1) = 3;
      newindx1(0) = 5;
      newindx1(1) = 4;
      hier::Box<NDIM> box3(newindx0,newindx1);
      double val_inbox1 = 0.0003;
      double val_inbox2 = 12345.0;
      double val_inbox3 = 21.0;
      val_not_inbox = 1.0;

      cvdata = tpatch->getPatchData(cdvindx[1]);
      for (pdat::CellIterator<NDIM> ci(cvdata->getBox()); ci && 
	      setToScalar_onBox_test_passed ; ci++) {
	 pdat::CellIndex<NDIM> cell_index = ci();
    
	 double value; 
	 if ( box1.contains(cell_index) ) { 
	    value = val_inbox1;
	 } else {
	    if ( box2.contains(cell_index) ) { 
	       value = val_inbox2;
	    } else {
	       if ( box3.contains(cell_index) ) { 
		  value = val_inbox3;
	       } else {
		  value = val_not_inbox;
	       }
	    }
	 }

	 if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
	    setToScalar_onBox_test_passed = false;
	 }
      }

      if (!setToScalar_onBox_test_passed) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #13: math::PatchCellDataOpsReal<NDIM>::setToScalar() on box\n"
		    << "Expected: cddata1 = 0.0003 in [(3,1),(3,2)]\n"
		    << "          cddata1 = 12345.0 in [(1,2),(1,2)]\n"
		    << "          cddata1 = 21.0 in [(5,3),(5,4)]\n"
		    << "          cddata1 = 1.0 everywhere else\n" << endl;
	 cdops_double.printData(cddata1,  tpatch->getBox(), tbox::plog);
      }
 
      // Test #14: math::PatchCellDataOpsReal<NDIM>::max() on box [(3,1),(7,4)]
      // Expected: lmax = 21.0 
      hier::Index<NDIM> indx2(4);
      indx2(0) = 7;
      double lmax = cdops_double.max(cddata1, hier::Box<NDIM>(indx0,indx2));
      if ( !tbox::MathUtilities<double>::equalEps(lmax,21.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #14: math::PatchCellDataOpsReal<NDIM>::max() on box [(3,1),(7,4)]\n"
		    << "Expected value = 21.0, Computed value = "
		    << lmax << endl;
      }
 
      // Test #15: math::PatchCellDataOpsReal<NDIM>::max() in box [(0,0),(9,4)]
      // Expected: lmax = 12345.0
      lmax = cdops_double.max(cddata1, tpatch->getBox());
      if ( !tbox::MathUtilities<double>::equalEps(lmax,12345.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #15: math::PatchCellDataOpsReal<NDIM>::max() in box [(0,0),(9,4)]\n"
		    << "Expected value = 12345.0, Computed value = "
		    << lmax << endl;
      }
 
// check axpy, axmy routines

      // set cddata0, cddata1, cddata2 
      cdops_double.setToScalar(cddata0, 0.0, cddata0->getGhostBox());
      cdops_double.setToScalar(cddata1, 1.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 2.0, cddata2->getGhostBox());

      // Test #16: math::PatchCellDataOpsReal<NDIM>::axpy() 
      // Expected: cddata0 = 0.5 * 1.0 + 2.0 = 2.5
      cdops_double.axpy(cddata0, 0.5, cddata1, cddata2, tpatch->getBox());
      double val_axpy = 2.5;
      if (!doubleDataSameAsValue(cdvindx[0],val_axpy,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #16: math::PatchCellDataOpsReal<NDIM>::axpy()\n"
		    << "Expected: cddata0 = " << val_axpy << endl;
	 cdops_double.printData(cddata0, cddata0->getBox(), tbox::plog);
      }
 
      // Test #17: math::PatchCellDataOpsReal<NDIM>::axmy() 
      // Expected: cddata0 = 1.5 * 2.0 - 1.0 = 2.0
      cdops_double.axmy(cddata0, 1.5, cddata2, cddata1, tpatch->getBox());
      double val_axmy = 2.0;
      if (!doubleDataSameAsValue(cdvindx[0],val_axmy,tpatch)) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #17: math::PatchCellDataOpsReal<NDIM>::axmy()\n"
		    << "Expected: cddata0 = " << val_axmy << endl;
	 cdops_double.printData(cddata0, cddata0->getBox(), tbox::plog);
      }
 
// Test the norm ops stuff

      // Test #18a: math::PatchCellDataOpsReal<NDIM>::sumControlVolumes() for cddata1
      // Expected: lsum = 0.5
      double lsum = cdops_double.sumControlVolumes(cddata1, weight, tpatch->getBox());
      if ( !tbox::MathUtilities<double>::equalEps(lsum, 0.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #18a: math::PatchCellDataOpsReal<NDIM>::sumControlVolumes() for cddata1\n"
		    << "Expected value = 0.5, Computed value = "
		    << lsum << endl;
      }
 
      // Test #18b: math::PatchCellDataOpsReal<NDIM>::sumControlVolumes() for cddata2
      // Expected: lsum = 0.5
      lsum = cdops_double.sumControlVolumes(cddata2, weight, tpatch->getBox());
      if ( !tbox::MathUtilities<double>::equalEps(lsum, 0.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #18b: math::PatchCellDataOpsReal<NDIM>::sumControlVolumes() for cddata2\n"
		    << "Expected value = 0.5, Computed value = "
		    << lsum << endl;
      }
 
      // Test #19a: math::PatchCellDataOpsReal<NDIM>::L1norm() for cddata1
      // Expected: l1norm = 0.5
      double l1norm = cdops_double.L1Norm(cddata1, tpatch->getBox(), weight);
      if ( !tbox::MathUtilities<double>::equalEps(l1norm, 0.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #19a: math::PatchCellDataOpsReal<NDIM>::L1norm() for cddata1\n"
		    << "Expected value = 0.5, Computed value = "
		    << l1norm << endl;
      }
 
      // Test #19b: math::PatchCellDataOpsReal<NDIM>::L1norm() for cddata2
      // Expected: l1norm = 1.0
      l1norm = cdops_double.L1Norm(cddata2, tpatch->getBox(), weight);
      if ( !tbox::MathUtilities<double>::equalEps(l1norm, 1.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #19b: math::PatchCellDataOpsReal<NDIM>::L1norm() for cddata2\n"
		    << "Expected value = 1.0, Computed value = "
		    << l1norm << endl;
      }

      // Test #20: math::PatchCellDataOpsReal<NDIM>::L2norm() for cddata2
      // Expected: l2norm = sqrt(2) = 1.4142135623731
      double l2norm = cdops_double.L2Norm(cddata2, tpatch->getBox(), weight);
      if ( !tbox::MathUtilities<double>::equalEps(l2norm, 1.4142135623731) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #20: math::PatchCellDataOpsReal<NDIM>::L2norm() for cddata2\n"
		    << "Expected value = sqrt(2) = 1.4142135623731, Computed value = "
		    << setprecision(12) << l2norm << endl;
      }

      // Reset cddata1 to 0.5
      cdops_double.setToScalar(cddata1, 0.5, tpatch->getBox());

      // Test #21: math::PatchCellDataOpsReal<NDIM>::weightedL2norm() for cddata2
      // Expected: wl2norm = sqrt(0.5) = 0.70710678118655
      double wl2norm = cdops_double.weightedL2Norm(cddata2,cddata1,tpatch->getBox(),weight);
      if ( !tbox::MathUtilities<double>::equalEps(wl2norm, 0.70710678118655) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #21: math::PatchCellDataOpsReal<NDIM>::weightedL2norm() for cddata2\n"
		    << "Expected value = sqrt(0.5) = 0.70710678118655, Computed value = "
		    << wl2norm << endl;
      }

      // Test #22: math::PatchCellDataOpsReal<NDIM>::RMSNorm() for cddata2
      // Expected: rmsnorm= L2-Norm/sqrt(control volume) = 2.0
      double rmsnorm = cdops_double.RMSNorm(cddata2,tpatch->getBox(),weight);
      if ( !tbox::MathUtilities<double>::equalEps(rmsnorm, 2.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #22: math::PatchCellDataOpsReal<NDIM>::RMSNorm() for cddata2\n"
		    << "Expected value = L2-Norm/sqrt(control volume) = 2.0, "
		    << "Computed value = " << rmsnorm << endl;
      }

      // Test #23: math::PatchCellDataOpsReal<NDIM>::weightedRMSNorm() for cddata2
      // Expected: wrmsnorm= Weighted L2-Norm/sqrt(control volume) = 1.0
      double wrmsnorm = cdops_double.weightedRMSNorm(cddata2,cddata1,tpatch->getBox(),weight);
      if ( !tbox::MathUtilities<double>::equalEps(wrmsnorm, 1.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #23: math::PatchCellDataOpsReal<NDIM>::weightedRMSNorm() for cddata2\n"
		    << "Expected value = Weighted L2-Norm/sqrt(control volume) = 1.0, "
		    << "Computed value = " << wrmsnorm << endl;
      }

      // Test #24: math::PatchCellDataOpsReal<NDIM>::maxNorm() for cddata2
      // Expected: maxnorm = 2.0
      double maxnorm = cdops_double.maxNorm(cddata2,tpatch->getBox(),weight);
      if ( !tbox::MathUtilities<double>::equalEps(maxnorm, 2.0) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #24: math::PatchCellDataOpsReal<NDIM>::maxNorm() for cddata2\n"
		    << "Expected value = 2.0, Computed value = " 
		    << maxnorm << endl;
      }

      // Reset cddata1 and cddata2
      cdops_double.setToScalar(cddata1, 5.0, cddata1->getGhostBox());
      cdops_double.setToScalar(cddata2, 3.0, cddata2->getGhostBox());

      // Test #25: math::PatchCellDataOpsReal<NDIM>::dotp() - (cddata1) * (cddata2)
      // Expected: dotp = 7.5
      double dotp = cdops_double.dot(cddata1, cddata2,tpatch->getBox(),weight);
      if ( !tbox::MathUtilities<double>::equalEps(dotp, 7.5) ) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #25: math::PatchCellDataOpsReal<NDIM>::dotp() - (cddata1) * (cddata2)\n"
		    << "Expected value = 7.5, Computed value = " 
		    << dotp<< endl;
      }
   
      // Test #26: Check state of hier::Patch<NDIM> before deallocating storage
      if (tpatch->getBox() != patch_box) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #26a: hier::Patch<NDIM> box incorrectly set\n"
		    << "Expected: d_box = " << patch_box  << "\n"
		    << "Set to: d_box = " << tpatch->getBox() << endl;
      }
      if (tpatch->getPatchNumber() != -1) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #26b: hier::Patch<NDIM> number incorrectly set\n"
		    << "Expected: d_patch_number = -1\n"
		    << "Set to: d_patch_number = " << tpatch->getPatchNumber() << endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; desc_id++) {

	 if (!tpatch->checkAllocated(desc_id)) {
	    num_failures++;
	    tbox::perr << "FAILED: - Test #26c.0: Descriptor index " << desc_id 
		       << " should be allocated but isn't!" << endl;
	 } else {

	    string patch_data_name = typeid(*tpatch->getPatchData(desc_id)).name();

	    hier::IntVector<NDIM> ghost_width = 
	       tpatch->getPatchData(desc_id)->getGhostCellWidth();
   
	    switch (desc_id) {
	       case 0:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.0.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(1)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.0.b: Ghost width incorrect\n"
				<< "Expected: (1,1)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
      
	       case 1:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.1.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(2)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.1.b: Ghost width incorrect\n"
				<< "Expected: (2,2)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 2:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.2.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(3)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.2.b: Ghost width incorrect\n"
				<< "Expected: (3,3)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 3:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, double >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.3.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, double >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(0)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.3.b: Ghost width incorrect\n"
				<< "Expected: (0,0)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 4:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, int >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.4.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, int >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(1)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.4.b: Ghost width incorrect\n"
				<< "Expected: (1,1)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	       case 5:
		  if (typeid(*tpatch->getPatchData(desc_id)) != 
		      typeid(pdat::CellData<NDIM, int >)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.5.a: hier::Patch<NDIM> Data name incorrect\n"
				<< "Expected: pdat::CellData<NDIM, int >\n"
				<< "Actual: " << patch_data_name << endl;
		  }
		  if (ghost_width != hier::IntVector<NDIM>(2)) {
		     num_failures++;
		     tbox::perr << "FAILED: - Test #26c.5.b: Ghost width incorrect\n"
				<< "Expected: (2,2)\n"
				<< "Actual: " << ghost_width << endl;
		  }
		  break;
   
	    };
	 }
      }

      // Deallocate all data on patch
      tpatch->deallocatePatchData(patch_components);

      // Test #27: Check state of hier::Patch<NDIM> after deallocating storage
      if (tpatch->getBox() != patch_box) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #27a: hier::Patch<NDIM> box incorrectly set\n"
		    << "Expected: d_box = " << patch_box  << "\n"
		    << "Set to: d_box = " << tpatch->getBox() << endl;
      }
      if (tpatch->getPatchNumber() != -1) {
	 num_failures++;
	 tbox::perr << "FAILED: - Test #27b: hier::Patch<NDIM> number incorrectly set\n"
		    << "Expected: d_patch_number = -1\n"
		    << "Set to: d_patch_number = " << tpatch->getPatchNumber() << endl;
      }
/* SAME ISSUE AS ABOVE FOR NUMBER OF COMPONENTS */
      for (desc_id = 0; desc_id < 6; desc_id++) {

	 if (tpatch->checkAllocated(desc_id)) {
	    num_failures++;
	    tbox::perr << "FAILED: - Test #27c: Descriptor index " << desc_id 
		       << " should be deallocated but isn't!" << endl;
	 } 
      }

      cwgt.setNull();
      cell_double_variable.setNull();
      cell_int_variable.setNull();

      if (num_failures == 0) {
	 tbox::pout << "\nPASSED:  cell patchtest" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(num_failures); 
}

/*
 * Returns true if all the data in the patch is equal to the specified
 * value.  Returns false otherwise.
 */
bool
doubleDataSameAsValue(int desc_id, double value,
		      tbox::Pointer<hier::Patch<NDIM> > patch)
{
   bool test_passed = true;

   tbox::Pointer< pdat::CellData<NDIM,double> > cvdata = patch->getPatchData(desc_id);

   for (pdat::CellIterator<NDIM> c(cvdata->getBox()); c && test_passed ; c++) {
      pdat::CellIndex<NDIM> cell_index = c();
      if ( !tbox::MathUtilities<double>::equalEps((*cvdata) (cell_index),value) ) {
         test_passed = false;
      }
   }

   return (test_passed);
}

