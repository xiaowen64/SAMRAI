//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/examples/dataops/indx_dataops.C $
// Package:     SAMRAI test
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1917 $
// Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
// Description: Main program to test index data operations
//

#include "SAMRAI_config.h"



// class holding information stored in index data
#include "SampleIndexData.h"

#include "tbox/SAMRAIManager.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIterator.h"
#include "IndexData.h"
#include "IndexVariable.h"
#include "IntVector.h"
#include "ProcessorMapping.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Patch.h"
#include "tbox/IOStream.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"
#include "VariableDatabase.h"
#include "VariableContext.h"

using namespace SAMRAI;

int main( int argc, char *argv[] ) {

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::SAMRAI_MPI::init(&argc, &argv);
      tbox::SAMRAIManager::startup();
// tbox::PIO::logOnlyNodeZero("indx_dataops.log");
      tbox::PIO::logAllNodes("indx_dataops.log");


/*
************************************************************************
*
*   Create a simple 2-level hierarchy to test.
*   (NOTE: it is setup to work on at most 2 processors)
*
************************************************************************
*/
      double lo[2] = {0.0, 0.0};
      double hi[2] = {1.0, 0.5};

      hier::Box<NDIM> coarse0(hier::Index<NDIM>(0,0), hier::Index<NDIM>(9,2));
      hier::Box<NDIM> coarse1(hier::Index<NDIM>(0,3), hier::Index<NDIM>(9,4));
      hier::Box<NDIM> fine0(hier::Index<NDIM>(4,4), hier::Index<NDIM>(7,7));
      hier::Box<NDIM> fine1(hier::Index<NDIM>(8,4), hier::Index<NDIM>(13,7));
      hier::IntVector<NDIM> ratio(2);

      hier::BoxArray<NDIM> coarse_domain(2);
      hier::BoxArray<NDIM> fine_domain(2);
      coarse_domain(0) = coarse0;
      coarse_domain(1) = coarse1;
      fine_domain(0) = fine0;
      fine_domain(1) = fine1;

      tbox::Pointer<geom::CartesianGridGeometry<NDIM> > geometry =
	 new geom::CartesianGridGeometry<NDIM>("CartesianGeometry",lo, hi, coarse_domain);

      tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy = 
	 new hier::PatchHierarchy<NDIM>("PatchHierarchy",geometry);

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

      /*
       * Create an IndexData<SampleIndexData> variable and register it with 
       * the variable database.
       */
      hier::VariableDatabase<NDIM>* variable_db = hier::VariableDatabase<NDIM>::getDatabase();
      tbox::Pointer< hier::VariableContext > cxt = variable_db->getContext("dummy");
      const hier::IntVector<NDIM> no_ghosts(0);

      tbox::Pointer< pdat::IndexVariable<NDIM, SampleIndexData> > data = 
	 new pdat::IndexVariable<NDIM,SampleIndexData>("sample");
      int data_id = variable_db->registerVariableAndContext(
	 data, cxt, no_ghosts);


/*
************************************************************************
*
*   Set index data.
*
************************************************************************
*/

      /*
       * Loop over hierarchy levels and set index data on cells of patches
       */
      int counter = 0;
      ostream& os = tbox::plog;
      for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
	 tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

	 // allocate "sample" data
	 level->allocatePatchData(data_id);
	 os << "\nLevel: " << level->getLevelNumber() << " ";

	 // loop over patches on level
	 for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
	    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(ip());     
	    os << "Patch: " << patch->getPatchNumber() << endl;

	    // access sample data from patch
	    tbox::Pointer< pdat::IndexData<NDIM, SampleIndexData> > sample = 
	       patch->getPatchData(data_id);

	    // iterate over cells of patch and invoke one "SampleIndexData" 
	    // instance on each cell (its possible to do more).
	    for (pdat::CellIterator<NDIM> ic(patch->getBox()); ic; ic++) {
	       SampleIndexData sd(ic());
	       sd.setInt(counter);
	       sample->appendItem(ic(), sd);
	       counter++;
	    }

	    // iterate over the "SampleIndexData" index data stored on the patch
	    // and dump the integer stored on it.
	    for (pdat::IndexData<NDIM,SampleIndexData>::Iterator id(*sample); id; id++) {
	       os << "   Index: " << id().getIndex() 
		  << "      SampleIndexData data: " << id().getInt()
		  << endl;
	    }
         
	 }
      }


      geometry.setNull();
      hierarchy.setNull();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAI_MPI::finalize();
 
   return(0); 
}


