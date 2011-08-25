/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Main program to test index data operations
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// class holding information stored in index data
#include "SampleIndexData.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[]) {

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
// tbox::PIO::logOnlyNodeZero("indx_dataops.log");
   tbox::PIO::logAllNodes("indx_dataops.log");

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

/*
 ************************************************************************
 *
 *   Create a simple 2-level hierarchy to test.
 *   (NOTE: it is setup to work on at most 2 processors)
 *
 ************************************************************************
 */
      double lo[2] = { 0.0, 0.0 };
      double hi[2] = { 1.0, 0.5 };

      hier::Box coarse0(hier::Index(0, 0), hier::Index(9, 2));
      hier::Box coarse1(hier::Index(0, 3), hier::Index(9, 4));
      hier::Box fine0(hier::Index(4, 4), hier::Index(7, 7));
      hier::Box fine1(hier::Index(8, 4), hier::Index(13, 7));
      hier::IntVector<NDIM> ratio(2);

      hier::BoxList coarse_domain;
      hier::BoxList fine_domain;
      coarse_domain.appendItem(coarse0);
      coarse_domain.appendItem(coarse1);
      fine_domain.appendItem(fine0);
      fine_domain.appendItem(fine1);

      tbox::Pointer<geom::CartesianGridGeometry> geometry =
         new geom::CartesianGridGeometry("CartesianGeometry",
            lo,
            hi,
            coarse_domain);

      tbox::Pointer<hier::PatchHierarchy> hierarchy =
         new hier::PatchHierarchy("PatchHierarchy", geometry);

      // Note: For these simple tests we allow at most 2 processors.
      tbox::SAMRAI_MPI mpi(SAMRAIManager::getSAMRAICommWorld());
      const int nproc = mpi.getSize();
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

      hierarchy->makeNewPatchLevel(0, hier::IntVector<NDIM>(
            1), coarse_domain, mapping0);
      hierarchy->makeNewPatchLevel(1, ratio, fine_domain, mapping1);

      /*
       * Create an IndexData<SampleIndexData> variable and register it with
       * the variable database.
       */
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      tbox::Pointer<hier::VariableContext> cxt = variable_db->getContext(
            "dummy");
      const hier::IntVector<NDIM> no_ghosts(0);

      tbox::Pointer<pdat::IndexVariable<NDIM, SampleIndexData,
                                        pdat::CellGeometry> > data =
         new pdat::IndexVariable<NDIM, SampleIndexData, pdat::CellGeometry>(
            "sample");
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
      std::ostream& os = tbox::plog;
      for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; ln--) {
         tbox::Pointer<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         // allocate "sample" data
         level->allocatePatchData(data_id);
         os << "\nLevel: " << level->getLevelNumber() << " ";

         // loop over patches on level
         for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
            tbox::Pointer<hier::Patch> patch = level->getPatch(ip());
            os << "Patch: " << patch->getLocalId() << std::endl;

            // access sample data from patch
            tbox::Pointer<pdat::IndexData<NDIM, SampleIndexData,
                                          pdat::CellGeometry> > sample =
               patch->getPatchData(data_id);

            // iterate over cells of patch and invoke one "SampleIndexData"
            // instance on each cell (its possible to do more).
            for (pdat::CellIterator ic(patch->getBox()); ic; ic++) {
               SampleIndexData sd(ic());
               sd.setInt(counter);
               sample->appendItem(ic(), sd);
               counter++;
            }

            // iterate over the "SampleIndexData" index data stored on the patch
            // and dump the integer stored on it.
            for (pdat::IndexData<NDIM, SampleIndexData,
                                 pdat::CellGeometry>::Iterator id(*sample);
                 id;
                 id++) {
               os << "   Index: " << id().getIndex()
                  << "      SampleIndexData data: " << id().getInt()
                  << std::endl;
            }

         }
      }

      geometry.setNull();
      hierarchy.setNull();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
