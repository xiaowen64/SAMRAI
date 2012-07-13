/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test index data operations
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// class holding information stored in index data
#include "SampleIndexData.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include "boost/shared_ptr.hpp"

using namespace SAMRAI;

int main(
   int argc,
   char* argv[]) {

   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI::MAX_DIM_VAL);
   const tbox::Dimension dim(d);

   if (dim != tbox::Dimension(2)) {
      TBOX_ERROR("This test code is completed only for 2D!!!");
   }

   const std::string log_fn = std::string("indx_dataops.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

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

      hier::Box coarse0(hier::Index(0, 0), hier::Index(9, 2), hier::BlockId(0));
      hier::Box coarse1(hier::Index(0, 3), hier::Index(9, 4), hier::BlockId(0));
      hier::Box fine0(hier::Index(4, 4), hier::Index(7, 7), hier::BlockId(0));
      hier::Box fine1(hier::Index(8, 4), hier::Index(13, 7), hier::BlockId(0));
      hier::IntVector ratio(dim, 2);

      coarse0.initialize(coarse0, hier::LocalId(0), 0);
      coarse1.initialize(coarse1, hier::LocalId(1), 0);
      fine0.initialize(fine0, hier::LocalId(0), 0);
      fine1.initialize(fine1, hier::LocalId(1), 0);

      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_domain;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_domain.pushBack(fine0);
      fine_domain.pushBack(fine1);

      boost::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      const int nproc = mpi.getSize();

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_domain.size();

      hier::BoxLevel layer0(hier::IntVector(dim, 1), geometry);
      hier::BoxLevel layer1(ratio, geometry);

      hier::BoxContainer::iterator coarse_itr(coarse_domain);
      for (int ib = 0; ib < n_coarse_boxes; ++ib, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0.getMPI().getRank()) {
               layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib),
                  layer0.getMPI().getRank()));
            }
         } else {
            layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_itr(fine_domain);
      for (int ib = 0; ib < n_fine_boxes; ++ib) {
         if (nproc > 1) {            if (ib == layer1.getMPI().getRank()) {
               layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib),
                  layer1.getMPI().getRank()));
            }
         } else {
            layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      /*
       * Create an IndexData<SampleIndexData> variable and register it with
       * the variable database.
       */
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      boost::shared_ptr<hier::VariableContext> cxt(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim, 0);

      boost::shared_ptr<pdat::IndexVariable<SampleIndexData,
                                            pdat::CellGeometry> > data(
         new pdat::IndexVariable<SampleIndexData, pdat::CellGeometry>(
            dim, "sample"));
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
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));

         // allocate "sample" data
         level->allocatePatchData(data_id);
         os << "\nLevel: " << level->getLevelNumber() << " ";

         // loop over patches on level
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            boost::shared_ptr<hier::Patch> patch(*ip);
            os << "Patch: " << patch->getLocalId() << std::endl;

            // access sample data from patch
            boost::shared_ptr<pdat::IndexData<SampleIndexData,
                              pdat::CellGeometry> > sample(
               patch->getPatchData(data_id),
               boost::detail::dynamic_cast_tag());

            // iterate over cells of patch and invoke one "SampleIndexData"
            // instance on each cell (its possible to do more).
            pdat::CellIterator icend(patch->getBox(), false);
            for (pdat::CellIterator ic(patch->getBox(), true);
                 ic != icend; ++ic) {
               SampleIndexData sd;
               sd.setInt(counter);
               sample->appendItem(*ic, sd);
               ++counter;
            }

            // iterate over the "SampleIndexData" index data stored on the patch
            // and dump the integer stored on it.
            int currData = counter - 1;
	    pdat::IndexData<SampleIndexData, pdat::CellGeometry>::iterator idend(*sample, false);
            for (pdat::IndexData<SampleIndexData,
                                 pdat::CellGeometry>::iterator id(*sample, true);
                 id != idend;
                 ++id) {
               os << "      SampleIndexData data: " << id->getInt()
                  << std::endl;
               if (id->getInt() != currData) {
                  ++num_failures;
                  tbox::perr
                  << "FAILED: - Index data set incorrectly" << std::endl;
               }
               --currData;
            }

         }
      }

      geometry.reset();
      hierarchy.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  indx dataops" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}
