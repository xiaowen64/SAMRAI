//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/multiblock/MultiblockCoarsenSchedule.C $
// Package:	SAMRAI multiblock
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2141 $
// Modified:	$LastChangedDate: 2008-04-23 08:36:33 -0700 (Wed, 23 Apr 2008) $
// Description:	Coarsening schedule for data transfer between AMR levels
//

#ifndef included_xfer_MultiblockCoarsenSchedule_C
#define included_xfer_MultiblockCoarsenSchedule_C
 
#include "MultiblockCoarsenSchedule.h"

#include "CoarsenCopyTransaction.h"
#include "tbox/ArenaManager.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"
#include "tbox/MathUtilities.h"


namespace SAMRAI {
    namespace xfer {

/*
 * ************************************************************************
 *                                                                        *
 * Create a coarsening schedule that transfers data from the source       *
 * patch data components of the fine level into the destination patch     *
 * data components of the coarse level.  If the coarsening operators      *
 * require data in ghost cells on the source level, then those ghost	  *
 * cells must be filled before this call.				  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
MultiblockCoarsenSchedule<DIM>::MultiblockCoarsenSchedule(
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > crse_level,
   tbox::Pointer<hier::MultiblockPatchLevel<DIM> > fine_level,
   const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes,
   tbox::Pointer<hier::MultiblockPatchHierarchy<DIM> > multiblock,
   MultiblockCoarsenPatchStrategy<DIM>* patch_strategy,
   MultiblockRefinePatchStrategy<DIM>* refine_strategy,
   bool fill_coarse_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!crse_level.isNull());
   TBOX_ASSERT(!fine_level.isNull());
   TBOX_ASSERT(!coarsen_classes.isNull());
#endif

   /*
    * Initial values; some will change in setup operations.
    */

   d_mblk_crse_level             = crse_level;
   d_mblk_fine_level             = fine_level;
   d_mblk_temp_crse_level.setNull();

   d_mblk_coarsen_patch_strategy = patch_strategy;
   d_mblk_refine_strategy = refine_strategy;

   d_multiblock_hier = multiblock;

   d_fill_coarse_data       = fill_coarse_data;

   d_schedule.setNull();

   d_mblk_refine_alg.setNull();
   d_mblk_refine_sched.setNull();

   /*
    * Compute ratio between fine and coarse levels.
    */
   hier::IntVector<DIM> fine = d_mblk_fine_level->getRatio();
   hier::IntVector<DIM> crse = d_mblk_crse_level->getRatio();
   int i;
   for (i = 0; i < DIM; i++) {
      if (fine(i) > 1) {
         d_ratio(i) = fine(i) / crse(i);
      } else {
         d_ratio(i) = 
            tbox::MathUtilities<int>::Abs( crse(i) / fine(i) );
      }
   }

   setCoarsenItems(coarsen_classes);
   initialCheckCoarsenClassItems();

   /*
    * Set up refine schedule to fill temporary coarse level data 
    * before coarsening operations, if needed.
    */

   if (d_fill_coarse_data) {
      setupRefineAlgorithm();
   }

   t_gen_sched = tbox::TimerManager::getManager()->
      getTimer("xfer::CoarsenSchedule<DIM>::generate_comm_schedule");

   t_gen_sched->start();
   generateSchedule();
   t_gen_sched->stop();
}

/*
 * ************************************************************************
 *                                                                        *
 * The destructor for the coarsen schedule class implicitly deallocates   *
 * all of the data associated with the communication schedule.            *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
MultiblockCoarsenSchedule<DIM>::~MultiblockCoarsenSchedule()
{
   if (d_coarsen_items) {
      for (int ici = 0; ici < d_number_coarsen_items; ici++) {
         d_coarsen_items[ici] =
            (typename xfer::CoarsenClasses<DIM>::Data*)NULL;
      }
   }
   delete [] d_coarsen_items;
}

/*
 * ************************************************************************
 *                                                                        *
 * Return const pointer to equivalence classes used in schedule.          *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
const tbox::Pointer< xfer::CoarsenClasses<DIM> >&
MultiblockCoarsenSchedule<DIM>::getEquivalenceClasses() const
{
   return(d_coarsen_classes);
}

/*
 * ************************************************************************
 *                                                                        *
 * Execute the stored communication schedule that copies data into the	  *
 * the destination patch data components of the destination level from	  *
 * the source patch data components of the source level.  The steps	  *
 * to the algorithm are as follows:					  *
 *									  *
 *	(1) Allocate the source space on the temporary patch level.	  *
 *	(2) Coarsen the data from the fine patch level to the temporary	  *
 *	    patch level (local operation).				  *
 *	(3) Copy data from the source space of the temporary patch	  *
 *	    level into the destination space of the destination patch	  *
 *	    level (requires interprocessor communication).		  *
 *	(4) Deallocate the source space on the temporary patch level.	  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM>
void MultiblockCoarsenSchedule<DIM>::coarsenData() const
{
   tbox::Pointer<tbox::Timer> t_coarsen_data = 
      tbox::TimerManager::getManager()->
      getTimer("mblk::MultiblockCoarsenSchedule<DIM>::coarsenData");
   t_coarsen_data->start();

   /*
    * Set the refine items for the copy transactions.  These items
    * will be shared by all such transaction objects in the communication
    * schedule.
    */
   xfer::CoarsenCopyTransaction<DIM>::setCoarsenItems(d_coarsen_items,
                                                      d_number_coarsen_items);

   /*
    * Allocate the source data space on the temporary patch level.
    * We do not know the current time, so set it to zero.  It should
    * not make a difference, since the copy routines should not check
    * whether the time markers match.
    */

   d_mblk_temp_crse_level->allocatePatchData(d_sources, 0.0,
      tbox::ArenaManager::getManager()->getScratchAllocator());

   if (d_fill_coarse_data) {
      d_mblk_refine_sched->fillData(0.0, false);
   }

   /*
    * Coarsen the data from the sources on the fine data level into the
    * sources on the temporary data level
    */

   coarsenSourceData(d_mblk_coarsen_patch_strategy);

   /*
    * Copy data from the source interiors of the temporary patch level
    * into the destination interiors of the destination patch level.
    */

   d_schedule->communicate();

   /*
    * Deallocate the source data in the temporary patch level.
    */

   d_mblk_temp_crse_level->deallocatePatchData(d_sources);

   /*
    * Unset the refine items for the copy transactions.  These items
    * are shared by all such transaction objects in the communication
    * schedule.
    */
   xfer::CoarsenCopyTransaction<DIM>::unsetCoarsenItems();

   t_coarsen_data->stop();
}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to set up local array of coarsen items.      *
 *                                                                       *
 * ***********************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::setCoarsenItems(
   const tbox::Pointer< xfer::CoarsenClasses<DIM> > coarsen_classes)
{

   d_coarsen_classes        = coarsen_classes;
   d_number_coarsen_items   = 0;
   d_coarsen_items          =
      (const typename xfer::CoarsenClasses<DIM>::Data**)NULL;

   const int num_coarsen_classes =
      d_coarsen_classes->getNumberOfEquivalenceClasses();

   /*
    * Determine total number of coarsen items and set state of 
    * d_sources component selector.
    */
   d_sources.clrAllFlags();
   d_number_coarsen_items = 0;

   int nc;
   for (nc = 0; nc < num_coarsen_classes; nc++) {
      for (typename
           tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator
           l(d_coarsen_classes->getIterator(nc)); l; l++) {
         d_sources.setFlag(l().d_src);
         d_number_coarsen_items++;
      }
   }

   /*
    * Allocate and initialize array of coarsen items.
    */

   d_coarsen_items =
      new const typename
         xfer::CoarsenClasses<DIM>::Data*[d_number_coarsen_items];

   int ircount = 0;
   for (nc = 0; nc < num_coarsen_classes; nc++) {
      for (typename
           tbox::List<typename xfer::CoarsenClasses<DIM>::Data>::Iterator
           l(d_coarsen_classes->getIterator(nc)); l; l++) {
         d_coarsen_items[ircount] = &(l());
         ircount++;
      }
   }

}

/*
 * ***********************************************************************
 *                                                                       *
 * Private utility function to check coarsen items in initial setup to   *
 * see whether source and destination patch data components have         *
 * sufficient ghost cell widths to satisfy the "ghost width to coarsen"  *
 * functionality described in the xfer::CoarsenAlgorithm<DIM> class header.   *
 * Specifically, the destination data must have a ghost cell width at    *
 * least as large as the ghost cell width to coarsen.  The source data   *
 * must have a ghost cell width at least as large as the ghost cell      *
 * width to coarsen refined to the source (finer) level index space.     *
 * Other checks are also performed here by calling the                   *
 * xfer::CoarsenClasses<DIM>::checkCoarsenItem() routine.                     *
 *                                                                       *
 * ***********************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::initialCheckCoarsenClassItems() const
{
   tbox::Pointer< hier::PatchDescriptor<DIM> > pd;
   pd.setNull();

   for (int nb = 0; nb < d_mblk_crse_level->getNumberOfBlocks(); nb++) {
      tbox::Pointer< hier::PatchLevel<DIM> > block_level =
         d_mblk_crse_level->getPatchLevelForBlock(nb);
      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }
      if (!block_level.isNull()) {
          pd = block_level->getPatchDescriptor();
          break;
      }
      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!pd.isNull());
#endif

   hier::IntVector<DIM> user_gcw(0);
   if (d_mblk_coarsen_patch_strategy) {
      user_gcw = d_mblk_coarsen_patch_strategy->getCoarsenOpStencilWidth();
   }

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {

      const typename xfer::CoarsenClasses<DIM>::Data* const crs_item =
         d_coarsen_items[ici];

#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_coarsen_classes->checkCoarsenItem(*crs_item, pd)) {
#endif

         const int dst_id = crs_item->d_dst;
         const int src_id = crs_item->d_src;

         tbox::Pointer< hier::PatchDataFactory<DIM> > dfact =
            pd->getPatchDataFactory(dst_id);
         tbox::Pointer< hier::PatchDataFactory<DIM> > sfact =
            pd->getPatchDataFactory(src_id);

         const hier::IntVector<DIM>& dst_gcw = dfact->getGhostCellWidth();
         const hier::IntVector<DIM>& src_gcw = sfact->getGhostCellWidth();

         if (crs_item->d_gcw_to_coarsen > dst_gcw) {
            TBOX_ERROR("Bad data given to xfer::CoarsenSchedule<DIM>...\n"
                       << "`ghost cell width to coarsen' specified in\n"
                       << "registration of `Destination' patch data " 
                       << pd->mapIndexToName(dst_id)
                       << " with xfer::CoarsenAlgorithm<DIM>\n"
                       << " is larger than ghost cell width of data \n"
                       << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
                       << "\n data ghost cell width = " << dst_gcw << std::endl);
         }

         if ( (crs_item->d_gcw_to_coarsen * d_ratio) > src_gcw ) { 
            TBOX_ERROR("Bad data given to xfer::CoarsenSchedule<DIM>...\n"
                       << "`Source' patch data " << pd->mapIndexToName(src_id)
                       << " has ghost cell width too small to support the\n"
                       << "`ghost cell width to coarsen' specified in"
                       << " registration with xfer::CoarsenAlgorithm<DIM>\n"
                       << "data ghost cell width = " << src_gcw
                       << "d_gcw_to_coarsen = " << crs_item->d_gcw_to_coarsen
                       << "\nratio between levels = " << d_ratio
                       << "\n Thus, data ghost width must be >= "
                       << (crs_item->d_gcw_to_coarsen * d_ratio) 
                       << std::endl);
         }

#ifdef DEBUG_CHECK_ASSERTIONS
      }
#endif

   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Generate the temporary coarse level by coarsening the fine patch	  *
 * level boxes.  Note that no patch data components are allocated until   *
 * they are needed during the coarsening operation.			  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::generateTemporaryLevel()
{

#ifdef DEBUG_CHECK_ASSERTIONS
   int i;
   for (i = 0; i < DIM; i++) {
      TBOX_ASSERT( d_ratio(i) != 0 );
   }
   if (DIM > 1) {
      for (i = 0; i < DIM; i++) {
         if ( d_ratio(i)*d_ratio((i+1)%DIM) < 0 ) {
            TBOX_ASSERT( (d_ratio(i) == 1) || (d_ratio((i+1)%DIM) == 1) );
         }
      }
   }
#endif

   tbox::Array< tbox::Pointer< hier::PatchLevel<DIM> > >
      temp_crse_array(d_mblk_fine_level->getNumberOfBlocks());

   for (int nb = 0; nb < d_mblk_fine_level->getNumberOfBlocks(); nb++) {
      temp_crse_array[nb].setNull();

      tbox::Pointer< hier::PatchLevel<DIM> > fine_block_level =
         d_mblk_fine_level->getPatchLevelForBlock(nb);

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }
      

      if (!fine_block_level.isNull()) {
         hier::BoxArray<DIM> boxes = fine_block_level->getBoxes();
         boxes.coarsen(d_ratio);

         temp_crse_array[nb] = new hier::PatchLevel<DIM>();
         temp_crse_array[nb]->setCoarsenedPatchLevel(fine_block_level, d_ratio);
         temp_crse_array[nb]->setLevelNumber(d_mblk_crse_level->getLevelNumber());
         temp_crse_array[nb]->setNextCoarserHierarchyLevelNumber(
                         d_mblk_crse_level->getLevelNumber());

      }

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }

   }

   d_mblk_temp_crse_level = new hier::MultiblockPatchLevel<DIM>(temp_crse_array);
   if (d_fill_coarse_data) {
      d_multiblock_hier->adjustMultiblockPatchLevelBoundaries(d_mblk_temp_crse_level);
   }
}

/*
 * ************************************************************************
 *                                                                        *
 * Set up refine algorithm to fill temporary coarse level before          *
 * performing coarsening operations, if necessary.                        *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::setupRefineAlgorithm()
{

   d_mblk_refine_alg.setNull();

   tbox::Pointer< xfer::RefineAlgorithm<DIM> > refine_alg = new xfer::RefineAlgorithm<DIM>();
   d_mblk_refine_alg = new MultiblockRefineAlgorithm<DIM>(refine_alg,
                                                            d_multiblock_hier);

   for (int ici = 0; ici < d_number_coarsen_items; ici++) {
      const int src_id = d_coarsen_items[ici]->d_src;
      d_mblk_refine_alg->registerRefine(src_id, src_id, src_id, NULL);
   }

}

/*
 * ************************************************************************
 *                                                                        *
 * Generate communication schedule that copies source patch data          *
 * from the temporary level into the destination patch data of the        *
 * destination (coarse) level.  The source and destination	          *
 * spaces may be the same.						  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::generateSchedule()
{

   /*
    * Set up coarsened version of fine level for temporary data storage.
    * Then, create refine algorithm if needed to fill temporary coarse
    * level before coarsen operations occur.
    */
   generateTemporaryLevel();

   if (d_fill_coarse_data) {
      d_mblk_refine_sched =
         d_mblk_refine_alg->createSchedule(d_mblk_temp_crse_level, d_mblk_crse_level,
                                            d_mblk_refine_strategy);
   }

   d_schedule = new tbox::Schedule();

   /*
    * Define local constants for easier access to level information
    */

   for (int nb = 0; nb < d_mblk_crse_level->getNumberOfBlocks(); nb++) {

      tbox::Pointer< hier::PatchLevel<DIM> > crse_block_level =
         d_mblk_crse_level->getPatchLevelForBlock(nb);

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }

      if (!crse_block_level.isNull()) {

         tbox::Pointer< hier::PatchLevel<DIM> > temp_crse_block_level =
            d_mblk_temp_crse_level->getPatchLevelForBlock(nb);

         if (!temp_crse_block_level.isNull()) { 
            const int dst_npatches =
               crse_block_level->getNumberOfPatches();
            const int src_npatches =
               temp_crse_block_level->getNumberOfPatches();

            const hier::ProcessorMapping& dst_mapping =
               crse_block_level->getProcessorMapping();
            const hier::ProcessorMapping& src_mapping =
               temp_crse_block_level->getProcessorMapping();

            const hier::BoxArray<DIM>& dst_boxes = crse_block_level->getBoxes();
            const hier::BoxArray<DIM>& src_boxes = temp_crse_block_level->getBoxes();

            tbox::Pointer< hier::PatchDescriptor<DIM> > pd =
               crse_block_level->getPatchDescriptor();

            const hier::IntVector<DIM> no_shift(0);
            tbox::Pointer< hier::PatchDataFactory<DIM> > src_pdf;
            tbox::Pointer< hier::PatchDataFactory<DIM> > dst_pdf;
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap;

            const int num_equiv_classes =
               d_coarsen_classes->getNumberOfEquivalenceClasses();

            /*
             * Loop over all of the patches on the source and destination levels.
             * Check to see whether source or destination is local to this processor.
             * If not, then skip over schedule generation operations.
             */

            for (int dp = 0; dp < dst_npatches; dp++) {

               for (int sp = 0; sp < src_npatches; sp++) {

                  if (dst_mapping.isMappingLocal(dp)
                      || src_mapping.isMappingLocal(sp)) {

                     int citem_count = 0;
                     for (int nc = 0; nc < num_equiv_classes; nc++) {

                        const typename xfer::CoarsenClasses<DIM>::Data&
                           rep_item =
                              d_coarsen_classes->getClassRepresentative(nc);

                        const int rep_item_dst_id = rep_item.d_dst; 
                        const int rep_item_src_id = rep_item.d_src;

                        /*
                         * Get the patch data factories and calculate the overlap.
                         * A single overlap is calculated for all of the items
                         * in an equivalence class.
                         */

                        src_pdf = pd->getPatchDataFactory(rep_item_src_id);
                        dst_pdf = pd->getPatchDataFactory(rep_item_dst_id);

                        hier::Box<DIM> src_mask(src_boxes[sp]);
                        if (rep_item.d_gcw_to_coarsen != hier::IntVector<DIM>(0)) {
                           src_mask.grow(hier::IntVector<DIM>::min(
                                          rep_item.d_gcw_to_coarsen,
                                          src_pdf->getGhostCellWidth()));
                        }

                        overlap = dst_pdf->getBoxGeometry(dst_boxes[dp])
                                         ->calculateOverlap(
                                        *src_pdf->getBoxGeometry(src_boxes[sp]),
                                        src_mask, true, no_shift);

#ifdef DEBUG_CHECK_ASSERTIONS
	                if (overlap.isNull()) {
		           TBOX_ERROR("Internal xfer::CoarsenSchedule<DIM> error..."
                                      << "\n Overlap is NULL for "
                                      << "\n src box = " << src_boxes[sp]
                                      << "\n dst box = " << dst_boxes[dp]
                                      << "\n src mask = " << src_mask << std::endl);
                        }
#endif

                        if (!overlap->isOverlapEmpty()) {
                           /*
                            * Iterate over the components of the coarsen descriptor
                            * list and generate appropriate schedule information.
                            */
                           for (typename
                                tbox::List<typename
                                xfer::CoarsenClasses<DIM>::Data>::Iterator
                                l(d_coarsen_classes->getIterator(nc));
                                l; l++) {

                              d_schedule->addTransaction(
                                 new xfer::CoarsenCopyTransaction<DIM>(
                                    crse_block_level,
                                    temp_crse_block_level,
                                    overlap,
                                    dp, sp,
                                    citem_count) );
                              citem_count++;

                           } // iterate over coarsen components in equivalence class

                        } else {
                           for (typename
                                tbox::List<typename
                                xfer::CoarsenClasses<DIM>::Data>::Iterator
                                l(d_coarsen_classes->getIterator(nc));
                                l; l++) {
                              citem_count++;
                           }
                        }// if overlap is not empty

                     }  // iterate over all coarsen equivalence classes

                  }  // if either source or destination is local

               } // loop over source patches

            } // loop over destination patches
         } // if temp_crse_level exists for this block
      } // if crse_level exists for this block

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }

   } // loob over blocks in multiblock domain
}

/*
 * ************************************************************************
 *                                                                        *
 * Coarsen data from the source space on the fine patch level into the	  *
 * source space on the coarse temporary patch level.			  *
 *                                                                        *
 * ************************************************************************
 */

template<int DIM> void MultiblockCoarsenSchedule<DIM>::coarsenSourceData(
   MultiblockCoarsenPatchStrategy<DIM>* patch_strategy) const
{

   for (int nb = 0; nb < d_mblk_fine_level->getNumberOfBlocks(); nb++) {

      tbox::Pointer< hier::PatchLevel<DIM> > fine_block_level = 
         d_mblk_fine_level->getPatchLevelForBlock(nb);

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->setCoarsenBlockNumber(nb);
      }

      if (!fine_block_level.isNull()) {

         tbox::Pointer< hier::PatchLevel<DIM> > temp_crse_block_level =
            d_mblk_temp_crse_level->getPatchLevelForBlock(nb);

         if (!temp_crse_block_level.isNull()) {

            /*
             * Loop over all local patches (fine and temp have the same mapping)
             */

            for (typename hier::PatchLevel<DIM>::Iterator p(fine_block_level);
                 p; p++) {
               tbox::Pointer< hier::Patch<DIM> > fine_patch =
                  fine_block_level->getPatch(p());
               tbox::Pointer< hier::Patch<DIM> > temp_patch = 
                  temp_crse_block_level->getPatch(p());

               const hier::Box<DIM>& box = temp_patch->getBox();

               /*
                * Coarsen the fine space onto the temporary coarse space
                */

               if (patch_strategy) {
                  patch_strategy->preprocessCoarsen(*temp_patch,
                                                    *fine_patch, box, d_ratio);
               }

               for (int ici = 0; ici < d_number_coarsen_items; ici++) {
                  const typename xfer::CoarsenClasses<DIM>::Data* const
                     crs_item = d_coarsen_items[ici];
                  if (!(crs_item->d_opcoarsen.isNull())) {
                     const int source_id = crs_item->d_src;
                     crs_item->d_opcoarsen->coarsen(*temp_patch, *fine_patch, 
                                                    source_id, source_id, 
                                                    box, d_ratio);
                  }
               }

               if (patch_strategy) {
                  patch_strategy->postprocessCoarsen(*temp_patch,
                                                     *fine_patch, box, d_ratio);
               }
            }
         }
      }

      if (d_mblk_coarsen_patch_strategy) {
         d_mblk_coarsen_patch_strategy->clearCoarsenBlockNumber();
      }

   }
}

}
}
#endif
