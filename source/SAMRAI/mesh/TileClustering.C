/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Asynchronous Berger-Rigoutsos algorithm wrapper
 *
 ************************************************************************/
#ifndef included_mesh_TileClustering_C
#define included_mesh_TileClustering_C

#include <stdlib.h>

#include "SAMRAI/mesh/TileClustering.h"

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {

/*
 ************************************************************************
 * Constructor stores parameters of options for ussing
 * the asynchronous Berger-Rigoutsos implementation.
 ************************************************************************
 */
TileClustering::TileClustering(
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_box_size(hier::IntVector(d_dim, 8)),
   d_coalesce_boxes(true),
   d_log_cluster_summary(false),
   d_barrier_and_time(false),
   d_print_steps(false)
{
   getFromInput(input_db);
   setTimers();
   d_oca.setTimerPrefix("hier::TileClustering");
}

TileClustering::~TileClustering()
{
}

void
TileClustering::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{
   if (input_db) {

      if (input_db->isInteger("box_size")) {
         input_db->getIntegerArray("box_size",
            &d_box_size[0],
            d_dim.getValue());
      }

      d_coalesce_boxes =
         input_db->getBoolWithDefault("coalesce_boxes",
            d_coalesce_boxes);

      d_barrier_and_time =
         input_db->getBoolWithDefault("DEV_barrier_and_time",
            d_barrier_and_time);

      d_log_cluster =
         input_db->getBoolWithDefault("DEV_log_cluster",
            d_log_cluster);

      d_log_cluster_summary =
         input_db->getBoolWithDefault("DEV_log_cluster_summary",
            d_log_cluster_summary);

      d_print_steps =
         input_db->getBoolWithDefault("DEV_print_steps",
            d_print_steps);
   }
}

/*
 ************************************************************************
 * Implement the BoxGeneratorStrategy interface method using
 * the asynchronous Berger-Rigoutsos implementation.
 ************************************************************************
 */
void
TileClustering::findBoxesContainingTags(
   boost::shared_ptr<hier::BoxLevel>& new_box_level,
   boost::shared_ptr<hier::Connector>& tag_to_new,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_data_index,
   const int tag_val,
   const hier::BoxContainer& bound_boxes,
   const hier::IntVector& min_box,
   const hier::IntVector& max_gcw)
{
   NULL_USE(min_box);
   NULL_USE(max_gcw);

   TBOX_ASSERT(!bound_boxes.isEmpty());
   TBOX_ASSERT_OBJDIM_EQUALITY4(
      *tag_level,
      *(bound_boxes.begin()),
      min_box,
      max_gcw);

   if (d_barrier_and_time) {
      t_find_boxes_containing_tags->barrierAndStart();
   }

   t_cluster_setup->start();

   const hier::IntVector &zero_vector = hier::IntVector::getZero(tag_level->getDim());

   for (hier::BoxContainer::const_iterator bb_itr = bound_boxes.begin();
        bb_itr != bound_boxes.end(); ++bb_itr) {
      if (bb_itr->empty()) {
         TBOX_ERROR("TileClustering: empty bounding box not allowed.");
      }
   }

   const hier::BoxLevel& tag_box_level = *tag_level->getBoxLevel();

   new_box_level.reset( new hier::BoxLevel(
                           tag_box_level.getRefinementRatio(),
                           tag_box_level.getGridGeometry(),
                           tag_box_level.getMPI() ) );

   tag_to_new.reset(new hier::Connector(*tag_level->getBoxLevel(),
                                        *new_box_level,
                                        zero_vector));
   hier::Connector* new_to_tag = new hier::Connector(
      *new_box_level,
      *tag_level->getBoxLevel(),
      zero_vector);
   tag_to_new->setTranspose(new_to_tag, true);

   t_cluster_setup->stop();

   t_cluster->start();

   /*
    * Generate new_box_level and Connectors
    */
   for ( size_t pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {

      hier::Patch &patch = *tag_level->getPatch(pi);
      const hier::Box &patch_box = patch.getBox();
      const hier::BlockId &block_id = patch_box.getBlockId();

      TBOX_ASSERT( bound_boxes.begin(block_id) != bound_boxes.end(block_id) );
      const hier::Box &bounding_box = *bound_boxes.begin(block_id);

      if ( patch.getBox().intersects(bounding_box) ) {

         boost::shared_ptr<pdat::CellData<int> > tag_data(
            patch.getPatchData(tag_data_index), boost::detail::dynamic_cast_tag());

         hier::BoxContainer tiles;
         int num_coarse_tags =
            findTilesContainingTags( tiles, *tag_data, tag_val );

         if (d_print_steps) {
            tbox::plog << "Tile Clustering generated " << tiles.size()
                       << " clusters from " << num_coarse_tags
                       << " in patch " << patch.getBox().getBoxId() << '\n';
         }

         for ( hier::BoxContainer::iterator bi=tiles.begin(); bi!=tiles.end(); ++bi ) {
            new_box_level->addBoxWithoutUpdate(*bi);
            new_to_tag->insertLocalNeighbor( patch_box, bi->getBoxId() );
            tag_to_new->insertLocalNeighbor( *bi, patch_box.getBoxId() );
         }

      } // Patch is in bounding box

   } // Loop through tag level

   new_box_level->finalize();

   t_cluster->stop();


   if ( d_coalesce_boxes ) {

      /*
       * Try to coalesce the boxes in new_box_level.
       */
      hier::BoxContainer new_boxes;
      if (!new_box_level->getBoxes().isEmpty()) {

         t_coalesce->start();

         hier::LocalId local_id(0);

         const int nblocks = new_box_level->getGridGeometry()->getNumberBlocks();

         for (int b = 0; b < nblocks; ++b) {
            hier::BlockId block_id(b);

            hier::BoxContainer block_boxes(new_box_level->getBoxes(), block_id);

            if (!block_boxes.isEmpty()) {
               block_boxes.unorder();
               block_boxes.coalesce();
               new_boxes.spliceBack(block_boxes);
            }
         }

         t_coalesce->stop();

      }

      if ( d_print_steps ) {
         tbox::plog << "TileClustering coalesced " << new_box_level->getLocalNumberOfBoxes()
                    << " new boxes into " << new_boxes.size() << "\n";
      }

      if ( new_boxes.size() != static_cast<int>(new_box_level->getLocalNumberOfBoxes()) ) {

         t_coalesce_adjustment->start();

         /*
          * Coalesce changed the new boxes, so rebuild new_box_level and
          * Connectors.
          */
         new_box_level->initialize( new_box_level->getRefinementRatio(),
                                    new_box_level->getGridGeometry(),
                                    new_box_level->getMPI() );
         tag_to_new.reset( new hier::Connector( tag_box_level,
                                                *new_box_level,
                                                zero_vector ) );
         new_to_tag = new hier::Connector( *new_box_level,
                                           tag_box_level,
                                           zero_vector );
         tag_to_new->setTranspose(new_to_tag, true);

         const hier::BoxContainer &tag_boxes = tag_box_level.getBoxes();
         tag_boxes.makeTree( tag_box_level.getGridGeometry().get() );

         /*
          * Assign ids to coalesced boxes, add to BoxLevel and add
          * new--->tag edges.
          */
         hier::LocalId last_used_id(-1);
         const int rank = new_box_level->getMPI().getRank();
         for ( hier::BoxContainer::iterator bi=new_boxes.begin();
               bi!=new_boxes.end(); ++bi ) {

            bi->setId(hier::BoxId(++last_used_id,rank));

            hier::BoxContainer tmp_overlap_boxes;
            tag_boxes.findOverlapBoxes(tmp_overlap_boxes,
                                       *bi,
                                       tag_box_level.getRefinementRatio() );

            new_box_level->addBox(*bi);
            new_to_tag->insertNeighbors( tmp_overlap_boxes, bi->getBoxId() );

         }
         new_box_level->finalize();

         /*
          * Add tag--->new edges.
          */
         new_boxes.makeTree( new_box_level->getGridGeometry().get() );
         for ( hier::BoxContainer::const_iterator bi=tag_boxes.begin();
               bi!=tag_boxes.end(); ++bi ) {

            hier::BoxContainer tmp_overlap_boxes;
            new_boxes.findOverlapBoxes(tmp_overlap_boxes,
                                       *bi,
                                       tag_box_level.getRefinementRatio() );

            tag_to_new->insertNeighbors( tmp_overlap_boxes, bi->getBoxId() );
         }

         t_coalesce_adjustment->stop();

      }

   }



   /*
    * Get some global parameters.  Do it before logging to prevent
    * the logging flag from having an undue side effect on performance.
    */
   if (d_barrier_and_time) {
      t_global_reductions->barrierAndStart();
   }

   t_cluster_wrapup->start();

   new_box_level->getGlobalNumberOfBoxes();
   new_box_level->getGlobalNumberOfCells();
   if (d_barrier_and_time) {
      t_global_reductions->barrierAndStop();
   }
#if 0
   for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
        bi != bound_boxes.end(); ++bi) {
      new_box_level->getGlobalBoundingBox(bi->getBlockId().getBlockValue());
   }
#endif

   if (d_log_cluster) {
      tbox::plog << "TileClustering cluster log:\n"
      << "\tNew box_level clustered by TileClustering:\n" << new_box_level->format("\t\t",
         2)
      << "\tTileClustering tag_to_new:\n" << tag_to_new->format("\t\t", 2)
      << "\tTileClustering new_to_tag:\n" << new_to_tag->format("\t\t", 2);
   }
   if (d_log_cluster_summary) {
      /*
       * Log summary of clustering.
       */
      tbox::plog << "TileClustering summary:\n";

      for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
           bi != bound_boxes.end(); ++bi) {
         const int bn = bi->getBlockId().getBlockValue();
         tbox::plog << "Block " << bn
                    << " initial bounding box = " << *bi << ", "
                    << bi->size() << " cells, "
                    << "final global bounding box = "
                    << new_box_level->getGlobalBoundingBox(bn)
                    << ", "
                    << new_box_level->getGlobalBoundingBox(bn).size()
                    << " cells.\n\t";
      }

      tbox::plog << "Final output has "
                 << new_box_level->getGlobalNumberOfCells()
                 << " global cells [" << new_box_level->getMinNumberOfCells()
                 << "-" << new_box_level->getMaxNumberOfCells() << "], "
                 << new_box_level->getGlobalNumberOfBoxes()
                 << " global mapped boxes [" << new_box_level->getMinNumberOfBoxes()
                 << "-" << new_box_level->getMaxNumberOfBoxes() << "]\n"
                 << "\tTileClustering new_level summary:\n" << new_box_level->format("\t\t",0)
                 << "\tTileClustering new_level statistics:\n" << new_box_level->formatStatistics("\t\t")
                 << "\tTileClustering new_to_tag summary:\n" << new_to_tag->format("\t\t",0)
                 << "\tTileClustering new_to_tag statistics:\n" << new_to_tag->formatStatistics("\t\t")
                 << "\tTileClustering tag_to_new summary:\n" << tag_to_new->format("\t\t",0)
                 << "\tTileClustering tag_to_new statistics:\n" << tag_to_new->formatStatistics("\t\t")
                 << "\n";
   }

   t_cluster_wrapup->stop();

   if (d_barrier_and_time) {
      t_find_boxes_containing_tags->barrierAndStop();
   }
}



/*
 ***********************************************************************
 ***********************************************************************
 */
int
TileClustering::findTilesContainingTags(
   hier::BoxContainer &tiles,
   const pdat::CellData<int> &tag_data,
   int tag_val)
{
   tiles.clear();
   tiles.unorder();

   hier::Box coarsened_box(tag_data.getBox());
   coarsened_box.coarsen(d_box_size);

   const int num_coarse_cells = coarsened_box.size();

   for ( int coarse_offset=0; coarse_offset<num_coarse_cells; ++coarse_offset ) {
      const pdat::CellIndex coarse_cell_index(coarsened_box.index(coarse_offset));

      /*
       * Set the tile extent to cover the coarse cell and intersect
       * with tag box to make it nest in the tag box.  If any part
       * extends outside local tag boxes, (1) the tile might
       * overlap with remote clusters, (2) its overlap with remote
       * tag boxes might not be detected and (3) it may extend
       * outside the tag level.  Tiles extending past non-local tag
       * boxes can appear if the tag level patch boundaries do not
       * coincide with the tile cuts.
       */
      hier::Box tile_box(coarse_cell_index,coarse_cell_index,coarsened_box.getBlockId());
      tile_box.refine(d_box_size);
      tile_box *= tag_data.getBox();

      /*
       * Loop through fine cells in tile_box.  If any is tagged,
       * tile_box will be used as a cluster.
       */
      pdat::CellIterator finecend(pdat::CellGeometry::end(tile_box));
      for ( pdat::CellIterator fineci(pdat::CellGeometry::begin(tile_box));
            fineci!=finecend; ++fineci ) {
         if ( tag_data(*fineci) == tag_val ) {
            /*
             * Make a cluster from tile_box.
             * Choose a LocalId that is independent of ordering so that
             * results are independent of multi-threading.
             */
            hier::LocalId local_id(coarsened_box.getLocalId()*1000000 +
                                   coarse_offset);

            tile_box.initialize( tile_box,
                                 local_id,
                                 coarsened_box.getOwnerRank() );
            tiles.pushBack(tile_box);

            break;
         }

      } // Loop through fine cells in the tile.

   } // Loop through coarse cells (tiles).

   const int num_coarse_tags = tiles.size();

   tiles.order();

   if ( d_coalesce_boxes && tiles.size() > 1 ) {
      hier::LocalId last_used_id = tiles.back().getLocalId();
      // Coalesce the tiles in this patch and assign ids if they changed.
      hier::BoxContainer unordered_tiles( tiles.begin(), tiles.end(), false );
      unordered_tiles.coalesce();
      if ( unordered_tiles.size() != num_coarse_tags ) {
         tiles.clear();
         tiles.order();
         for ( hier::BoxContainer::iterator bi=unordered_tiles.begin(); bi!=unordered_tiles.end(); ++bi ) {
            bi->setLocalId(++last_used_id);
            tiles.insert(*bi);
         }
      }
   }

   return num_coarse_tags;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
boost::shared_ptr<pdat::CellData<int> >
TileClustering::makeCoarsenedTagData(const pdat::CellData<int> &tag_data,
                                           int tag_val) const
{
   hier::Box coarsened_box(tag_data.getBox());
   coarsened_box.coarsen(d_box_size);

   boost::shared_ptr<pdat::CellData<int> > coarsened_tag_data(
      new pdat::CellData<int>(coarsened_box,
                              1,
                              hier::IntVector::getZero(tag_data.getDim())));
   coarsened_tag_data->fill(0, 0);

   size_t tag_count = 0;
   size_t coarse_tag_count = 0;
   pdat::CellIterator finecend(pdat::CellGeometry::end(tag_data.getBox()));
   for ( pdat::CellIterator fineci(pdat::CellGeometry::begin(tag_data.getBox()));
         fineci!=finecend; ++fineci ) {

      if ( tag_data(*fineci) == tag_val ) {
         pdat::CellIndex coarseci = pdat::CellIndex( *fineci / d_box_size );

         coarse_tag_count += ( (*coarsened_tag_data)(coarseci) != tag_val );
         ++tag_count;

         (*coarsened_tag_data)(coarseci) = tag_val;
      }

   }
   if (d_print_steps) {
      tbox::plog << "TileClustering coarsened box " << tag_data.getBox()
                 << " (" << tag_count << " tags) to " << coarsened_box
                 << " (" << coarse_tag_count << " tags).\n";
   }

   return coarsened_tag_data;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
TileClustering::setTimers()
{
   t_find_boxes_containing_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()");
   t_cluster = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()_cluster");
   t_coalesce = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()_coalesce");
   t_coalesce_adjustment = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()_coalesce_adjustment");
   t_cluster_setup = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()_setup");
   t_cluster_wrapup = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::findBoxesContainingTags()_wrapup");
   t_global_reductions = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::global_reductions");
}

}
}
#endif
