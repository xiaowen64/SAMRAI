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
#include "SAMRAI/tbox/OpenMPUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {

const std::string TileClustering::s_default_timer_prefix("mesh::TileClustering");
std::map<std::string, TileClustering::TimerStruct> TileClustering::s_static_timers;

int TileClustering::s_primary_mpi_tag = 1234;
int TileClustering::s_secondary_mpi_tag = 1235;
int TileClustering::s_first_data_length = 1000;

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
   d_allow_remote_tile_extent(false),
   d_coalesce_boxes(true),
   d_log_cluster_summary(false),
   d_log_cluster(false),
   d_barrier_and_time(false),
   d_print_steps(false)
{
   TBOX_omp_init_lock(&l_outputs);
   TBOX_omp_init_lock(&l_interm);
   getFromInput(input_db);
   setTimerPrefix(s_default_timer_prefix);
   d_oca.setTimerPrefix(s_default_timer_prefix);
}

TileClustering::~TileClustering()
{
   TBOX_omp_destroy_lock(&l_outputs);
   TBOX_omp_destroy_lock(&l_interm);
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

      d_allow_remote_tile_extent =
         input_db->getBoolWithDefault("allow_remote_tile_extent",
            d_allow_remote_tile_extent);

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
      d_object_timers->t_find_boxes_containing_tags->barrierAndStart();
   }

   d_object_timers->t_cluster_setup->start();

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

   tag_box_level.getBoxes().makeTree(tag_box_level.getGridGeometry().get());

   d_object_timers->t_cluster_setup->stop();

   int tiles_have_remote_extent = 0;

   if ( d_allow_remote_tile_extent ) {

      clusterWholeTiles(
         *new_box_level,
         tag_to_new,
         tiles_have_remote_extent,
         tag_level,
         bound_boxes,
         tag_data_index,
         tag_val);

      new_box_level->getMPI().AllReduce( &tiles_have_remote_extent, 1, MPI_MAX );

      if ( tiles_have_remote_extent ) {
         detectSemilocalEdges( tag_to_new );
         /*
          * Remove duplicated new tiles.  For each set of coinciding tiles,
          * determine the process with the greatest tag overlap and keep only
          * the copy from that process.  Discard the others.
          */
         removeDuplicateTiles( *new_box_level, *tag_to_new );
      }
   }

   else {

      clusterWithinProcessBoundaries(
         *new_box_level,
         *tag_to_new,
         tag_level,
         bound_boxes,
         tag_data_index,
         tag_val);

   }


   if ( d_coalesce_boxes ) {
      coalesceClusters(*new_box_level, tag_to_new, tiles_have_remote_extent);
   }



   /*
    * Get some global parameters.  Do it before logging to prevent the
    * logging flag from having a side-effect on performance timers.
    */

   d_object_timers->t_cluster_wrapup->barrierAndStart();

   if (d_barrier_and_time) {
      d_object_timers->t_global_reductions->start();
   }
   new_box_level->getGlobalNumberOfBoxes();
   new_box_level->getGlobalNumberOfCells();
   if (d_barrier_and_time) {
      d_object_timers->t_global_reductions->barrierAndStop();
   }

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

   d_object_timers->t_cluster_wrapup->stop();

   if (d_barrier_and_time) {
      d_object_timers->t_find_boxes_containing_tags->barrierAndStop();
   }
}



/*
***********************************************************************
* Cluster tags into tiles, but limit tiles to lie within local process
* boundaries, so that no tiles extend past process boundaries.
***********************************************************************
*/
void
TileClustering::clusterWithinProcessBoundaries(
   hier::BoxLevel &new_box_level,
   hier::Connector &tag_to_tile,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const hier::BoxContainer& bound_boxes,
   int tag_data_index,
   int tag_val)
{
   d_object_timers->t_cluster->start();

   // Determine max number of tiles any local patch can generate.
   int max_tiles_for_any_patch = 0;
   for ( int pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {
      hier::Box coarsened_box = tag_level->getPatch(pi)->getBox();
      coarsened_box.coarsen(d_box_size);
      hier::IntVector number_tiles = coarsened_box.numberCells();
      number_tiles *= 3; // Possible merging of smaller tiles on either side of it.
      max_tiles_for_any_patch = tbox::MathUtilities<int>::Max(
         max_tiles_for_any_patch, number_tiles.getProduct() );
   }

   hier::Connector &tiles_to_tag = tag_to_tile.getTranspose();

   /*
    * Generate new_box_level and Connectors
    */
#pragma omp parallel if ( tag_level->getLocalNumberOfPatches() > 4*omp_get_max_threads() )
#pragma omp for schedule(dynamic)
   for ( int pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {

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
            findTilesContainingTags( tiles, *tag_data, tag_val,
                                     pi*max_tiles_for_any_patch );

         if (d_print_steps) {
            tbox::plog << "Tile Clustering generated " << tiles.size()
                       << " clusters from " << num_coarse_tags
                       << " in patch " << patch.getBox().getBoxId() << '\n';
         }

         TBOX_omp_set_lock(&l_outputs);
         for ( hier::BoxContainer::iterator bi=tiles.begin(); bi!=tiles.end(); ++bi ) {
            new_box_level.addBoxWithoutUpdate(*bi);
            tiles_to_tag.insertLocalNeighbor( patch_box, bi->getBoxId() );
            tag_to_tile.insertLocalNeighbor( *bi, patch_box.getBoxId() );
         }
         TBOX_omp_unset_lock(&l_outputs);

      } // Patch is in bounding box

   } // Loop through tag level

   new_box_level.finalize();

   d_object_timers->t_cluster->stop();

   return;
}



/*
***********************************************************************
* Cluster tags into whole tiles.  The tiles are not cut up, even where
* they cross process boundaries or level boundaries.
*
* This requires tag<==>tag to have a width of at least the tile size,
* but it doesn't require any communication.
*
* Any tile with a local tag will be added locally.  If tile crosses
* patch boundaries, this method does not detect overlaps with remote
* tag boxes.  If the tile has tags on multiple patches, the tile will
* be duplicated (with different BoxIds).  Missing edges and tile
* duplication must be corrected by a postprocessing step after this
* method.
*
* This method is local.
***********************************************************************
*/
void
TileClustering::clusterWholeTiles(
   hier::BoxLevel &tile_box_level,
   boost::shared_ptr<hier::Connector> &tag_to_tile,
   int &local_tiles_have_remote_extent,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const hier::BoxContainer& bound_boxes,
   int tag_data_index,
   int tag_val)
{
   d_object_timers->t_cluster->start();

   const hier::BoxLevel &tag_box_level = *tag_level->getBoxLevel();
   const hier::Connector &tag_to_tag = tag_box_level.findConnector(
      tag_box_level, d_box_size, hier::CONNECTOR_IMPLICIT_CREATION_RULE, true);

   hier::BoxContainer visible_tag_boxes;
   tag_to_tag.getLocalNeighbors(visible_tag_boxes);

   hier::Connector &tiles_to_tag = tag_to_tile->getTranspose();

   hier::LocalId last_used_local_id(-1);

   /*
    * Generate tile_box_level and Connectors
    */

   local_tiles_have_remote_extent = 0;

   for ( int pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {

      hier::Patch &patch = *tag_level->getPatch(pi);
      const hier::Box &patch_box = patch.getBox();
      const hier::BlockId &block_id = patch_box.getBlockId();

      TBOX_ASSERT( bound_boxes.begin(block_id) != bound_boxes.end(block_id) );
      const hier::Box &bounding_box = *bound_boxes.begin(block_id);

      if ( !patch.getBox().intersects(bounding_box) ) { continue; }

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch.getPatchData(tag_data_index), boost::detail::dynamic_cast_tag());

      boost::shared_ptr<pdat::CellData<int> > coarsened_tag_data =
         makeCoarsenedTagData(*tag_data, tag_val);

      const hier::Box &coarsened_tag_box = coarsened_tag_data->getBox();
      const int num_coarse_cells = coarsened_tag_box.size();

      for ( int coarse_offset=0; coarse_offset<num_coarse_cells; ++coarse_offset ) {

         const pdat::CellIndex coarse_cell_index(coarsened_tag_box.index(coarse_offset));

         if ( (*coarsened_tag_data)(coarse_cell_index) == tag_val ) {

            hier::Box whole_tile(coarse_cell_index,coarse_cell_index,
                                 coarsened_tag_box.getBlockId(),
                                 ++last_used_local_id,
                                 tag_box_level.getMPI().getRank());
            whole_tile.refine(d_box_size);

            hier::BoxContainer overlapping_tag_boxes(whole_tile);
            if ( !patch_box.contains(whole_tile) ) {
               visible_tag_boxes.findOverlapBoxes( overlapping_tag_boxes, whole_tile );
            }

            tile_box_level.addBox(whole_tile);
            for ( hier::BoxContainer::iterator bi=overlapping_tag_boxes.begin();
                  bi!=overlapping_tag_boxes.end(); ++bi ) {
               if ( bi->getOwnerRank() == tag_box_level.getMPI().getRank() ) {
                  tag_to_tile->insertLocalNeighbor( whole_tile, bi->getBoxId() );
               }
               tiles_to_tag.insertLocalNeighbor( *bi, whole_tile.getBoxId() );
            }

            if ( !local_tiles_have_remote_extent ) {
               for ( hier::BoxContainer::const_iterator bi=overlapping_tag_boxes.begin();
                     bi!=overlapping_tag_boxes.end(); ++ bi ) {
                  if ( bi->getOwnerRank() != patch_box.getOwnerRank() ) {
                     local_tiles_have_remote_extent = true;
                     break;
                  }
               }
            }

         } // coarse_cell_index has local tag.

      } // Loop through coarsened tag cells

   } // Loop through tag level

   tile_box_level.finalize();

   d_object_timers->t_cluster->stop();

   return;
}



/*
***********************************************************************
* The method clusterWholeTiles may generate duplicate tiles when a
* tile crosses any tag box boundary.  This methods removes the
* duplicates and keep only the one on the process with the greatest
* overlap.
***********************************************************************
*/
void
TileClustering::removeDuplicateTiles(
   hier::BoxLevel &tile_box_level,
   hier::Connector &tag_to_tile)
{

   const hier::BoxLevel &tag_box_level = tag_to_tile.getBase();

   hier::Connector &tiles_to_tag = tag_to_tile.getTranspose();


   /*
    * Get tiles_crossing_patch_boundaries.  These are
    * - local tiles with multiple tag neighbors, and
    * - remote tiles visible locally
    */

   hier::BoxContainer visible_tiles;
   tag_to_tile.getLocalNeighbors(visible_tiles);
   hier::BoxContainer tiles_crossing_patch_boundaries;

   for ( hier::BoxContainer::const_iterator bi=visible_tiles.begin();
         bi!=visible_tiles.end(); ++bi ) {
      if ( bi->getOwnerRank() != tiles_to_tag.getMPI().getRank() ||
           tiles_to_tag.numLocalNeighbors(bi->getBoxId()) > 1 ) {
         tiles_crossing_patch_boundaries.pushBack(*bi);
      }
   }

   visible_tiles.clear();


   // Chosen tiles among all the duplicate tiles.
   std::vector<hier::Box> chosen_tiles;
   // Map from a duplicate tile to (index of) the chosen tile.
   std::map<hier::BoxId,size_t> changes;

   /*
    * Look for duplicate_tiles and choose one from each group of
    * duplicates.  Chose the first among the duplicate tiles.  Because
    * duplicate_tiles are sorted, we are arbitrarily choosing the
    * first in sorted order.  An alternative is to choose one from the
    * process with most overlap.
    */
   while ( !tiles_crossing_patch_boundaries.isEmpty() ) {

      hier::BoxContainer duplicate_tiles(tiles_crossing_patch_boundaries.front(), true);
      tiles_crossing_patch_boundaries.popFront();

      /*
       * If the O(N) search for duplicates is too slow, it can be
       * replaced by ordering the tiles by the lower corner and doing
       * an O(lg N) search.
       */
      for ( hier::BoxContainer::iterator bi=tiles_crossing_patch_boundaries.begin();
            bi!=tiles_crossing_patch_boundaries.end(); /* incremented in loop */ ) {
         if ( bi->lower() == duplicate_tiles.front().lower() ) {
            TBOX_ASSERT( bi->upper() == duplicate_tiles.front().upper() );
            duplicate_tiles.insert(*bi);
            changes[bi->getBoxId()] = chosen_tiles.size();
            tiles_crossing_patch_boundaries.erase(bi++);
         }
         else {
            ++bi;
         }
      }

      chosen_tiles.push_back( *duplicate_tiles.begin() );

   }

   /*
    * Change tile_box_level and Connectors based on the change map.
    */
   hier::BoxLevel tmp_tiles_box_level(
      tile_box_level.getRefinementRatio(),
      tile_box_level.getGridGeometry(),
      tile_box_level.getMPI() );
   hier::Connector tmp_tag_to_tile( tag_box_level,
                                     tmp_tiles_box_level,
                                     tag_to_tile.getConnectorWidth() );
   hier::Connector tmp_tiles_to_tag( tmp_tiles_box_level,
                                     tag_box_level,
                                     tiles_to_tag.getConnectorWidth() );

   for ( hier::BoxContainer::const_iterator bi=tile_box_level.getBoxes().begin();
         bi!=tile_box_level.getBoxes().end(); ++bi ) {

      const hier::Box &possibly_duplicated_tile(*bi);

      std::map<hier::BoxId,size_t>::const_iterator chosen_box_itr = changes.find(possibly_duplicated_tile.getBoxId());

      const hier::Box &unique_tile =
         chosen_box_itr == changes.end() ? possibly_duplicated_tile : chosen_tiles[chosen_box_itr->second];

      if ( unique_tile.getOwnerRank() == tmp_tiles_box_level.getMPI().getRank() ) {
         tmp_tiles_box_level.addBoxWithoutUpdate(unique_tile);
         hier::Connector::ConstNeighborhoodIterator neighborhood =
            tiles_to_tag.find(bi->getBoxId());
         for ( hier::Connector::ConstNeighborIterator na=tiles_to_tag.begin(neighborhood);
               na!=tiles_to_tag.end(neighborhood); ++na ) {
            tmp_tiles_to_tag.insertLocalNeighbor( *na, unique_tile.getBoxId() );
         }
      }

   }

   for ( hier::Connector::ConstNeighborhoodIterator ni=tag_to_tile.begin();
         ni!=tag_to_tile.end(); ++ni ) {

      const hier::BoxId &tag_box_id = *ni;
      const hier::Box &tag_box = *tag_box_level.getBoxStrict(*ni);

      for ( hier::Connector::ConstNeighborIterator na=tag_to_tile.begin(ni);
            na!=tag_to_tile.end(ni); ++na ) {

         std::map<hier::BoxId,size_t>::const_iterator chosen_box_itr = changes.find(na->getBoxId());

         if ( chosen_box_itr != changes.end() ) {
            tmp_tag_to_tile.insertLocalNeighbor(chosen_tiles[chosen_box_itr->second], *ni);
         }
         else {
            tmp_tag_to_tile.insertLocalNeighbor(*na, *ni);
         }

      }
   }

   tile_box_level.finalize();

   hier::BoxLevel::swap(tile_box_level, tmp_tiles_box_level);
   tag_to_tile = tmp_tag_to_tile;
   tiles_to_tag = tmp_tiles_to_tag;

{
   // There should be no overlaps.
   hier::BoxContainer visible_tiles;
   tag_to_tile.getLocalNeighbors(visible_tiles);
   visible_tiles.makeTree( tile_box_level.getGridGeometry().get() );
   for ( hier::BoxContainer::const_iterator bi=visible_tiles.begin();
         bi!=visible_tiles.end(); ++bi ) {
      hier::BoxContainer overlaps;
      visible_tiles.findOverlapBoxes( overlaps, *bi );
      assert( overlaps.size() == 1 );
      assert( overlaps.front().isIdEqual(*bi) );
      assert( overlaps.front().isSpatiallyEqual(*bi) );
   }
}

   return;
}



/*
***********************************************************************
***********************************************************************
*/
void
TileClustering::detectSemilocalEdges(
   boost::shared_ptr<hier::Connector> &tag_to_tile )
{
   const hier::BoxLevel &tag_box_level = tag_to_tile->getBase();
   const hier::Connector &tag_to_tag = tag_box_level.findConnector(
      tag_box_level, d_box_size, hier::CONNECTOR_IMPLICIT_CREATION_RULE, true);

   /*
    * Bridge tag<==>tag<==>new to get the complete tag<==>new.
    * Currently, tag boxes don't know about any remote new boxes that
    * may overlap them.
    *
    * Note: The bridge is convenient but overkill.  We can get same
    * information with a lighter weight communication and no
    * communication at all when no tiles cross process boundaries.
    */
   d_oca.bridge( tag_to_tile,
                 tag_to_tag,
                 hier::Connector(*tag_to_tile) /* verify that copying is really needed */,
                 hier::IntVector::getZero(d_dim),
                 true /* compute transpose */ );

   return;
}



/*
***********************************************************************
* Cluster tags into tiles, but limit tiles to lie within the tag
* level, so that no tiles extend past the tag level.
***********************************************************************
*/
void
TileClustering::clusterWithinLevelBoundaries(
   hier::BoxLevel &new_box_level,
   hier::Connector &tag_to_tile,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const hier::BoxContainer& bound_boxes,
   int tag_data_index,
   int tag_val)
{
   TBOX_ERROR("TileClustering::clusterWithinLevelBoundaries is an incomplete WIP.");

   d_object_timers->t_cluster->start();

   // Determine max number of tiles any local patch can generate.
   int max_tiles_for_any_patch = 0;
   for ( int pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {
      hier::Box coarsened_box = tag_level->getPatch(pi)->getBox();
      coarsened_box.coarsen(d_box_size);
      hier::IntVector number_tiles = coarsened_box.numberCells();
      number_tiles *= 3; // Possible merging of smaller tiles on either side of it.
      max_tiles_for_any_patch = tbox::MathUtilities<int>::Max(
         max_tiles_for_any_patch, number_tiles.getProduct() );
   }

   hier::Connector &tiles_to_tag = tag_to_tile.getTranspose();

   const hier::BoxLevel &tag_box_level = tiles_to_tag.getHead();
   const hier::BoxContainer &local_tag_boxes = tag_box_level.getBoxes();
   local_tag_boxes.makeTree( tag_box_level.getGridGeometry().get() );

   const hier::Connector &tag_to_tag = tag_box_level.findConnector(
      tag_box_level, d_box_size, hier::CONNECTOR_IMPLICIT_CREATION_RULE, true);

   hier::BoxContainer visible_tag_boxes;
   tag_to_tag.getLocalNeighbors(visible_tag_boxes);

   /*
    * Collected data for communication.
    */
   std::set<int> incoming_ranks;
   std::map<int,boost::shared_ptr<tbox::MessageStream> > outgoing_streams;

   hier::LocalId last_used_local_id(-1);

   /*
    * Generate new_box_level and Connectors
    */
   for ( int pi=0; pi<tag_level->getLocalNumberOfPatches(); ++pi ) {

      hier::Patch &patch = *tag_level->getPatch(pi);
      const hier::Box &patch_box = patch.getBox();
      const hier::BlockId &block_id = patch_box.getBlockId();

      TBOX_ASSERT( bound_boxes.begin(block_id) != bound_boxes.end(block_id) );
      const hier::Box &bounding_box = *bound_boxes.begin(block_id);

      if ( !patch.getBox().intersects(bounding_box) ) { continue; }

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch.getPatchData(tag_data_index), boost::detail::dynamic_cast_tag());

      boost::shared_ptr<pdat::CellData<int> > coarsened_tag_data =
         makeCoarsenedTagData(*tag_data, tag_val);

      const hier::Box &coarsened_tag_box = coarsened_tag_data->getBox();
      const int num_coarse_cells = coarsened_tag_box.size();

      for ( int coarse_offset=0; coarse_offset<num_coarse_cells; ++coarse_offset ) {

         const pdat::CellIndex coarse_cell_index(coarsened_tag_box.index(coarse_offset));
         hier::Box whole_tile(coarse_cell_index,coarse_cell_index,
                              coarsened_tag_box.getBlockId());
         whole_tile.refine(d_box_size);

         // valid_tile_parts are parts of whole tile inside the tag level.
         hier::BoxContainer valid_tile_parts(whole_tile);
         valid_tile_parts.intersectBoxes(visible_tag_boxes);
         valid_tile_parts.coalesce();

         // local_tile_parts are parts of whole tile inside the local tag boxes.
         hier::BoxContainer local_tile_parts(valid_tile_parts);
         local_tile_parts.intersectBoxes(local_tag_boxes);
         local_tile_parts.coalesce();
         TBOX_ASSERT( !local_tile_parts.isEmpty() );

         if ( local_tile_parts.size() == 1 &&
              whole_tile.isSpatiallyEqual(local_tile_parts.front()) ) {
            // Tile is completely local.  Add the whole tile.
            const hier::LocalId local_id = ++last_used_local_id;
            whole_tile.initialize( whole_tile,
                                   local_id,
                                   coarsened_tag_box.getOwnerRank() );
            new_box_level.addBox(whole_tile);
            tag_to_tile.insertLocalNeighbor( whole_tile, patch_box.getBoxId() );
            tiles_to_tag.insertLocalNeighbor( patch_box, whole_tile.getBoxId() );
         }
         else {
            /*
             * Tile crosses local process boundary.
             * Process with most overlaps with this tile should own it.
             */

            hier::BoxContainer overlapping_tag_boxes;
            visible_tag_boxes.findOverlapBoxes( overlapping_tag_boxes, whole_tile );

            // process_overlap[r] is how much process r overlaps the tile.
            std::map<int,int> process_overlap;
            for ( hier::BoxContainer::iterator bi=overlapping_tag_boxes.begin();
                  bi!=overlapping_tag_boxes.end(); ++bi ) {
               process_overlap[bi->getOwnerRank()] += (*bi * whole_tile).size();
            }

            // owner_overlap corresponds to process with greatest overlap.
            std::map<int,int>::const_iterator owner_overlap = process_overlap.begin();
            for ( std::map<int,int>::const_iterator mi=process_overlap.begin();
                  mi!=process_overlap.end(); ++mi ) {
               if ( mi->second > owner_overlap->second ||
                    ( mi->second == owner_overlap->second &&
                      mi->first < owner_overlap->first ) ) {
                  owner_overlap = mi;
               }
            }

            // If tile is tagged and locally owned, add all valid parts.
            if ( (*coarsened_tag_data)(coarse_cell_index) == tag_val &&
                 owner_overlap->first == tag_box_level.getMPI().getRank() ) {
               for ( hier::BoxContainer::const_iterator bi=valid_tile_parts.begin();
                     bi!=valid_tile_parts.end(); ++bi ) {
                  const hier::LocalId local_id = ++last_used_local_id;
                  const hier::Box valid_tile( *bi,
                                              local_id,
                                              coarsened_tag_box.getOwnerRank() );
                  new_box_level.addBox(valid_tile);
                  tag_to_tile.insertLocalNeighbor( valid_tile, patch_box.getBoxId() );
                  tiles_to_tag.insertLocalNeighbor( patch_box, valid_tile.getBoxId() );
               }
            }

            if ( owner_overlap->first == tag_box_level.getMPI().getRank() ) {
               // Set up receives from processes that have overlaps with tile.
               for ( std::map<int,int>::const_iterator mi=process_overlap.begin();
                     mi!=process_overlap.end(); ++mi ) {
                  incoming_ranks.insert(mi->first);
               }
            }
            else {
               /*
                * Set up sends to processes that have overlaps with
                * tile.  The message contains information other
                * ranks need to set up edges to this tile: the
                * coarse_cell_index, whether local process has a
                * tag in it, and local tag boxes overlapping it.
                */
               tbox::MessageStream tmpstream;
               tmpstream.pack( &coarse_cell_index(0), coarse_cell_index.getDim().getValue() );
               tmpstream << static_cast<int>((*coarsened_tag_data)(coarse_cell_index) == tag_val);
               tmpstream << static_cast<int>(overlapping_tag_boxes.size());
               for ( hier::BoxContainer::iterator bi=overlapping_tag_boxes.begin();
                     bi!=overlapping_tag_boxes.end(); ++bi ) {
                  bi->putToMessageStream(tmpstream);
               }
               for ( std::map<int,int>::const_iterator mi=process_overlap.begin();
                     mi!=process_overlap.end(); ++mi ) {
                  boost::shared_ptr<tbox::MessageStream> &mstream = outgoing_streams[mi->first];
                  if ( !mstream ) {
                     mstream.reset(new tbox::MessageStream);
                  }
                  mstream->pack(tmpstream);
               }
            }
         }

      } // Loop through coarsened tag cells

   } // Loop through tag level


   /*
    * Send and receive data about tiles that cross process boundaries.
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<char> *comm_recvs = new tbox::AsyncCommPeer<char>[incoming_ranks.size()];
   tbox::AsyncCommPeer<char> *comm_sends = new tbox::AsyncCommPeer<char>[outgoing_streams.size()];
   size_t counter = 0;
   for ( std::set<int>::const_iterator si=incoming_ranks.begin();
         si!=incoming_ranks.end(); ++si ) {
      comm_recvs[counter].initialize(&comm_stage);
      comm_recvs[counter].setMPI(tag_box_level.getMPI());
      comm_recvs[counter].setPeerRank(*si);
      comm_recvs[counter].setMPITag(s_primary_mpi_tag, s_secondary_mpi_tag);
      comm_recvs[counter].limitFirstDataLength(s_first_data_length);
      comm_recvs[counter].beginRecv();
      if ( comm_recvs[counter].isDone() ) {
         comm_recvs[counter].pushToCompletionQueue();
      }
      ++counter;
   }
   counter = 0;
   for ( std::map<int,boost::shared_ptr<tbox::MessageStream> >::const_iterator mi=outgoing_streams.begin();
         mi!=outgoing_streams.end(); ++mi ) {
      comm_sends[counter].initialize(&comm_stage);
      comm_sends[counter].setMPI(tag_box_level.getMPI());
      comm_sends[counter].setPeerRank(mi->first);
      comm_sends[counter].setMPITag(s_primary_mpi_tag, s_secondary_mpi_tag);
      comm_sends[counter].limitFirstDataLength(s_first_data_length);
      comm_sends[counter].beginSend(
         static_cast<const char*>(mi->second->getBufferStart()),
         mi->second->getCurrentSize());
      if ( comm_sends[counter].isDone() ) {
         comm_sends[counter].pushToCompletionQueue();
      }
      ++counter;
   }

   while ( comm_stage.numberOfCompletedMembers() > 0 ||
           comm_stage.advanceSome() ) {

      tbox::AsyncCommPeer<char>* completed_comm =
         CPP_CAST<tbox::AsyncCommPeer<char> *>(comm_stage.popCompletionQueue());

      TBOX_ASSERT(completed_comm != 0);
      TBOX_ASSERT(completed_comm->isDone());
      if ( completed_comm->isReceiver() ) {

         const int sender = completed_comm->getPeerRank();

         tbox::MessageStream incoming_stream(
            completed_comm->getRecvSize() * sizeof(char),
            tbox::MessageStream::Read,
            completed_comm->getRecvData(),
            false /* don't use deep copy */ );

         while ( incoming_stream.canCopyOut(1) ) {

            pdat::CellIndex coarse_cell_index(d_dim);
            incoming_stream.unpack( &coarse_cell_index(0), coarse_cell_index.getDim().getValue() );

            int has_remote_tag = -1, num_remote_overlapping_tag_boxes = -1;
            incoming_stream >> has_remote_tag >> num_remote_overlapping_tag_boxes;

            for ( int i=0; i<num_remote_overlapping_tag_boxes; ++i ) {
               hier::Box remote_overlapping_tag_box(d_dim);
               remote_overlapping_tag_box.getFromMessageStream(incoming_stream);
            }

            if ( has_remote_tag ) {
            }
         }

         completed_comm->clearRecvData();
      } else {
         // No further action required for completed send.
      }
   }

   new_box_level.finalize();

   d_object_timers->t_cluster->stop();

   return;
}


/*
 ***********************************************************************
 ***********************************************************************
 */
int
TileClustering::findTilesContainingTags(
   hier::BoxContainer &tiles,
   const pdat::CellData<int> &tag_data,
   int tag_val,
                                        int first_tile_index)
{
   tiles.clear();
   tiles.unorder();

   hier::Box coarsened_box(tag_data.getBox());
   coarsened_box.coarsen(d_box_size);

   const int num_coarse_cells = coarsened_box.size();

#pragma omp parallel
#pragma omp for schedule(dynamic)
   for ( int coarse_offset=0; coarse_offset<num_coarse_cells; ++coarse_offset ) {
      // if (coarse_offset==0) tbox::plog << "Inner loop has " << TBOX_omp_get_num_threads() << " threads over " << num_coarse_cells << " coarse cells." << std::endl;
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
            hier::LocalId local_id(first_tile_index + coarse_offset);
            if ( local_id < hier::LocalId::getZero() ) {
               TBOX_ERROR("TileClustering code cannot compute a valid non-zero\n"
                          <<"LocalId for a tile.\n");
            }

            tile_box.initialize( tile_box,
                                 local_id,
                                 coarsened_box.getOwnerRank() );
            TBOX_omp_set_lock(&l_interm);
            tiles.pushBack(tile_box);
            TBOX_omp_unset_lock(&l_interm);

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

   size_t coarse_tag_count = 0;

   const int num_coarse_cells = coarsened_box.size();
#pragma omp parallel
#pragma omp for schedule(dynamic)
   for ( int offset=0; offset<num_coarse_cells; ++offset ) {
      const pdat::CellIndex coarse_cell_index(coarsened_box.index(offset));

      hier::Box fine_cells_box(coarse_cell_index,coarse_cell_index,coarsened_box.getBlockId());
      fine_cells_box.refine(d_box_size);
      fine_cells_box *= tag_data.getBox();

      pdat::CellIterator finecend(pdat::CellGeometry::end(fine_cells_box));
      for ( pdat::CellIterator fineci(pdat::CellGeometry::begin(fine_cells_box));
            fineci!=finecend; ++fineci ) {
         if ( tag_data(*fineci) == tag_val ) {
            (*coarsened_tag_data)(coarse_cell_index) = tag_val;
            ++coarse_tag_count;
            break;
         }
      }
   }
   if (d_print_steps) {
      tbox::plog << "TileClustering coarsened box " << tag_data.getBox()
                 << " to " << coarsened_box
                 << " (" << coarse_tag_count << " tags).\n";
   }

   return coarsened_tag_data;
}


/*
 ***********************************************************************
 * This method currently does no communication, assuming tiles don't
 * cross process boundaries on the tag level.  Now, with whole tiles,
 * this should be modified so it does the appropriate communication.
 ***********************************************************************
 */
void
TileClustering::coalesceClusters(
   hier::BoxLevel &tile_box_level,
   boost::shared_ptr<hier::Connector> &tag_to_tile,
   int tiles_have_remote_extent)
{
   /*
    * Try to coalesce the boxes in tile_box_level.
    */
   std::vector<hier::Box> box_vector;
   if (!tile_box_level.getBoxes().isEmpty()) {

      d_object_timers->t_coalesce->start();

      hier::LocalId local_id(0);

      const int nblocks = tile_box_level.getGridGeometry()->getNumberBlocks();

#pragma omp parallel
#pragma omp for schedule(dynamic)
      for (int b = 0; b < nblocks; ++b) {
         hier::BlockId block_id(b);

         hier::BoxContainer block_boxes(tile_box_level.getBoxes(), block_id);

         if (!block_boxes.isEmpty()) {
            block_boxes.unorder();
            block_boxes.coalesce();
            TBOX_omp_set_lock(&l_outputs);
            box_vector.insert(box_vector.end(), block_boxes.begin(), block_boxes.end() );
            TBOX_omp_unset_lock(&l_outputs);
         }
      }

      d_object_timers->t_coalesce->stop();

   }

   if ( d_print_steps ) {
      tbox::plog << "TileClustering coalesced " << tile_box_level.getLocalNumberOfBoxes()
                 << " tiles into " << box_vector.size() << "\n";
   }

   if ( box_vector.size() != tile_box_level.getLocalNumberOfBoxes() ) {

      d_object_timers->t_coalesce_adjustment->start();

      /*
       * Coalesce changed the tiles, so rebuild tile_box_level and
       * Connectors.
       */
      const hier::IntVector &zero_vector = hier::IntVector::getZero(d_dim);
      tile_box_level.initialize( tile_box_level.getRefinementRatio(),
                                tile_box_level.getGridGeometry(),
                                tile_box_level.getMPI() );
      tag_to_tile.reset( new hier::Connector( tag_to_tile->getBase(),
                                             tile_box_level,
                                             zero_vector ) );
      hier::Connector *tiles_to_tag = new hier::Connector( tile_box_level,
                                                         tag_to_tile->getBase(),
                                                         zero_vector );
      tag_to_tile->setTranspose(tiles_to_tag, true);

      const hier::BoxContainer &tag_boxes = tag_to_tile->getBase().getBoxes();
      tag_boxes.makeTree( tag_to_tile->getBase().getGridGeometry().get() );

      /*
       * Assign ids to coalesced boxes, add to BoxLevel and add
       * tile--->tag edges.
       */
      const int rank = tile_box_level.getMPI().getRank();
#pragma omp parallel
#pragma omp for schedule(dynamic)
      for ( size_t i=0; i<box_vector.size(); ++i ) {

         box_vector[i].setId(hier::BoxId(hier::LocalId(static_cast<int>(i)),rank));

         hier::BoxContainer tmp_overlap_boxes;
         tag_boxes.findOverlapBoxes(tmp_overlap_boxes,
                                    box_vector[i],
                                    tag_to_tile->getBase().getRefinementRatio() );

         TBOX_omp_set_lock(&l_outputs);
         tile_box_level.addBox(box_vector[i]);
         tiles_to_tag->insertNeighbors( tmp_overlap_boxes, box_vector[i].getBoxId() );
         TBOX_omp_unset_lock(&l_outputs);

      }
      tile_box_level.finalize();

      /*
       * Add tag--->tile edges.
       */
      hier::BoxContainer tiles;
      for ( size_t i=0; i<box_vector.size(); ++i ) tiles.pushBack(box_vector[i]);
      tiles.makeTree( tile_box_level.getGridGeometry().get() );
      std::vector<hier::Box> real_box_vector, periodic_image_box_vector;
      tag_boxes.separatePeriodicImages( real_box_vector, periodic_image_box_vector );

      for ( size_t ib=0; ib<real_box_vector.size(); ++ib ) {

         hier::BoxContainer tmp_overlap_boxes;
         tiles.findOverlapBoxes(tmp_overlap_boxes,
                                real_box_vector[ib],
                                tag_to_tile->getBase().getRefinementRatio() );

         TBOX_omp_set_lock(&l_outputs);
         tag_to_tile->insertNeighbors( tmp_overlap_boxes,
                                       real_box_vector[ib].getBoxId() );
         TBOX_omp_unset_lock(&l_outputs);
      }

      d_object_timers->t_coalesce_adjustment->stop();

   }

   return;
}


/*
 ***********************************************************************
 ***********************************************************************
 */
void
TileClustering::setTimerPrefix(
   const std::string& timer_prefix)
{
   std::map<std::string, TimerStruct>::iterator ti(
      s_static_timers.find(timer_prefix));

   if (ti != s_static_timers.end()) {
      d_object_timers = &(ti->second);
   } else {

      d_object_timers = &s_static_timers[timer_prefix];

      tbox::TimerManager *tm = tbox::TimerManager::getManager();

      d_object_timers->t_find_boxes_containing_tags = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()");
      d_object_timers->t_cluster = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()_cluster");
      d_object_timers->t_coalesce = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()_coalesce");
      d_object_timers->t_coalesce_adjustment = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()_coalesce_adjustment");
      d_object_timers->t_cluster_setup = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()_setup");
      d_object_timers->t_cluster_wrapup = tm->getTimer(
         timer_prefix + "::findBoxesContainingTags()_wrapup");
      d_object_timers->t_global_reductions = tm->getTimer(
         timer_prefix + "::global_reductions");

   }

}

}
}
#endif
