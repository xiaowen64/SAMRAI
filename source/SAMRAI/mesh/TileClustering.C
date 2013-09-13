/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2013 Lawrence Livermore National Security, LLC
 * Description:   Tile clustering algorithm.
 *
 ************************************************************************/
#ifndef included_mesh_TileClustering_C
#define included_mesh_TileClustering_C

#include <stdlib.h>

#include "SAMRAI/mesh/TileClustering.h"

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
 ************************************************************************
 */
TileClustering::TileClustering(
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_box_size(hier::IntVector(d_dim, 8)),
   d_allow_remote_tile_extent(true),
   d_coalesce_boxes(true),
   d_debug_checks(false),
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
   d_mca.setTimerPrefix(s_default_timer_prefix);
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

      d_debug_checks =
         input_db->getBoolWithDefault("DEV_debug_checks",
            d_debug_checks);
   }
}

/*
 ************************************************************************
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

      shearTilesAtBlockBoundaries( *new_box_level, *tag_to_new );

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
      tbox::plog << "TileClustering summary:\n"
                 << "\tClustered with box_size = " << d_box_size << '\n';

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

   if ( d_debug_checks ) {

      tag_to_new->assertConsistencyWithBase();
      tag_to_new->assertConsistencyWithHead();
      tag_to_new->assertOverlapCorrectness();
      tag_to_new->getTranspose().assertConsistencyWithBase();
      tag_to_new->getTranspose().assertConsistencyWithHead();
      tag_to_new->getTranspose().assertOverlapCorrectness();
      tag_to_new->assertTransposeCorrectness(tag_to_new->getTranspose());

      // There should be no overlaps.
      hier::BoxContainer visible_tiles(true);
      tag_to_new->getLocalNeighbors(visible_tiles);
      visible_tiles.makeTree( tag_to_new->getBase().getGridGeometry().get() );
      for ( hier::BoxContainer::const_iterator bi=visible_tiles.begin();
            bi!=visible_tiles.end(); ++bi ) {
         const hier::Box &tile = *bi;
         hier::BoxContainer overlaps;
         visible_tiles.findOverlapBoxes( overlaps, tile );
         TBOX_ASSERT( overlaps.size() == 1 );
         TBOX_ASSERT( overlaps.front().isIdEqual(tile) );
         TBOX_ASSERT( overlaps.front().isSpatiallyEqual(tile) );
      }

   }


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

   hier::Connector &tile_to_tag = tag_to_tile.getTranspose();

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
            tile_to_tag.insertLocalNeighbor( patch_box, bi->getBoxId() );
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
* be duplicated (with different BoxIds).  tile--->tag will be complete
* (using info from tag--->tag), but tag--->tile will have missing
* edges.  Missing edges and tile duplication must be corrected by a
* postprocessing step after this method.  See detectSemilocalEdges()
* and removeDuplicateTiles().
*
* This method does no communication.
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
   visible_tag_boxes.makeTree(tag_box_level.getGridGeometry().get());

   hier::Connector &tile_to_tag = tag_to_tile->getTranspose();

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
                                 patch_box.getBlockId(),
                                 ++last_used_local_id,
                                 patch_box.getOwnerRank());
            whole_tile.refine(d_box_size);

            tile_box_level.addBox(whole_tile);

            /*
             * Compare whole_tile with overlapping_tag_boxes to determine
             * edges and remote extents.
             */
            hier::BoxContainer overlapping_tag_boxes;
            visible_tag_boxes.findOverlapBoxes( overlapping_tag_boxes,
                                                whole_tile,
                                                tag_box_level.getRefinementRatio() );

            for ( hier::BoxContainer::iterator bi=overlapping_tag_boxes.begin();
                  bi!=overlapping_tag_boxes.end(); ++bi ) {

               tile_to_tag.insertLocalNeighbor( *bi, whole_tile.getBoxId() );
               if ( bi->getOwnerRank() == whole_tile.getOwnerRank() ) {
                  tag_to_tile->insertLocalNeighbor( whole_tile, bi->getBoxId() );
               }

               local_tiles_have_remote_extent |= bi->getOwnerRank() != whole_tile.getOwnerRank();
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
* Methods clusterWholeTiles() and detectSemilocalEdges(), preceding
* this one in execution order, may generate duplicate tiles when a
* tile crosses any tag box boundary.  This methods removes the
* duplicates.
*
* This method does no communication.
***********************************************************************
*/
void
TileClustering::removeDuplicateTiles(
   hier::BoxLevel &tile_box_level,
   hier::Connector &tag_to_tile)
{

   hier::Connector &tile_to_tag = tag_to_tile.getTranspose();

   /*
    * Get tiles_crossing_patch_boundaries.  These are
    * - local tiles with multiple tag neighbors, and
    * - remote tiles visible locally (found by detectSemilocalEdges)
    *
    * The latter may not have multiple local tag neighbors.
    */

   hier::BoxContainer visible_tiles(true);
   tag_to_tile.getLocalNeighbors(visible_tiles);

   hier::BoxContainer tiles_crossing_patch_boundaries;
   for ( hier::BoxContainer::const_iterator ti=visible_tiles.begin();
         ti!=visible_tiles.end(); ++ti ) {
      const hier::Box &tile(*ti);
      if ( tile.getOwnerRank() != tile_to_tag.getMPI().getRank() ||
           tile_to_tag.numLocalNeighbors(tile.getBoxId()) > 1 ) {
         tiles_crossing_patch_boundaries.pushBack(tile);
      }
   }
   visible_tiles.clear(); // No longer needed.


   // Chosen tiles among all the duplicate tiles.
   std::vector<hier::Box> chosen_tiles;
   // Map from a duplicate tile to (index of) the chosen tile.
   std::map<hier::BoxId,size_t> changes;

   /*
    * Look for similar_tiles (tiles with same extents) and choose one
    * from each group of similars.  Chose the first in the group.
    * Because similar_tiles are sorted, we are arbitrarily choosing
    * the first in sorted order.  An alternative is to choose one from
    * the process with most overlap.
    */
   while ( !tiles_crossing_patch_boundaries.isEmpty() ) {

      hier::BoxContainer similar_tiles(tiles_crossing_patch_boundaries.front(), true);
      tiles_crossing_patch_boundaries.popFront();

      /*
       * Search for tiles with same extent as the first one in
       * similar_tiles.  If the O(N) search for duplicates is too
       * slow, it can be replaced by ordering the tiles by the lower
       * corner and doing an O(lg N) search.
       */
      for ( hier::BoxContainer::iterator bi=tiles_crossing_patch_boundaries.begin();
            bi!=tiles_crossing_patch_boundaries.end(); /* incremented in loop */ ) {
         if ( bi->lower() == similar_tiles.front().lower() ) {
            TBOX_ASSERT( bi->upper() == similar_tiles.front().upper() );
            similar_tiles.insert(*bi);
            changes[bi->getBoxId()] = chosen_tiles.size();
            tiles_crossing_patch_boundaries.erase(bi++);
         }
         else {
            ++bi;
         }
      }

      if ( similar_tiles.size() > 1 ) {
         chosen_tiles.push_back( *similar_tiles.begin() );
      }

   }

   /*
    * Change tile_box_level and Connectors based on the change map.
    */

   for ( hier::BoxContainer::const_iterator tile_itr=tile_box_level.getBoxes().begin();
         tile_itr!=tile_box_level.getBoxes().end(); /* incremented in loop */ ) {

      const hier::Box &possibly_duplicated_tile(*tile_itr);

      std::map<hier::BoxId,size_t>::const_iterator chosen_box_itr =
         changes.find(possibly_duplicated_tile.getBoxId());

      if ( chosen_box_itr != changes.end() &&
           !chosen_tiles[chosen_box_itr->second].isIdEqual(possibly_duplicated_tile) ) {

         const hier::Box &unique_tile = chosen_tiles[chosen_box_itr->second];


         // Add unique_tile if it's local.
         if ( unique_tile.getOwnerRank() == tile_box_level.getMPI().getRank() ) {
            tile_box_level.addBoxWithoutUpdate(unique_tile);
            hier::Connector::ConstNeighborhoodIterator neighborhood =
               tile_to_tag.find(possibly_duplicated_tile.getBoxId());
            for ( hier::Connector::ConstNeighborIterator na=tile_to_tag.begin(neighborhood);
                  na!=tile_to_tag.end(neighborhood); ++na ) {
               tile_to_tag.insertLocalNeighbor( *na, unique_tile.getBoxId() );
            }
         }

         // Remove duplicated tile.
         tile_to_tag.eraseLocalNeighborhood(tile_itr->getBoxId());
         tile_box_level.eraseBoxWithoutUpdate(*(tile_itr++));

      }
      else {
         ++tile_itr;
      }

   }
   tile_box_level.finalize();
   tag_to_tile.setHead(tile_box_level, true);
   tile_to_tag.setBase(tile_box_level, true);
   tile_box_level.deallocateGlobalizedVersion();

   for ( hier::Connector::ConstNeighborhoodIterator ni=tag_to_tile.begin();
         ni!=tag_to_tile.end(); ++ni ) {

      for ( hier::Connector::ConstNeighborIterator na=tag_to_tile.begin(ni);
            na!=tag_to_tile.end(ni); /* incremented in loop */ ) {

         std::map<hier::BoxId,size_t>::const_iterator chosen_box_itr =
            changes.find(na->getBoxId());

         if ( chosen_box_itr != changes.end() ) {
            tag_to_tile.insertLocalNeighbor(chosen_tiles[chosen_box_itr->second], *ni);
            tag_to_tile.eraseNeighbor(*(na++), *ni);
         }
         else {
            ++na;
         }

      }
   }

   return;
}



/*
***********************************************************************
* Detect semilocal edges missing from the outputs of
* clusterWholeTiles().
*
* Methods clusterWholeTiles(), preceding this one in execution order,
* may generate tile extending past local tag boxes.  It doesn't
* generate any semilocal edges because it is a completely local
* algorithm.  This method generates the missing semilocal edges.
*
* On entry, tile_to_tag must be complete, but tag_to_tile may be
* missing semilocal edges.  On exit, both would be complete overlap
* Connectors.
*
* This method does a bridge communication.
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
    * Note: Bridging is convenient but overkill.  We can get same
    * information with much lighter weight communication and no
    * communication at all where no tiles cross process boundaries.
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
* Shear tiles at tag level block boundaries.
* Remove interblock neighbors.
* Fix Connectors.
***********************************************************************
*/
void
TileClustering::shearTilesAtBlockBoundaries(
   hier::BoxLevel &tile_box_level,
   hier::Connector &tag_to_tile)
{
   const hier::BoxContainer &tiles = tile_box_level.getBoxes();
   const hier::BoxLevel &tag_box_level = tag_to_tile.getBase();
   hier::Connector &tile_to_tag = tag_to_tile.getTranspose();

   const boost::shared_ptr<const hier::BaseGridGeometry> &grid_geom =
      tag_box_level.getGridGeometry();

   hier::LocalId last_used_id = tile_box_level.getLastLocalId();

   std::map<hier::BlockId,hier::BoxContainer> domain_blocks;

   // Map from changed tiles to their replacements.
   std::map<hier::LocalId,hier::BoxContainer> changes;

   hier::BoxLevel sheared_tile_box_level( tile_box_level.getRefinementRatio(),
                                          tile_box_level.getGridGeometry() );
   hier::MappingConnector tile_to_sheared( tile_box_level,
                                           sheared_tile_box_level,
                                           hier::IntVector::getZero(d_dim) );

   for ( hier::BoxContainer::const_iterator ti=tiles.begin();
         ti!=tiles.end(); ++ti ) {

      const hier::Box tile = *ti;

      hier::BoxContainer &block_boxes = domain_blocks[tile.getBlockId()];
      if ( block_boxes.isEmpty() ) {
         grid_geom->computePhysicalDomain( block_boxes,
                                           tile_box_level.getRefinementRatio(),
                                           tile.getBlockId() );
         block_boxes.makeTree();
      }

      hier::BoxContainer &inside_block = changes[tile.getLocalId()];
      inside_block.pushFront(tile);
      inside_block.intersectBoxes(block_boxes);

      if ( inside_block.front().isSpatiallyEqual(tile) ) {
         // No change to this tile.
         TBOX_ASSERT( inside_block.size() == 1 );
         TBOX_ASSERT( inside_block.front().isIdEqual(tile) );
         changes.erase(tile.getLocalId());
         sheared_tile_box_level.addBoxWithoutUpdate(tile);
      }
      else {
         inside_block.coalesce();
         hier::BoxContainer tag_neighbors;
         tile_to_tag.getNeighborBoxes( tile.getBoxId(), tag_neighbors );
         for ( hier::BoxContainer::iterator ii=inside_block.begin();
               ii!=inside_block.end(); ++ii ) {
            ii->setLocalId(++last_used_id);
            tile_box_level.addBoxWithoutUpdate(*ii);
            tile_to_sheared.insertNeighbors( inside_block, tile.getBoxId() );
         }
      }

   }

   sheared_tile_box_level.finalize();
   tile_box_level.deallocateGlobalizedVersion();

   d_mca.modify( tag_to_tile,
                 tile_to_sheared,
                 &tile_box_level,
                 &sheared_tile_box_level );

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
 * Coalesce tile clusters and update tag<==>tile.
 *
 * This method uses a modify operation to update tag<==>tile, but if
 * tag<==>tile is local (tiles have no remote extent), the modify
 * operation should properly degenerate to a local (non-communicating)
 * operation.
 ***********************************************************************
 */
void
TileClustering::coalesceClusters(
   hier::BoxLevel &tile_box_level,
   boost::shared_ptr<hier::Connector> &tag_to_tile,
   int tiles_have_remote_extent)
{
   /*
    * If tiles have no remote extent, we can make significant
    * optimizations.
    */
   if ( !tiles_have_remote_extent ) {
      coalesceClusters(tile_box_level, tag_to_tile);
      return;
   }

   /*
    * Coalesce the boxes and give coalesced boxes unique ids.
    */
   const hier::BoxContainer &pre_boxes = tile_box_level.getBoxes();
   hier::BoxContainer post_boxes(pre_boxes.begin(), pre_boxes.end(), false);
   d_object_timers->t_coalesce->start();
   post_boxes.coalesce();
   hier::LocalId last_used_id(tile_box_level.getLastLocalId());
   for ( hier::BoxContainer::iterator bi=post_boxes.begin();
         bi!=post_boxes.end(); ++bi ) {
      bi->setId(hier::BoxId(++last_used_id, tile_box_level.getMPI().getRank()));
   }
   d_object_timers->t_coalesce->stop();


   if ( d_print_steps ) {
      tbox::plog << "TileClustering coalesced "
                 << tile_box_level.getLocalNumberOfBoxes()
                 << " tiles into " << post_boxes.size() << "\n";
   }


   d_object_timers->t_coalesce_adjustment->start();


   /*
    * Build a map that represents the changes from pre- to post-coalesce.
    */
   const hier::IntVector &zero_vector = hier::IntVector::getZero(d_dim);
   hier::BoxLevel tmp_tile_box_level(
      tile_box_level.getRefinementRatio(),
      tile_box_level.getGridGeometry(),
      tile_box_level.getMPI() );

   hier::MappingConnector pre_to_post(
      tile_box_level,
      tmp_tile_box_level,
      zero_vector );

   pre_boxes.makeTree( tile_box_level.getGridGeometry().get() );

   for ( hier::BoxContainer::const_iterator post_itr=post_boxes.begin();
         post_itr!=post_boxes.end(); ++post_itr ) {

      hier::BoxContainer tmp_overlap_boxes;
      pre_boxes.findOverlapBoxes( tmp_overlap_boxes, *post_itr,
                                  tile_box_level.getRefinementRatio() );

      TBOX_ASSERT( !tmp_overlap_boxes.isEmpty() );
      if ( tmp_overlap_boxes.size() == 1 ) {
         // pre- and post-box are the same.  No mapping edge.
         TBOX_ASSERT( tmp_overlap_boxes.front().isSpatiallyEqual(*post_itr) );
         tmp_tile_box_level.addBoxWithoutUpdate(tmp_overlap_boxes.front());
      }
      else {
         // Add coalesced box and edges to it.
         tmp_tile_box_level.addBoxWithoutUpdate(*post_itr);
         for ( hier::BoxContainer::const_iterator pre_itr=tmp_overlap_boxes.begin();
               pre_itr!=tmp_overlap_boxes.end(); ++pre_itr ) {
            TBOX_ASSERT( post_itr->getOwnerRank() == pre_itr->getOwnerRank() );
            pre_to_post.insertLocalNeighbor(*post_itr, pre_itr->getBoxId());
         }
      }

   }

   tmp_tile_box_level.finalize();
   TBOX_ASSERT( pre_to_post.isLocal() );


   /*
    * Apply the modifications.
    */
   if ( d_debug_checks ) {
      d_mca.setSanityCheckMethodPreconditions(true);
      d_mca.setSanityCheckMethodPostconditions(true);
   }
   d_mca.modify( *tag_to_tile,
                 pre_to_post,
                 &tile_box_level,
                 &tmp_tile_box_level );
   d_mca.setSanityCheckMethodPreconditions(false);
   d_mca.setSanityCheckMethodPostconditions(false);

   d_object_timers->t_coalesce_adjustment->stop();

   return;
}


/*
 ***********************************************************************
 * This method does no communication but requires that tiles don't
 * cross process boundaries on the tag level.
 *
 ***********************************************************************
 */
void
TileClustering::coalesceClusters(
   hier::BoxLevel &tile_box_level,
   boost::shared_ptr<hier::Connector> &tag_to_tile)
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

   tile_box_level.deallocateGlobalizedVersion();

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
      hier::Connector *tile_to_tag = new hier::Connector( tile_box_level,
                                                          tag_to_tile->getBase(),
                                                          zero_vector );
      tag_to_tile->setTranspose(tile_to_tag, true);

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
         tile_to_tag->insertNeighbors( tmp_overlap_boxes, box_vector[i].getBoxId() );
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
