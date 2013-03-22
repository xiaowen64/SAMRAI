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
   d_log_cluster_summary(false),
   d_barrier_and_time(false)
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
   const double efficiency_tol,
   const double combine_tol,
   const hier::IntVector& max_gcw)
{
   NULL_USE(efficiency_tol);
   NULL_USE(combine_tol);
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

   if (d_barrier_and_time) {
      t_cluster->start();
   }

   const hier::BoxContainer &local_tag_boxes = tag_box_level.getBoxes();
   local_tag_boxes.makeTree();

   // Generate new_box_level and Connectors
   for ( hier::PatchLevel::iterator pi=tag_level->begin();
         pi!=tag_level->end(); ++pi ) {

      const hier::Box &patch_box = pi->getBox();
      const hier::BlockId &block_id = patch_box.getBlockId();

      TBOX_ASSERT( bound_boxes.begin(block_id) != bound_boxes.end(block_id) );
      const hier::Box &bounding_box = *bound_boxes.begin(block_id);

      if ( pi->getBox().intersects(bounding_box) ) {

         boost::shared_ptr<pdat::CellData<int> > tag_data(
            pi->getPatchData(tag_data_index), boost::detail::dynamic_cast_tag());

         boost::shared_ptr<pdat::CellData<int> > coarsened_tag_data =
            makeCoarsenedTagData(*tag_data, tag_val);

         hier::Box coarsened_box = coarsened_tag_data->getBox();

         pdat::CellIterator ccend(pdat::CellGeometry::end(coarsened_box));
         for ( pdat::CellData<int>::iterator cci(pdat::CellGeometry::begin(coarsened_box));
               cci!=ccend; ++cci ) {

            pdat::CellIndex cindex = *cci;

            if ( (*coarsened_tag_data)(cindex) == tag_val ) {

               /*
                * new_box must nest inside local boxes.  If any part
                * extends outside local_tag_boxes, its overlap with
                * remote tag boxes might not be detected.  Remove
                * those parts and insert the rest into new_box_level.
                * A new_box overlapping non-local tag boxes can appear
                * if the some tag level patch boundaries do not
                * coincide with the tile cuts.
                */

               hier::Box tmp_new_box(cindex, cindex, coarsened_box.getBlockId());
               tmp_new_box.refine(d_box_size);

               hier::BoxContainer tmp_new_boxes(tmp_new_box);
               tmp_new_boxes.intersectBoxes(local_tag_boxes);
               tmp_new_boxes.coalesce();

               for ( hier::BoxContainer::const_iterator bi=tmp_new_boxes.begin();
                     bi!=tmp_new_boxes.end(); ++bi ) {
                  hier::Box new_box( *bi,
                                     new_box_level->getLastLocalId()+1,
                                     new_box_level->getMPI().getRank() );

                  new_box_level->addBoxWithoutUpdate(new_box);
                  hier::BoxContainer::const_iterator new_box_iter =
                     new_box_level->getBoxStrict(new_box);

                  new_to_tag->insertLocalNeighbor( patch_box, new_box_iter->getBoxId() );
                  tag_to_new->insertLocalNeighbor( *new_box_iter, patch_box.getBoxId() );
               }

            }

         }

      }

   }

   new_box_level->finalize();


   if ( d_coalesce_boxes ) {

      /*
       * Try to coalesce the boxes in new_box_level.
       */
      hier::BoxContainer new_boxes(new_box_level->getBoxes());
      new_boxes.unorder();
      new_boxes.coalesce();
      if ( new_boxes.size() != size_t(new_box_level->getLocalNumberOfBoxes()) ) {

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
         int rank = new_box_level->getMPI().getRank();
         for ( hier::BoxContainer::iterator bi=new_boxes.begin();
               bi!=new_boxes.end(); ++bi ) {
            bi->setId(hier::BoxId(++last_used_id,rank));
            new_box_level->addBox(*bi);

            hier::BoxContainer tmp_overlap_boxes;
            tag_boxes.findOverlapBoxes(tmp_overlap_boxes, *bi);

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
            new_boxes.findOverlapBoxes(tmp_overlap_boxes, *bi);

            tag_to_new->insertNeighbors( tmp_overlap_boxes, bi->getBoxId() );
         }

      }

   }

   if (d_barrier_and_time) {
      t_cluster->barrierAndStop();
   }



   /*
    * Get some global parameters.  Do it before logging to prevent
    * the logging flag from having an undue side effect on performance.
    */
   if (d_barrier_and_time) {
      t_global_reductions->barrierAndStart();
   }
   new_box_level->getGlobalNumberOfBoxes();
   new_box_level->getGlobalNumberOfCells();
   if (d_barrier_and_time) {
      t_global_reductions->barrierAndStop();
   }
   for (hier::BoxContainer::const_iterator bi = bound_boxes.begin();
        bi != bound_boxes.end(); ++bi) {
      new_box_level->getGlobalBoundingBox(bi->getBlockId().getBlockValue());
   }

   if (d_log_cluster) {
      tbox::plog << "TileClustering cluster log:\n"
      << "\tNew box_level clustered by TileClustering:\n" << new_box_level->format("",
         2)
      << "\tTileClustering tag_to_new:\n" << tag_to_new->format("", 2)
      << "\tTileClustering new_to_tag:\n" << new_to_tag->format("", 2);
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

   if (d_barrier_and_time) {
      t_find_boxes_containing_tags->barrierAndStop();
   }
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

   pdat::CellIterator finecend(pdat::CellGeometry::end(tag_data.getBox()));
   for ( pdat::CellIterator fineci(pdat::CellGeometry::begin(tag_data.getBox()));
         fineci!=finecend; ++fineci ) {

      if ( tag_data(*fineci) == tag_val ) {
         pdat::CellIndex coarseci =
            pdat::CellIndex( *fineci / d_box_size );
         (*coarsened_tag_data)(coarseci) = tag_val;
      }

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
      getTimer("mesh::TileClustering::cluster");
   t_global_reductions = tbox::TimerManager::getManager()->
      getTimer("mesh::TileClustering::global_reductions");
}

}
}
#endif
