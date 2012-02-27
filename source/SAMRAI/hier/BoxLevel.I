/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of boxes in a box_level of a distributed box graph.
 *
 ************************************************************************/
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/BoxLevelHandle.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"

#include <boost/make_shared.hpp>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
const tbox::Dimension& BoxLevel::getDim() const
{
   return d_ratio.getDim();
}

SAMRAI_INLINE_KEYWORD
bool BoxLevel::isInitialized() const
{
   return d_ratio(0) != 0;
}

SAMRAI_INLINE_KEYWORD
BoxLevel::ParallelState
BoxLevel::getParallelState() const
{
   return d_parallel_state;
}

/*
 ***********************************************************************
 * Detach this object from the handle it has been using.
 *
 * Postcondition: Object that cached the handle would no longer
 * be able to access this BoxLevel from the handle.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::detachMyHandle() const
{
   if (d_handle) {
      d_handle->detachMyBoxLevel();
      d_handle.reset();
   }
}

SAMRAI_INLINE_KEYWORD
const BoxContainer& BoxLevel::getBoxes() const
{
   return d_boxes;
}

SAMRAI_INLINE_KEYWORD
const BoxContainer& BoxLevel::getGlobalBoxes() const
{
   return d_global_boxes;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::clearPersistentOverlapConnectors()
{
   if (d_persistent_overlap_connectors != NULL) {
      d_persistent_overlap_connectors->clear();
   }
}

/*
 ****************************************************************************
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::deallocateGlobalizedVersion() const
{
   if (d_globalized_version != NULL) {
      TBOX_ASSERT(d_globalized_version->getParallelState() == GLOBALIZED);
      delete d_globalized_version;
      d_globalized_version = NULL;
   }
}

/*
 ****************************************************************************
 * Clear certain data that would become inconsistent when a box changes.
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::clearForBoxChanges(bool isInvalid)
{
   deallocateGlobalizedVersion();
   clearPersistentOverlapConnectors();
   if (isInvalid) {
      /*
       * Box removal can lead on inconsistent Connectors holding on to
       * handle, so detach the handle.  Box addition does NOT lead to
       * such inconsistencies, so we can leave the handle alone in
       * those cases.
       */
      detachMyHandle();
   }
}

SAMRAI_INLINE_KEYWORD
const IntVector& BoxLevel::getRefinementRatio() const
{
   return d_ratio;
}

SAMRAI_INLINE_KEYWORD
size_t BoxLevel::getLocalNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   return d_local_number_of_boxes;
}

SAMRAI_INLINE_KEYWORD
size_t BoxLevel::getLocalNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   return d_local_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getGlobalNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_boxes;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getMaxNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_max_number_of_boxes;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getMinNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_min_number_of_boxes;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getGlobalNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_global_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getMaxNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_max_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
int BoxLevel::getMinNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   cacheGlobalReducedData();
   return d_min_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
const Box& BoxLevel::getGlobalBoundingBox(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_bounding_box[block_num];
}

SAMRAI_INLINE_KEYWORD
const Box& BoxLevel::getLocalBoundingBox(int block_num) const
{
   return d_local_bounding_box[block_num];
}

SAMRAI_INLINE_KEYWORD
const IntVector& BoxLevel::getLocalMaxBoxSize(int block_num) const
{
   return d_local_max_box_size[block_num];
}

SAMRAI_INLINE_KEYWORD
const IntVector& BoxLevel::getGlobalMaxBoxSize(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_max_box_size[block_num];
}

SAMRAI_INLINE_KEYWORD
const IntVector& BoxLevel::getLocalMinBoxSize(int block_num) const
{
   return d_local_min_box_size[block_num];
}

SAMRAI_INLINE_KEYWORD
const IntVector& BoxLevel::getGlobalMinBoxSize(int block_num) const
{
   cacheGlobalReducedData();
   return d_global_min_box_size[block_num];
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxLevel::hasBox(
   const GlobalId& global_id,
   const PeriodicId& periodic_id) const
{
   const Box box(
      getDim(),
      global_id,
      periodic_id);
   return hasBox(box);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxLevel::hasBox(
   const BoxId& box_id) const
{
   const Box box(getDim(), box_id);
   return hasBox(box);
}

SAMRAI_INLINE_KEYWORD
const tbox::SAMRAI_MPI& BoxLevel::getMPI() const
{
   return d_mpi;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::addBoxWithoutUpdate(
   const Box& box)
{
   if (d_parallel_state == GLOBALIZED) {
      d_global_boxes.insert(box);
   }
   d_boxes.insert(box);
   return;
}

/*
 ****************************************************************************
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::eraseBoxWithoutUpdate(
   const Box& box)
{
   d_boxes.erase(box);
   return;
}

/*
 ****************************************************************************
 * TODO: This method puts finer in an inconsistent state by not
 * updating the attributes depenent on what has been changed.  This
 * method seems to be an initializer, but it is not clear from the
 * name or documentation if that is so.  The same goes for coarsenBoxes.
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::refineBoxes(
   BoxLevel& finer,
   const IntVector& ratio,
   const IntVector& final_ratio) const
{
   finer.d_boxes = d_boxes;
   finer.d_boxes.refine(ratio);
   finer.d_ratio = final_ratio;
   return;
}

/*
 ****************************************************************************
 ****************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::coarsenBoxes(
   BoxLevel& coarser,
   const IntVector& ratio,
   const IntVector& final_ratio) const
{
   coarser.d_boxes = d_boxes;
   coarser.d_boxes.coarsen(ratio);
   coarser.d_ratio = final_ratio;
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator BoxLevel::getBox(
   const Box& box) const
{
   if (box.getOwnerRank() == d_mpi.getRank()) {
      return d_boxes.find(box);
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "BoxLevel::getBox: cannot get remote box "
            << box << " without being in globalized state.");
      }
#endif
      return d_global_boxes.find(box);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator BoxLevel::getBox(
   const BoxId& box_id) const
{

   const Box box(getDim(), box_id);
   return getBox(box);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<const GridGeometry>&
BoxLevel::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const boost::shared_ptr<BoxLevelHandle>& BoxLevel::getBoxLevelHandle() const
{
   if (!d_handle) {
      /*
       * No handle yet.  Generate one attached to this object.
       */
      d_handle.reset(new BoxLevelHandle(this));
   }
   if (d_handle->d_box_level != this) {
      /*
       * Sanity check: The handle for this object should be attached
       * to this object.
       */
      TBOX_ERROR("Library error in BoxLevelHandle::getBoxLevel");
   }
   return d_handle;
}


/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::initializeCallback()
{
   t_initialize_private = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::initializePrivate()");
   t_acquire_remote_boxes = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::acquireRemoteBoxes()");
   t_cache_global_reduced_data = tbox::TimerManager::getManager()->
      getTimer("hier::BoxLevel::cacheGlobalReducedData()");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void BoxLevel::finalizeCallback()
{
   t_initialize_private.reset();
   t_acquire_remote_boxes.reset();
   t_cache_global_reduced_data.reset();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
