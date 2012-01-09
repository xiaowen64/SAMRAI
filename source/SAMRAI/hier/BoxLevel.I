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
#include "SAMRAI/hier/BoxContainerIterator.h"
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"

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

SAMRAI_INLINE_KEYWORD
bool BoxLevel::operator == (
   const BoxLevel& r) const
{
   if (this == &r) return true;

   if (d_ratio != r.d_ratio) return false;

   if (d_mpi != r.d_mpi) return false;

   if (getBoxes() != r.getBoxes()) return false;

   return true;
}

SAMRAI_INLINE_KEYWORD
bool BoxLevel::operator != (
   const BoxLevel& r) const
{
   if (this == &r) return false;

   if (d_ratio != r.d_ratio) return true;

   if (d_mpi != r.d_mpi) return true;

   if (getBoxes() != r.getBoxes()) return true;

   return false;
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

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxLevel& BoxLevel::operator = (
   const BoxLevel& rhs)
{
   if (&rhs != this) {
      /*
       * Protect this block from assignment to self because it is
       * inefficient and it removes d_boxes data before resetting it.
       */

      deallocateGlobalizedVersion();
      clearPersistentOverlapConnectors();
      detachMyHandle();

      d_grid_geometry = rhs.d_grid_geometry;

      d_parallel_state = rhs.d_parallel_state;
      d_mpi = rhs.d_mpi;
      d_ratio = rhs.d_ratio;

      d_local_number_of_cells = rhs.d_local_number_of_cells;
      d_local_number_of_boxes = rhs.d_local_number_of_boxes;
      d_global_number_of_cells = rhs.d_global_number_of_cells;
      d_global_number_of_boxes = rhs.d_global_number_of_boxes;

      d_local_max_box_size = rhs.d_local_max_box_size;
      d_global_max_box_size = rhs.d_global_max_box_size;
      d_local_min_box_size = rhs.d_local_min_box_size;
      d_global_min_box_size = rhs.d_global_min_box_size;

      d_local_bounding_box = rhs.d_local_bounding_box;
      d_local_bounding_box_up_to_date = rhs.d_local_bounding_box_up_to_date;
      d_global_bounding_box = rhs.d_global_bounding_box;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;

      d_boxes = rhs.d_boxes;
      d_global_boxes = rhs.d_global_boxes;
      d_grid_geometry = rhs.d_grid_geometry;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
LocalId BoxLevel::getFirstLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxContainer& boxes = getBoxes();
   if (boxes.isEmpty()) {
      return s_negative_one_local_id;
   }
   BoxContainer::ConstIterator ni = boxes.begin();
   while (ni->isPeriodicImage()) {
      TBOX_ASSERT(ni != boxes.end());   // There should be a real box!
      ++ni;
   }
   return ni->getLocalId();
}

SAMRAI_INLINE_KEYWORD
LocalId BoxLevel::getLastLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxContainer& boxes = getBoxes();
   if (boxes.isEmpty()) {
      return s_negative_one_local_id;
   }
   LocalId last_local_id(0);
   for (BoxContainer::ConstIterator ni = boxes.begin();
        ni != boxes.end(); ++ni) {
      if (last_local_id < ni->getLocalId()) {
         last_local_id = ni->getLocalId();
      }
   }
   return last_local_id;
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
size_t BoxLevel::getLocalNumberOfBoxes(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_mpi.getRank()) {
      TBOX_ERROR(
         "Non-local boxes are not available in DISTRIBUTED mode.");
   }
   TBOX_ASSERT(rank >= 0 && rank < d_mpi.getSize());
#endif

   if (rank == d_mpi.getRank()) {
      return d_local_number_of_boxes;
   } else {
      size_t count = 0;
      BoxContainerSingleOwnerIterator mbi(d_global_boxes, rank);
      for ( ; mbi.isValid(); mbi++) {
         if (!(*mbi).isPeriodicImage()) {
            ++count;
         }
      }

      return count;
   }
}

SAMRAI_INLINE_KEYWORD
size_t BoxLevel::getLocalNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   return d_local_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
size_t BoxLevel::getLocalNumberOfCells(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_mpi.getRank()) {
      TBOX_ERROR(
         "Non-local boxes are not available in DISTRIBUTED mode.");
   }
   TBOX_ASSERT(rank >= 0 && rank < d_mpi.getSize());
#endif

   if (rank == d_mpi.getRank()) {
      return d_local_number_of_cells;
   } else {
      size_t count = 0;
      BoxContainerSingleOwnerIterator mbi(d_global_boxes, rank);
      for ( ; mbi.isValid(); mbi++) {
         if (!(*mbi).isPeriodicImage()) {
            count += (*mbi).size();
         }
      }

      return count;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool BoxLevel::hasBox(
   const Box& box) const
{
   if (box.getOwnerRank() == d_mpi.getRank()) {

      BoxContainer::ConstIterator ni = d_boxes.find(box);
      return ni != d_boxes.end();

   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state == DISTRIBUTED) {
         TBOX_ERROR("BoxLevel: Cannot check on remote Box "
            << box << " while in DISTRIBUTED mode.\n"
            << "See BoxLevel::setParallelState().");
      }
#endif
      BoxContainer::ConstIterator ni = d_global_boxes.find(box);
      return ni != d_global_boxes.end();
   }
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
BoxContainer::ConstIterator
BoxLevel::getBox(
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
BoxContainer::ConstIterator
BoxLevel::getBox(
   const BoxId& box_id) const
{

   const Box box(getDim(),
                        box_id);
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
BoxContainer::ConstIterator
BoxLevel::getBoxStrict(
   const Box& box) const
{
   if (box.getOwnerRank() == d_mpi.getRank()) {
      BoxContainer::ConstIterator ni = d_boxes.find(box);
      if (ni == d_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level.");
      }

      return ni;
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "BoxLevel::getBox: cannot get remote box "
            << box << " without being in globalized state.");
      }
#endif
      BoxContainer::ConstIterator ni = d_global_boxes.find(box);
      if (ni == d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level.");
      }

      return ni;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxContainer::ConstIterator
BoxLevel::getBoxStrict(
   const BoxId& box_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (box_id.getOwnerRank() != d_mpi.getRank() && d_parallel_state != GLOBALIZED) {
      TBOX_ERROR(
         "BoxLevel::getBoxStrict: cannot get remote box " << box_id
                                                                 <<
         " without being in globalized state.");
   }
#endif

   Box box(getDim(),
                  box_id);
   if (box.getOwnerRank() == d_mpi.getRank()) {
      BoxContainer::ConstIterator ni = d_boxes.find(box);
      if (ni == d_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level.");
      }
      return ni;
   } else {
      BoxContainer::ConstIterator ni = d_global_boxes.find(box);
      if (ni == d_global_boxes.end()) {
         TBOX_ERROR(
            "BoxContainer::getBoxStrict: requested box "
            << box
            << " does not exist in the box_level.");
      }
      return ni;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Pointer<const GridGeometry>&
BoxLevel::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Pointer<BoxLevelHandle>& BoxLevel::
getBoxLevelHandle() const
{
   if (!d_handle) {
      /*
       * No handle yet.  Generate one attached to this object.
       */
      tbox::Pointer<BoxLevelHandle> ptr(new BoxLevelHandle(this));
      d_handle = ptr;
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

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
