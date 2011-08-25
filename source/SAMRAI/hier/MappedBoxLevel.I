/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of mapped_boxes in a mapped_box_level of a distributed box graph. 
 *
 ************************************************************************/
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/MappedBoxLevelHandle.h"
#include "SAMRAI/hier/BoxSetSingleOwnerIterator.h"

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
const tbox::Dimension& MappedBoxLevel::getDim() const
{
   return d_ratio.getDim();
}

SAMRAI_INLINE_KEYWORD
bool MappedBoxLevel::isInitialized() const
{
   return d_ratio(0) != 0;
}

SAMRAI_INLINE_KEYWORD
MappedBoxLevel::ParallelState
MappedBoxLevel::getParallelState() const
{
   return d_parallel_state;
}

/*
 ***********************************************************************
 * Detach this object from the handle it has been using.
 *
 * Postcondition: Object that cached the handle would no longer
 * be able to access this MappedBoxLevel from the handle.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void MappedBoxLevel::detachMyHandle() const
{
   if (!d_handle.isNull()) {
      d_handle->detachMyMappedBoxLevel();
      d_handle.setNull();
   }
}

SAMRAI_INLINE_KEYWORD
const BoxSet& MappedBoxLevel::getMappedBoxes() const
{
   return d_mapped_boxes;
}

SAMRAI_INLINE_KEYWORD
const BoxSet& MappedBoxLevel::getGlobalMappedBoxes() const
{
   return d_global_mapped_boxes;
}

SAMRAI_INLINE_KEYWORD
bool MappedBoxLevel::operator == (
   const MappedBoxLevel& r) const
{
   if (this == &r) return true;

   if (d_ratio != r.d_ratio) return false;

   if (d_mpi != r.d_mpi) return false;

   if (getMappedBoxes() != r.getMappedBoxes()) return false;

   return true;
}

SAMRAI_INLINE_KEYWORD
bool MappedBoxLevel::operator != (
   const MappedBoxLevel& r) const
{
   if (this == &r) return false;

   if (d_ratio != r.d_ratio) return true;

   if (d_mpi != r.d_mpi) return true;

   if (getMappedBoxes() != r.getMappedBoxes()) return true;

   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void MappedBoxLevel::clearPersistentOverlapConnectors()
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
void MappedBoxLevel::deallocateGlobalizedVersion() const
{
   if ( d_globalized_version != NULL ) {
      TBOX_ASSERT( d_globalized_version->getParallelState() == GLOBALIZED );
      delete d_globalized_version;
      d_globalized_version = NULL;
   }
   return;
}

/*
****************************************************************************
Clear certain data that would become inconsistent when a box changes.
****************************************************************************
*/
SAMRAI_INLINE_KEYWORD
void MappedBoxLevel::clearForBoxChanges( bool isInvalid )
{
   deallocateGlobalizedVersion();
   clearPersistentOverlapConnectors();
   if ( isInvalid ) {
      /*
       * Box removal can lead on inconsistent Connectors holding on to
       * handle, so detach the handle.  Box addition does NOT lead to
       * such inconsistencies, so we can leave the handle alone in
       * those cases.
       */
      detachMyHandle();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
MappedBoxLevel& MappedBoxLevel::operator = (
   const MappedBoxLevel& rhs)
{
   if (&rhs != this) {
      /*
       * Protect this block from assignment to self because it is
       * inefficient and it removes d_mapped_boxes data before resetting it.
       */

      deallocateGlobalizedVersion();
      clearPersistentOverlapConnectors();
      detachMyHandle();

      d_grid_geometry = rhs.d_grid_geometry;

      d_parallel_state = rhs.d_parallel_state;
      d_mpi = rhs.d_mpi;
      d_nproc = rhs.d_nproc;
      d_rank = rhs.d_rank;
      d_ratio = rhs.d_ratio;

      d_local_number_of_cells = rhs.d_local_number_of_cells;
      d_local_number_of_mapped_boxes = rhs.d_local_number_of_mapped_boxes;
      d_global_number_of_cells = rhs.d_global_number_of_cells;
      d_global_number_of_mapped_boxes = rhs.d_global_number_of_mapped_boxes;

      d_local_max_box_size = rhs.d_local_max_box_size;
      d_global_max_box_size = rhs.d_global_max_box_size;
      d_local_min_box_size = rhs.d_local_min_box_size;
      d_global_min_box_size = rhs.d_global_min_box_size;

      d_local_bounding_box = rhs.d_local_bounding_box;
      d_local_bounding_box_up_to_date = rhs.d_local_bounding_box_up_to_date;
      d_global_bounding_box = rhs.d_global_bounding_box;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;

      d_mapped_boxes = rhs.d_mapped_boxes;
      d_global_mapped_boxes = rhs.d_global_mapped_boxes;
      d_grid_geometry = rhs.d_grid_geometry;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
LocalId MappedBoxLevel::getFirstLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxSet& mapped_boxes = getMappedBoxes();
   if (mapped_boxes.empty()) {
      return s_negative_one_local_id;
   }
   BoxSet::const_iterator ni = mapped_boxes.begin();
   while (ni->isPeriodicImage()) {
      TBOX_ASSERT(ni != mapped_boxes.end());   // There should be a real mapped_box!
      ++ni;
   }
   return ni->getLocalId();
}

SAMRAI_INLINE_KEYWORD
LocalId MappedBoxLevel::getLastLocalId() const
{
   TBOX_ASSERT(isInitialized());

   const BoxSet& mapped_boxes = getMappedBoxes();
   if (mapped_boxes.empty()) {
      return s_negative_one_local_id;
   }
   LocalId last_local_id(0);
   for ( BoxSet::const_reverse_iterator ni = mapped_boxes.rbegin();
         ni!=mapped_boxes.rend(); ++ni ) {
      if ( last_local_id < ni->getLocalId() ) {
         last_local_id = ni->getLocalId();
      }
   }
   return last_local_id;
}

SAMRAI_INLINE_KEYWORD
const IntVector& MappedBoxLevel::getRefinementRatio() const
{
   return d_ratio;
}

SAMRAI_INLINE_KEYWORD
size_t MappedBoxLevel::getLocalNumberOfBoxes() const
{
   TBOX_ASSERT(isInitialized());

   return d_local_number_of_mapped_boxes;
}

SAMRAI_INLINE_KEYWORD
size_t MappedBoxLevel::getLocalNumberOfBoxes(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_rank) {
      TBOX_ERROR(
         "Non-local mapped_boxes are not available in DISTRIBUTED mode.");
   }
   TBOX_ASSERT(rank >= 0 && rank < d_nproc);
#endif

   if (rank == d_rank) {
      return d_local_number_of_mapped_boxes;
   } else {
      size_t count = 0;
      BoxSetSingleOwnerIterator mbi(d_global_mapped_boxes, rank);
      for ( ; mbi.isValid(); mbi++) {
         if (!(*mbi).isPeriodicImage()) {
            ++count;
         }
      }

      return count;
   }
}

SAMRAI_INLINE_KEYWORD
size_t MappedBoxLevel::getLocalNumberOfCells() const
{
   TBOX_ASSERT(isInitialized());

   return d_local_number_of_cells;
}

SAMRAI_INLINE_KEYWORD
size_t MappedBoxLevel::getLocalNumberOfCells(
   int rank) const
{
   TBOX_ASSERT(isInitialized());
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == DISTRIBUTED && rank != d_rank) {
      TBOX_ERROR(
         "Non-local mapped_boxes are not available in DISTRIBUTED mode.");
   }
   TBOX_ASSERT(rank >= 0 && rank < d_nproc);
#endif

   if (rank == d_rank) {
      return d_local_number_of_cells;
   } else {
      size_t count = 0;
      BoxSetSingleOwnerIterator mbi(d_global_mapped_boxes, rank);
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
bool MappedBoxLevel::hasMappedBox(
   const Box& mapped_box) const
{
   if (mapped_box.getOwnerRank() == d_rank) {

      BoxSet::const_iterator ni = d_mapped_boxes.find(mapped_box);
      return ni != d_mapped_boxes.end();
 
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state == DISTRIBUTED) {
         TBOX_ERROR("MappedBoxLevel: Cannot check on remote Box "
            << mapped_box << " while in DISTRIBUTED mode.\n"
            << "See MappedBoxLevel::setParallelState().");
      }
#endif
      BoxSet::const_iterator ni = d_global_mapped_boxes.find(mapped_box);
      return ni != d_global_mapped_boxes.end();
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool MappedBoxLevel::hasMappedBox(
   const GlobalId& global_id,
   const BlockId& block_id,
   const PeriodicId &periodic_id) const
{
   const Box mapped_box(
      getDim(),
      global_id,
      block_id,
      periodic_id);
   return hasMappedBox(mapped_box);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool MappedBoxLevel::hasMappedBox(
   const BoxId& mapped_box_id) const
{
   const Box mapped_box(getDim(), mapped_box_id);
   return hasMappedBox(mapped_box);
}

SAMRAI_INLINE_KEYWORD
const tbox::SAMRAI_MPI& MappedBoxLevel::getMPI() const
{
   return d_mpi;
}

SAMRAI_INLINE_KEYWORD
int MappedBoxLevel::getRank() const
{
   return d_rank;
}

SAMRAI_INLINE_KEYWORD
int MappedBoxLevel::getNproc() const
{
   return d_nproc;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxSet::const_iterator
MappedBoxLevel::getMappedBox(
   const Box& mapped_box) const
{
   BoxSet::const_iterator ni;
   if (mapped_box.getOwnerRank() == d_rank) {
      ni = d_mapped_boxes.find(mapped_box);
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "MappedBoxLevel::getMappedBox: cannot get remote mapped_box "
            << mapped_box << " without being in globalized state.");
      }
#endif
      ni = d_global_mapped_boxes.find(mapped_box);
   }
   return ni;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxSet::const_iterator
MappedBoxLevel::getMappedBox(
   const BoxId& mapped_box_id) const
{

   const Box mapped_box(getDim(),
                              mapped_box_id);
   BoxSet::const_iterator ni;
   if (mapped_box.getOwnerRank() == d_rank) {
      ni = d_mapped_boxes.find(mapped_box);
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "MappedBoxLevel::getMappedBox: cannot get remote mapped_box "
            << mapped_box << " without being in globalized state.");
      }
#endif
      ni = d_global_mapped_boxes.find(mapped_box);
   }
   return ni;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxSet::const_iterator
MappedBoxLevel::getMappedBoxStrict(
   const Box& mapped_box) const
{
   BoxSet::const_iterator ni;
   if (mapped_box.getOwnerRank() == d_rank) {
      ni = d_mapped_boxes.find(mapped_box);
      if (ni == d_mapped_boxes.end()) {
         TBOX_ERROR(
            "BoxSet::getMappedBoxStrict: requested mapped_box "
            << mapped_box
            <<
            " does not exist in the mapped_box_level.");
      }
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parallel_state != GLOBALIZED) {
         TBOX_ERROR(
            "MappedBoxLevel::getMappedBox: cannot get remote mapped_box "
            << mapped_box << " without being in globalized state.");
      }
#endif
      ni = d_global_mapped_boxes.find(mapped_box);
      if (ni == d_global_mapped_boxes.end()) {
         TBOX_ERROR(
            "BoxSet::getMappedBoxStrict: requested mapped_box "
            << mapped_box
            <<
            " does not exist in the mapped_box_level.");
      }
   }
   return ni;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
BoxSet::const_iterator
MappedBoxLevel::getMappedBoxStrict(
   const BoxId& mapped_box_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (mapped_box_id.getOwnerRank() != d_rank && d_parallel_state != GLOBALIZED) {
      TBOX_ERROR(
         "MappedBoxLevel::getMappedBoxStrict: cannot get remote mapped_box " << mapped_box_id
         <<
         " without being in globalized state.");
   }
#endif

   Box mapped_box(getDim(),
                        mapped_box_id);
   BoxSet::const_iterator ni;
   if (mapped_box.getOwnerRank() == d_rank) {
      ni = d_mapped_boxes.find(mapped_box);
      if (ni == d_mapped_boxes.end()) {
         TBOX_ERROR(
            "BoxSet::getMappedBoxStrict: requested mapped_box "
            << mapped_box
            <<
            " does not exist in the mapped_box_level.");
      }
   } else {
      ni = d_global_mapped_boxes.find(mapped_box);
      if (ni == d_global_mapped_boxes.end()) {
         TBOX_ERROR(
            "BoxSet::getMappedBoxStrict: requested mapped_box "
            << mapped_box
            <<
            " does not exist in the mapped_box_level.");
      }
   }
   return ni;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::ConstPointer<GridGeometry>&
MappedBoxLevel::getGridGeometry() const
{
   return d_grid_geometry;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Pointer<MappedBoxLevelHandle>& MappedBoxLevel::
getMappedBoxLevelHandle() const
{
   if (d_handle.isNull()) {
      /*
       * No handle yet.  Generate one attached to this object.
       */
      tbox::Pointer<MappedBoxLevelHandle> ptr(new MappedBoxLevelHandle(this));
      d_handle = ptr;
   }
   if (d_handle->d_mapped_box_level != this) {
      /*
       * Sanity check: The handle for this object should be attached
       * to this object.
       */
      TBOX_ERROR("Library error in MappedBoxLevelHandle::getMappedBoxLevel");
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
