/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Set of relationships incident from a mapped_box_level of a distributed box graph.
 *
 ************************************************************************/
#include "SAMRAI/hier/Connector.h"

#include "SAMRAI/hier/BoxLevelHandle.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
bool Connector::isInitialized() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * The base and head mapped_box_level handles should have
    * been initialized when this object was initialized.
    * They should remain initialized as long as they are
    * referenced by a Connector.  Because this check may
    * occur when the base and head BoxLevels are in
    * transition, we don't require that the BoxLevels
    * themselves be initialized.
    */
   if (d_nproc != BAD_INT) {
      TBOX_ASSERT(!d_base_handle.isNull());
      TBOX_ASSERT(!d_head_handle.isNull());
   }
#endif
   return d_nproc != BAD_INT;
}

SAMRAI_INLINE_KEYWORD
BoxLevel::ParallelState Connector::getParallelState() const
{
   return d_parallel_state;
}

SAMRAI_INLINE_KEYWORD
const tbox::SAMRAI_MPI& Connector::getMPI() const
{
   TBOX_ASSERT(isInitialized());

   return d_base_handle->getBoxLevel().getMPI();
}

SAMRAI_INLINE_KEYWORD
int Connector::getRank() const
{
   return d_rank;
}

SAMRAI_INLINE_KEYWORD
int Connector::getNproc() const
{
   return d_nproc;
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& Connector::getConnectorWidth() const
{
   TBOX_ASSERT(isInitialized());
   return d_base_width;
}

SAMRAI_INLINE_KEYWORD
const hier::BoxLevel& Connector::getBase() const
{
   TBOX_ASSERT(isInitialized());
   return d_base_handle->getBoxLevel();
}

SAMRAI_INLINE_KEYWORD
const hier::BoxLevel& Connector::getHead() const
{
   TBOX_ASSERT(isInitialized());
   return d_head_handle->getBoxLevel();
}

SAMRAI_INLINE_KEYWORD
const hier::IntVector& Connector::getRatio() const
{
   TBOX_ASSERT(isInitialized());
   return d_ratio;
}

SAMRAI_INLINE_KEYWORD
bool Connector::ratioIsExact() const
{
   TBOX_ASSERT(isInitialized());
   return d_ratio_is_exact;
}

SAMRAI_INLINE_KEYWORD
bool Connector::getHeadCoarserFlag() const
{
   TBOX_ASSERT(isInitialized());
   return d_head_coarser;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::hasNeighborSet(
   const BoxId& mapped_box_id) const
{
   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_rank);
   }
   const NeighborhoodSet& relationships = d_parallel_state == BoxLevel::DISTRIBUTED ?
      d_relationships : d_global_relationships;
   BoxId non_per_id(mapped_box_id.getGlobalId(),
                    mapped_box_id.getBlockId(),
                    PeriodicId::zero());
   NeighborhoodSet::const_iterator ei = relationships.find(non_per_id);
   return ei != relationships.end();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const Connector::NeighborSet& Connector::getNeighborSet(
   const BoxId& mapped_box_id) const
{
   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_rank);
   }
   const NeighborhoodSet& relationships = d_parallel_state == BoxLevel::DISTRIBUTED ?
      d_relationships : d_global_relationships;
   BoxId non_per_id(mapped_box_id.getGlobalId(),
                    mapped_box_id.getBlockId(),
                    PeriodicId::zero());
   NeighborhoodSet::const_iterator ei = relationships.find(non_per_id);
   if (ei == relationships.end()) {
      TBOX_ERROR("Connector::getNeighborSet: No neighbor set exists for\n"
         << "mapped_box " << mapped_box_id << ".\n");
   }
   return ei->second;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const Connector& Connector::operator = (
   const Connector& rhs)
{
   if (this != &rhs) {
      d_base_handle = rhs.d_base_handle;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;
      d_global_number_of_neighbor_sets = rhs.d_global_number_of_neighbor_sets;
      d_head_handle = rhs.d_head_handle;
      d_global_number_of_relationships = rhs.d_global_number_of_relationships;
      d_relationships = rhs.d_relationships;
      d_global_relationships = rhs.d_global_relationships;
      d_base_width = rhs.d_base_width;
      d_ratio = rhs.d_ratio;
      d_head_coarser = rhs.d_head_coarser;
      d_parallel_state = rhs.d_parallel_state;
      d_nproc = rhs.d_nproc;
      d_rank = rhs.d_rank;
      d_connector_type = rhs.d_connector_type;
   }
   return *this;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::operator == (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return true;
   }
   // Note: two uninitialized Connectors always compare equal.
   if (!isInitialized() && !rhs.isInitialized()) {
      return true;
   }
   if (!isInitialized() && rhs.isInitialized()) {
      return false;
   }
   if (isInitialized() && !rhs.isInitialized()) {
      return false;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return false;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return false;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return false;
   }
   if (getNeighborhoodSets() != rhs.getNeighborhoodSets()) {
      return false;
   }

   return true;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::operator != (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return false;
   }
   // Note: two uninitialized Connectors always compare equal.
   if (!isInitialized() && !rhs.isInitialized()) {
      return false;
   }
   if (!isInitialized() && rhs.isInitialized()) {
      return true;
   }
   if (isInitialized() && !rhs.isInitialized()) {
      return true;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return true;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return true;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return true;
   }
   if (getNeighborhoodSets() != rhs.getNeighborhoodSets()) {
      return true;
   }

   return false;
}

}
}
