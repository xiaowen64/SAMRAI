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
   return !d_base_handle.isNull();
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
      TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mpi.getRank());
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
   return find(mapped_box_id)->second;
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
      d_mpi = rhs.d_mpi;
      d_base_width = rhs.d_base_width;
      d_ratio = rhs.d_ratio;
      d_head_coarser = rhs.d_head_coarser;
      d_parallel_state = rhs.d_parallel_state;
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
   if (d_relationships != rhs.d_relationships) {
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
   if (d_relationships != rhs.d_relationships) {
      return true;
   }

   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::coarsenLocalNeighbors(
   Connector& coarser,
   const IntVector& ratio) const
{
   if (coarser != *this) {
      coarser.d_relationships.clear();
   }
   d_relationships.coarsenNeighbors(coarser.d_relationships, ratio);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::refineLocalNeighbors(
   Connector& finer,
   const IntVector& ratio) const
{
   if (finer != *this) {
      finer.d_relationships.clear();
   }
   d_relationships.refineNeighbors(finer.d_relationships, ratio);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::growLocalNeighbors(
   Connector& grown,
   const IntVector& growth) const
{
   if (grown != *this) {
      grown.d_relationships.clear();
   }
   d_relationships.growNeighbors(grown.d_relationships, growth);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::insertLocalNeighbor(
   const Box& neighbor,
   const BoxId& box_id)
{
   d_relationships[box_id].insert(neighbor);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::makeEmptyLocalNeighborhood(
      const BoxId& box_id)
{
   d_relationships[box_id];
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::eraseLocalNeighborhood(
      const BoxId& box_id)
{
   d_relationships.erase(box_id);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::hasLocalNeighbor(
   const BoxId& box_id,
   const Box& neighbor) const
{
   NeighborhoodSet::const_iterator itr = d_relationships.find(box_id);
   if (itr == d_relationships.end()) {
      return false;
   }
   else if ((*itr).second.find(neighbor) == (*itr).second.orderedEnd()) {
      return false;
   }
   else {
      return true;
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::isEmptyNeighborhood(
   const BoxId& box_id) const
{
   TBOX_ASSERT(hasNeighborSet(box_id));

   const NeighborhoodSet& relationships =
      d_parallel_state == BoxLevel::DISTRIBUTED ?
      d_relationships : d_global_relationships;

   return relationships.find(box_id)->second.isEmpty();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
int Connector::numLocalNeighbors(
   const BoxId& box_id) const
{
   TBOX_ASSERT(hasNeighborSet(box_id));

   BoxId non_per_id(box_id.getGlobalId(),
                    box_id.getBlockId(),
                    PeriodicId::zero());
   return d_relationships.find(non_per_id)->second.size();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalNeighbors(
   NeighborSet& neighbors) const
{
   d_relationships.getNeighbors(neighbors);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalNeighbors(
   std::map<BlockId, BoxList>& neighbors) const
{
   d_relationships.getNeighbors(neighbors);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalOwners(std::set<int>& owners) const
{
   d_relationships.getOwners(owners);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::clearLocalNeighborhoods()
{
   d_relationships.clear();
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborhoodIterator Connector::begin() const
{
   return d_relationships.begin();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborhoodIterator Connector::end() const
{
   return d_relationships.end();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborIterator Connector::begin(
   ConstNeighborhoodIterator& nbrhd) const
{
   TBOX_ASSERT(nbrhd != end());
   return nbrhd->second.orderedBegin();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborIterator Connector::end(
   ConstNeighborhoodIterator& nbrhd) const
{
   TBOX_ASSERT(nbrhd != end());
   return nbrhd->second.orderedEnd();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborhoodIterator Connector::findLocal(
   const BoxId& box_id) const
{
   BoxId non_per_id(box_id.getGlobalId(),
                    box_id.getBlockId(),
                    PeriodicId::zero());
   return d_relationships.find(non_per_id);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborhoodIterator Connector::find(
   const BoxId& box_id) const
{
   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      TBOX_ASSERT(box_id.getOwnerRank() == d_mpi.getRank());
   }
   const NeighborhoodSet& relationships = d_parallel_state == BoxLevel::DISTRIBUTED ?
      d_relationships : d_global_relationships;
   BoxId non_per_id(box_id.getGlobalId(),
                    box_id.getBlockId(),
                    PeriodicId::zero());
   NeighborhoodSet::const_iterator ei = relationships.find(non_per_id);
   if (ei == relationships.end()) {
      TBOX_ERROR("Connector::find: No neighbor set exists for\n"
         << "box " << box_id << ".\n");
   }
   return ei;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::localNeighborhoodsEqual(const Connector& other) const
{
   return d_relationships == other.d_relationships;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::neighborhoodEqual(
   const BoxId& box_id,
   const Connector& other) const
{
   return getNeighborSet(box_id) == other.getNeighborSet(box_id);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::putNeighborhoodsToDatabase(tbox::Database& database) const
{
   d_relationships.putToDatabase(database);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getNeighborhoodsFromDatabase(tbox::Database& database)
{
   d_relationships.getFromDatabase(database);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::writeNeighborhoodsToErrorStream(
   const std::string& border) const
{
   tbox::perr << d_relationships.format(border, 2);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::writeNeighborhoodToErrorStream(
   const BoxId& box_id,
   const std::string& border) const
{
   tbox::perr << getNeighborSet(box_id).format(border);
   return;
}

}
}
