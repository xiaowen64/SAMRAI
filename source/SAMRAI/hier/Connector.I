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
bool Connector::isFinalized() const
{
   return d_finalized;
}

SAMRAI_INLINE_KEYWORD
BoxLevel::ParallelState Connector::getParallelState() const
{
   return d_parallel_state;
}

SAMRAI_INLINE_KEYWORD
const tbox::SAMRAI_MPI& Connector::getMPI() const
{
   TBOX_ASSERT(isFinalized());
   return d_base_handle->getBoxLevel().getMPI();
}

SAMRAI_INLINE_KEYWORD
const IntVector& Connector::getConnectorWidth() const
{
   TBOX_ASSERT(isFinalized());
   return d_base_width;
}

SAMRAI_INLINE_KEYWORD
const BoxLevel& Connector::getBase() const
{
   TBOX_ASSERT(isFinalized());
   return d_base_handle->getBoxLevel();
}

SAMRAI_INLINE_KEYWORD
const BoxLevel& Connector::getHead() const
{
   TBOX_ASSERT(isFinalized());
   return d_head_handle->getBoxLevel();
}

SAMRAI_INLINE_KEYWORD
const IntVector& Connector::getRatio() const
{
   TBOX_ASSERT(isFinalized());
   return d_ratio;
}

SAMRAI_INLINE_KEYWORD
bool Connector::ratioIsExact() const
{
   TBOX_ASSERT(isFinalized());
   return d_ratio_is_exact;
}

SAMRAI_INLINE_KEYWORD
bool Connector::getHeadCoarserFlag() const
{
   TBOX_ASSERT(isFinalized());
   return d_head_coarser;
}

SAMRAI_INLINE_KEYWORD
bool Connector::isLocal() const
{
   return d_relationships.isLocal(getMPI().getRank());
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::hasNeighborSet(
   const BoxId& mapped_box_id) const
{
   const BoxNeighborhoodCollection& relationships =
      getRelations(mapped_box_id);
   BoxId non_per_id(mapped_box_id.getGlobalId(),
                    PeriodicId::zero());
   ConstNeighborhoodIterator ei = relationships.find(non_per_id);
   return ei != relationships.end();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
const BoxNeighborhoodCollection& Connector::getRelations(
   const BoxId& mapped_box_id) const
{
   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      TBOX_ASSERT(mapped_box_id.getOwnerRank() == d_mpi.getRank());
   }
   const BoxNeighborhoodCollection& relationships =
      d_parallel_state == BoxLevel::DISTRIBUTED ?
      d_relationships : d_global_relationships;
   return relationships;
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
      d_finalized = rhs.d_finalized;
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
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return true;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return false;
   }
   if (isFinalized() && !rhs.isFinalized()) {
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
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return false;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return true;
   }
   if (isFinalized() && !rhs.isFinalized()) {
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
   const IntVector& ratio)
{
   d_relationships.coarsenNeighbors(ratio);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::refineLocalNeighbors(
   const IntVector& ratio)
{
   d_relationships.refineNeighbors(ratio);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::growLocalNeighbors(
   const IntVector& growth)
{
   d_relationships.growNeighbors(growth);
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
   TBOX_ASSERT( box_id.getOwnerRank() == d_mpi.getRank() );
   d_relationships.insert(box_id, neighbor);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::insertLocalNeighbor(
   const Box& neighbor,
   NeighborhoodIterator& base_box_itr)
{
   TBOX_ASSERT( base_box_itr->getOwnerRank() == d_mpi.getRank() );
   d_relationships.insert(base_box_itr, neighbor);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::NeighborhoodIterator Connector::makeEmptyLocalNeighborhood(
   const BoxId& box_id)
{
   TBOX_ASSERT( box_id.getOwnerRank() == d_mpi.getRank() );
   return d_relationships.insert(box_id).first;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::eraseLocalNeighborhood(
   const BoxId& box_id)
{
   TBOX_ASSERT( box_id.getOwnerRank() == d_mpi.getRank() );
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
   TBOX_ASSERT( box_id.getOwnerRank() == d_mpi.getRank() );
   return d_relationships.hasNeighbor(box_id, neighbor);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
bool Connector::isEmptyNeighborhood(
   const BoxId& box_id) const
{
   return getRelations(box_id).emptyBoxNeighborhood(box_id);
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
                    PeriodicId::zero());
   return d_relationships.numNeighbors(non_per_id);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalNeighbors(
   BoxContainer& neighbors) const
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
   std::map<BlockId, BoxContainer>& neighbors) const
{
   d_relationships.getNeighbors(neighbors);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalOwners(
   std::set<int>& owners) const
{
   d_relationships.getOwners(owners);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getLocalOwners(
   ConstNeighborhoodIterator& base_boxes_itr,
   std::set<int>& owners) const
{
   d_relationships.getOwners(base_boxes_itr, owners);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::clearNeighborhoods()
{
   d_relationships.clear();
   d_global_relationships.clear();
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
Connector::NeighborhoodIterator Connector::begin()
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
Connector::NeighborhoodIterator Connector::end()
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
   return nbrhd.d_collection->begin(nbrhd);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::NeighborIterator Connector::begin(
   NeighborhoodIterator& nbrhd)
{
   BoxNeighborhoodCollection* tmp =
      const_cast<BoxNeighborhoodCollection*>(nbrhd.d_collection);
   return tmp->begin(nbrhd);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::ConstNeighborIterator Connector::end(
   ConstNeighborhoodIterator& nbrhd) const
{
   return nbrhd.d_collection->end(nbrhd);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::NeighborIterator Connector::end(
   NeighborhoodIterator& nbrhd)
{
   BoxNeighborhoodCollection* tmp =
      const_cast<BoxNeighborhoodCollection*>(nbrhd.d_collection);
   return tmp->end(nbrhd);
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
                    PeriodicId::zero());
   return d_relationships.find(non_per_id);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
Connector::NeighborhoodIterator Connector::findLocal(
   const BoxId& box_id)
{
   BoxId non_per_id(box_id.getGlobalId(),
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
   const BoxNeighborhoodCollection& relationships = getRelations(box_id);
   BoxId non_per_id(box_id.getGlobalId(),
                    PeriodicId::zero());
   ConstNeighborhoodIterator ei = relationships.find(non_per_id);
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
bool Connector::localNeighborhoodsEqual(
   const Connector& other) const
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
   const BoxNeighborhoodCollection& relationships = getRelations(box_id);
   const BoxNeighborhoodCollection& other_relationships =
      other.getRelations(box_id);
   BoxId non_per_id(box_id.getGlobalId(),
                    PeriodicId::zero());
   return relationships.neighborhoodEqual(box_id, other_relationships);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::putNeighborhoodsToDatabase(
   tbox::Database& database) const
{
   d_relationships.putToDatabase(database);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void Connector::getNeighborhoodsFromDatabase(
   tbox::Database& database)
{
   d_relationships.getFromDatabase(database);
   return;
}

}
}
