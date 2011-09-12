/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   A class describing the adjacency of Boxes.
 *
 ************************************************************************/

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::BoxNeighborhoodCollection(
   int rank) :
   d_rank(rank)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::~BoxNeighborhoodCollection()
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::begin()
{
   return Iterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator
BoxNeighborhoodCollection::begin() const
{
   return ConstIterator(*this);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::end()
{
   return Iterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator BoxNeighborhoodCollection::end() const
{
   return ConstIterator(*this, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::begin(
   const Box& root)
{
   Iterator itr(find(root));
   return begin(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::begin(
   const Box& root) const
{
   return begin(find(root));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::begin(
   const BoxId& root_id)
{
   Iterator itr(find(root_id));
   return begin(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::begin(
   const BoxId& root_id) const
{
   return begin(find(root_id));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::begin(
   BoxNeighborhoodCollection::Iterator& root_itr)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoolCollection::begin root_itr points to end");
   }

   return NeighborIterator(*this, root_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::begin(
   const BoxNeighborhoodCollection::ConstIterator& root_itr) const
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoolCollection::begin root_itr points to end");
   }

   return ConstNeighborIterator(*this, root_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::end(
   const Box& root)
{
   Iterator itr(find(root));
   return end(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::end(
   const Box& root) const
{
   return end(find(root));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::end(
   const BoxId& root_id)
{
   Iterator itr(find(root_id));
   return end(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::end(
   const BoxId& root_id) const
{
   return end(find(root_id));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::end(
   BoxNeighborhoodCollection::Iterator& root_itr)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::end root_itr points to end");
   }

   return NeighborIterator(*this, root_itr, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::end(
   const BoxNeighborhoodCollection::ConstIterator& root_itr) const
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::end root_itr points to end");
   }

   return ConstNeighborIterator(*this, root_itr, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator BoxNeighborhoodCollection::find(
   const BoxId& root_id) const
{
   return find(Box(tbox::Dimension::getMaxDimension(), root_id));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::find(
   const BoxId& root_id)
{
   return find(Box(tbox::Dimension::getMaxDimension(), root_id));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::empty() const
{
   return begin() == end();
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::emptyBoxNeighborhood(
   const Box& root) const
{
   return emptyBoxNeighborhood(find(root));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::emptyBoxNeighborhood(
   const BoxId& root_id) const
{
   return emptyBoxNeighborhood(find(root_id));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::emptyBoxNeighborhood(
   const ConstIterator& root_itr) const
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::emptyBoxNeighborhood root_itr points to end");
   }

   return root_itr.d_itr->second.size() == 0;
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numNeighbors(
   const Box& root) const
{
   return numNeighbors(find(root));
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numNeighbors(
   const BoxId& root_id) const
{
   return numNeighbors(find(root_id));
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numNeighbors(
   const ConstIterator& root_itr) const
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::numNeighbors root_itr points to end");
   }

   return root_itr.d_itr->second.size();
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::hasPeriodicNeighborhood(
   const Box& root) const
{
   return hasPeriodicNeighborhood(find(root));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::hasPeriodicNeighborhood(
   const BoxId& root_id) const
{
   return hasPeriodicNeighborhood(find(root_id));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::isRoot(
   const Box& root) const
{
   return find(root) != end();
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::isRoot(
   const BoxId& root_id) const
{
   return find(root_id) != end();
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::insert(
   const BoxId& root_id,
   const Box& new_nbr)
{
   Iterator itr(find(root_id));
   insert(itr, new_nbr);
   return itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::insert(
   const BoxId& root_id,
   const BoxContainer& new_nbrs)
{
   Iterator itr(find(root_id));
   insert(itr, new_nbrs);
   return itr;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const Box& root,
   const Box& nbr)
{
   Iterator itr(find(root));
   erase(itr, nbr);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& root_id,
   const Box& nbr)
{
   Iterator itr(find(root_id));
   erase(itr, nbr);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   Iterator& root_itr,
   const Box& nbr)
{
   if (root_itr == end()) {
      TBOX_ERROR("hier::BoxNeighborhoodCollection::erase root_itr points to end");
   }

   root_itr.d_itr->second.erase(&nbr);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const Box& root,
   const BoxContainer& nbrs)
{
   Iterator itr(find(root));
   erase(itr, nbrs);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& root_id,
   const BoxContainer& nbrs)
{
   Iterator itr(find(root_id));
   erase(itr, nbrs);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const Box& root,
   bool erase_root)
{
   Iterator itr(find(root));
   erase(itr, erase_root);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& root_id,
   bool erase_root)
{
   Iterator itr(find(root_id));
   erase(itr, erase_root);
   return;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   AdjListItr itr) :
   d_collection(&nbrhds),
   d_itr(itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::Iterator(
   const Iterator& other) :
   d_collection(other.d_collection),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator&
BoxNeighborhoodCollection::Iterator::operator = (
   const Iterator& rhs)
{
   d_collection = rhs.d_collection;
   d_itr = rhs.d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::~Iterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxNeighborhoodCollection::Iterator::operator * ()
{
   return d_itr->first;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::Iterator::operator == (
   const Iterator& rhs) const
{
   return d_collection == rhs.d_collection && d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::Iterator::operator != (
   const Iterator& rhs) const
{
   return d_collection != rhs.d_collection || d_itr != rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   AdjListConstItr itr) :
   d_collection(&nbrhds),
   d_itr(itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const ConstIterator& other) :
   d_collection(other.d_collection),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator&
BoxNeighborhoodCollection::ConstIterator::operator = (
   const ConstIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_itr = rhs.d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::~ConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxNeighborhoodCollection::ConstIterator::operator * () const
{
   return d_itr->first;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstIterator::operator == (
   const ConstIterator& rhs) const
{
   return d_collection == rhs.d_collection && d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstIterator::operator != (
   const ConstIterator& rhs) const
{
   return d_collection != rhs.d_collection || d_itr != rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator::NeighborIterator(
   BoxNeighborhoodCollection& nbrhds,
   Iterator& root,
   bool from_start) :
   d_collection(&nbrhds),
   d_root(&(root.d_itr->first)),
   d_itr(from_start ? root.d_itr->second.begin() :
                      root.d_itr->second.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator::NeighborIterator(
   const NeighborIterator& other) :
   d_collection(other.d_collection),
   d_root(other.d_root),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator&
BoxNeighborhoodCollection::NeighborIterator::operator = (
   const NeighborIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_root = rhs.d_root;
   d_itr = rhs.d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator::~NeighborIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box& BoxNeighborhoodCollection::NeighborIterator::operator * () const
{
   return *(*d_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator&
BoxNeighborhoodCollection::NeighborIterator::operator ++ ()
{
   ++d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::NeighborIterator::operator == (
   const NeighborIterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_root == rhs.d_root &&
          d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::NeighborIterator::operator != (
   const NeighborIterator& rhs) const
{
   return d_collection != rhs.d_collection ||
          d_root != rhs.d_root ||
          d_itr != rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const BoxNeighborhoodCollection& nbrhds,
   const ConstIterator& root,
   bool from_start) :
   d_collection(&nbrhds),
   d_root(&(root.d_itr->first)),
   d_itr(from_start ? root.d_itr->second.begin() :
                      root.d_itr->second.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const ConstNeighborIterator& other) :
   d_collection(other.d_collection),
   d_root(other.d_root),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator&
BoxNeighborhoodCollection::ConstNeighborIterator::operator = (
   const ConstNeighborIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_root = rhs.d_root;
   d_itr = rhs.d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::~ConstNeighborIterator()
{
}

SAMRAI_INLINE_KEYWORD
const Box&
BoxNeighborhoodCollection::ConstNeighborIterator::operator * () const
{
   return *(*d_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator&
BoxNeighborhoodCollection::ConstNeighborIterator::operator ++ ()
{
   ++d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstNeighborIterator::operator == (
   const ConstNeighborIterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_root == rhs.d_root &&
          d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstNeighborIterator::operator != (
   const ConstNeighborIterator& rhs) const
{
   return d_collection != rhs.d_collection ||
          d_root != rhs.d_root ||
          d_itr != rhs.d_itr;
}

}
}
