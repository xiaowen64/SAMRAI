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
BoxNeighborhoodCollection::BoxNeighborhoodCollection()
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
   const BoxId& base_box_id)
{
   Iterator itr(find(base_box_id));
   return begin(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::begin(
   const BoxId& base_box_id) const
{
   return begin(find(base_box_id));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::begin(
   BoxNeighborhoodCollection::Iterator& base_box_itr)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return NeighborIterator(base_box_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::begin(
   const BoxNeighborhoodCollection::ConstIterator& base_box_itr) const
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return ConstNeighborIterator(base_box_itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::end(
   const BoxId& base_box_id)
{
   Iterator itr(find(base_box_id));
   return end(itr);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::end(
   const BoxId& base_box_id) const
{
   return end(find(base_box_id));
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator BoxNeighborhoodCollection::end(
   BoxNeighborhoodCollection::Iterator& base_box_itr)
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return NeighborIterator(base_box_itr, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator
BoxNeighborhoodCollection::end(
   const BoxNeighborhoodCollection::ConstIterator& base_box_itr) const
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return ConstNeighborIterator(base_box_itr, false);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator BoxNeighborhoodCollection::find(
   const BoxId& base_box_id) const
{
   BaseBoxPoolItr base_boxes_itr = d_base_boxes.find(base_box_id);
   if (base_boxes_itr == d_base_boxes.end()) {
      return end();
   }
   else {
      return ConstIterator(*this, d_adj_list.find(&(*base_boxes_itr)));
   }
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::find(
   const BoxId& base_box_id)
{
   BaseBoxPoolItr base_boxes_itr = d_base_boxes.find(base_box_id);
   if (base_boxes_itr == d_base_boxes.end()) {
      return end();
   }
   else {
      return Iterator(*this, d_adj_list.find(&(*base_boxes_itr)));
   }
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::empty() const
{
   return d_base_boxes.empty();
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numBoxNeighborhoods() const
{
   return static_cast<int>(d_base_boxes.size());
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::emptyBoxNeighborhood(
   const BoxId& base_box_id) const
{
   return emptyBoxNeighborhood(find(base_box_id));
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::emptyBoxNeighborhood(
   const ConstIterator& base_box_itr) const
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return base_box_itr.d_itr->second.empty();
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numNeighbors(
   const BoxId& base_box_id) const
{
   return numNeighbors(find(base_box_id));
}

SAMRAI_INLINE_KEYWORD
int BoxNeighborhoodCollection::numNeighbors(
   const ConstIterator& base_box_itr) const
{
   TBOX_ASSERT(base_box_itr.d_collection == this);
   TBOX_ASSERT(base_box_itr != end());

   return static_cast<int>(base_box_itr.d_itr->second.size());
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::hasNeighbor(
   const BoxId& base_box_id,
   const Box& nbr) const
{
   return hasNeighbor(find(base_box_id), nbr);
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::isBaseBox(
   const BoxId& base_box_id) const
{
   return find(base_box_id) != end();
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::insert(
   const BoxId& base_box_id,
   const Box& new_nbr)
{
   Iterator base_boxes_itr = insert(base_box_id).first;
   insert(base_boxes_itr, new_nbr);
   return base_boxes_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator BoxNeighborhoodCollection::insert(
   const BoxId& base_box_id,
   const BoxContainer& new_nbrs)
{
   Iterator base_boxes_itr = insert(base_box_id).first;
   insert(base_boxes_itr, new_nbrs);
   return base_boxes_itr;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& base_box_id,
   const Box& nbr)
{
   Iterator itr(find(base_box_id));
   erase(itr, nbr);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& base_box_id,
   const BoxContainer& nbrs)
{
   Iterator itr(find(base_box_id));
   erase(itr, nbrs);
   return;
}

SAMRAI_INLINE_KEYWORD
void BoxNeighborhoodCollection::erase(
   const BoxId& base_box_id)
{
   Iterator itr(find(base_box_id));
   erase(itr);
   return;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   bool from_start) :
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
                      nbrhds.d_adj_list.end()),
   d_base_boxes_itr(from_start ? nbrhds.d_base_boxes.begin() :
                                 nbrhds.d_base_boxes.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::Iterator(
   BoxNeighborhoodCollection& nbrhds,
   AdjListItr itr) :
   d_collection(&nbrhds),
   d_itr(itr),
   d_base_boxes_itr(nbrhds.d_base_boxes.find(*(itr->first)))
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::Iterator(
   const Iterator& other) :
   d_collection(other.d_collection),
   d_itr(other.d_itr),
   d_base_boxes_itr(other.d_base_boxes_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator&
BoxNeighborhoodCollection::Iterator::operator = (
   const Iterator& rhs)
{
   d_collection = rhs.d_collection;
   d_itr = rhs.d_itr;
   d_base_boxes_itr = rhs.d_base_boxes_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator::~Iterator()
{
}

SAMRAI_INLINE_KEYWORD
const BoxId& BoxNeighborhoodCollection::Iterator::operator * () const
{
   return *(d_itr->first);
}

SAMRAI_INLINE_KEYWORD
const BoxId* BoxNeighborhoodCollection::Iterator::operator -> () const
{
   return d_itr->first;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::Iterator&
BoxNeighborhoodCollection::Iterator::operator ++ ()
{
   // Go to the next base Box.
   if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
      ++d_base_boxes_itr;
      ++d_itr;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::Iterator::operator == (
   const Iterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_itr == rhs.d_itr &&
          d_base_boxes_itr == rhs.d_base_boxes_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::Iterator::operator != (
   const Iterator& rhs) const
{
   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   bool from_start) :
   d_collection(&nbrhds),
   d_itr(from_start ? nbrhds.d_adj_list.begin() :
                      nbrhds.d_adj_list.end()),
   d_base_boxes_itr(from_start ? nbrhds.d_base_boxes.begin() :
                                 nbrhds.d_base_boxes.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const BoxNeighborhoodCollection& nbrhds,
   AdjListConstItr itr) :
   d_collection(&nbrhds),
   d_itr(itr),
   d_base_boxes_itr(nbrhds.d_base_boxes.find(*(itr->first)))
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const ConstIterator& other) :
   d_collection(other.d_collection),
   d_itr(other.d_itr),
   d_base_boxes_itr(other.d_base_boxes_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::ConstIterator(
   const Iterator& other) :
   d_collection(other.d_collection),
   d_itr(other.d_itr),
   d_base_boxes_itr(other.d_base_boxes_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator&
BoxNeighborhoodCollection::ConstIterator::operator = (
   const ConstIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_itr = rhs.d_itr;
   d_base_boxes_itr = rhs.d_base_boxes_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator&
BoxNeighborhoodCollection::ConstIterator::operator = (
   const Iterator& rhs)
{
   d_collection = rhs.d_collection;
   d_itr = rhs.d_itr;
   d_base_boxes_itr = rhs.d_base_boxes_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator::~ConstIterator()
{
}

SAMRAI_INLINE_KEYWORD
const BoxId& BoxNeighborhoodCollection::ConstIterator::operator * () const
{
   return *(d_itr->first);
}

SAMRAI_INLINE_KEYWORD
const BoxId* BoxNeighborhoodCollection::ConstIterator::operator -> () const
{
   return d_itr->first;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstIterator&
BoxNeighborhoodCollection::ConstIterator::operator ++ ()
{
   // Go to the next base Box.
   if (d_base_boxes_itr != d_collection->d_base_boxes.end()) {
      ++d_base_boxes_itr;
      ++d_itr;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstIterator::operator == (
   const ConstIterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_itr == rhs.d_itr &&
          d_base_boxes_itr == rhs.d_base_boxes_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstIterator::operator != (
   const ConstIterator& rhs) const
{
   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator::NeighborIterator(
   Iterator& base_box_itr,
   bool from_start) :
   d_collection(base_box_itr.d_collection),
   d_base_box(base_box_itr.d_itr->first),
   d_itr(from_start ? base_box_itr.d_itr->second.begin() :
                      base_box_itr.d_itr->second.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator::NeighborIterator(
   const NeighborIterator& other) :
   d_collection(other.d_collection),
   d_base_box(other.d_base_box),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator&
BoxNeighborhoodCollection::NeighborIterator::operator = (
   const NeighborIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_base_box = rhs.d_base_box;
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
const Box* BoxNeighborhoodCollection::NeighborIterator::operator -> () const
{
   return *d_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::NeighborIterator&
BoxNeighborhoodCollection::NeighborIterator::operator ++ ()
{
   if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
      ++d_itr;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::NeighborIterator::operator == (
   const NeighborIterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_base_box == rhs.d_base_box &&
          d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::NeighborIterator::operator != (
   const NeighborIterator& rhs) const
{
   return !(*this == rhs);
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const ConstIterator& base_box_itr,
   bool from_start) :
   d_collection(base_box_itr.d_collection),
   d_base_box(base_box_itr.d_itr->first),
   d_itr(from_start ? base_box_itr.d_itr->second.begin() :
                      base_box_itr.d_itr->second.end())
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const ConstNeighborIterator& other) :
   d_collection(other.d_collection),
   d_base_box(other.d_base_box),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator::ConstNeighborIterator(
   const NeighborIterator& other) :
   d_collection(other.d_collection),
   d_base_box(other.d_base_box),
   d_itr(other.d_itr)
{
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator&
BoxNeighborhoodCollection::ConstNeighborIterator::operator = (
   const ConstNeighborIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_base_box = rhs.d_base_box;
   d_itr = rhs.d_itr;
   return *this;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator&
BoxNeighborhoodCollection::ConstNeighborIterator::operator = (
   const NeighborIterator& rhs)
{
   d_collection = rhs.d_collection;
   d_base_box = rhs.d_base_box;
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
const Box*
BoxNeighborhoodCollection::ConstNeighborIterator::operator -> () const
{
   return *d_itr;
}

SAMRAI_INLINE_KEYWORD
BoxNeighborhoodCollection::ConstNeighborIterator&
BoxNeighborhoodCollection::ConstNeighborIterator::operator ++ ()
{
   if (d_itr != d_collection->d_adj_list.find(d_base_box)->second.end()) {
      ++d_itr;
   }
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstNeighborIterator::operator == (
   const ConstNeighborIterator& rhs) const
{
   return d_collection == rhs.d_collection &&
          d_base_box == rhs.d_base_box &&
          d_itr == rhs.d_itr;
}

SAMRAI_INLINE_KEYWORD
bool BoxNeighborhoodCollection::ConstNeighborIterator::operator != (
   const ConstNeighborIterator& rhs) const
{
   return !(*this == rhs);
}

}
}
