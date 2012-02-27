/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
NeighborhoodSet&
NeighborhoodSet::operator = (
   const NeighborhoodSet& r)
{
   d_map = r.d_map;
   return *this;
}

SAMRAI_INLINE_KEYWORD
bool
NeighborhoodSet::operator == (
   const NeighborhoodSet& r) const
{
   return (d_map == r.d_map);
}

SAMRAI_INLINE_KEYWORD
bool
NeighborhoodSet::operator != (
   const NeighborhoodSet& r) const
{
   return (d_map != r.d_map);
}

SAMRAI_INLINE_KEYWORD
void
NeighborhoodSet::swap(
   NeighborhoodSet& r)
{
   d_map.swap(r.d_map);
}

/*
 ***********************************************************************
 * Removes periodic neighbors in a NeighborhoodSet.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
NeighborhoodSet::removePeriodicNeighbors()
{
   for (iterator ei = begin(); ei != end(); ++ei) {
      ei->second.removePeriodicImageBoxes();
   }
}

/*
 ***********************************************************************
 * Print out a NeighborhoodSet according to settings in the Outputter.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
std::ostream&
operator << (
   std::ostream& s,
   const NeighborhoodSet::Outputter& format)
{
   format.d_set.recursivePrint(s, format.d_border, format.d_detail_depth);
   return s;
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the NeighborhoodSet to a stream.
 ***********************************************************************
 */
SAMRAI_INLINE_KEYWORD
NeighborhoodSet::Outputter
NeighborhoodSet::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::iterator
NeighborhoodSet::begin() {
   return d_map.begin();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::iterator
NeighborhoodSet::end() {
   return d_map.end();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::const_iterator
NeighborhoodSet::begin() const {
   return d_map.begin();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::const_iterator
NeighborhoodSet::end() const {
   return d_map.end();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::reverse_iterator
NeighborhoodSet::rbegin() {
   return d_map.rbegin();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::reverse_iterator
NeighborhoodSet::rend() {
   return d_map.rend();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::iterator
NeighborhoodSet::insert(
   iterator i,
   const value_type& v) {
   return d_map.insert(i, v);
}

SAMRAI_INLINE_KEYWORD
void
NeighborhoodSet::erase(
   iterator i) {
   d_map.erase(i);
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::size_type
NeighborhoodSet::erase(
   const key_type& k) {
   TBOX_ASSERT(k.getPeriodicId().getPeriodicValue() == 0);
   return d_map.erase(k);
}

SAMRAI_INLINE_KEYWORD
void
NeighborhoodSet::erase(
   iterator first,
   iterator last) {
   d_map.erase(first, last);
}

SAMRAI_INLINE_KEYWORD
size_t
NeighborhoodSet::size() const {
   return d_map.size();
}

SAMRAI_INLINE_KEYWORD
bool
NeighborhoodSet::empty() const {
   return d_map.empty();
}

SAMRAI_INLINE_KEYWORD
void
NeighborhoodSet::clear() {
   d_map.clear();
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::NeighborSet&
NeighborhoodSet::operator [] (
   const key_type& k) {
   TBOX_ASSERT(k.getPeriodicId().getPeriodicValue() == 0);
   return d_map[k];
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::iterator
NeighborhoodSet::find(
   const key_type& k) {
   TBOX_ASSERT(k.getPeriodicId().getPeriodicValue() == 0);
   return d_map.find(k);
}

SAMRAI_INLINE_KEYWORD
NeighborhoodSet::const_iterator
NeighborhoodSet::find(
   const key_type& k) const {
   TBOX_ASSERT(k.getPeriodicId().getPeriodicValue() == 0);
   return d_map.find(k);
}

}
}
