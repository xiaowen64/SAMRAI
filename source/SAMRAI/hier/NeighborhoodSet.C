/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#ifndef included_hier_NeighborhoodSet_C
#define included_hier_NeighborhoodSet_C

#include "SAMRAI/hier/NeighborhoodSet.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"

#include "SAMRAI/tbox/Utilities.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/NeighborhoodSet.I"
#endif

namespace SAMRAI {
namespace hier {

const int NeighborhoodSet::HIER_EDGE_SET_VERSION = 0;

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 ***********************************************************************
 * Find all the Boxes corresponding to a given rank as a range in
 * the container.
 ***********************************************************************
 */
NeighborhoodSet::Range
NeighborhoodSet::findRanksRange(
   int rank)
{
   Range rval(d_map.end(), d_map.end());
   if (!d_map.empty()) {
      /*
       * start_key is an iterator pointing to the first Box in
       * the range.  stop_key is an iterator pointing just past the
       * last Box in the range.
       */
      BoxId start_key(LocalId::getZero(), rank);
      rval.first = d_map.lower_bound(start_key);
      BoxId stop_key(LocalId::getZero(), rank + 1);
      rval.second = d_map.lower_bound(stop_key);
   }
   return rval;
}

/*
 ***********************************************************************
 * Find all the Boxes corresponding to a given rank as a range in
 * the container.
 ***********************************************************************
 */
NeighborhoodSet::ConstRange
NeighborhoodSet::findRanksRange(
   int rank) const
{
   ConstRange rval(d_map.end(), d_map.end());
   if (!d_map.empty()) {
      /*
       * start_key is an iterator pointing to the first Box in
       * the range.  stop_key is an iterator pointing just past the
       * last Box in the range.
       */
      BoxId start_key(LocalId::getZero(), rank);
      rval.first = d_map.lower_bound(start_key);
      BoxId stop_key(LocalId::getZero(), rank + 1);
      rval.second = d_map.lower_bound(stop_key);
   }
   return rval;
}

/*
 ***********************************************************************
 * Coarsen the Neighbors in a NeighborhoodSet.
 ***********************************************************************
 */
void
NeighborhoodSet::coarsenNeighbors(
   NeighborhoodSet& output_edges,
   const IntVector& ratio) const
{
   for (const_iterator ei = begin(); ei != end(); ++ei) {
      const BoxId& mapped_box_id = (*ei).first;
      const NeighborSet& orig_nabrs = (*ei).second;
      orig_nabrs.coarsen(output_edges[mapped_box_id], ratio);
   }
}

/*
 ***********************************************************************
 * Refine the Neighbors in a NeighborhoodSet.
 ***********************************************************************
 */
void
NeighborhoodSet::refineNeighbors(
   NeighborhoodSet& output_edges,
   const IntVector& ratio) const
{
   for (const_iterator ei = begin(); ei != end(); ++ei) {
      const BoxId& mapped_box_id = (*ei).first;
      const NeighborSet& orig_nabrs = (*ei).second;
      orig_nabrs.refine(output_edges[mapped_box_id], ratio);
   }
}

/*
 ***********************************************************************
 * Grow the Neighbors in a NeighborhoodSet.
 ***********************************************************************
 */
void
NeighborhoodSet::growNeighbors(
   NeighborhoodSet& output_edges,
   const IntVector& growth) const
{
   for (const_iterator ei = begin(); ei != end(); ++ei) {
      const BoxId& mapped_box_id = (*ei).first;
      const NeighborSet& orig_nabrs = (*ei).second;
      orig_nabrs.grow(output_edges[mapped_box_id], growth);
   }
}

/*
 ***********************************************************************
 * Removes periodic neighbors in a NeighborhoodSet.
 ***********************************************************************
 */
void
NeighborhoodSet::removePeriodicNeighbors()
{
   for (iterator ei = begin(); ei != end(); ++ei) {
      ei->second.removePeriodicImageBoxes();
   }
}

/*!
 ***********************************************************************
 * Insert all neighbors from a NeighborhoodSet into a single NeighborSet.
 ***********************************************************************
 */
void
NeighborhoodSet::getNeighbors(
   NeighborSet& all_nabrs) const
{
   for (const_iterator ei = begin(); ei != end(); ++ei) {
      const NeighborSet& nabrs = (*ei).second;
      all_nabrs.insert(nabrs.begin(), nabrs.end());
   }
}

/*!
 ***********************************************************************
 * Insert all neighbors from a NeighborhoodSet into a single NeighborSet.
 ***********************************************************************
 */
void
NeighborhoodSet::getNeighbors(
   BoxList& all_nabrs) const
{
   NeighborSet tmp_nabrs;
   getNeighbors(tmp_nabrs);
   for (BoxSet::const_iterator ei = tmp_nabrs.begin();
        ei != tmp_nabrs.end(); ++ei) {
      all_nabrs.appendItem(*ei);
   }
}

/*!
 ***********************************************************************
 * Insert all neighbors from a NeighborhoodSet into a single NeighborSet.
 ***********************************************************************
 */
void
NeighborhoodSet::getNeighbors(
   BoxList& all_nabrs,
   const BlockId& block_id) const
{
   NeighborSet tmp_nabrs;
   getNeighbors(tmp_nabrs);
   for (BoxSetSingleBlockIterator ei(tmp_nabrs, block_id);
        ei.isValid(); ++ei) {
      all_nabrs.appendItem(*ei);
   }
}

/*!
 ***********************************************************************
 * Insert all neighbors from a NeighborhoodSet into a single NeighborSet.
 ***********************************************************************
 */
void
NeighborhoodSet::getNeighbors(
   std::map<BlockId, BoxList>& all_nabrs) const
{
   NeighborSet tmp_nabrs;
   getNeighbors(tmp_nabrs);
   for (BoxSet::const_iterator ei = tmp_nabrs.begin();
        ei != tmp_nabrs.end(); ++ei) {
      all_nabrs[ei->getBlockId()].appendItem(*ei);
   }
}

/*!
 ***********************************************************************
 * Insert all owners of neighbors from a NeighborhoodSet into a single
 * set container.
 ***********************************************************************
 */
void
NeighborhoodSet::getOwners(
   std::set<int>& owners) const
{
   for (const_iterator i_nabrhood = begin(); i_nabrhood != end(); ++i_nabrhood) {
      const value_type& nabrhood = *i_nabrhood;
      const NeighborSet& nabr = nabrhood.second;
      nabr.getOwners(owners);
   }
}

/*
 ***********************************************************************
 * Write the NeighborhoodSet to a database.
 ***********************************************************************
 */

void NeighborhoodSet::putToDatabase(
   tbox::Database& database) const
{
   // This appears to be used in the RedistributedRestartUtility.
   database.putBool("d_is_edge_set", true);

   database.putInteger(
      "HIER_EDGE_SET_VERSION", HIER_EDGE_SET_VERSION);
   database.putInteger("number_of_sets", size());

   if (!empty()) {

      std::vector<int> block_ids;
      std::vector<int> owners;
      std::vector<int> local_indices;
      std::vector<int> periodic_ids;
      block_ids.reserve(size());
      owners.reserve(size());
      local_indices.reserve(size());
      periodic_ids.reserve(size());
      for (const_iterator ei = begin(); ei != end(); ++ei) {
         block_ids.push_back(ei->first.getBlockId().getBlockValue());
         owners.push_back(ei->first.getOwnerRank());
         local_indices.push_back(ei->first.getLocalId().getValue());
         periodic_ids.push_back(ei->first.getPeriodicId().getPeriodicValue());
      }

      database.putIntegerArray("block_ids", &block_ids[0], size());
      database.putIntegerArray("owners", &owners[0], size());
      database.putIntegerArray("local_indices", &local_indices[0], size());
      database.putIntegerArray("periodic_ids", &periodic_ids[0], size());

      const std::string set_db_string("set_for_local_id_");

      for (const_iterator ei = begin(); ei != end(); ++ei) {
         const BoxId& mbid = ei->first;
         const NeighborSet& mapped_boxes = ei->second;
         const std::string set_name =
            set_db_string
            + tbox::Utilities::blockToString(mbid.getBlockId().getBlockValue())
            + tbox::Utilities::processorToString(mbid.getOwnerRank())
            + tbox::Utilities::patchToString(mbid.getLocalId().getValue())
            + tbox::Utilities::intToString(mbid.getPeriodicId().getPeriodicValue());
         mapped_boxes.putToDatabase(*database.putDatabase(set_name));
      }

   }
}

/*
 ***********************************************************************
 * Read the NeighborhoodSet from a database.
 ***********************************************************************
 */

void NeighborhoodSet::getFromDatabase(
   tbox::Database& database)
{
   const unsigned int number_of_sets = database.getInteger("number_of_sets");
   if (number_of_sets > 0) {

      std::vector<int> block_ids(number_of_sets);
      std::vector<int> owners(number_of_sets);
      std::vector<int> local_indices(number_of_sets);
      std::vector<int> periodic_ids(number_of_sets);
      database.getIntegerArray("block_ids", &block_ids[0], number_of_sets);
      database.getIntegerArray("owners", &owners[0], number_of_sets);
      database.getIntegerArray("local_indices", &local_indices[0], number_of_sets);
      database.getIntegerArray("periodic_ids", &periodic_ids[0], number_of_sets);

      const std::string set_db_string("set_for_local_id_");

      std::pair<BoxId, NeighborSet> tmp_pair;
      for (size_t i = 0; i < number_of_sets; ++i) {
         tmp_pair.first.initialize(LocalId(local_indices[i]),
            owners[i],
            BlockId(block_ids[i]),
            PeriodicId(periodic_ids[i]));
         const std::string set_name =
            set_db_string
            + tbox::Utilities::blockToString(tmp_pair.first.getBlockId().getBlockValue())
            + tbox::Utilities::processorToString(tmp_pair.first.getOwnerRank())
            + tbox::Utilities::patchToString(tmp_pair.first.getLocalId().getValue())
            + tbox::Utilities::intToString(tmp_pair.first.getPeriodicId().getPeriodicValue());
         iterator mi = insert(end(), tmp_pair);
         mi->second.getFromDatabase(*database.getDatabase(set_name));
      }

   }
}

/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 * Print out global bounding box only if it has been computed already.
 ***********************************************************************
 */

void NeighborhoodSet::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   const std::string indented_border = border + "  ";
   co << border << "  " << size() << " neigborhoods:\n";
   for (NeighborhoodSet::const_iterator ei = begin(); ei != end(); ++ei) {
      const BoxId& mbid = ei->first;
      const NeighborSet& nabrs = ei->second;
      co << border << "  " << mbid << "\n";
      co << border << "    Neighbors (" << nabrs.size() << "):\n";
      if (detail_depth > 1) {
         nabrs.recursivePrint(co, indented_border, detail_depth - 1);
      }
   }
}

/*
 ***********************************************************************
 * Construct a NeighborhoodSet Outputter with formatting parameters.
 ***********************************************************************
 */

NeighborhoodSet::Outputter::Outputter(
   const NeighborhoodSet& neighborhood_set,
   const std::string& border,
   int detail_depth):
   d_set(neighborhood_set),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a NeighborhoodSet according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
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

NeighborhoodSet::Outputter NeighborhoodSet::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
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

#endif
