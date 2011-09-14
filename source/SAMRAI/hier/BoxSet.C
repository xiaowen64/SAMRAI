/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std
 *
 ************************************************************************/
#ifndef included_hier_BoxSet_C
#define included_hier_BoxSet_C

#include "SAMRAI/hier/BoxSet.h"

#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/BoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/BoxSetSingleOwnerIterator.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/BoxSet.I"
#endif

namespace SAMRAI {
namespace hier {

const int BoxSet::HIER_BOX_SET_VERSION = 0;

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

/*
 ***********************************************************************
 ***********************************************************************
 */
bool BoxSet::isLocallyEqual(
   const BoxSet& other,
   int rank) const
{
   if (this == &other) return true;

   bool is_equal = true;
   BoxSetSingleOwnerIterator this_iter(*this, rank);
   BoxSetSingleOwnerIterator other_iter(other, rank);

   while (this_iter.isValid() && other_iter.isValid()) {
      if (!(*this_iter).isIdEqual((*other_iter))) {
         is_equal = false;
         break;
      }
      this_iter++;
      other_iter++;
   }

   if (is_equal && (this_iter.isValid() != other_iter.isValid())) {
      // local count in this is different from local count in other.
      is_equal = false;
   }

   return is_equal;

}

/*
 ***********************************************************************
 * Write the BoxSet to a database.
 ***********************************************************************
 */
void BoxSet::putToDatabase(
   tbox::Database& database) const
{
   database.putInteger(
      "HIER_BOX_SET_VERSION", HIER_BOX_SET_VERSION);

   const int mbs_size = size();
   database.putInteger("mapped_box_set_size", mbs_size);
   if (mbs_size > 0) {

      std::vector<int> local_ids;
      std::vector<int> ranks;
      std::vector<int> block_ids;
      std::vector<int> periodic_ids;
      local_ids.reserve(mbs_size);
      ranks.reserve(mbs_size);
      block_ids.reserve(mbs_size);
      periodic_ids.reserve(mbs_size);

      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      int counter = -1;
      for (BoxSet::const_iterator ni = begin();
           ni != end(); ++ni) {
         local_ids.push_back(ni->getLocalId().getValue());
         ranks.push_back(ni->getOwnerRank());
         block_ids.push_back(ni->getBlockId().getBlockValue());
         periodic_ids.push_back(ni->getPeriodicId().getPeriodicValue());
         db_box_array[++counter] = *ni;
      }

      database.putIntegerArray(
         "local_indices", &local_ids[0], mbs_size);
      database.putIntegerArray(
         "ranks", &ranks[0], mbs_size);
      database.putIntegerArray(
         "block_ids", &block_ids[0], mbs_size);
      database.putIntegerArray(
         "periodic_ids", &periodic_ids[0], mbs_size);
      database.putDatabaseBoxArray(
         "boxes", &db_box_array[0], mbs_size);
   }
}

/*
 ***********************************************************************
 * Read the BoxSet from a database.
 ***********************************************************************
 */
void BoxSet::getFromDatabase(
   tbox::Database& database)
{
   const unsigned int mbs_size = database.getInteger("mapped_box_set_size");
   if (mbs_size > 0) {
      std::vector<int> local_ids(mbs_size);
      std::vector<int> ranks(mbs_size);
      std::vector<int> block_ids(mbs_size);
      std::vector<int> periodic_ids(mbs_size);
      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      database.getIntegerArray(
         "local_indices", &local_ids[0], mbs_size);
      database.getIntegerArray(
         "ranks", &ranks[0], mbs_size);
      database.getIntegerArray(
         "block_ids", &block_ids[0], mbs_size);
      database.getIntegerArray(
         "periodic_ids", &periodic_ids[0], mbs_size);
      database.getDatabaseBoxArray(
         "boxes", &db_box_array[0], mbs_size);

      for (unsigned int i = 0; i < mbs_size; ++i) {
         Box box(db_box_array[i]);
         Box mapped_box(
            box,
            LocalId(local_ids[i]),
            ranks[i],
            BlockId(block_ids[i]),
            PeriodicId(periodic_ids[i]));
         insert(end(), mapped_box);
      }
   }
}

/*
 ***********************************************************************
 * Construct the BoxList consisting of the Boxes in this BoxSet
 * in the requested block.
 ***********************************************************************
 */
tbox::Pointer<BoxList>
BoxSet::getSingleBlockBoxList(
   const tbox::Dimension& dim,
   const BlockId& which_block) const
{
   BoxSetSingleBlockIterator itr(*this, which_block);
   BoxList* boxes_in_block = new BoxList(dim);
   while (itr.isValid()) {
      const Box& mapped_box = *itr;
      TBOX_ASSERT(dim == mapped_box.getDim());
      boxes_in_block->pushBack(mapped_box);
      ++itr;
   }
   return tbox::Pointer<BoxList>(boxes_in_block);
}

/*
 ***********************************************************************
 * Refine the boxes of a BoxSet.
 ***********************************************************************
 */
void BoxSet::refine(
   BoxSet& output_mapped_boxes,
   const IntVector& ratio) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.refine(ratio);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      BoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.refine(ratio);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Coarsen the boxes of a BoxSet.
 ***********************************************************************
 */
void BoxSet::coarsen(
   BoxSet& output_mapped_boxes,
   const IntVector& ratio) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.coarsen(ratio);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      BoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.coarsen(ratio);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Grow the boxes of a BoxSet.
 ***********************************************************************
 */
void BoxSet::grow(
   BoxSet& output_mapped_boxes,
   const IntVector& growth) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.grow(growth);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      BoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         Box n = *na;
         n.grow(growth);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Remove periodic image Boxes from a BoxSet.
 ***********************************************************************
 */
void BoxSet::removePeriodicImageBoxes(
   BoxSet& output_mapped_boxes) const
{
   iterator hint = output_mapped_boxes.begin();
   for (const_iterator na = begin(); na != end(); ++na) {
      const Box& n = *na;
      if (!n.isPeriodicImage()) {
         hint = output_mapped_boxes.insert(hint, n);
      }
   }
}

/*
 ***********************************************************************
 * Unshift periodic image Boxes from a BoxSet.
 ***********************************************************************
 */
void BoxSet::unshiftPeriodicImageBoxes(
   BoxSet& output_mapped_boxes,
   const IntVector& refinement_ratio) const
{
   iterator hint = output_mapped_boxes.begin();

   if (!empty()) {
      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      for (const_iterator na = begin(); na != end(); ++na) {
         if (na->isPeriodicImage()) {
            const Box unshifted_mapped_box(
               *na, zero_shift_number, refinement_ratio);
            hint = output_mapped_boxes.insert(hint, unshifted_mapped_box);
         } else {
            hint = output_mapped_boxes.insert(hint, *na);
         }
      }
   }
}

/*
 ***********************************************************************
 * Remove from a BoxList its intersection with a BoxSet.
 ***********************************************************************
 */
void BoxSet::removeBoxListIntersections(
   BoxList& boxes) const
{
   for (const_iterator na = begin(); na != end(); ++na) {
      const Box& nabr = *na;
      boxes.removeIntersections(nabr);
   }
}

/*
 ***********************************************************************
 * Insert Box owners into a single set container.
 ***********************************************************************
 */
void BoxSet::getOwners(
   std::set<int>& owners) const
{
   for (const_iterator i_nabr = begin(); i_nabr != end(); ++i_nabr) {
      const int owner = (*i_nabr).getOwnerRank();
      owners.insert(owner);
   }
}

/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 ***********************************************************************
 */
void BoxSet::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   NULL_USE(detail_depth);
   for (const_iterator bi = begin(); bi != end(); ++bi) {
      Box mapped_box = *bi;
      co << border << "    "
         << mapped_box << "   "
         << mapped_box.numberCells() << '\n';
   }
}

/*
 ***********************************************************************
 * Construct a BoxSet Outputter with formatting parameters.
 ***********************************************************************
 */

BoxSet::Outputter::Outputter(
   const BoxSet& mapped_box_set,
   const std::string& border,
   int detail_depth):
   d_set(mapped_box_set),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a BoxSet according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
   std::ostream& s,
   const BoxSet::Outputter& format)
{
   format.d_set.recursivePrint(s, format.d_border, format.d_detail_depth);
   return s;
}

/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxSet to a stream.
 ***********************************************************************
 */

BoxSet::Outputter BoxSet::format(
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
