/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Extension of a std 
 *
 ************************************************************************/
#ifndef included_hier_MappedBoxSet_C
#define included_hier_MappedBoxSet_C

#include "SAMRAI/hier/MappedBoxSet.h"

#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/hier/MappedBoxSetSingleBlockIterator.h"
#include "SAMRAI/hier/MappedBoxSetSingleOwnerIterator.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"

#ifndef SAMRAI_INLINE
#include "SAMRAI/hier/MappedBoxSet.I"
#endif

namespace SAMRAI {
namespace hier {

const int MappedBoxSet::HIER_MAPPED_BOX_SET_VERSION = 0;

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
bool MappedBoxSet::isLocallyEqual(
   const MappedBoxSet& other,
   int rank) const
{
   if (this == &other) return true;

   bool is_equal = true;
   MappedBoxSetSingleOwnerIterator this_iter(*this, rank);
   MappedBoxSetSingleOwnerIterator other_iter(other, rank);

   while (this_iter.isValid() && other_iter.isValid()) {
      if ((*this_iter) != (*other_iter)) {
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
 * Write the MappedBoxSet to a database.
 ***********************************************************************
 */
void MappedBoxSet::putToDatabase(
   tbox::Database& database) const
{
   database.putInteger(
      "HIER_MAPPED_BOX_SET_VERSION", HIER_MAPPED_BOX_SET_VERSION);

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
      for (MappedBoxSet::const_iterator ni = begin();
           ni != end(); ++ni) {
         local_ids.push_back(ni->getLocalId().getValue());
         ranks.push_back(ni->getOwnerRank());
         block_ids.push_back(ni->getBlockId().getBlockValue());
         periodic_ids.push_back(ni->getPeriodicId().getPeriodicValue());
         db_box_array[++counter] = ni->getBox();
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
 * Read the MappedBoxSet from a database.
 ***********************************************************************
 */
void MappedBoxSet::getFromDatabase(
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
         MappedBox mapped_box(
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
 * Construct the BoxList consisting of the Boxes in this MappedBoxSet
 * in the requested block.
 ***********************************************************************
 */
tbox::Pointer<BoxList>
MappedBoxSet::getSingleBlockBoxList(
   const tbox::Dimension& dim,
   const BlockId& which_block) const
{
   MappedBoxSetSingleBlockIterator itr(*this, which_block);
   BoxList* boxes_in_block = new BoxList(dim);
   while (itr.isValid()) {
      const MappedBox& mapped_box = *itr;
      TBOX_ASSERT(dim == mapped_box.getDim());
      boxes_in_block->appendItem(mapped_box.getBox());
      ++itr;
   }
   return tbox::Pointer<BoxList>(boxes_in_block);
}

/*
 ***********************************************************************
 * Refine the boxes of a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::refine(
   MappedBoxSet& output_mapped_boxes,
   const IntVector& ratio) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().refine(ratio);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      MappedBoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().refine(ratio);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Coarsen the boxes of a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::coarsen(
   MappedBoxSet& output_mapped_boxes,
   const IntVector& ratio) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().coarsen(ratio);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      MappedBoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().coarsen(ratio);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Grow the boxes of a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::grow(
   MappedBoxSet& output_mapped_boxes,
   const IntVector& growth) const
{
   if (this != &output_mapped_boxes) {
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().grow(growth);
         output_mapped_boxes.insert(output_mapped_boxes.end(), n);
      }
   } else {
      MappedBoxSet tmp_mapped_boxes;
      for (const_iterator na = begin(); na != end(); ++na) {
         MappedBox n = *na;
         n.getBox().grow(growth);
         tmp_mapped_boxes.insert(tmp_mapped_boxes.end(), n);
      }
      output_mapped_boxes.swap(tmp_mapped_boxes);
   }
}

/*
 ***********************************************************************
 * Remove periodic image MappedBoxes from a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::removePeriodicImageMappedBoxes(
   MappedBoxSet& output_mapped_boxes) const
{
   iterator hint = output_mapped_boxes.begin();
   for (const_iterator na = begin(); na != end(); ++na) {
      const MappedBox& n = *na;
      if (!n.isPeriodicImage()) {
         hint = output_mapped_boxes.insert(hint, n);
      }
   }
}

/*
 ***********************************************************************
 * Unshift periodic image MappedBoxes from a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::unshiftPeriodicImageMappedBoxes(
   MappedBoxSet& output_mapped_boxes,
   const IntVector& refinement_ratio) const
{
   iterator hint = output_mapped_boxes.begin();

   if (!empty()) {
      const MappedBox& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                     first_element.getDim())->
                                  getZeroShiftNumber());

      for (const_iterator na = begin(); na != end(); ++na) {
         if (na->isPeriodicImage()) {
            const MappedBox unshifted_mapped_box(
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
 * Remove from a BoxList its intersection with a MappedBoxSet.
 ***********************************************************************
 */
void MappedBoxSet::removeBoxListIntersections(
   BoxList& boxes) const
{
   for (const_iterator na = begin(); na != end(); ++na) {
      const MappedBox& nabr = *na;
      boxes.removeIntersections(nabr.getBox());
   }
}

/*
 ***********************************************************************
 * Insert MappedBox owners into a single set container.
 ***********************************************************************
 */
void MappedBoxSet::getOwners(
   std::set<int>& owners) const
{
   for (const_iterator i_nabr = begin(); i_nabr != end(); ++i_nabr) {
      const int owner = (*i_nabr).getOwnerRank();
      owners.insert(owner);
   }
}

/*
 ***********************************************************************
 * Convert a MappedBoxSet to a MappedBoxList.
 ***********************************************************************
 */
void MappedBoxSet::convertToBoxList(
   BoxList& box_list) const
{
   for (const_iterator ni = begin(); ni != end(); ++ni) {
      box_list.appendItem(ni->getBox());
   }
}

/*
 ***********************************************************************
 * Avoid communication in this method.  It is often used for debugging.
 ***********************************************************************
 */
void MappedBoxSet::recursivePrint(
   std::ostream& co,
   const std::string& border,
   int detail_depth) const
{
   NULL_USE(detail_depth);
   for (const_iterator bi = begin(); bi != end(); ++bi) {
      MappedBox mapped_box = *bi;
      co << border << "    "
         << mapped_box << "   "
         << mapped_box.getBox().numberCells() << '\n';
   }
}


/*
 ***********************************************************************
 * Construct a MappedBoxSet Outputter with formatting parameters.
 ***********************************************************************
 */

MappedBoxSet::Outputter::Outputter(
   const MappedBoxSet &mapped_box_set,
   const std::string& border,
   int detail_depth )
   : d_set(mapped_box_set),
     d_border(border),
     d_detail_depth(detail_depth)
{
   return;
}


/*
 ***********************************************************************
 * Print out a MappedBoxSet according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream& operator << (
   std::ostream& s,
   const MappedBoxSet::Outputter& format)
{
   format.d_set.recursivePrint( s, format.d_border, format.d_detail_depth );
   return s;
}


/*
 ***********************************************************************
 * Return a Outputter that can dump the MappedBoxSet to a stream.
 ***********************************************************************
 */

MappedBoxSet::Outputter MappedBoxSet::format(
   const std::string& border,
   int detail_depth ) const
{
   return Outputter( *this, border, detail_depth);
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
