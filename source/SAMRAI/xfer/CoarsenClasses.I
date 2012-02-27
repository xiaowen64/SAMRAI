/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes.
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
int
CoarsenClasses::getNumberOfCoarsenItems() const
{
   return d_num_coarsen_items;
}

SAMRAI_INLINE_KEYWORD
CoarsenClasses::Data&
CoarsenClasses::getCoarsenItem(
   const int coarsen_items_array_id)
{
   return d_coarsen_classes_data_items[coarsen_items_array_id];
}

SAMRAI_INLINE_KEYWORD
int
CoarsenClasses::getNumberOfEquivalenceClasses() const
{
   return d_equivalence_class_indices.size();
}

SAMRAI_INLINE_KEYWORD
int
CoarsenClasses::getCoarsenItemArraySize() const
{
   return d_coarsen_classes_data_items.size();
}

/*
 *************************************************************************
 *
 * Return representative item for given equivalence class (first in list)
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const CoarsenClasses::Data&
CoarsenClasses::getClassRepresentative(
   int equiv_class_index) const
{
   TBOX_ASSERT((equiv_class_index >= 0) &&
      (equiv_class_index < getNumberOfEquivalenceClasses()));
   return d_coarsen_classes_data_items[
             d_equivalence_class_indices[equiv_class_index].getFirstItem()];
}

/*
 *************************************************************************
 *
 * Return iterator for list of coarsen items for given equivalence class
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
tbox::List<int>::Iterator
CoarsenClasses::getIterator(
   int equiv_class_index)
{
   TBOX_ASSERT((equiv_class_index >= 0) &&
      (equiv_class_index < getNumberOfEquivalenceClasses()));
   return tbox::List<int>::
          Iterator(d_equivalence_class_indices[equiv_class_index]);
}

/*
 *************************************************************************
 *
 * Increase the data items array to the specified size.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CoarsenClasses::increaseCoarsenItemArraySize(
   const int size,
   const tbox::Dimension& dim)
{
   if (size > d_coarsen_classes_data_items.size()) {
      d_coarsen_classes_data_items.resizeArray(size, Data(dim));
   }
}

}
}
