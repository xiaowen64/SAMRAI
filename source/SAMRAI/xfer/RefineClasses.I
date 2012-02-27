/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence
 *                classes.
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
int
RefineClasses::getNumberOfRefineItems() const
{
   return d_num_refine_items;
}

SAMRAI_INLINE_KEYWORD
RefineClasses::Data&
RefineClasses::getRefineItem(
   const int refine_items_array_id)
{
   return d_refine_classes_data_items[refine_items_array_id];
}

SAMRAI_INLINE_KEYWORD
int
RefineClasses::getNumberOfEquivalenceClasses() const
{
   return d_equivalence_class_indices.size();
}

SAMRAI_INLINE_KEYWORD
int
RefineClasses::getRefineItemArraySize() const
{
   return d_refine_classes_data_items.size();
}

/*
 *************************************************************************
 *
 * Return representative item for given equivalence class (first in list)
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const RefineClasses::Data&
RefineClasses::getClassRepresentative(
   int equiv_class_index) const
{
   TBOX_ASSERT((equiv_class_index >= 0) &&
      (equiv_class_index < getNumberOfEquivalenceClasses()));
   return d_refine_classes_data_items[
             d_equivalence_class_indices[equiv_class_index].getFirstItem()];
}

/*
 *************************************************************************
 *
 * Return iterator for list of refine items for given equivalence class
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
tbox::List<int>::Iterator
RefineClasses::getIterator(
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
RefineClasses::increaseRefineItemArraySize(
   const int size)
{
   if (size > d_refine_classes_data_items.size()) {
      d_refine_classes_data_items.resizeArray(size);
   }
}

}
}
