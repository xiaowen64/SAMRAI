/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing refinement data in equivalence classes. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
int RefineClasses::getNumberOfRefineItems() const
{
   return d_num_refine_items;
}

SAMRAI_INLINE_KEYWORD
RefineClasses::Data& RefineClasses::getRefineItem(
   const int refine_items_array_id)
{
   return d_refine_classes_data_items[refine_items_array_id];
}

SAMRAI_INLINE_KEYWORD
int RefineClasses::getNumberOfEquivalenceClasses() const
{
   return d_equivalence_class_indices.size();
}

SAMRAI_INLINE_KEYWORD
int RefineClasses::getRefineItemArraySize() const
{
   return d_refine_classes_data_items.size();
}

}
}
