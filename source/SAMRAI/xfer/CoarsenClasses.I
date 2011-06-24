/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Simple structure for managing coarsening data in equivalence classes. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
int CoarsenClasses::getNumberOfCoarsenItems() const
{
   return d_num_coarsen_items;
}

SAMRAI_INLINE_KEYWORD
CoarsenClasses::Data& CoarsenClasses::getCoarsenItem(
   const int coarsen_items_array_id)
{
   return d_coarsen_classes_data_items[coarsen_items_array_id];
}

SAMRAI_INLINE_KEYWORD
int CoarsenClasses::getNumberOfEquivalenceClasses() const
{
   return d_equivalence_class_ids.size();
}

SAMRAI_INLINE_KEYWORD
int CoarsenClasses::getCoarsenItemArraySize() const
{
   return d_coarsen_classes_data_items.size();
}

}
}
