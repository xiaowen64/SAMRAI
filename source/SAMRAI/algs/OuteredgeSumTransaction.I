/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for summing outeredge data
 *
 ************************************************************************/

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
void
OuteredgeSumTransaction::setRefineItems(
   const xfer::RefineClasses::Data** refine_items,
   int num_refine_items)
{
   TBOX_ASSERT(refine_items != (const xfer::RefineClasses::Data **)NULL);
   TBOX_ASSERT(num_refine_items >= 0);

   s_refine_items = refine_items;
   s_num_refine_items = num_refine_items;
}

SAMRAI_INLINE_KEYWORD
void
OuteredgeSumTransaction::unsetRefineItems()
{
   s_refine_items = (const xfer::RefineClasses::Data **)NULL;
   s_num_refine_items = 0;
}

}
}
