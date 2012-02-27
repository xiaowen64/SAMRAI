/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data
 *                refining
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Initialization, set/unset functions for static array of refine items.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
RefineCopyTransaction::setRefineItems(
   const RefineClasses::Data** refine_items,
   int num_refine_items)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(refine_items != (const RefineClasses::Data **)NULL);
   TBOX_ASSERT(num_refine_items >= 0);
#endif
   s_refine_items = refine_items;
   s_num_refine_items = num_refine_items;
}

SAMRAI_INLINE_KEYWORD
void
RefineCopyTransaction::unsetRefineItems()
{
   s_refine_items = (const RefineClasses::Data **)NULL;
   s_num_refine_items = 0;
}

}
}
