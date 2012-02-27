/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for time interpolation during data
 *                refining.
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

SAMRAI_INLINE_KEYWORD
void
RefineTimeTransaction::setTransactionTime(
   const double time)
{
   s_time = time;
}

SAMRAI_INLINE_KEYWORD
void
RefineTimeTransaction::setRefineItems(
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
RefineTimeTransaction::unsetRefineItems()
{
   s_refine_items = (const RefineClasses::Data **)NULL;
   s_num_refine_items = 0;
}

}
}
