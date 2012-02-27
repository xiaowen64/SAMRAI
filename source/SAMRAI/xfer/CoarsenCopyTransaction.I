/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data
 *                coarsening
 *
 ************************************************************************/

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Initialization, set/unset functions for static array of coarsen items.
 *
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
void
CoarsenCopyTransaction::setCoarsenItems(
   const CoarsenClasses::Data** coarsen_items,
   int num_coarsen_items)
{
   TBOX_ASSERT(coarsen_items != (const CoarsenClasses::Data **)NULL);
   TBOX_ASSERT(num_coarsen_items >= 0);

   s_coarsen_items = coarsen_items;
   s_num_coarsen_items = num_coarsen_items;
}

SAMRAI_INLINE_KEYWORD
void
CoarsenCopyTransaction::unsetCoarsenItems()
{
   s_coarsen_items = (const CoarsenClasses::Data **)NULL;
   s_num_coarsen_items = 0;
}

}
}
