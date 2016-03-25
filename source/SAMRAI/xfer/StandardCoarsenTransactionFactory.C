/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory for create standard copy transactions
 *                for coarsen schedules.
 *
 ************************************************************************/

#ifndef included_xfer_StandardCoarsenTransactionFactory_C
#define included_xfer_StandardCoarsenTransactionFactory_C

#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"

#include "SAMRAI/xfer/CoarsenCopyTransaction.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

StandardCoarsenTransactionFactory::StandardCoarsenTransactionFactory()
{
}

StandardCoarsenTransactionFactory::~StandardCoarsenTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Set/unset information for transactions managed by this factory class.
 *
 *************************************************************************
 */

void StandardCoarsenTransactionFactory::setCoarsenItems(
   const CoarsenClasses::Data** coarsen_items,
   int num_coarsen_items)
{
   xfer::CoarsenCopyTransaction::setCoarsenItems(coarsen_items,
      num_coarsen_items);
   d_coarsen_items = coarsen_items;
   d_num_coarsen_items = num_coarsen_items;
}

void StandardCoarsenTransactionFactory::unsetCoarsenItems()
{
   xfer::CoarsenCopyTransaction::unsetCoarsenItems();
   d_coarsen_items = (const xfer::CoarsenClasses::Data **)NULL;
   d_num_coarsen_items = 0;
}

/*
 *************************************************************************
 *
 * Allocate appropriate transaction object.
 *
 *************************************************************************
 */

tbox::Pointer<tbox::Transaction>
StandardCoarsenTransactionFactory::allocate(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   tbox::Pointer<hier::BoxOverlap> overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   int citem_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);

   CoarsenCopyTransaction* transaction =
      new CoarsenCopyTransaction(dst_level, src_level,
         overlap,
         dst_mapped_box,
         src_mapped_box,
         citem_id);
   return tbox::Pointer<tbox::Transaction>(transaction);
}

}
}
#endif
