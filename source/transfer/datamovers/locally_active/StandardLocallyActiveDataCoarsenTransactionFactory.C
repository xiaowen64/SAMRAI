//
// File:	StandardLocallyActiveDataCoarsenTransactionFactory.C
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 684 $
// Modified:	$Date: 2005-10-21 14:59:38 -0700 (Fri, 21 Oct 2005) $
// Description:	Concrete factory for create standard copy transactions
//              for locally-active data coarsen schedules.
//

#ifndef included_xfer_StandardLocallyActiveDataCoarsenTransactionFactory_C
#define included_xfer_StandardLocallyActiveDataCoarsenTransactionFactory_C

#include "StandardLocallyActiveDataCoarsenTransactionFactory.h"

#include "CoarsenCopyTransaction.h"
#include "tbox/ArenaManager.h"

namespace SAMRAI {
    namespace xfer {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   * 
*                                                                       *
*************************************************************************
*/

template<int DIM> 
StandardLocallyActiveDataCoarsenTransactionFactory<DIM>::
   StandardLocallyActiveDataCoarsenTransactionFactory()
{
}

template<int DIM> 
StandardLocallyActiveDataCoarsenTransactionFactory<DIM>::
   ~StandardLocallyActiveDataCoarsenTransactionFactory()
{
}

/*
*************************************************************************
*                                                                       *
* Set/unset information for transactions managed by this factory class. *
*                                                                       *
*************************************************************************
*/

template<int DIM> void 
StandardLocallyActiveDataCoarsenTransactionFactory<DIM>::setCoarsenItems(
   const typename CoarsenClasses<DIM>::Data** coarsen_items,
   int num_coarsen_items)
{
   xfer::CoarsenCopyTransaction<DIM>::setCoarsenItems(coarsen_items,
                                                      num_coarsen_items);
   d_coarsen_items = coarsen_items;
   d_num_coarsen_items = num_coarsen_items;
}

template<int DIM> void 
StandardLocallyActiveDataCoarsenTransactionFactory<DIM>::unsetCoarsenItems()
{
   xfer::CoarsenCopyTransaction<DIM>::unsetCoarsenItems();
   d_coarsen_items = (const typename xfer::CoarsenClasses<DIM>::Data**)NULL;
   d_num_coarsen_items = 0;
}

/*
*************************************************************************
*                                                                       *
* Allocate appropriate transaction object.                              *
*                                                                       *
*************************************************************************
*/

template<int DIM>
tbox::Pointer<tbox::Transaction>
StandardLocallyActiveDataCoarsenTransactionFactory<DIM>::allocate(
   tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer< hier::PatchLevel<DIM> > src_level,
   tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
   int dst_patch_id,
   int src_patch_id,
   int citem_id,
   tbox::Pointer<tbox::Arena> pool) const
{

   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   CoarsenCopyTransaction<DIM>* transaction =
      new (pool) CoarsenCopyTransaction<DIM>(dst_level, src_level,
                                             overlap,
                                             dst_patch_id,
                                             src_patch_id,
                                             citem_id);
   return(tbox::Pointer<tbox::Transaction>(transaction, pool));
}

}
}
#endif
