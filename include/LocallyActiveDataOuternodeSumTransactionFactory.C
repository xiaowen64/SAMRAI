//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/femutils/locally_active/LocallyActiveDataOuternodeSumTransactionFactory.C $
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Concrete factory for creating outernode sum transaction objects
//              for locally-active data refine schedules.
//

#ifndef included_algs_LocallyActiveDataOuternodeSumTransactionFactory_C
#define included_algs_LocallyActiveDataOuternodeSumTransactionFactory_C

#include "LocallyActiveDataOuternodeSumTransactionFactory.h"
#include "OuternodeSumTransaction.h"
#include "tbox/ArenaManager.h"

namespace SAMRAI {
    namespace algs {

/*
*************************************************************************
*                                                                       *
* Default constructor and destructor.                                   *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
LocallyActiveDataOuternodeSumTransactionFactory<DIM>::
LocallyActiveDataOuternodeSumTransactionFactory()
{
}

template<int DIM> 
LocallyActiveDataOuternodeSumTransactionFactory<DIM>::
~LocallyActiveDataOuternodeSumTransactionFactory()
{
}

/*
*************************************************************************
*                                                                       *
* Set/unset information for transactions managed by this factory class. *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataOuternodeSumTransactionFactory<DIM>::setRefineItems(
   const typename xfer::RefineClasses<DIM>::Data** refine_items,
   int num_refine_items)
{
   algs::OuternodeSumTransaction<DIM>::setRefineItems(refine_items,
                                                      num_refine_items);
   d_refine_items = refine_items;
   d_number_refine_items = num_refine_items;
}

template<int DIM>
void LocallyActiveDataOuternodeSumTransactionFactory<DIM>::unsetRefineItems()
{
   algs::OuternodeSumTransaction<DIM>::unsetRefineItems();
   d_refine_items = (const typename xfer::RefineClasses<DIM>::Data**)NULL;
   d_number_refine_items = 0;
}

/*
*************************************************************************
*                                                                       *
* Allocate outernode sum transaction object.                            *
*                                                                       *
*************************************************************************
*/

template<int DIM> 
tbox::Pointer<tbox::Transaction>
LocallyActiveDataOuternodeSumTransactionFactory<DIM>::allocate(
   tbox::Pointer<hier::PatchLevel<DIM> > dst_level,
   tbox::Pointer<hier::PatchLevel<DIM> > src_level,
   tbox::Pointer<hier::BoxOverlap<DIM> > overlap,
   int dst_patch_id,
   int src_patch_id,
   int ritem_id,
   const hier::Box<DIM>& box,
   bool use_time_interpolation, 
   tbox::Pointer<tbox::Arena> pool ) const
{
   (void) box;
   (void) use_time_interpolation;

   if (pool.isNull()) {
      pool = tbox::ArenaManager::getManager()->getStandardAllocator();
   }

   OuternodeSumTransaction<DIM> *transaction =
      new (pool) OuternodeSumTransaction<DIM>(dst_level,
                                              src_level,
                                              overlap,
                                              dst_patch_id,
                                              src_patch_id,
                                              ritem_id);
   return(tbox::Pointer<tbox::Transaction>(transaction, pool));
}

/*
*************************************************************************
*                                                                       *
* Initialize (to 0.0) scratch storage for sum transactions.             *
*                                                                       *
*************************************************************************
*/

template<int DIM>
void LocallyActiveDataOuternodeSumTransactionFactory<DIM>::preprocessScratchSpace(
   tbox::Pointer< hier::PatchLevel<DIM> > level,
   double fill_time,
   const hier::LocallyActiveDataPatchLevelManager<DIM>& preprocess_mgr) const
{
   (void) fill_time;
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!level.isNull());
   TBOX_ASSERT(d_refine_items != (const typename xfer::RefineClasses<DIM>::Data**)NULL);
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++ ) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      for (int iri = 0; iri < d_number_refine_items; iri++) {
         const int data_id = d_refine_items[iri]->d_scratch;
         if ( preprocess_mgr.getPatchDataActive( hier::PatchDataId(data_id), 
                                                 hier::PatchNumber(ip()) ) ) {
            tbox::Pointer< pdat::OuternodeData<DIM,double> > onode_data =
               patch->getPatchData(data_id);
            onode_data->fillAll(0.0);
         }
      }
   }
}

}
}
#endif

