//
// File:        OuternodeSumTransactionFactory.C
// Package:     SAMRAI algorithms
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 696 $
// Modified:    $Date: 2005-11-03 12:27:01 -0800 (Thu, 03 Nov 2005) $
// Description: Factory for creating outernode sum transaction objects
//

#ifndef included_algs_OuternodeSumTransactionFactory_C
#define included_algs_OuternodeSumTransactionFactory_C

#include "OuternodeSumTransactionFactory.h"
#include "OuternodeSumTransaction.h"
#include "tbox/ArenaManager.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

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
OuternodeSumTransactionFactory<DIM>::OuternodeSumTransactionFactory()
{
}

template<int DIM> 
OuternodeSumTransactionFactory<DIM>::~OuternodeSumTransactionFactory()
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
void OuternodeSumTransactionFactory<DIM>::setRefineItems(
   const typename xfer::RefineClasses<DIM>::Data** refine_items,
   int num_refine_items)
{
   algs::OuternodeSumTransaction<DIM>::setRefineItems(refine_items,
                                                      num_refine_items);
   d_refine_items = refine_items;
   d_number_refine_items = num_refine_items;
}

template<int DIM>
void OuternodeSumTransactionFactory<DIM>::unsetRefineItems()
{
   OuternodeSumTransaction<DIM>::unsetRefineItems();
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
OuternodeSumTransactionFactory<DIM>::allocate(
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
void OuternodeSumTransactionFactory<DIM>::preprocessScratchSpace(
   tbox::Pointer<hier::PatchLevel<DIM> > level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   (void) fill_time;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!level.isNull());
#endif

   for (typename hier::PatchLevel<DIM>::Iterator ip(level); ip; ip++ ) {
      tbox::Pointer< hier::Patch<DIM> > patch = level->getPatch(ip());

      const int ncomponents = preprocess_vector.getSize();
      for (int n=0; n<ncomponents; ++n ) {
         if ( preprocess_vector.isSet(n) ) {
            tbox::Pointer<pdat::OuternodeData<DIM, double> > onode_data = 
               patch->getPatchData(n);
            onode_data->fillAll(0.0);
         }
      }

   }
}

}
}
#endif

