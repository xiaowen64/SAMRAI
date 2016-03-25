//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/standard/CoarsenTransactionFactory.h $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Interface for factory objects that create transactions for
//              coarsen schedules. 
//

#ifndef included_xfer_CoarsenTransactionFactory
#define included_xfer_CoarsenTransactionFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxOverlap.h"
#include "PatchLevel.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"
#include "tbox/DescribedClass.h"
#include "tbox/Transaction.h"
#include "CoarsenClasses.h"


namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Abstract base class defining the interface for all concrete transaction 
 * factory objects that generate data transaction objects used with a CoarsenSchedule<DIM> 
 * object.  A concrete subclass will allocate new transaction objects.  This class 
 * is an example of the ``Abstract Factory'' method described in the Design Patterns 
 * book by Gamma, et al.
 *
 * To add a new type of Transaction object MyCoarsenTransaction:
 *
 * -# Implement a concrete CoarsenTransactionFactory<DIM> object as a subclass
 *       that is derived from this CoarsenTransactionFactory<DIM> base class.
 *       Implement the abstract virtual functions as appropriate for the 
 *       concrete subclass; in particular, the allocate() function must return 
 *       a new instance of the desired transaction object.
 * -# The type of the transaction allocated by the concrete factory is a 
 *       Transaction<DIM>.  Thus, the new transaction object must be derived
 *       from the Transaction<DIM> base class and implement the abstract 
 *       virtual functions declared by the base class.
 * 
 * @see tbox::Transaction
 */

template<int DIM> 
class CoarsenTransactionFactory : public tbox::DescribedClass
{
public:
   /*!
    * @brief Default constructor.
    */
   CoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~CoarsenTransactionFactory();

   /*!
    * @brief Pure virtual function to set the array of CoarsenClass::Data items 
    * associated with the coarsen schedule.  Typical concrete transactions used by
    * the schedule use this information to communicate data.  This operation
    * is called by the coarsen schedule during the execution of the 
    * CoarsenSchedule<DIM>::fillData() routine before data communication 
    * operations begin.
    */
   virtual void setCoarsenItems(
      const typename CoarsenClasses<DIM>::Data** coarsen_items,
      int num_coarsen_items) = 0;

   /*!
    * @brief Pure virtual function to clear the array of CoarsenClass::Data items 
    * associated with the coarsen schedule.  This operation is called by the 
    * coarsen schedule after data communication operations are complete.
    */
   virtual void unsetCoarsenItems() = 0;

   /*!
    * @brief Pure virtual function to allocate a concrete coarsen transaction object.
    * This routine is called by the coarsen schedule during construction of the 
    * schedule.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param citem_id       Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    * @param pool           Optional pointer to memory pool from which the
    *                       refine transaction may be allocated. Default is null.
    */
   virtual tbox::Pointer<tbox::Transaction>
   allocate(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
            tbox::Pointer< hier::PatchLevel<DIM> > src_level,
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
            int dst_patch_id,
            int src_patch_id,
            int citem_id,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const = 0;

private:
   CoarsenTransactionFactory(const CoarsenTransactionFactory<DIM>&); // not implemented
   void operator=(const CoarsenTransactionFactory<DIM>&);	  // not implemented

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CoarsenTransactionFactory.C"
#endif
