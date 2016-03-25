//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/datamovers/locally_active/LocallyActiveDataRefineTransactionFactory.h $
// Package:	SAMRAI transfer
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Interface for factory objects that create transactions for
//              locally-active data refine schedules. 
//

#ifndef included_xfer_LocallyActiveDataRefineTransactionFactory
#define included_xfer_LocallyActiveDataRefineTransactionFactory

#include "SAMRAI_config.h"
#include "Box.h"
#include "BoxOverlap.h"
#include "PatchLevel.h"
#include "LocallyActiveDataPatchLevelManager.h"
#include "tbox/Arena.h"
#include "tbox/Pointer.h"
#include "tbox/DescribedClass.h"
#include "tbox/Transaction.h"
#include "RefineClasses.h"


namespace SAMRAI {
    namespace xfer {

/*!
 * @brief Abstract base class defining the interface for all concrete transaction 
 * factory objects that generate data transaction objects used with a 
 * LocallyActiveDataRefineSchedule<DIM> object.  A concrete subclass will allocate 
 * new transaction objects.  This class is an example of the ``Abstract Factory'' 
 * method described in the Design Patterns book by Gamma, et al.
 *
 * To add a new type of Transaction object MyRefineTransaction:
 *
 * -# Implement a concrete LocallyActiveDataRefineTransactionFactory<DIM> object as 
 *       a subclass of this LocallyActiveDataRefineTransactionFactory<DIM> base class.
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
class LocallyActiveDataRefineTransactionFactory : public tbox::DescribedClass
{
public:
   /*!
    * @brief Default constructor.
    */
   LocallyActiveDataRefineTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~LocallyActiveDataRefineTransactionFactory();

   /*!
    * @brief Pure virtual function to set the array of RefineClass::Data items 
    * associated with the refine schedule.  Typical concrete transactions used by
    * the schedule use this information to communicate data.  This operation
    * is called by the refine schedule during the execution of the 
    * LocallyActiveDataRefineSchedule<DIM>::fillData() routine before 
    * data communication operations begin.
    */
   virtual void setRefineItems(
      const typename RefineClasses<DIM>::Data** refine_items,
      int num_refine_items) = 0;

   /*!
    * @brief Pure virtual function to clear the array of RefineClass::Data items 
    * associated with the refine schedule.  This operation is called by the 
    * refine schedule after data communication operations are complete.
    */
   virtual void unsetRefineItems() = 0;

   /*!
    * @brief Pure virtual function to allocate a concrete refine transaction object.
    * This routine is called by the refine schedule during construction of the 
    * schedule.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param ritem_id       Integer index of RefineClass::Data item associated
    *                       with transaction.
    * @param box            Optional const reference to box defining region of 
    *                       refine transaction.  Default is an empty box.
    * @param use_time_interpolation  Optional boolean flag indicating whether the
    *                       refine transaction involves time interpolation.
    *                       Default is false.
    * @param pool           Optional pointer to memory pool from which the
    *                       refine transaction may be allocated. Default is null.
    */
   virtual tbox::Pointer<tbox::Transaction>
   allocate(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
            tbox::Pointer< hier::PatchLevel<DIM> > src_level,
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
            int dst_patch_id,
            int src_patch_id,
            int ritem_id,
            const hier::Box<DIM>& box = hier::Box<DIM>(),
            bool use_time_interpolation = false,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const = 0;

   /*!
    * @brief Virtual function to set simulation time for transaction objects.
    * This operation is called by the refine schedule during the execution of 
    * the LocallyActiveDataRefineSchedule<DIM>::fillData() routine before data 
    * communication operations begin.  This function is optional for the concrete 
    * transaction factory object.  The default implementation is a no-op.
    */
   virtual void setTransactionTime(double fill_time);

   /*!
    * @brief Virtual function allowing transaction factory to preprocess scratch 
    * space data before transactactions use it if they need to.  This function is 
    * optional for the concrete transaction factory object.  
    * The default implementation is a no-op.
    *
    * @param level        tbox::Pointer to patch level holding scratch data.
    * @param fill_time    Double value of simulation time corresponding to 
    *                     LocallyActiveDataRefineSchedule<DIM> operations.
    * @param preprocess_mgr Const reference to LocallyActiveDataPatchLevelManager<DIM>
    *                     objects that indicates patch data array indices of which scratch 
    *                     patch data objects on each patch to preprocess.
    */
   virtual void preprocessScratchSpace(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      double fill_time,
      const hier::LocallyActiveDataPatchLevelManager<DIM>& preprocess_mgr) const;

private:
   // The following two functions are not implemented
   LocallyActiveDataRefineTransactionFactory(
      const LocallyActiveDataRefineTransactionFactory<DIM>&);
   void operator=(const LocallyActiveDataRefineTransactionFactory<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "LocallyActiveDataRefineTransactionFactory.C"
#endif
