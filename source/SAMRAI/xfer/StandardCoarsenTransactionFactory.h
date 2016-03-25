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

#ifndef included_xfer_StandardCoarsenTransactionFactory
#define included_xfer_StandardCoarsenTransactionFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Transaction.h"
#include "SAMRAI/xfer/CoarsenClasses.h"
#include "SAMRAI/xfer/CoarsenTransactionFactory.h"

namespace SAMRAI {
namespace xfer {

/*!
 * @brief Concrete subclass of CoarsenTransactionFactory base class that allocates
 * CoarsenCopyTransaction objects for a CoarsenSchedule object.
 *
 * @see xfer::CoarsenCopyTransaction
 */

class StandardCoarsenTransactionFactory:public CoarsenTransactionFactory
{
public:
   /*!
    * @brief Default constructor.
    */
   StandardCoarsenTransactionFactory();

   /*!
    * @brief Virtual destructor.
    */
   virtual ~StandardCoarsenTransactionFactory();

   /*!
    * @brief Set the array of CoarsenClass::Data items used by the transactions.
    */
   void
   setCoarsenItems(
      const CoarsenClasses::Data ** coarsen_items,
      int num_coarsen_items);

   /*!
    * @brief Clear the array of CoarsenClass::Data items used by the transactions.
    */
   void
   unsetCoarsenItems();

   /*!
    * @brief Allocate a CoarsenCopyTransaction object.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param citem_id       Integer index of CoarsenClass::Data item associated
    *                       with transaction.
    */
   tbox::Pointer<tbox::Transaction>
   allocate(
      tbox::Pointer<hier::PatchLevel> dst_level,
      tbox::Pointer<hier::PatchLevel> src_level,
      tbox::Pointer<hier::BoxOverlap> overlap,
      const hier::Box& dst_mapped_box,
      const hier::Box& src_mapped_box,
      int citem_id) const;

private:
   // The following two functions are not implemented
   StandardCoarsenTransactionFactory(
      const StandardCoarsenTransactionFactory&);
   void
   operator = (
      const StandardCoarsenTransactionFactory&);

   const xfer::CoarsenClasses::Data** d_coarsen_items;
   int d_num_coarsen_items;

};

}
}
#endif
