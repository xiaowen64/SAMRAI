//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/algorithm/femutils/standard/OuteredgeSumTransactionFactory.h $
// Package:	SAMRAI algorithms
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1808 $
// Modified:	$LastChangedDate: 2007-12-19 16:38:32 -0800 (Wed, 19 Dec 2007) $
// Description:	Factory for creating outeredge sum transaction objects
//

#ifndef included_algs_OuteredgeSumTransactionFactory
#define included_algs_OuteredgeSumTransactionFactory

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Arena
#include "tbox/Arena.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_hier_ComponentSelector
#include "ComponentSelector.h"
#endif
#ifndef included_hier_PatchLevel
#include "PatchLevel.h"
#endif
#ifndef included_xfer_RefineClasses
#include "RefineClasses.h"
#endif
#ifndef included_xfer_RefineTransactionFactory
#include "RefineTransactionFactory.h"
#endif

#ifndef NULL
#define NULL (0)
#endif

namespace SAMRAI {
    namespace algs {

/*!
 * @brief Concrete subclass of the xfer::RefineTransactionFactory<DIM> base class that
 * allocates outeredge sum transaction objects for a xfer::RefineSchedule<DIM> object.
 *
 * @see xfer::RefineTransactionFactory
 * @see xfer::OuteredgeSumTransaction
 */

template<int DIM> 
class OuteredgeSumTransactionFactory : public xfer::RefineTransactionFactory<DIM>
{
public:
   /*!
    * @brief Default constructor.
    */
   OuteredgeSumTransactionFactory();

   /*!
    * @brief Virtual destructor for base class.
    */
   virtual ~OuteredgeSumTransactionFactory();

   /*!
    * @brief Set the array of xfer::RefineClass<DIM>::Data items used by the transactions.
    */
   void setRefineItems(const typename xfer::RefineClasses<DIM>::Data** refine_items,
                       int num_refine_items);

   /*!
    * @brief Clear the array of xfer::RefineClass<DIM>::Data items used by the transactions.
    */
   void unsetRefineItems();

   /*!
    * @brief Allocate an OuteredgeSumTransaction<DIM> object.
    *
    * @param dst_level      tbox::Pointer to destination patch level.
    * @param src_level      tbox::Pointer to source patch level.
    * @param overlap        tbox::Pointer to overlap region between patches.
    * @param dst_patch_id   Integer index of destination patch in destination
    *                       patch level.
    * @param src_patch_id   Integer index of source patch in source patch level.
    * @param ritem_id       Integer index of xfer::RefineClass<DIM>::Data item 
    *                       associated with transaction.
    * @param box            Optional const reference to box defining region of
    *                       refine transaction.  Default is an empty box.
    * @param use_time_interpolation  Optional boolean flag indicating whether the
    *                       refine transaction involves time interpolation.
    *                       Default is false.
    * @param pool           Optional pointer to memory pool from which the
    *                       refine transaction may be allocated. Default is null.
    */
   tbox::Pointer<tbox::Transaction>
   allocate(tbox::Pointer< hier::PatchLevel<DIM> > dst_level,
            tbox::Pointer< hier::PatchLevel<DIM> > src_level,
            tbox::Pointer< hier::BoxOverlap<DIM> > overlap,
            int dst_patch_id,
            int src_patch_id,
            int ritem_id,
            const hier::Box<DIM>& box = hier::Box<DIM>(),
            bool use_time_interpolation = false,
            tbox::Pointer<tbox::Arena> pool = (tbox::Arena *) NULL) const;

   /*!
    * @brief Function to initialize scratch space data for the sum transactions
    * (patch data components indicated by the component selector) to zero.
    *
    * @param level        tbox::Pointer to patch level holding scratch data.
    * @param fill_time    Double value of simulation time at which preprocess
    *                     operation is called.
    * @param preprocess_vector Const reference to hier::ComponentSelector indicating
    *                     patch data array indices of scratch patch data objects
    *                     to preprocess.
    */
   void preprocessScratchSpace(
      tbox::Pointer< hier::PatchLevel<DIM> > level,
      double fill_time,
      const hier::ComponentSelector& preprocess_vector) const;

private:
   // The following two functions are not implemented
   OuteredgeSumTransactionFactory<DIM>(const OuteredgeSumTransactionFactory<DIM>&);
   void operator=(const OuteredgeSumTransactionFactory<DIM>&);

   const typename xfer::RefineClasses<DIM>::Data** d_refine_items;
   int d_number_refine_items;

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OuteredgeSumTransactionFactory.C"
#endif



