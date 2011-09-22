/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Factory for creating outeredge sum transaction objects
 *
 ************************************************************************/

#ifndef included_algs_OuteredgeSumTransactionFactory_C
#define included_algs_OuteredgeSumTransactionFactory_C

#include "SAMRAI/algs/OuteredgeSumTransactionFactory.h"

#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/algs/OuteredgeSumTransaction.h"

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

OuteredgeSumTransactionFactory::OuteredgeSumTransactionFactory()
{
}

OuteredgeSumTransactionFactory::~OuteredgeSumTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Set/unset information for transactions managed by this factory class.
 *
 *************************************************************************
 */

void OuteredgeSumTransactionFactory::setRefineItems(
   const xfer::RefineClasses::Data** refine_items,
   int num_refine_items)
{
   algs::OuteredgeSumTransaction::setRefineItems(refine_items,
      num_refine_items);
   d_refine_items = refine_items;
   d_number_refine_items = num_refine_items;
}

void OuteredgeSumTransactionFactory::unsetRefineItems()
{
   algs::OuteredgeSumTransaction::unsetRefineItems();
   d_refine_items = (const xfer::RefineClasses::Data **)NULL;
   d_number_refine_items = 0;
}

/*
 *************************************************************************
 *
 * Allocate outeredge sum transaction object.
 *
 *************************************************************************
 */

tbox::Pointer<tbox::Transaction>
OuteredgeSumTransactionFactory::allocate(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   tbox::Pointer<hier::BoxOverlap> overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   int ritem_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   NULL_USE(box);
   NULL_USE(use_time_interpolation);

   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level, *src_level, dst_node, src_node);

   OuteredgeSumTransaction* transaction =
      new OuteredgeSumTransaction(dst_level,
         src_level,
         overlap,
         dst_node,
         src_node,
         ritem_id);
   return tbox::Pointer<tbox::Transaction>(transaction);
}

tbox::Pointer<tbox::Transaction>
OuteredgeSumTransactionFactory::allocate(
   tbox::Pointer<hier::PatchLevel> dst_level,
   tbox::Pointer<hier::PatchLevel> src_level,
   tbox::Pointer<hier::BoxOverlap> overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   int ritem_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level, *src_level, dst_node, src_node);

   return allocate(dst_level,
      src_level,
      overlap,
      dst_node,
      src_node,
      ritem_id,
      hier::Box::getEmptyBox(dst_level->getDim()),
      false);
}

/*
 *************************************************************************
 *
 * Initialize (to 0.0) scratch storage for sum transactions.
 *
 *************************************************************************
 */

void OuteredgeSumTransactionFactory::preprocessScratchSpace(
   tbox::Pointer<hier::PatchLevel> level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(fill_time);
   TBOX_ASSERT(!level.isNull());

   for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
      tbox::Pointer<hier::Patch> patch = *ip;

      const int ncomponents = preprocess_vector.getSize();
      for (int n = 0; n < ncomponents; ++n) {
         if (preprocess_vector.isSet(n)) {
            tbox::Pointer<pdat::OuteredgeData<double> > oedge_data =
               patch->getPatchData(n);
            oedge_data->fillAll(0.0);
         }
      }

   }
}

}
}
#endif
