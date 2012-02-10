/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface for factory objects that create transactions for
 *                refine schedules.
 *
 ************************************************************************/

#ifndef included_xfer_RefineTransactionFactory_C
#define included_xfer_RefineTransactionFactory_C

#include "SAMRAI/xfer/RefineTransactionFactory.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

RefineTransactionFactory::RefineTransactionFactory()
{
}

RefineTransactionFactory::~RefineTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Invoke the allocate method with default values for box and
 * use_time_interpolation. Can't use optional argument as a Dimension
 * value is needed.
 *
 *************************************************************************
 */
boost::shared_ptr<tbox::Transaction>
RefineTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   int ritem_id) const {
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);

   return allocate(
      dst_level,
      src_level,
      overlap,
      dst_mapped_box,
      src_mapped_box,
      ritem_id,
      hier::Box::getEmptyBox(src_level->getDim()),
      false);
}

/*
 *************************************************************************
 *
 * Default no-op implementations of optional virtual functions.
 *
 *************************************************************************
 */

void RefineTransactionFactory::setTransactionTime(
   double fill_time)
{
   NULL_USE(fill_time);
}

void RefineTransactionFactory::preprocessScratchSpace(
   const boost::shared_ptr<hier::PatchLevel>& level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(level);
   NULL_USE(fill_time);
   NULL_USE(preprocess_vector);
}

}
}
#endif
