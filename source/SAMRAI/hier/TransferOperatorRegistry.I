/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Singleton registry for all tranfer operators.
 *
 ************************************************************************/

namespace SAMRAI {
namespace hier {

SAMRAI_INLINE_KEYWORD
void
TransferOperatorRegistry::addTimeInterpolateOperator(
   const boost::shared_ptr<TimeInterpolateOperator>& time_op)
{
   d_time_operators.addItem(time_op);
}

/*
 *************************************************************************
 * Set the mininum value to be returned by
 * getMaxTransferOpStencilWidth().
 *************************************************************************
 */

SAMRAI_INLINE_KEYWORD
void
TransferOperatorRegistry::setMinTransferOpStencilWidth(
   const IntVector& min_width)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_min_stencil_width, min_width);
   d_min_stencil_width = min_width;
}

/*
 *************************************************************************
 * Get the dimension of the GridGeometry object holding this singleton.
 *************************************************************************
 */
SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
TransferOperatorRegistry::getDim() const
{
   return d_dim;
}

}
}
