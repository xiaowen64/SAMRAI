/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic integration scheme.
 *
 ************************************************************************/

namespace SAMRAI {
namespace algs {

SAMRAI_INLINE_KEYWORD
boost::shared_ptr<hier::VariableContext>
HyperbolicPatchStrategy::getDataContext()
const
{
   return d_data_context;
}

SAMRAI_INLINE_KEYWORD
void
HyperbolicPatchStrategy::setDataContext(
   const boost::shared_ptr<hier::VariableContext>& context)
{
   d_data_context = context;
}

SAMRAI_INLINE_KEYWORD
void
HyperbolicPatchStrategy::clearDataContext()
{
   d_data_context.reset();
}

SAMRAI_INLINE_KEYWORD
const tbox::Dimension&
HyperbolicPatchStrategy::getDim() const
{
   return d_dim;
}

}
}
